"""Leakage-conscious, single-label multiclass research training.

Run as python -m src.robust_variant_training --help. See the companion research
guide for schema, label-availability semantics, and limits of internal validation.
Only explicitly allowlisted columns become features. IDs only define split groups.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.feature_extraction import DictVectorizer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, classification_report, confusion_matrix
from sklearn.model_selection import ParameterSampler, StratifiedGroupKFold, TimeSeriesSplit
from sklearn.preprocessing import LabelEncoder, MaxAbsScaler
from sklearn.utils.class_weight import compute_sample_weight


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metrics(y, pred, classes):
    report = classification_report(y, pred, labels=classes, output_dict=True, zero_division=0)
    return {
        "n": len(y), "accuracy": float(accuracy_score(y, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "macro_precision": report["macro avg"]["precision"],
        "macro_recall": report["macro avg"]["recall"],
        "macro_f1": report["macro avg"]["f1-score"],
        "classification_report": report, "class_order": list(classes),
        "confusion_matrix": confusion_matrix(y, pred, labels=classes).tolist(),
    }


def objective(y, pred):
    # Lexicographic selection: macro-F1 breaks balanced-accuracy ties only.
    m = metrics(y, pred, np.unique(y))
    return m["balanced_accuracy"], m["macro_f1"]


def validate_schema(frame, schema):
    features = schema.get("numeric", []) + schema.get("categorical", []) + schema.get("tokens", [])
    reserved = {schema["target"], schema["group"], schema.get("time"), schema.get("label_time")}
    if not features or len(features) != len(set(features)) or reserved.intersection(features):
        raise ValueError("Use nonempty, distinct feature lists excluding target, IDs and timestamps.")
    required = set(features) | {schema["target"], schema["group"]}
    if schema.get("time"):
        if not schema.get("label_time"):
            raise ValueError("Temporal training requires label_time (when the outcome became available).")
        required.update([schema["time"], schema["label_time"]])
    absent = required - set(frame.columns)
    if absent:
        raise ValueError(f"Missing columns: {sorted(absent)}")
    if frame[[schema["target"], schema["group"]]].isna().any().any():
        raise ValueError("Missing targets/groups are not allowed.")
    if frame[schema["target"]].nunique() < 2:
        raise ValueError("At least two mutually exclusive target classes are required.")
    for col in (schema["target"], schema["group"]):
        if frame[col].astype(str).str.strip().eq("").any():
            raise ValueError(f"Empty {col} values are not allowed.")
    if schema.get("time"):
        for col in (schema["time"], schema["label_time"]):
            frame[col] = pd.to_datetime(frame[col], utc=True, errors="raise")
            if frame[col].isna().any():
                raise ValueError(f"Missing dates in {col}")
        if (frame[schema["label_time"]] < frame[schema["time"]]).any():
            raise ValueError("Outcome availability cannot precede the prediction timestamp.")


def add_history(frame, history, schema):
    """Attach external, strictly earlier observations; never shift the target.

    history columns: group column, available_at, observed_label, evidence_count.
    Each row is a cumulative snapshot, unique per entity and availability time.
    The caller must supply genuine historical records, not backdated current data.
    """
    if not schema.get("time"):
        raise ValueError("History requires prediction timestamps.")
    group, time = schema["group"], schema["time"]
    history = history.copy()
    history["available_at"] = pd.to_datetime(history["available_at"], utc=True, errors="raise")
    if history[[group, "available_at", "observed_label", "evidence_count"]].isna().any().any():
        raise ValueError("History fields must be nonmissing.")
    if history.duplicated([group, "available_at"]).any():
        raise ValueError("Ambiguous duplicate history snapshots.")
    history["evidence_count"] = pd.to_numeric(history["evidence_count"], errors="raise")
    if not np.isfinite(history["evidence_count"]).all() or (history["evidence_count"] < 0).any():
        raise ValueError("Evidence counts must be finite and nonnegative.")
    history = history.sort_values([group, "available_at"])
    changed = history["observed_label"].ne(history.groupby(group)["observed_label"].shift())
    history["last_change"] = history["available_at"].where(changed).groupby(history[group]).ffill()
    history["history_evidence_delta"] = history.groupby(group)["evidence_count"].diff().fillna(0)
    history = history.rename(columns={"observed_label": "history_previous_label", "evidence_count": "history_evidence_count"})
    generated = ["history_previous_label", "history_evidence_count", "history_evidence_delta", "history_days_since_change"]
    if set(generated).intersection(frame.columns):
        raise ValueError("Input already contains generated history columns.")
    left = frame.copy()
    left["_original_order"] = np.arange(len(left))
    merged = pd.merge_asof(left.sort_values(time), history.sort_values("available_at"),
                           left_on=time, right_on="available_at", by=group,
                           direction="backward", allow_exact_matches=False)
    merged["history_days_since_change"] = (merged[time] - merged["last_change"]).dt.total_seconds() / 86400
    return merged.sort_values("_original_order").drop(columns=["_original_order", "available_at", "last_change"]).reset_index(drop=True)


def records(frame, schema):
    """Sparse features, explicit missingness; no fitting outside training folds."""
    output = []
    for row in frame.to_dict("records"):
        item = {}
        for col in schema.get("numeric", []):
            value = pd.to_numeric(row.get(col), errors="coerce")
            missing = pd.isna(value) or not np.isfinite(value)
            item[col] = 0.0 if missing else float(value)
            item[col + "__missing"] = float(missing)
        for col in schema.get("categorical", []):
            value = row.get(col)
            item[col] = "__MISSING__" if pd.isna(value) else str(value)
        for col in schema.get("tokens", []):
            value = row.get(col)
            terms = [] if pd.isna(value) else str(value).split("|")
            for token in sorted({term.strip() for term in terms if term.strip()}):
                item[col + "=" + token] = 1.0
        output.append(item)
    return output


class VariantEnsemble(ClassifierMixin, BaseEstimator):
    """Persisted predictor includes preprocessing and multiclass decision offsets."""

    def __init__(self, schema, params=None, seed=42, backend="sklearn", sampler="none"):
        self.schema = schema
        self.params = params
        self.seed = seed
        self.backend = backend
        self.sampler = sampler

    def fit(self, X, y):
        params = self.params or {}
        self.encoder_ = LabelEncoder().fit(y)
        self.classes_ = self.encoder_.classes_
        labels = self.encoder_.transform(y)
        self.vectorizer_ = DictVectorizer(dtype=np.float64)
        self.scaler_ = MaxAbsScaler()
        matrix = self.scaler_.fit_transform(self.vectorizer_.fit_transform(records(X, self.schema)))
        if self.sampler != "none":
            if self.schema.get("categorical") or self.schema.get("tokens"):
                raise ValueError("SMOTE variants are restricted to numeric schemas; use weights for HPO/categorical inputs.")
            from imblearn.over_sampling import SMOTE, BorderlineSMOTE, ADASYN
            n = int(np.bincount(labels).min())
            if n < 2:
                raise ValueError("SMOTE requires at least two examples per class in every training fold.")
            cls = {"smote": SMOTE, "borderline": BorderlineSMOTE, "adasyn": ADASYN}[self.sampler]
            kwargs = {"n_neighbors" if self.sampler == "adasyn" else "k_neighbors": min(5, n - 1)}
            matrix, labels = cls(random_state=self.seed, **kwargs).fit_resample(matrix, labels)
        weights = compute_sample_weight("balanced", labels) ** params.get("weight_power", 1.0)
        weights /= weights.mean()
        trees = params.get("n_estimators", 100)
        common = dict(n_estimators=trees, max_depth=params.get("max_depth", 12),
                      min_samples_leaf=params.get("min_samples_leaf", 3),
                      max_features="sqrt", n_jobs=1, random_state=self.seed)
        models = [RandomForestClassifier(**common), ExtraTreesClassifier(**common)]
        if self.backend == "lightgbm":
            from lightgbm import LGBMClassifier
            models[1] = LGBMClassifier(n_estimators=trees, num_leaves=15, learning_rate=.05,
                                       reg_lambda=5, verbosity=-1, n_jobs=1, random_state=self.seed)
        elif self.backend == "xgboost":
            from xgboost import XGBClassifier
            models[1] = XGBClassifier(n_estimators=trees, max_depth=4, learning_rate=.05,
                                      reg_lambda=5, n_jobs=1, random_state=self.seed)
        elif self.backend == "catboost":
            from catboost import CatBoostClassifier
            models[1] = CatBoostClassifier(iterations=trees, depth=5, learning_rate=.05,
                                           l2_leaf_reg=5, thread_count=1, random_seed=self.seed,
                                           verbose=False, allow_writing_files=False)
        models.append(LogisticRegression(C=params.get("C", 1.0), solver="lbfgs", max_iter=1500))
        self.models_ = [model.fit(matrix, labels, sample_weight=weights) for model in models]
        self.offsets_ = np.zeros(len(self.classes_))
        return self

    def predict_proba(self, X):
        matrix = self.scaler_.transform(self.vectorizer_.transform(records(X, self.schema)))
        return np.mean([model.predict_proba(matrix) for model in self.models_], axis=0)

    def predict(self, X):
        scores = np.log(np.clip(self.predict_proba(X), 1e-12, 1)) + self.offsets_
        return self.classes_[scores.argmax(axis=1)]


def checked_split(frame, train, valid, schema):
    """Purge future/unavailable labels and overlapping loci before fitting."""
    group = schema["group"]
    train = np.asarray(train)
    valid = np.asarray(valid)
    train = train[~frame.iloc[train][group].isin(frame.iloc[valid][group]).to_numpy()]
    if schema.get("time"):
        cutoff = frame.iloc[valid][schema["time"]].min()
        train = train[(frame.iloc[train][schema["label_time"]] < cutoff).to_numpy()]
        if len(train) and frame.iloc[train][schema["time"]].max() >= cutoff:
            raise ValueError("Temporal boundary is not strictly chronological.")
    classes = set(frame[schema["target"]])
    if set(frame.iloc[train][schema["target"]]) != classes or set(frame.iloc[valid][schema["target"]]) != classes:
        raise ValueError("A split lacks classes after purging. Use more data, fewer folds or different prespecified time windows.")
    return train, valid


def folds(frame, schema, count, seed):
    if schema.get("time"):
        times = np.array(sorted(frame[schema["time"]].unique()))
        for a, b in TimeSeriesSplit(n_splits=count).split(times):
            yield checked_split(frame, np.flatnonzero(frame[schema["time"]].isin(times[a])),
                                np.flatnonzero(frame[schema["time"]].isin(times[b])), schema)
    else:
        split = StratifiedGroupKFold(n_splits=count, shuffle=True, random_state=seed)
        for a, b in split.split(frame, frame[schema["target"]], frame[schema["group"]]):
            yield checked_split(frame, a, b, schema)


def partitions(frame, schema, seed):
    if schema.get("time"):
        times = sorted(frame[schema["time"]].unique())
        if len(times) < 10:
            raise ValueError("Temporal evaluation requires at least ten distinct prediction times.")
        cut = times[int(len(times) * .8)]
        dev, test = checked_split(frame, np.flatnonzero(frame[schema["time"]] < cut),
                                  np.flatnonzero(frame[schema["time"]] >= cut), schema)
    else:
        dev, test = next(folds(frame, schema, 5, seed))
    development = frame.iloc[dev].reset_index(drop=True)
    if schema.get("time"):
        times = sorted(development[schema["time"]].unique())
        cut = times[int(len(times) * .8)]
        train, tune = checked_split(development, np.flatnonzero(development[schema["time"]] < cut),
                                    np.flatnonzero(development[schema["time"]] >= cut), schema)
    else:
        train, tune = next(folds(development, schema, 5, seed + 1))
    return dev[train], dev[tune], test


def tune_offsets(model, X, y):
    """Coordinate search for multiclass log-probability offsets on tuning data only.

    These are decision scores, not calibrated clinical probabilities. The final
    model is NOT refit after this step, so the selected score scale is preserved.
    """
    logp = np.log(np.clip(model.predict_proba(X), 1e-12, 1))
    offsets = np.zeros(len(model.classes_))
    best = objective(y, model.classes_[logp.argmax(axis=1)])
    for _ in range(3):
        for cls in range(1, len(offsets)):
            for value in np.linspace(-2, 2, 17):
                candidate = offsets.copy()
                candidate[cls] = value
                score = objective(y, model.classes_[(logp + candidate).argmax(axis=1)])
                if score > best:
                    best, offsets = score, candidate
    model.offsets_ = offsets
    return offsets.tolist()


def run(args):
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)  # No silent reuse of a test experiment.
    schema = json.loads(Path(args.schema).read_text(encoding="utf-8"))
    frame = pd.read_csv(args.data, dtype={schema["target"]: str, schema["group"]: str})
    validate_schema(frame, schema)
    if args.history:
        history = pd.read_csv(args.history, dtype={schema["group"]: str})
        frame = add_history(frame, history, schema)
        schema.setdefault("categorical", []).append("history_previous_label")
        schema.setdefault("numeric", []).extend(["history_evidence_count", "history_evidence_delta", "history_days_since_change"])
    train, tune, test = partitions(frame, schema, args.seed)
    target = schema["target"]
    X = frame.iloc[train].reset_index(drop=True)
    y = X[target].to_numpy()
    classes = np.unique(frame[target])
    splits = list(folds(X, schema, args.folds, args.seed))
    membership = pd.DataFrame({"row": np.arange(len(frame)), "partition": "purged"})
    for name, ids in [("train", train), ("decision_tuning", tune), ("test", test)]:
        membership.loc[ids, "partition"] = name
    membership.to_csv(out / "partitions.csv", index=False)
    audit = {"schema": schema, "seed": args.seed, "backend": args.backend, "sampler": args.sampler,
             "evaluation_status": args.evaluation_status,
             "data_sha256": file_hash(args.data),
             "history_sha256": file_hash(args.history) if args.history else None,
             "support": {name: frame.iloc[ids][target].value_counts().to_dict()
                         for name, ids in [("train", train), ("decision_tuning", tune), ("test", test)]},
             "versions": {pkg: importlib.metadata.version(pkg) for pkg in ["numpy", "pandas", "scikit-learn", "joblib"]}}
    for package in ([args.backend] if args.backend != "sklearn" else []) + (["imbalanced-learn"] if args.sampler != "none" else []):
        audit["versions"][package] = importlib.metadata.version(package)
    save_json(out / "audit.json", audit)
    space = {"weight_power": [.5, 1.0, 1.5], "max_depth": [8, 16, None],
             "min_samples_leaf": [2, 5, 10], "C": [.1, 1., 10.], "n_estimators": [args.trees]}
    candidates = [{"weight_power": 0., "max_depth": 12, "min_samples_leaf": 3,
                   "C": 1., "n_estimators": args.trees}]
    candidates += list(ParameterSampler(space, n_iter=args.trials, random_state=args.seed))
    experiments, best, best_params = [], (-1., -1.), None
    for number, params in enumerate(candidates):
        reports, truth, predictions = [], [], []
        for fold, (a, b) in enumerate(splits):
            model = VariantEnsemble(schema, params, args.seed, args.backend, args.sampler).fit(X.iloc[a], y[a])
            pred = model.predict(X.iloc[b])
            reports.append(metrics(y[b], pred, classes))
            truth.extend(y[b]); predictions.extend(pred)
        pooled = metrics(truth, predictions, classes)
        score = (float(np.mean([r["balanced_accuracy"] for r in reports])),
                 float(np.mean([r["macro_f1"] for r in reports])))
        experiments.append({"candidate": number, "params": params, "folds": reports,
                            "pooled_cv": pooled, "mean_fold_objective": list(score)})
        save_json(out / "development.json", experiments)
        print(f"Candidate {number}: CV balanced accuracy={score[0]:.4f}, macro-F1={score[1]:.4f}", flush=True)
        if score > best:
            best, best_params = score, params
    model = VariantEnsemble(schema, best_params, args.seed, args.backend, args.sampler).fit(X, y)
    tuning = frame.iloc[tune]
    before = metrics(tuning[target], model.predict(tuning), classes)
    offsets = tune_offsets(model, tuning, tuning[target].to_numpy())
    after = metrics(tuning[target], model.predict(tuning), classes)
    save_json(out / "decision_tuning.json", {"before": before, "after": after,
                                             "caution": "Optimization-set estimates, not unbiased performance."})
    # Permute original allowlisted columns on tuning data, never test data.
    feature_columns = schema.get("numeric", []) + schema.get("categorical", []) + schema.get("tokens", [])
    importance = permutation_importance(model, tuning[feature_columns], tuning[target],
                                        scoring="balanced_accuracy", n_repeats=3,
                                        random_state=args.seed, n_jobs=1)
    pd.DataFrame({"feature": feature_columns, "mean": importance.importances_mean,
                  "std": importance.importances_std}).sort_values("mean", ascending=False).to_csv(out / "feature_importance.csv", index=False)
    joblib.dump(model, out / "model.joblib", compress=3)
    save_json(out / "selection_frozen.json", {"params": best_params, "offsets": offsets,
                                               "cv_objective": list(best), "refit_after_tuning": False})
    heldout = frame.iloc[test]
    probabilities = model.predict_proba(heldout)
    final = {"evaluation_status": args.evaluation_status,
             "default_argmax": metrics(heldout[target], classes[probabilities.argmax(axis=1)], classes),
             "selected_policy": metrics(heldout[target], model.predict(heldout), classes)}
    save_json(out / "test_metrics.json", final)
    print(json.dumps(final, indent=2), flush=True)
    return model, final


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, help="CSV or CSV.gz with one mutually exclusive label per row")
    parser.add_argument("--schema", required=True, help="JSON allowlist; see research/robust_variant_schema.json")
    parser.add_argument("--output", required=True, help="New directory; existing paths are rejected")
    parser.add_argument("--history", help="Independent timestamped snapshot CSV")
    parser.add_argument("--evaluation-status", required=True, choices=["exploratory_reused_source", "new_prespecified_holdout"])
    parser.add_argument("--backend", choices=["sklearn", "lightgbm", "xgboost", "catboost"], default="sklearn")
    parser.add_argument("--sampler", choices=["none", "smote", "borderline", "adasyn"], default="none")
    parser.add_argument("--trials", type=int, default=12)
    parser.add_argument("--folds", type=int, default=3)
    parser.add_argument("--trees", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if args.trials < 1 or args.folds < 2 or args.trees < 1:
        parser.error("Require trials >= 1, folds >= 2, trees >= 1")
    run(args)


if __name__ == "__main__":
    # Persist the importable class, never __main__.VariantEnsemble in joblib.
    from src.robust_variant_training import main as entrypoint
    entrypoint()
