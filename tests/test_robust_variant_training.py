"""Regression tests for split, history, serialization and target-leakage guards."""
import tempfile
import unittest
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from src.robust_variant_training import (
    VariantEnsemble, add_history, checked_split, folds, partitions,
    records, tune_offsets, validate_schema,
)


class RobustTrainingTests(unittest.TestCase):
    def setUp(self):
        self.schema = {"target": "label", "group": "locus", "numeric": ["score"],
                       "categorical": ["gene"], "tokens": ["hpo"]}
        self.frame = pd.DataFrame({"label": ["a", "b", "c"] * 80,
                                   "locus": [str(i // 2) for i in range(240)],
                                   "score": np.arange(240) % 7,
                                   "gene": ["G1", "G2"] * 120,
                                   "hpo": ["HP:1|HP:2"] * 240})

    def test_group_disjoint_partitions_and_folds(self):
        ids = partitions(self.frame, self.schema, 42)
        for i in range(3):
            for j in range(i):
                self.assertFalse(set(self.frame.iloc[ids[i]].locus) & set(self.frame.iloc[ids[j]].locus))
        for a, b in folds(self.frame, self.schema, 3, 42):
            self.assertFalse(set(self.frame.iloc[a].locus) & set(self.frame.iloc[b].locus))

    def test_target_and_identifier_never_features(self):
        first = records(self.frame, self.schema)
        changed = self.frame.assign(label="secret", locus="secret", unused_leak="secret")
        self.assertEqual(first, records(changed, self.schema))
        bad = dict(self.schema, numeric=["label"])
        with self.assertRaises(ValueError):
            validate_schema(self.frame, bad)

    def test_unseen_categories_and_serialized_decisions(self):
        model = VariantEnsemble(self.schema, {"n_estimators": 4, "max_depth": 3}).fit(self.frame, self.frame.label)
        tune_offsets(model, self.frame, self.frame.label.to_numpy())
        unseen = self.frame.iloc[:3].assign(gene="unseen", hpo="HP:unknown")
        pred = model.predict(unseen)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "model.joblib"
            joblib.dump(model, path)
            np.testing.assert_array_equal(pred, joblib.load(path).predict(unseen))
        self.assertTrue(np.allclose(model.predict_proba(unseen).sum(axis=1), 1))

    def test_temporal_boundaries_and_delayed_outcomes(self):
        frame = self.frame.copy()
        frame["time"] = pd.to_datetime([f"2020-01-{i // 12 + 1:02d}" for i in range(240)], utc=True)
        frame["available"] = frame.time
        schema = dict(self.schema, time="time", label_time="available")
        validate_schema(frame, schema)
        train, tune, test = partitions(frame, schema, 42)
        self.assertLess(frame.iloc[train].available.max(), frame.iloc[tune].time.min())
        self.assertLess(frame.iloc[tune].available.max(), frame.iloc[test].time.min())
        for a, b in folds(frame, schema, 3, 42):
            self.assertLess(frame.iloc[a].available.max(), frame.iloc[b].time.min())
        frame.loc[0, "available"] = pd.Timestamp("2030-01-01", tz="UTC")
        a, _ = checked_split(frame, np.arange(120), np.arange(120, 240), schema)
        self.assertNotIn(0, a)

    def test_strictly_past_history_and_future_invariance(self):
        schema = dict(self.schema, time="time")
        query = pd.DataFrame({"locus": ["1"] * 3, "time": pd.to_datetime(["2020-01-01", "2020-01-02", "2020-01-03"], utc=True)})
        history = pd.DataFrame({"locus": ["1"] * 3, "available_at": ["2020-01-01", "2020-01-02", "2030-01-01"],
                                "observed_label": ["a", "b", "c"], "evidence_count": [1, 3, 100]})
        result = add_history(query, history, schema)
        self.assertTrue(pd.isna(result.loc[0, "history_previous_label"]))
        self.assertEqual(result.loc[1, "history_previous_label"], "a")
        self.assertEqual(result.loc[2, "history_evidence_delta"], 2)
        history.loc[2, "observed_label"] = "future_modified"
        pd.testing.assert_frame_equal(result, add_history(query, history, schema))

    def test_smote_rejects_categorical_schema(self):
        with self.assertRaisesRegex(ValueError, "numeric schemas"):
            VariantEnsemble(self.schema, sampler="smote").fit(self.frame, self.frame.label)

    def test_missing_class_fails(self):
        with self.assertRaisesRegex(ValueError, "lacks classes"):
            checked_split(self.frame, np.arange(60), np.array([60]), self.schema)


if __name__ == "__main__":
    unittest.main()
