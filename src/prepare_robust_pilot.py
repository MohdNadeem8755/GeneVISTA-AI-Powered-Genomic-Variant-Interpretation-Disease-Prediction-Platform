"""Export a reproducible exploratory pilot from the legacy TRAIN partition only.

Never treat the resulting holdout as a new independent ClinVar benchmark.
The cache is a trusted local joblib artifact created by this repository.
"""
import argparse
from pathlib import Path

import joblib
import pandas as pd

from src.robust_variant_training import save_json
from src.variant_features import variant_features


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", default="data/processed/variant_prediction_v1/cohort.joblib")
    parser.add_argument("--output", required=True)
    parser.add_argument("--rows", type=int, default=6000)
    args = parser.parse_args()
    if args.rows < 100:
        parser.error("Use at least 100 rows")
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    frame, _, _ = joblib.load(args.cache)
    pool = frame[frame.split.eq("train")]
    selected = pool.sample(n=min(args.rows, len(pool)), random_state=42).reset_index(drop=True)
    # Reuse only sequence-derived numeric features; gene/type stay categorical.
    numeric_rows = [{key: value for key, value in variant_features(row).items() if "=" not in key}
                    for row in selected.to_dict("records")]
    numeric = pd.DataFrame(numeric_rows)
    exported = pd.concat([selected[["label", "site", "gene_symbol", "variant_type"]], numeric], axis=1)
    exported.to_csv(out / "cohort.csv.gz", index=False)
    save_json(out / "schema.json", {"target": "label", "group": "site", "numeric": list(numeric.columns),
                                    "categorical": ["gene_symbol", "variant_type"], "tokens": []})
    save_json(out / "provenance.json", {"source": str(Path(args.cache).resolve()), "source_partition": "train",
                                        "rows": len(exported), "status": "exploratory_reused_source", "seed": 42})


if __name__ == "__main__":
    main()
