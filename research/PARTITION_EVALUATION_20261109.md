# Frozen ensemble partition evaluation — November 9, 2026

User requested training, validation and testing scores on the dashboard.
Scored the existing `biological_retuning_fresh_20261031/selected.joblib` without
fitting, changing weights, choosing thresholds, or deploying it for predictions.
Model SHA256: `4bcb109a50aa7ff2175d588c1353b878396cffb859cdb1631ddbd45a86300d19`.

| Partition | Records | Accuracy | Balanced accuracy | Macro F1 |
| --- | ---: | ---: | ---: | ---: |
| Fitting | 295,067 | 81.64% | 65.55% | 0.5781 |
| Early stopping | 32,595 | 81.33% | 64.43% | 0.5715 |
| Validation | 73,770 | 81.51% | 64.98% | 0.5732 |
| Internal test | 200,134 | 81.50% | 65.35% | 0.5755 |

Training and stopping rows reconstruct the original sampled train-locus split.
Validation reconstructs the original selection sample. The internal test uses
the existing annotation-matched test extension, excluding legacy v2 sampled
loci, as for the active v4 model. This historical test has been inspected during
earlier project experiments; it is not an untouched later/external evaluation.
No new clinical-performance claim is supported. Calibration rows were not used.

Features use the original frozen source files and October 29 feature cutoff.
Development features come from the original cache; test features are generated
using the same enrichment functions, column schema and missing-value handling.
The original development-only training loader still rejects test data.

Verification: source hashes and annotation provenance checked, zero cross-split
locus overlap, finite normalized probabilities, exact original validation metrics
reproduced, expected partition sizes checked, selected model hash unchanged.
Ten relevant tests passed, including endpoint rejection of mismatched model
reports. Production frontend build passed.

Runner: `src/evaluate_dashboard_partitions.py`. The runner refuses to overwrite
an existing completed evaluation. Results and protocol reside beside the model:
`reports/biological_retuning_fresh_20261031/partition_evaluation.json` and
`partition_evaluation_protocol.json`. The evaluation JSON includes all per-class
metrics and confusion matrices, and is published atomically only after checks.

Dashboard now displays partition metrics for this single saved ensemble and
explains cleaning, feature preparation, class weights, early stopping and model
selection. It explicitly distinguishes original cohort exclusions from the later
label audit, which flagged changed/conflicting labels without correcting them.
The separately accessible active-model report continues to describe v4; individual
variant inference still uses that deployed model.
