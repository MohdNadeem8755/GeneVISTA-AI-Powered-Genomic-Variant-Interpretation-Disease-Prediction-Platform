# Minority-class error investigation — 2026-10-30

Analyzed the unchanged active `annotation_prediction_v4` on **370,562** matched
validation records. This validation set was previously used for model selection;
results are descriptive development findings, not a fresh independent benchmark.
No test labels were accessed, and no labels, thresholds or models were changed.

## Where the errors occur

Counts below use the original cohort label as reference, not clinical ground truth.
“Opposite group” means benign/likely benign versus pathogenic/likely pathogenic.

| Reference class | Support | Recall | Predicted paired class | Predicted VUS | Opposite group |
| --- | ---: | ---: | ---: | ---: | ---: |
| Benign | 5,094 | 77.40% | 405 | 714 | 32 |
| Likely benign | 19,892 | 67.87% | 1,632 | 4,510 | 250 |
| Likely pathogenic | 11,136 | 64.57% | 2,858 | 1,052 | 36 |
| Pathogenic | 12,263 | 64.41% | 3,704 | 624 | 37 |

Of **15,854 errors** in these four classes, 8,599 are within-pair swaps, 6,900
go to VUS, and 355 go to the opposite group. The two problems need different
responses: benign-side errors mostly become VUS; pathogenic-side errors mostly
confuse the strength of the pathogenic classification.

## Main findings

**Population-frequency missingness is a priority.** Benign recall is 8.70%
(28/322) when joint gnomAD AF is absent, versus 82.04% (3,915/4,772) when it is
present. Missing AF is not evidence of rarity or zero frequency. This association
may also reflect different genes, variants, sequencing coverage and ascertainment;
the analysis does not establish a causal effect of missingness. Inspect reference
matching and population coverage before treating this as a modeling issue.

**Review tiers behave differently by class.** Benign recall is 95.99% for
multiple-submitters/no-conflicts records (2,441 cases), versus 53.80% for
single-submitter records (2,108 cases). For Likely benign the pattern reverses:
58.14% versus 69.83%. Expert subsets are small and were weaker in Step 1.
Consequently, simply filtering to one review tier does not resolve all boundaries.
Review status remains a diagnostic variable, not an added target-proxy feature.

**Label drift explains only a small identified subset.** Of the 15,854 errors,
605 have changed labels in the October 5 snapshot and 328 are newly conflicting.
These counts overlap. Reconcile them first, but doing so cannot explain away all
errors. “Unchanged” means no change flagged by the Step 1 audit, not a verified
correct label. Combined labels must not be forced into a particular exact class.

**The balancing policy trades overall accuracy for minority recall.** On these
same validation records:

| Decision policy | Accuracy | Balanced accuracy | Macro F1 |
| --- | ---: | ---: | ---: |
| Active class-balanced rule | 76.60% | 70.45% | 0.5272 |
| Highest raw probability, reference only | 91.29% | 57.16% | 0.6274 |

Relative to highest probability, the active rule corrects 496 Benign, 8,757 Likely
benign, 3,302 Likely pathogenic and 593 Pathogenic predictions, while harming
0, 137, 0 and 19 respectively. It also turns 67,453 previously correct VUS
predictions into errors. Neither policy reaches the requested targets. These
results do not authorize switching policies based on overall accuracy alone.

**Gene-level failure differs.** Among gene/class cells with at least 100 validation
records, FBN1 Likely pathogenic is a useful review target: 45/105 correct, with
55 predicted Pathogenic. In contrast, TTN Likely pathogenic is 299/307 correct.
These are exploratory, selected examples, not evidence that a gene-specific model
will necessarily improve performance. The full table includes all genes and
support counts to avoid hiding small-sample failures.

**New feature coverage is uneven.** The Step 2 cache has matched functional
assays for only 52 erroneous records in the four focus classes, but ClinGen
context for 11,154. These features were not used by the current model. Their
presence does not establish predictive benefit; gene context can also encourage
gene memorization and needs gene-held-out evaluation.

## Recommended Step 4 experiments

1. Audit allele matching and coverage for the 322 Benign records missing gnomAD
   AF; preserve missingness instead of imputing zero frequency.
2. Resolve the flagged drift/conflict subset through evidence review without
   replacing labels based on model predictions.
3. On fixed development folds, compare the original features with Step 2
   enrichment, separately ablating population, ClinGen, splice and assay features.
4. Compare a five-class model against a hierarchical model that first separates
   benign-like/pathogenic-like/VUS, then predicts the exact class within each
   pair. Report final five-class metrics so the task remains unchanged.
5. Tune class costs and decisions on development folds, reporting the VUS recall
   cost explicitly. Do not inflate overall accuracy by sacrificing minority classes.
6. Evaluate a frozen choice on an untouched future dataset; the existing validation
   and historical test have already been inspected.

These are proposed experiments, not completed training or promises of 95%.

## Artifacts and reproduction

`reports/minority_errors_20261030/report.json` contains full classification
reports, confusion matrices, policy effects and model SHA-256.

`priority_examples.csv` is a bounded review queue, stratified by reference class
and error type and ordered by conflict/drift flags then selected-class probability.
It is not a random sample or a replacement training set.

Additional outputs: `class_errors.csv`, `review_tier_errors.csv`,
`label_drift_errors.csv`, `gene_errors.csv`, per-feature coverage/error tables,
`feature_error_profiles.csv`, and `validation_predictions.csv.gz`.

```powershell
$env:OMP_NUM_THREADS='4'
$env:OPENBLAS_NUM_THREADS='4'
.venv\Scripts\python.exe -m src.investigate_class_errors --output reports/minority_errors_NEW
.venv\Scripts\python.exe -m unittest tests.test_class_errors -v
```

Inputs are opened read-only and existing output directories are rejected.
Error categories have a regression test. Output totals and confusion matrices
were checked against the complete prediction export.
