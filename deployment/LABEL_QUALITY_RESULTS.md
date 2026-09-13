# Label-quality and gene-specialist feasibility audit

The audit reads training and validation only from the existing annotated 10% position sample. Expert-panel/practice-guideline labels and multiple-submitter/no-conflict labels are identified by exact review-status strings. Combined labels are not converted into a single class.

No gene met the exploratory specialist feasibility screen of at least 30 training and 10 validation examples in every one of the five classes. Even BRCA2 had only 13 training Likely pathogenic and 13 training Benign records. This is a sample-size limitation, not evidence that gene-specialist models cannot work.

## Controlled validation experiment

All candidates were compared on identical existing validation rows. These are exploratory validation results, not a new independent test.

| Candidate | All matched validation accuracy | All matched macro F1 | High-quality subset macro F1 |
|---|---:|---:|---:|
| Existing annotation candidate | 91.05% | 0.638 | 0.717 |
| High-quality training only | 90.39% | 0.599 | 0.724 |
| High-quality training, partial class balancing | 88.52% | 0.595 | 0.731 |

The new candidates were not promoted. Training only on cleaner labels slightly improves macro F1 on a similarly curated subset but worsens performance over the full matched validation set. No 99% result was achieved. Calibration and previously evaluated test partitions were not read.

Run `python -m src.audit_annotation_labels`, then `python -m src.compare_annotation_quality` for the local experiment. Existing comparison output is protected against accidental overwrite. Detailed reports and candidate models stay under reports/annotation_v3 on D: and are ignored by Git.

Next: expand the cached ClinVar cohort for candidate genes beyond the 10% sample, preserving position groups and excluding previously used evaluation loci from any claimed fresh test. Check per-class support again before specialist training. No new external download is required for that step. Source-predictor training overlap remains a separate concern for independent evaluation.
