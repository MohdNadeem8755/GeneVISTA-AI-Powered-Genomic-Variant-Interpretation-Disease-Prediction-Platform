# Expanded cohort: five-class gene specialist comparison

The full cached GRCh38 SNV cohort was imported using verified dbNSFP5.4a annotations. Eight genes met the exploratory minimum of 30 training and 10 validation records per class: BRCA1, BRCA2, DNAH11, KMT2D, MECP2, MUTYH, SYNGAP1 and TSC2. Only BRCA1 also met the high-quality-label screen.

`python -m src.train_gene_specialists` trains two gene-aware pooled controls and two specialist configurations per gene. Pooled and specialist candidates use the same expanded training records and evaluate the same validation records. Pooled controls receive gene identity. Constant/all-missing columns are removed using training data only; remaining missing values stay missing. A regression test protects this behavior after a local histogram-binning failure on a constant column.

| Candidate | Validation accuracy | Balanced accuracy | Macro F1 |
|---|---:|---:|---:|
| Previous global annotation candidate | 77.35% | 43.60% | 0.456 |
| Expanded pooled model | 81.16% | 56.91% | 0.597 |
| Expanded pooled, partial class balancing | 79.45% | 64.95% | 0.636 |
| Gene specialists | 80.90% | 57.94% | 0.610 |
| Gene specialists, partial class balancing | 78.50% | 60.78% | 0.608 |

These results concern only the eight support-selected genes and are exploratory validation results. They are not overall catalog accuracy and cannot be directly compared with the existing 90.69% internal test result. Previous-global versus expanded-candidate improvement includes additional training data. Specialists do not outperform the strongest pooled control across the reported metrics. No 99% accuracy was achieved and no candidate was promoted.

No calibration or test rows were read. Eligibility used validation label counts, previously evaluated test loci remain used, and upstream source predictor overlap is unresolved. A fresh independent assessment remains necessary before a generalization claim. Candidate models are uncalibrated research artifacts.

Reports and artifacts remain local under reports/annotation_v3/specialists_variable_features; the expanded database remains under data/processed/annotation_v3/full_cohort.sqlite. Current deployed/local production predictions are unchanged. The next evidence-backed research direction is improved biological/functional evidence and a prospectively defined evaluation cohort, rather than repeatedly tuning to the old test.
