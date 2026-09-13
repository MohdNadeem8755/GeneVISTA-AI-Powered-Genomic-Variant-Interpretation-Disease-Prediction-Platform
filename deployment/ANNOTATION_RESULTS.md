# dbNSFP annotation experiment: September 2026

The local download is stored in `D:\GeneVISTA\data\incoming`: dbNSFP5.4a_grch38.gz, its .md5 and .tbi. The verified MD5 is d2ad53e1a974208fea6050117be8fee1. Import scanned 91,571,221 rows and matched 246,793 unique alleles in the sampled cohort. The raw download is not a deployment artifact.

## Reproducible training

Run Python modules from the repository root. Set TEMP/TMP to D:\GeneVISTA\.runtime and constrain OMP_NUM_THREADS to 2. `src.annotation_pipeline` prepares and imports the cohort. `python -m src.train_annotation_candidates` compares eight configurations on the full matched training/validation partitions. `python -m src.finalize_annotation_model` calibrates the selected frozen estimator and evaluates the locked test once; its exclusive evaluation marker prevents accidental repeated test selection.

The selected model uses 15 leaves, 250 boosting iterations, no class reweighting, annotation features and numeric sequence features. Training has 163,219 matched records, validation 36,922, calibration 24,544 and test 22,105. The comparison selected by validation macro F1, not by test performance.

## Measured result

| Metric | Locked matched test |
|---|---:|
| Accuracy | 90.69% |
| Balanced accuracy | 54.10% |
| Macro F1 | 0.610 |
| Log loss | 0.2973 |

99% accuracy was not achieved. High overall accuracy does not imply strong performance on every class. These are same-snapshot, position-separated internal results, not external clinical validation. Supervised source predictors may have used overlapping ClinVar training data. A separate population/conservation ablation excluded supervised pathogenicity scores; upstream overlap remains a limitation of the selected model. Labels and review metadata are not input features. Missing values remain missing.

## Local application

The additional annotation panel is served by `/api/variants/{id}/annotation-prediction`. A covered example in the demo catalog is variant 14094. The panel displays one class, five probabilities, cohort membership and actual test metrics. Existing sequence predictions and disease annotation rankings remain separate. Disease scores are not patient disease risks.

Local artifacts are `models/annotation_prediction_v2/predictor.joblib` and `model_card.json`. The optional read-only lookup defaults to `data/processed/annotation_v2/download_verified_cohort.sqlite`; override with GENEVISTA_ANNOTATION_DB. The request must match the cached variant name and gene, and lookup provenance must match model provenance. Missing or ambiguous matches return an explicit reason.

Raw data, lookup and experimental model artifacts are kept locally and excluded from Git. Thus the new annotation panel requires local artifacts; an existing Render deployment without them reports that the optional model is not configured. The existing sequence model remains usable there. Shipping annotation resources needs a deliberate deployment package appropriate for the source terms and host resources.
