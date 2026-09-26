# Protein language-model and CatBoost experiment

This experiment uses cached model scores, not a chat model generating diagnoses or synthetic labels. ESM1b is a protein language model; MutFormer applies BERT-style modelling to protein context. popEVE adds evolutionary/population modelling. No paid API or new transformer download is required. These scores do not provide patient disease probabilities.

Source: [dbNSFP column definitions](https://dist.genos.us/release/dbNSFP5.4a_variant.columns.txt). Only an explicit feature allowlist is imported. Clinical classification, review status, variant IDs, and disease-label columns are excluded.

The compressed download is rechecked against its original MD5. An exact GRCh38 chromosome, position, reference, and alternate allele join imports features into a separate database. The active annotation database is unchanged. The extraction reads training and validation cohort keys only, never test labels. Individual missing values stay missing. Duplicate source rows are pooled explicitly.

ESM1b, popEVE, and MutFormer summaries preserve minimum, median, and maximum scores; the lower direction of ESM1b/popEVE is not accidentally reversed. These are allele-level pooled summaries, not a claim of transcript-specific matching. MPC/VEST4 use MANE Select entries only when score and transcript vectors align exactly. This does not establish that the chosen transcript is appropriate for every clinical condition.

CatBoost compares existing annotation/sequence inputs against inputs augmented with these protein features, using the same locus-sampled training and validation records. Class weighting uses training labels only. Candidate selection maximizes the smaller of accuracy and balanced accuracy; macro F1 breaks ties. Checkpoints and candidate metrics are recorded. CatBoost categorical processing only learns from its training pool. Neither the protein-feature importer nor training loader can request test labels.

The already-inspected internal test is not reused as an untouched benchmark. Development validation has been used repeatedly and is optimistic for model selection. Upstream predictor training overlap remains unresolved. A genuine generalization claim requires new, independently curated evaluation records. The current deployed model is retained unless a separately evaluated replacement is deliberately installed.

## Run on D:

Install `research/requirements-protein-lm.txt` in the existing project environment. Set TEMP and TMP to D:/GeneVISTA/.runtime; use four OpenMP/BLAS threads.

```powershell
.venv/Scripts/python.exe -m src.protein_lm_features
.venv/Scripts/python.exe -m src.train_protein_catboost --scope baseline --output reports/catboost_NEW_BASELINE
.venv/Scripts/python.exe -m src.train_protein_catboost --scope protein_lm --output reports/catboost_NEW_PROTEIN
```

The importer refuses to overwrite an existing feature database. Each training run requires a new output directory. Reports contain same-cohort comparisons; model checkpoints use CatBoost's native format. The optional training libraries are not added to web-server deployment requirements.

References: [CatBoost fit API](https://catboost.ai/docs/en/concepts/python-reference_catboostclassifier_fit), [ESM1b paper](https://doi.org/10.1038/s41588-023-01465-0), [MutFormer paper](https://doi.org/10.1016/j.xinn.2023.100487).

The UI exposes experimental validation metrics under Latest training experiments. They refresh independently of the current model test metrics.

## Pilot comparison

All rows below use the same 163,620 training and 74,330 validation records.

| Features | Training class-weight exponent | Accuracy | Balanced accuracy |
|---|---:|---:|---:|
| baseline | 0.5 | 89.23% | 67.01% |
| baseline | 1.0 | 75.71% | 71.47% |
| protein_lm | 0.5 | 89.36% | 67.36% |
| protein_lm | 1.0 | 76.14% | 71.76% |

Feature extraction retained 2,003,946 exact development-cohort alleles: ESM1b scores for 1,949,576, popEVE for 1,554,734, and MutFormer for 1,856,836. Missing scores remain missing.

## Completed full-data candidate

The protein-language-model feature candidate trained on 1,633,384 records and was evaluated on the same 74,330 development-validation records as the incumbent. At 500 trees, accuracy is **77.4707%**, balanced accuracy **74.1590%**, and macro F1 **0.54518**. On this identical validation cohort the incumbent scored 76.5909% accuracy and 69.9006% balanced accuracy: improvements of 0.8799 and 4.2585 percentage points respectively. Neither metric reached 90%.

These are reused development-validation results, not independent test or clinical-validation results. No calibration or test labels were read in this experiment; the production model remains unchanged. The saved native model passed probability round-trip verification. Offline inference was checked on three cached validation variants and an absent variant; absent or ambiguous alleles do not receive inferred predictions.

Artifacts: `reports/catboost_protein_v6_full/`. The UI's Latest training experiments panel reads the completed metrics. To inspect a cached research prediction:

```powershell
.venv/Scripts/python.exe -m src.protein_candidate --variant 17
```

Class estimates are uncalibrated weighted model outputs, not disease occurrence probabilities. The inference helper deliberately requires exact alleles and matching source provenance.
