# Balanced five-class annotation experiment

This is separate from notebooks 04–06, which model future VUS reclassification and disease-annotation rankings.

## Reproduce

Run from D:\GeneVISTA using the project environment:

```powershell
$env:TEMP='D:\GeneVISTA\.runtime'
$env:TMP=$env:TEMP
$env:OMP_NUM_THREADS='4'
$env:OPENBLAS_NUM_THREADS='4'
.venv\Scripts\python.exe -m src.balanced_annotation --output reports/balanced_annotation_NEW_RUN
```

An existing run directory cannot be silently reused. Repeated experiments after viewing test results must not describe the same test extension as untouched.

## Protocol

- Read verified dbNSFP 5.4a GRCh38 annotations already cached in annotation_v3/full_cohort.sqlite.
- Audit chromosome/position separation. Exclude review status, clinical labels, disease annotations and variant IDs from features.
- Use population frequency, conservation, sequence-change features and existing pretrained AlphaMissense scores. Gene vocabulary contains the 128 most frequent training genes; unseen genes remain missing.
- Compare five prespecified configurations on deterministic locus-sampled training/validation data. Include a population/conservation ablation excluding pathogenicity predictors.
- Select by validation balanced accuracy, with macro F1 and accuracy as tie-breakers. Refit the winning configuration on all matched training records.
- Calibrate on the separate calibration split. Select a decision-policy exponent from a fixed validation grid. Freeze the artifact before reading test labels.
- Evaluate the new test extension, excluding all loci in the legacy v2 sampled cohort. Compare the incumbent and candidate on exactly the same records.
- Upstream score training overlap is unresolved. This remains internal annotation prediction, not independent clinical or unseen-gene validation.

## Decision policy and UI

The predicted class maximizes `calibrated_probability / training_class_prior ** exponent`. This targets more even class recall and can lower overall accuracy. Probabilities displayed in the UI remain the original calibrated class estimates. They are not renamed or rescaled decision scores, and are not patient disease probabilities. The UI explains when the selected class differs from the largest probability.

The optional v4 artifact lives under models/annotation_prediction_v4. If absent, the API retains v2. Its default lookup is the full local v3 annotation database. GENEVISTA_ANNOTATION_DB can override the lookup, but provenance and exact allele/details matching remain required. The existing sequence and disease-ranking models remain separate.

## Evidence files

The experiment folder contains audit.json, development.json, selection_frozen.json, test_started.json, model_card.json and the saved predictor. These record coverage, class-wise results, calibration, provenance and the accuracy/balanced-accuracy trade-off. Do not compare accuracy across different cohorts as if they were the same benchmark.
