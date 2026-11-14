# Annotated ensemble retraining — 2026-10-29

Retrained two CatBoost classifiers (depths 6/8, class-weight powers 0.5/1)
and a weighted histogram gradient boosting classifier on 326,951 annotated
variants. Used 300 boosting iterations per learner. Features include population
frequency, conservation, sequence changes and existing pathogenicity annotations.
Audited zero locus overlap across the existing partitions.

Calibrated each learner on a separate sampled calibration partition. Selected
four soft-voting mixtures and five class-prior decision corrections on 74,330
validation rows by balanced accuracy, then macro-F1. The winning mixture weights
were 2:1:1 and prior exponent 1. The same validation partition had previously
been used for research, so selection estimates are exploratory.

After freezing selection, evaluated the previously inspected internal test:

| Model | Accuracy | Balanced accuracy | Macro F1 |
| --- | ---: | ---: | ---: |
| New ensemble | 71.71% | 70.40% | 0.4902 |
| Active model (UI, rounded) | 76.60% | 70.9% | 0.527 |

The new ensemble did not improve on the active model on validation or test.
It was **not installed**. Neither metric exceeded 95%. More learners alone did
not resolve the limited separation between the five recorded annotation classes.
Independent functional features and carefully audited minority labels remain
necessary directions; do not reuse this test to claim independent improvements.

Full exact comparison and candidate classification report/confusion matrix:
`reports/annotation_ensemble_20261029/results.json`. Candidate artifact:
`reports/annotation_ensemble_20261029/candidate.joblib`. These research artifacts
are not production bundles. Script: `src/retrain_annotation_ensemble.py`.

Reproduce with a new output directory:

```powershell
$env:OMP_NUM_THREADS='4'
$env:OPENBLAS_NUM_THREADS='4'
.venv\Scripts\python.exe -m src.retrain_annotation_ensemble --output reports/annotation_ensemble_NEW --modulus 5 --iterations 300
```

UI now displays only the active annotation accuracy section, removes old-model
and training-comparison accuracy panels, and directs navigation to that section.
It explains that negative BRCA1 functional scores are valid. The cached c.190T>G
functional and RNA scores match the displayed -1.951 and -5.421. The assay's
non-functional cutoff is -1.328, consistent with the
[original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6181777/).

Validation: frontend production build passed; nine assay, active-model,
decision-policy and ensemble class-order tests passed. Local HTTP response was 200.
