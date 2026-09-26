# ANN experiment: 23 September 2026

Three class-weighted MLP candidates were compared on 163,620 training and 74,330 validation records. The selected 128/64-unit ANN was then trained on all 1,633,384 matched training records. The selected full-data checkpoint was epoch 11; validation patience stopped training. No test or calibration records were read. Test accuracy is therefore not established for this ANN.

| Model | Same validation accuracy | Same validation balanced accuracy |
|---|---:|---:|
| Active v4 | 76.59% | 69.90% |
| Full-data ANN | 75.62% | 70.30% |

The ANN improved balanced accuracy by 0.40 percentage points but reduced overall accuracy by 0.97 points. Neither reached both metrics above 90%. It remains a research candidate; the current model is unchanged. Validation was reused for early stopping and selection, so it is not an independent performance estimate. Weighted ANN outputs are not calibrated patient-risk probabilities.

Numeric imputation, scaling, missing indicators, and gene one-hot encoding are fit on training rows only. Unknown genes are handled without fitting on validation. Chromosome/position overlap is checked. Labels, review status, disease annotations, and IDs are excluded from inputs. Cached AlphaMissense scores already provide a pretrained deep-learning feature, with unresolved upstream training overlap.

Reproduce from D:/GeneVISTA:

```powershell
.venv/Scripts/python.exe -m src.ann_annotation --output reports/ann_NEW_PILOT
.venv/Scripts/python.exe -m src.ann_annotation --output reports/ann_NEW_FULL --candidate ann128_balanced --train-modulus 1
```

Use a new output directory each time. Configure TEMP/TMP to D:/GeneVISTA/.runtime and BLAS/OpenMP threads to four. Models and detailed reports are saved in the run folders. A larger ANN alone did not resolve the class-confusion problem; further work needs additional independently sourced, transcript-aligned evidence and a fresh external evaluation.

## UI

The main result and quality panel identify the active annotation model. Original sequence predictions and evaluation are collapsed under comparison sections. Per-variant class estimates are explicitly distinguished from model accuracy. The API returns evaluation metadata from the same loaded artifact used for predictions. The interface refreshes saved results every minute, on focus, or manually; it does not manufacture real-time accuracy. Browser automation was unavailable due to the local sandbox helper, so verification used the production build, API checks, and regression tests.
