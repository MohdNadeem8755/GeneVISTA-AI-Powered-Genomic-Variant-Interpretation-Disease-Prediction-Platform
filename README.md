# GeneVISTA

A research workspace for ClinVar classification history, ClinGen gene–disease evidence, exploratory disease-annotation ranking, and an interactive DNA → gene → variant → protein lesson.

The current model produces sigmoid-calibrated estimates for five exact ClinVar annotation classes and independent MONDO disease-annotation matches from gene and parsed sequence/protein-change features. It does **not** estimate patient disease risk or diagnose disease. Recorded classifications, model predictions, and gene-level curation evidence are displayed separately. The earlier VUS reclassification experiment remains available under a collapsed historical section.

## Run locally on D:

Use the existing Python installation at `D:\AncondaAPP\python.exe`. Install `backend/requirements.txt` into a Python environment on D: if needed. Run from `D:\GeneVISTA`:

```powershell
$env:TEMP='D:\GeneVISTA\.runtime'
$env:TMP=$env:TEMP
$env:GENEVISTA_AUTH_ENABLED='false'
& D:\AncondaAPP\python.exe -m uvicorn backend.app.main:app --host 127.0.0.1 --port 8000
```

In `frontend`, run `pnpm install --frozen-lockfile --store-dir D:/GeneVISTA/.runtime/pnpm-store` and `pnpm dev`. After `pnpm build`, the backend serves the complete frontend at port 8000 with same-origin API requests. Local development can also use `VITE_API_BASE_URL`.

## Notebook 06 and model

Open `notebooks/06_disease_ranking.ipynb` with the D: Python kernel. It reuses `data/processed/disease_ranking_v1/prepared.joblib` when present. The SGD input is converted to CSR with checked 32-bit indices. Training selects disease categories using training support only and saves the model, model card, and temporal evaluation under `models/disease_ranking_v1`.

The UI ranks eligible VUS records with genes observed during training. Non-VUS records and unseen genes receive explicit eligibility messages. Performance charts compare average precision with prevalence and show event support. The 2025–2026 evaluation is exploratory because that period was previously examined.

## Login

Local mode can run without authentication. Production enables a single administrator-provisioned research account. Generate a salted password hash with `python deployment/password_hash.py`; set `GENEVISTA_USERNAME` and `GENEVISTA_PASSWORD_HASH` privately. Never commit passwords or environment secrets. Sessions use random, HTTP-only cookies, expire after eight hours, and are revoked on logout. Restarting the single server worker signs users out. Login attempts are limited. This is a shared workspace account, not public registration or patient-record storage.

For local login testing set `GENEVISTA_AUTH_ENABLED=true` and `GENEVISTA_COOKIE_SECURE=false`; production requires HTTPS and secure cookies.

## Deploy to Render

See [deployment/README.md](deployment/README.md). The Docker build serves the frontend and API together on a Free service. A bundled, clearly labeled demo database needs no external download or paid disk. The Blueprint configures authentication and health checks. Full local data and caches stay out of Git; the small demo database, trained model, and prepared ClinGen evidence are versioned.

## Verification

```powershell
& D:\AncondaAPP\python.exe -m unittest discover -s tests -v
```

Build the frontend with `pnpm build`. Runtime checks use the cached SQLite database. Notebook 06 contains executed outputs and a saved-model reload check.

## Five-class research predictor

`src/train_variant_prediction.py` trains from the cached September 2026 clean ClinVar snapshot. It uses approximately 10% of eligible genomic loci with exact five-class labels. Training, sigmoid calibration, and internal testing use disjoint chromosome/start/stop loci (70/15/15); example loci 17660 and 37565 are reserved for testing. This is an internal snapshot experiment, not prospective or external validation.

Inputs are only gene symbol, variant type, and parsed nucleotide/protein changes from the variant name. Classification, review status, submitter counts, disease names, and IDs are not model features. Training annotations still reflect ClinVar reporting and ascertainment biases.

The reused internal test contains 59,112 variants: accuracy 83.50%, balanced accuracy 53.83%, macro F1 0.531. Recall for Likely pathogenic is 6.71% and for Benign is 6.65%. The requested 99% accuracy has not been established. The model is not clinically reliable for distinguishing all five classes. The UI preserves disagreements between recorded evidence and predictions.

Disease percentages estimate annotation matches, conditional on a record having a MONDO annotation. Up to 30 training-supported disease categories are modeled; inference candidates require observed training support for an input gene. The percentages are independent and do not sum to 100%. They can reflect conditions being investigated, including for benign variants, and are not penetrance or patient disease risk.

Artifacts and detailed metrics are in `models/variant_prediction_v1`. Canonical disease names are resolved from the EMBL-EBI MONDO ontology service and saved locally. Training can reuse `data/processed/variant_prediction_v1/cohort.joblib`. The full local database and training cache remain on D:; Docker includes the small saved predictor and the labeled demo database.

To reproduce the nonlinear comparison, first run `python src/train_variant_prediction.py`, then `python src/improve_variant_prediction.py`. The second script compares a linear baseline and two histogram gradient-boosting candidates on a locus-disjoint inner validation split, selects by macro F1, refits on training data, and calibrates separately. It uses 256 training-selected features to bound memory. The original test was previously inspected; reported final scores are reused internal estimates, not external validation. Experiment measurements are saved in `models/variant_prediction_v1/training_experiments.json`.
