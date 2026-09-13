# GeneVISTA

A research workspace for ClinVar classification history, ClinGen gene–disease evidence, exploratory disease-annotation ranking, and an interactive DNA → gene → variant → protein lesson.

The disease model ranks future record-level MONDO annotation events when a starting VUS becomes pathogenic/likely pathogenic. It does **not** estimate patient disease risk, diagnose a disease, or provide calibrated clinical-class probabilities. Current classifications and curated gene-level evidence remain separate from model outputs.

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

See [deployment/README.md](deployment/README.md). The Docker build serves the frontend and API together. The Blueprint configures authentication, health checks, and a persistent database disk. Full data and caches stay out of Git. The small trained disease model and prepared ClinGen evidence are versioned.

## Verification

```powershell
& D:\AncondaAPP\python.exe -m unittest discover -s tests -v
```

Build the frontend with `pnpm build`. Runtime checks use the cached SQLite database. Notebook 06 contains executed outputs and a saved-model reload check.
