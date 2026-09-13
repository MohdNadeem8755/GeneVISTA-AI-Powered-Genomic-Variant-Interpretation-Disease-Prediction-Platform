# GeneVISTA frontend

React, Vite, and CSS. The workspace includes search and classification timelines, ClinGen disease evidence, VUS disease-annotation rankings, temporal performance charts, an animated DNA logo and molecular lesson, and login/logout.

From `D:\GeneVISTA\frontend`:

```powershell
pnpm install --frozen-lockfile --store-dir D:/GeneVISTA/.runtime/pnpm-store
pnpm dev
```

Run the backend from `D:\GeneVISTA` with `D:\AncondaAPP\python.exe -m uvicorn backend.app.main:app --host 127.0.0.1 --port 8000`. The development frontend is at `http://127.0.0.1:5173`.

`pnpm build` writes `frontend/dist`. The backend serves this production build and its API on the same origin. Render uses this arrangement; no production API URL is needed. `VITE_API_BASE_URL` supports an explicit API origin for other deployments, which also require corresponding CORS and cookie settings.

Local development defaults to unauthenticated research mode. See the root README for enabling the workspace account. Authentication, session expiry, and logout are enforced by the backend.

Ranking scores are uncalibrated and apply only to eligible VUS records. They are not disease-risk probabilities or clinical diagnoses. The charts display average precision, event prevalence, and positive support from saved temporal evaluation.

Local caches and draft directories are ignored by Git. Keep local package caches on D: using the store option above.
