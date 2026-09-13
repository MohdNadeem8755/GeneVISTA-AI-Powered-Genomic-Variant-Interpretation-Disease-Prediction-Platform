# Render deployment

1. Push this repository, then create a Render Blueprint from `render.yaml`. It uses a paid Starter service and a 5 GB persistent disk; review the current plan cost before creating it.
2. Provision the full `data/app/genevista.sqlite` database separately. GitHub cannot store this 1.4 GB file in normal Git. Upload it to storage you control and provide a direct HTTPS download URL as `GENEVISTA_DB_URL`, or copy it to `/var/data/genevista.sqlite` on the Render disk before starting. A temporary signed URL is acceptable for first provisioning; keep a reproducible copy available for disk replacement.
3. Calculate the local file checksum with `Get-FileHash D:\GeneVISTA\data\app\genevista.sqlite -Algorithm SHA256` and set `GENEVISTA_DB_SHA256`. Downloads are checked before activation. Existing data is verified against the configured checksum at startup.
4. Run `D:\AncondaAPP\python.exe deployment/password_hash.py`. Set the resulting hash and a chosen username in `GENEVISTA_PASSWORD_HASH` and `GENEVISTA_USERNAME`. Keep authentication enabled and secure cookies enabled on Render.
5. Deploy. Startup validates the database and model before serving. Open the Render URL, sign in, search BRCA1, select a VUS, and verify its ranking and the evaluation charts. `/health` exposes readiness without credentials; the data APIs require a session.

The frontend and API share an origin. No production API URL or cross-site cookie setup is required. One server worker is intentional because sessions and rate limits are stored in memory. Restarts invalidate sessions. The persistent disk stores only the read-only scientific database. Increase the service memory if the measured workload needs it.

To change the database version, replace the disk file and update its checksum together. A changed URL alone does not overwrite an existing database.

Reference: https://render.com/docs/blueprint-spec

Current cached database SHA-256: `369094a1b332590821a98489eb2fe6158572c2eb0550183560b157945ebf082f`. Recompute this if the data changes.
