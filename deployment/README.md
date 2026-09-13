# Free Render demo

The default Docker image includes a small read-only database extracted from the cached full dataset. No paid disk, external database upload, download URL, or checksum environment variable is required. The full local dataset remains at `D:\GeneVISTA\data\app\genevista.sqlite`.

## Render setup

Use the GeneVISTA GitHub repository, branch `main`, runtime **Docker**, empty Root Directory, Dockerfile `./Dockerfile`, and the **Free** instance. Set the health check path to `/health` under Advanced.

Set these environment variables:

| Name | Value |
| --- | --- |
| GENEVISTA_USERNAME | Your chosen workspace login name |
| GENEVISTA_PASSWORD_HASH | Output of the local password hash utility |
| GENEVISTA_AUTH_ENABLED | true |
| GENEVISTA_COOKIE_SECURE | true |
| GENEVISTA_DB | /app/data/demo/genevista.sqlite |

Remove any old `GENEVISTA_DB_URL` and `GENEVISTA_DB_SHA256` variables. Replace an old `/var/data/genevista.sqlite` path with the demo path above. Leave Docker Command empty. Do not add a persistent disk.

Run `D:\AncondaAPP\python.exe D:\GeneVISTA\deployment\password_hash.py` to generate the salted login hash. Paste the hash into Render, never the original password. Use the original password to log in to the deployed app. Credentials remain private; the repository does not include a default password.

## Demo scope

The builder selects up to 50 variants per gene and grouped classification for 20 named genes, plus the existing example IDs. Their stored variant details, gene links, and six-snapshot timelines are copied exactly. This is a selected demo, not a representative sample. The UI labels it and reports only demo record counts. The disease model and temporal evaluation remain the originally trained full-cohort results.

Rebuild locally with `D:\AncondaAPP\python.exe deployment/build_demo.py`. The source database is attached in read-only mode and never modified. Commit the resulting small `data/demo/genevista.sqlite` file to update the deployment sample.

## Free service behavior

Render can spin down idle free services. Opening the app after inactivity may require waiting for startup and using Retry. Sessions are kept in server memory, so a restart signs users out. The database is part of the image and is restored with each deployment. See https://render.com/docs/free.

A live Render deployment still requires connecting the repository and entering your private login variables. Local startup checks do not verify Render's Docker build or service availability.
