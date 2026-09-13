FROM node:22-bookworm-slim AS web
WORKDIR /build
COPY frontend/package.json frontend/pnpm-lock.yaml ./
RUN corepack enable && corepack prepare pnpm@10.11.0 --activate && pnpm install --frozen-lockfile
COPY frontend/ ./
RUN rm -f .npmrc
RUN pnpm build

FROM python:3.14-slim
WORKDIR /app
ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 GENEVISTA_DB=/app/data/demo/genevista.sqlite GENEVISTA_AUTH_ENABLED=true GENEVISTA_COOKIE_SECURE=true OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
COPY backend/requirements.txt backend/requirements.txt
RUN pip install --no-cache-dir -r backend/requirements.txt
COPY backend/ backend/
COPY src/ src/
COPY models/disease_ranking_v1/ models/disease_ranking_v1/
COPY data/processed/clingen_gene_evidence_mapped.csv.gz data/processed/clingen_gene_evidence_mapped.csv.gz
COPY data/demo/genevista.sqlite data/demo/genevista.sqlite
COPY deployment/ deployment/
COPY --from=web /build/dist frontend/dist
CMD ["python", "deployment/start.py"]
