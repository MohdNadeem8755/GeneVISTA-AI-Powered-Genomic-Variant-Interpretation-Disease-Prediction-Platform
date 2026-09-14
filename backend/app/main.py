"""Run from project root: python -m uvicorn backend.app.main:app --host 127.0.0.1 --port 8000"""

import logging
import os
import sqlite3
from pathlib import Path

from fastapi import FastAPI, HTTPException, Path as PathParameter, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from .database import Database
from . import auth
from .ranking import rank, evaluation
from .prediction import predict, bundle as prediction_bundle
from .annotation_prediction import predict as annotation_predict
from .tp53_evidence import evidence as tp53_evidence
from fastapi.staticfiles import StaticFiles

project_root = Path(__file__).resolve().parents[2]
database = Database(os.environ.get("GENEVISTA_DB", project_root / "data" / "app" / "genevista.sqlite"))
app = FastAPI(
    title="GeneVISTA API",
    version="0.1.0",
    description="Research five-class predictions, ranked disease annotation matches, and observed variant evidence.",
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=True,
    allow_methods=["GET", "POST"],
    allow_headers=["Accept", "Content-Type", "X-GeneVISTA-Request"],
)


app.include_router(auth.router)

def prediction_support(result):
    return {'clinical_class_probability':result['status']=='available',
            'disease_annotation_match_probability':result.get('diseases',{}).get('status')=='available',
            'disease_probability':False,
            'reason':'Model class and annotation-match estimates are research outputs; patient disease risk is not estimated.'}


@app.middleware("http")
async def workspace_access(request: Request, call_next):
    if request.method == 'POST' and request.url.path.startswith('/api/auth/'):
        if request.headers.get('x-genevista-request') != '1':
            return JSONResponse(status_code=403,content={'detail':'Missing request header.'})
    if auth.enabled() and request.url.path.startswith('/api/') and not request.url.path.startswith('/api/auth/'):
        if not auth.session(request):
            return JSONResponse(status_code=401,content={'detail':'Please sign in.'})
    response = await call_next(request)
    response.headers['X-Content-Type-Options'] = 'nosniff'
    response.headers['Referrer-Policy'] = 'strict-origin-when-cross-origin'
    if request.url.path.startswith('/api/'):
        response.headers['Cache-Control'] = 'no-store'
    return response

@app.exception_handler(sqlite3.Error)
async def database_error(request: Request, error: sqlite3.Error):
    logging.exception("Database access failed", exc_info=error)
    return JSONResponse(status_code=503, content={"detail": "The application database is unavailable. Check the server log and database path."})


@app.get("/api")
def root():
    return {"name": "GeneVISTA API", "docs": "/docs", "health": "/health"}


@app.get("/health")
def health():
    try:
        stats = database.stats()
    except (ValueError, KeyError) as error:
        raise HTTPException(status_code=503, detail="Database metadata is incompatible.") from error
    return {"status": "ok", "database": "ready", "latest_snapshot": stats["latest_snapshot"]}


@app.get("/api/stats")
def stats():
    return database.stats()


@app.get("/api/variants/search")
def search(
    q: str = Query(..., min_length=1, max_length=100, description="Exact gene symbol or numeric VariationID"),
    limit: int = Query(20, ge=1, le=100),
    after_id: int = Query(0, ge=0, le=9223372036854775807),
):
    try:
        result = database.search(q, limit, after_id)
        for item in result['items']:
            item['prediction'] = predict(item.get('latest_details'), item['variation_id'], include_diseases=False)
            item['prediction_support'] = prediction_support(item['prediction'])
        return result
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@app.get("/api/variants/{variation_id}/timeline")
def timeline(variation_id: int = PathParameter(..., ge=1, le=9223372036854775807)):
    result = database.timeline(variation_id)
    if result is None:
        raise HTTPException(status_code=404, detail="Variant not found in the prepared timelines.")
    return result


@app.get("/api/variants/{variation_id}")
def detail(variation_id: int = PathParameter(..., ge=1, le=9223372036854775807)):
    result = database.detail(variation_id)
    if result is None:
        raise HTTPException(status_code=404, detail="Variant not found in the prepared timelines.")
    result['prediction'] = predict(result.get('latest_details'), variation_id)
    result['prediction_support'] = prediction_support(result['prediction'])
    return result


@app.get('/api/model/evaluation')
def model_evaluation():
    try:
        return evaluation()
    except FileNotFoundError as error:
        raise HTTPException(503,'Disease ranking model is not installed.') from error

@app.get('/api/variants/{variation_id}/ranking')
def disease_ranking(variation_id: int = PathParameter(...,ge=1,le=9223372036854775807)):
    record = database.detail(variation_id)
    if record is None: raise HTTPException(404,'Variant not found.')
    return rank(record['latest_details'])

@app.get('/api/model/prediction-evaluation')
def prediction_evaluation():
    try:
        return prediction_bundle()['metadata']
    except FileNotFoundError as error:
        raise HTTPException(503,'The research prediction model is not installed.') from error

@app.get('/api/variants/{variation_id}/annotation-prediction')
def annotation_prediction(variation_id: int):
    record = database.detail(variation_id)
    if record is None:
        raise HTTPException(404, 'Variant not found.')
    return annotation_predict(record.get('latest_details'), variation_id)

@app.get('/api/model/annotation-evaluation')
def annotation_evaluation():
    import json
    card = project_root / 'models/annotation_prediction_v2/model_card.json'
    if not card.is_file():
        return {'status': 'not_configured'}
    return json.loads(card.read_text(encoding='utf-8'))

@app.get('/api/variants/{variation_id}/tp53-evidence')
def variant_tp53_evidence(variation_id: int = PathParameter(..., ge=1, le=9223372036854775807)):
    record = database.detail(variation_id)
    if record is None:
        raise HTTPException(404, 'Variant not found.')
    return tp53_evidence(record.get('latest_details'))

frontend_dist = project_root / 'frontend/dist'
if frontend_dist.exists():
    app.mount('/', StaticFiles(directory=frontend_dist,html=True),name='frontend')
