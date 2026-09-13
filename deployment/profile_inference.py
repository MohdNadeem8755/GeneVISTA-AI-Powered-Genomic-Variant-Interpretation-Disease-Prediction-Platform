"""Measure a fresh local app process; Windows measurements are not Render guarantees."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

def worker():
    import psutil
    from fastapi.testclient import TestClient
    from backend.app.main import app
    proc=psutil.Process()
    timings=[]
    with TestClient(app) as client:
        for path in ['/health','/api/stats','/api/model/prediction-evaluation',
                     '/api/model/evaluation','/api/variants/17660','/api/variants/37565',
                     '/api/variants/search?q=BRCA1&limit=100',
                     '/api/variants/search?q=TP53&limit=100']:
            start=time.perf_counter(); response=client.get(path)
            if response.status_code!=200: raise RuntimeError(f'{path}: {response.status_code}')
            timings.append({'path':path,'seconds':round(time.perf_counter()-start,4)})
    memory=proc.memory_info()
    return {'rss_mib':round(memory.rss/2**20,2),
        'peak_working_set_mib':round(getattr(memory,'peak_wset',memory.rss)/2**20,2),
        'requests':timings,'platform':sys.platform,
        'scope':'Fresh process, sequential local requests including 100-result gene searches. Does not establish Linux/container peak memory, cold-build requirements, concurrency capacity or hosting availability.'}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--worker',action='store_true')
    p.add_argument('--database',type=Path,default=ROOT/'data/demo/genevista.sqlite')
    p.add_argument('--memory-budget-mib',type=float,default=512)
    p.add_argument('--output',type=Path,default=ROOT/'reports/deployment/inference_profile.json')
    args=p.parse_args()
    if args.worker:
        print(json.dumps(worker())); return
    runtime=ROOT/'.runtime'; runtime.mkdir(exist_ok=True)
    env=dict(os.environ,GENEVISTA_DB=str(args.database.resolve()),GENEVISTA_AUTH_ENABLED='false',
        TEMP=str(runtime),TMP=str(runtime),PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    started=time.perf_counter()
    run=subprocess.run([sys.executable,__file__,'--worker'],env=env,cwd=ROOT,capture_output=True,text=True,check=True)
    report=json.loads(run.stdout)
    report.update(database=args.database.name,database_bytes=args.database.stat().st_size,
        predictor_bytes=(ROOT/'models/variant_prediction_v1/predictor.joblib').stat().st_size,
        total_startup_and_requests_seconds=round(time.perf_counter()-started,3),
        requested_memory_budget_mib=args.memory_budget_mib,
        below_80_percent_of_requested_budget=report['peak_working_set_mib']<args.memory_budget_mib*.8)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
