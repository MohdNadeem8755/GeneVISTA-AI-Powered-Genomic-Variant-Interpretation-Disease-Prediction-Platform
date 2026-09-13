"""Validate the bundled demo or provision an optional full database, then serve."""
import hashlib
import os
from pathlib import Path
import re
import sqlite3
from contextlib import closing
import sys
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

def provision():
    target=Path(os.environ.get('GENEVISTA_DB', ROOT/'data/app/genevista.sqlite'))
    expected=os.environ.get('GENEVISTA_DB_SHA256','').lower()
    if not target.exists():
        url=os.environ.get('GENEVISTA_DB_URL','')
        if not url.startswith('https://') or not re.fullmatch('[0-9a-f]{64}',expected):
            raise RuntimeError('Set an HTTPS GENEVISTA_DB_URL and GENEVISTA_DB_SHA256, or provision the database manually.')
        target.parent.mkdir(parents=True,exist_ok=True)
        temporary=target.with_suffix('.download')
        digest=hashlib.sha256()
        try:
            with urllib.request.urlopen(url,timeout=120) as source, temporary.open('wb') as output:
                if not source.url.startswith('https://'): raise RuntimeError('Database redirect must use HTTPS.')
                for chunk in iter(lambda:source.read(1024*1024),b''):
                    output.write(chunk);digest.update(chunk)
            if digest.hexdigest()!=expected: raise RuntimeError('Database checksum mismatch.')
            temporary.replace(target)
        finally:
            temporary.unlink(missing_ok=True)
    if expected:
        with target.open('rb') as source:
            if hashlib.file_digest(source,'sha256').hexdigest()!=expected:
                raise RuntimeError('Installed database checksum differs from configured version.')
    with closing(sqlite3.connect(target.as_uri()+'?mode=ro',uri=True)) as connection:
        if connection.execute('PRAGMA quick_check').fetchone()[0]!='ok':
            raise RuntimeError('Database integrity check failed.')
    from backend.app.database import Database
    Database(target).stats()
    from backend.app.ranking import bundle
    bundle()
    from backend.app.auth import enabled
    if enabled():
        if not os.environ.get('GENEVISTA_USERNAME') or not re.fullmatch('[0-9a-f]{32}:[0-9a-f]{128}',os.environ.get('GENEVISTA_PASSWORD_HASH','')):
            raise RuntimeError('Configure the workspace username and password hash before starting.')

if __name__=='__main__':
    provision()
    os.chdir(ROOT)
    os.execv(sys.executable,[sys.executable,'-m','uvicorn','backend.app.main:app','--host','0.0.0.0','--port',os.environ.get('PORT','8000'),'--workers','1'])
