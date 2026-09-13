"""Single research-workspace account with expiring, server-side sessions."""
import hashlib
import hmac
import os
import secrets
import time
from collections import OrderedDict
from threading import Lock
from fastapi import APIRouter, HTTPException, Request, Response
from pydantic import BaseModel, Field

router = APIRouter(prefix='/api/auth')
_sessions = {}
_attempts = OrderedDict()
_lock = Lock()
COOKIE = 'genevista_session'
TTL = 8 * 60 * 60

def enabled():
    return os.environ.get('GENEVISTA_AUTH_ENABLED','false').lower() == 'true'

def session(request):
    token = request.cookies.get(COOKIE,'')
    key = hashlib.sha256(token.encode()).hexdigest()
    with _lock:
        record = _sessions.get(key)
        if record and record[1] > time.time(): return record[0]
        _sessions.pop(key,None)
    return None

class Credentials(BaseModel):
    username: str = Field(min_length=1,max_length=128)
    password: str = Field(min_length=1,max_length=256)

@router.get('/session')
def current(request: Request):
    return {'enabled':enabled(),'username':session(request) if enabled() else 'Local researcher'}

@router.post('/login')
def login(body: Credentials, request: Request, response: Response):
    if not enabled(): raise HTTPException(409,'Authentication is disabled in local mode.')
    expected = os.environ.get('GENEVISTA_PASSWORD_HASH','')
    username = os.environ.get('GENEVISTA_USERNAME','')
    if not expected or not username: raise HTTPException(503,'Workspace login has not been configured.')
    address = request.client.host if request.client else 'unknown'
    now = time.time()
    with _lock:
        attempts = [t for t in _attempts.get(address,[]) if t > now-300]
        if len(attempts) >= 10: raise HTTPException(429,'Too many attempts. Try again in five minutes.')
        _attempts[address] = attempts + [now]
        _attempts.move_to_end(address)
        if len(_attempts)>10000: _attempts.popitem(last=False)
    try:
        salt, digest = expected.split(':')
        actual = hashlib.scrypt(body.password.encode(),salt=bytes.fromhex(salt),n=16384,r=8,p=1).hex()
        valid = hmac.compare_digest(actual,digest) & hmac.compare_digest(body.username.encode(),username.encode())
    except (ValueError,TypeError):
        raise HTTPException(503,'Workspace credentials are not configured correctly.')
    if not valid: raise HTTPException(401,'Invalid username or password.')
    token = secrets.token_urlsafe(32)
    with _lock:
        for key, value in list(_sessions.items()):
            if value[1] <= now: del _sessions[key]
        if len(_sessions)>=10000: raise HTTPException(503,'Session capacity reached.')
        _sessions[hashlib.sha256(token.encode()).hexdigest()] = (username,now+TTL)
    response.set_cookie(COOKIE,token,max_age=TTL,httponly=True,
        secure=os.environ.get('GENEVISTA_COOKIE_SECURE','true').lower()=='true',samesite='strict',path='/')
    return {'username':username}

@router.post('/logout')
def logout(request: Request, response: Response):
    with _lock:
        _sessions.pop(hashlib.sha256(request.cookies.get(COOKIE,'').encode()).hexdigest(),None)
    response.delete_cookie(COOKIE,path='/')
    return {'status':'signed_out'}
