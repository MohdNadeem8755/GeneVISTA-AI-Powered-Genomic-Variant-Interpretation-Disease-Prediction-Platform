import React, {useEffect,useState} from 'react';
import {DnaLogo} from './DnaAnimation.jsx';
import {authApi} from './api.js';
export default function Login({children}) {
 const [session,setSession]=useState(null),[error,setError]=useState(''),[busy,setBusy]=useState(false);
 async function check(){setError('');try{setSession(await authApi('session'));}catch(e){setError(e.message);}}
 useEffect(()=>{check();const expire=()=>setSession({enabled:true,username:null});window.addEventListener('genevista-auth-expired',expire);return()=>window.removeEventListener('genevista-auth-expired',expire);},[]);
 async function submit(e){e.preventDefault();setBusy(true);setError('');const data=new FormData(e.currentTarget);try{const user=await authApi('login',{username:data.get('username'),password:data.get('password')});setSession({enabled:true,...user});}catch(e){setError(e.message);}finally{setBusy(false);}}
 async function logout(){setBusy(true);try{await authApi('logout',{});setSession({enabled:true,username:null});}catch(e){setError(e.message);}finally{setBusy(false);}}
 if(session && (!session.enabled || session.username)) return <>{session.enabled&&<div className="account-bar"><span>{session.username}</span><button disabled={busy} onClick={logout}>Sign out</button>{error&&<span role="alert">{error}</span>}</div>}{children}</>;
 return <div className="login-page"><section className="login-card"><DnaLogo motion={true}/><p className="eyebrow">GENEVISTA RESEARCH WORKSPACE</p><h1>Welcome back.</h1><p>Sign in to explore variant evidence and disease annotation rankings.</p>{!session?<><p role="status">{error||'Connecting to the workspace…'}</p>{error&&<button onClick={check}>Retry connection</button>}</>:<form onSubmit={submit}><label>Username<input name="username" autoComplete="username" required maxLength={128}/></label><label>Password<input name="password" type="password" autoComplete="current-password" required maxLength={256}/></label>{error&&<p role="alert">{error}</p>}<button disabled={busy}>{busy?'Signing in…':'Sign in'}</button></form>}<small>Access is provided by your workspace administrator.</small></section></div>;
}
