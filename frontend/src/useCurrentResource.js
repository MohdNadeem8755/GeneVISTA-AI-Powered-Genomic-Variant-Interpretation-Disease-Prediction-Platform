import {useCallback, useEffect, useState} from 'react';
import {fetchApi, errorMessage} from './api.js';

// Refresh saved model results, never invent a live accuracy from predictions.
export default function useCurrentResource(path) {
 const [state,setState]=useState({data:null,error:'',checkedAt:null});
 const [attempt,setAttempt]=useState(0);
 const refresh=useCallback(()=>setAttempt(value=>value+1),[]);
 useEffect(()=>{
  if(!path){setState({data:null,error:'',checkedAt:null});return;}
  let controller;
  let disposed=false;
  const load=()=>{
   controller?.abort();controller=new AbortController();
   const current=controller;
   fetchApi(path,current.signal).then(data=>{
    if(!disposed&&!current.signal.aborted)setState({data,error:'',checkedAt:new Date()});
   }).catch(error=>{
    if(!disposed&&!current.signal.aborted)setState({data:null,error:errorMessage(error),checkedAt:null});
   });
  };
  setState({data:null,error:'',checkedAt:null});load();
  const interval=window.setInterval(()=>{if(!document.hidden)load();},60000);
  window.addEventListener('focus',load);
  return()=>{disposed=true;controller?.abort();clearInterval(interval);window.removeEventListener('focus',load);};
 },[path,attempt]);
 return {...state,refresh};
}
