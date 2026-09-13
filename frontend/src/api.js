const API_BASE = (import.meta.env.VITE_API_BASE_URL || (import.meta.env.DEV ? 'http://127.0.0.1:8000' : '')).replace(/\/$/, '');
export async function fetchApi(path, signal) {
  const response = await fetch(API_BASE + path, {
    credentials: 'include',
    signal: signal ? AbortSignal.any([signal, AbortSignal.timeout(20000)]) : AbortSignal.timeout(20000),
  });
  if (!response.ok) {
    if(response.status===401){window.dispatchEvent(new Event('genevista-auth-expired'));throw new Error('Your session has expired. Please sign in.');}
    if (response.status === 404) throw new Error('This variant was not found in the prepared dataset.');
    throw new Error(`The service returned HTTP ${response.status}. Please try again.`);
  }
  return response.json();
}
export function errorMessage(error) {
  if (error.name === 'TimeoutError') return 'The request took too long. Please try again.';
  if (error instanceof TypeError) return 'Cannot reach the data service. Please retry shortly.';
  return error.message || 'Something went wrong. Please try again.';
}
export const stateLabels = {
  vus: 'Uncertain significance', benign: 'Benign / likely benign',
  pathogenic: 'Pathogenic / likely pathogenic', conflicting: 'Conflicting',
  other: 'Other', absent_grch38: 'Absent in GRCh38', present_excluded: 'Present, excluded',
};

export async function authApi(action,body){
 const response=await fetch(API_BASE+'/api/auth/'+action,{method:body?'POST':'GET',credentials:'include',headers:body?{'Content-Type':'application/json','X-GeneVISTA-Request':'1'}:{},body:body?JSON.stringify(body):undefined,signal:AbortSignal.timeout(20000)});
 const data=await response.json();if(!response.ok)throw new Error(typeof data.detail==='string'?data.detail:'Unable to sign in.');return data;
}
