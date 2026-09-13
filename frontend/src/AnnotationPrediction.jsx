import React, {useEffect, useState} from 'react';
import {percentage} from './PredictionResult.jsx';
export default function AnnotationPrediction({selected}) {
 const id=selected?.data?.variation_id;
 const [state,setState]=useState(null);
 useEffect(()=>{
  if(!id){setState(null);return;}
  const controller=new AbortController();setState(null);
  fetch(`/api/variants/${id}/annotation-prediction`,{signal:controller.signal})
   .then(async response=>{if(!response.ok)throw new Error('Annotation prediction could not be loaded.');return response.json();})
   .then(setState).catch(error=>{if(error.name!=='AbortError')setState({reason:error.message});});
  return ()=>controller.abort();
 },[id]);
 if(!id)return null;
 return <section className="panel" aria-label="Annotation model prediction">
  <p className="eyebrow">DBNSFP 5.4A · ADDITIONAL RESEARCH MODEL</p>
  <h2>Annotation-based predicted result</h2>
  {!state?<p role="status">Loading annotation evidence…</p>:state.status!=='available'?<p>{state.reason}</p>:<>
   <h3>{state.predicted_class==='Uncertain significance'?'VUS — Uncertain significance':state.predicted_class}</h3>
   <p>{percentage(state.top_probability)} model class estimate · {state.allele}</p>
   <div className="class-spectrum">{state.probabilities.map(row=><div key={row.class_name}><strong>{row.class_name}</strong><p>{percentage(row.probability)}</p><progress value={row.probability} max="1" aria-label={`${row.class_name} probability`}/></div>)}</div>
   {state.disagrees_with_recorded&&<p className="prediction-warning">This estimate differs from the recorded classification. Review the original evidence.</p>}
   <p>Test accuracy: {percentage(state.evaluation.accuracy)} · Balanced accuracy: {percentage(state.evaluation.balanced_accuracy)} · Macro F1: {state.evaluation.macro_f1.toFixed(3)} · {state.test_count.toLocaleString()} matched test variants.</p>
   <p>Model cohort: {state.cohort_membership}. Training-record predictions are not independent validation.</p>
   <details><summary>Coverage and interpretation</summary><p>{state.scope}</p><p>{state.caution}</p><p>Class estimates are not the probability of developing a disease. This model is separate from the sequence-based disease annotation rankings.</p></details>
  </>}
 </section>;
}
