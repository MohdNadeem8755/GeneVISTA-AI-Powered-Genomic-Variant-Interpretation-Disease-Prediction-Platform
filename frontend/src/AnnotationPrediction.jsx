import React from 'react';
import useCurrentResource from './useCurrentResource.js';
import {percentage} from './PredictionResult.jsx';
export default function AnnotationPrediction({selected}) {
 const id=selected?.data?.variation_id;
 const {data:state,error,refresh,checkedAt}=useCurrentResource(id?`/api/variants/${id}/annotation-prediction`:null);
 if(!id)return null;
 return <section className="panel" aria-label="Annotation model prediction">
  <p className="eyebrow">CURRENT ANNOTATION MODEL · RESEARCH</p>
  <h2>Current model prediction</h2>
  <button className="secondary" onClick={refresh}>Refresh prediction</button>
  {error?<p role="alert">{error}</p>:!state?<p role="status">Loading annotation evidence…</p>:state.status!=='available'?<p>{state.reason}</p>:<>
   <p className="footnote">Model: {state.model_id} · Checked {checkedAt?.toLocaleTimeString()}</p>
   <h3>{state.predicted_class==='Uncertain significance'?'VUS — Uncertain significance':state.predicted_class}</h3>
   <p>{percentage(state.top_probability)} class probability estimate (not model accuracy) · {state.allele}</p>
   {state.decision_policy&&<p>{state.decision_policy}{state.most_probable_class!==state.predicted_class&&` Highest probability: ${state.most_probable_class}.`}</p>}
   <div className="class-spectrum">{state.probabilities.map(row=><div key={row.class_name}><strong>{row.class_name}</strong><p>{percentage(row.probability)}</p><progress value={row.probability} max="1" aria-label={`${row.class_name} probability`}/></div>)}</div>
   {state.disagrees_with_recorded&&<p className="prediction-warning">This estimate differs from the recorded classification. Review the original evidence.</p>}
   <p>Saved test accuracy: {percentage(state.evaluation.accuracy)} · Balanced accuracy: {percentage(state.evaluation.balanced_accuracy)} · Macro F1: {state.evaluation.macro_f1.toFixed(3)} · {state.test_count.toLocaleString()} matched test variants.</p>
   <p>Model cohort: {state.cohort_membership}. Training-record predictions are not independent validation.</p>
   <details><summary>Coverage and interpretation</summary><p>{state.scope}</p><p>{state.caution}</p><p>Class estimates are not the probability of developing a disease. This model is separate from the sequence-based disease annotation rankings.</p></details>
  </>}
 </section>;
}
