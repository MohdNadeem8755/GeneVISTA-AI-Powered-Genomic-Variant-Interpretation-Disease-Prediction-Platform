import React, {useEffect,useState} from 'react';
import {percentage} from './PredictionResult.jsx';
export default function AnnotationEvaluation(){
 const [data,setData]=useState(null),[error,setError]=useState('');
 useEffect(()=>{const c=new AbortController();fetch('/api/model/annotation-evaluation',{signal:c.signal}).then(r=>{if(!r.ok)throw new Error('Could not load annotation evaluation.');return r.json();}).then(setData).catch(e=>{if(e.name!=='AbortError')setError(e.message);});return()=>c.abort();},[]);
 if(data?.status==='not_configured')return null;
 return <section className="panel model-charts" id="annotation-model-quality">
 <p className="eyebrow">NEW ANNOTATION MODEL · DBNSFP 5.4A</p><h2>Annotation model quality</h2>
 {error?<p role="alert">{error}</p>:!data?<p role="status">Loading annotation evaluation…</p>:<>
 <p>Used for the additional annotation-based prediction when an exact cached allele match is available. The sequence model below has its own evaluation. These test cohorts differ, so their accuracy figures are not a direct comparison.</p>
 <div className="evaluation-stats"><div><small>Accuracy</small><strong>{(data.test.accuracy*100).toFixed(2)}%</strong></div><div><small>Balanced accuracy</small><strong>{percentage(data.test.balanced_accuracy)}</strong></div><div><small>Macro F1</small><strong>{data.test.macro_f1.toFixed(3)}</strong></div><div><small>Matched test records</small><strong>{data.test_count.toLocaleString()}</strong></div></div>
 <div className="table-scroll"><table><thead><tr><th>Class</th><th>Precision</th><th>Recall</th><th>Test support</th></tr></thead><tbody>{['Pathogenic','Likely pathogenic','Uncertain significance','Likely benign','Benign'].map(label=>{const row=data.test.per_class[label];return <tr key={label}><td>{label}</td><td>{percentage(row.precision)}</td><td>{percentage(row.recall)}</td><td>{row.support.toLocaleString()}</td></tr>;})}</tbody></table></div>
 <details><summary>Evaluation scope and limitations</summary><p>{data.scope}</p><p>{data.upstream_overlap}</p><p>{data.calibration_count.toLocaleString()} separate calibration records. Test labels were not used to select the model. 99% accuracy was not achieved.</p></details>
 </>}
 </section>;
}
