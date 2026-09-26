import React from 'react';
import useCurrentResource from './useCurrentResource.js';
import {percentage} from './PredictionResult.jsx';
export default function AnnotationEvaluation(){
 const {data,error,refresh,checkedAt}=useCurrentResource('/api/model/annotation-evaluation');
 if(data?.status==='not_configured')return null;
 return <section className="panel model-charts" id="annotation-model-quality">
 <p className="eyebrow">ACTIVE MODEL · LATEST SAVED EVALUATION</p><h2>Current model accuracy</h2>
 <button className="secondary" onClick={refresh}>Refresh metrics</button>
 {error?<p role="alert">{error}</p>:!data?<p role="status">Loading annotation evaluation…</p>:<>
 <p>These metrics belong to the active annotation model. They refresh every minute and when you return to this page. Accuracy changes after a newly evaluated model is installed, not after each prediction.</p>
 <p className="footnote">Model: {data.model_id} · Checked {checkedAt?.toLocaleTimeString()}</p>
 {data.incumbent_same_test&&<details><summary>Previous model comparison</summary><p>On the same test records, the previous model scored {percentage(data.incumbent_same_test.accuracy)} accuracy and {percentage(data.incumbent_same_test.balanced_accuracy)} balanced accuracy.</p></details>}
 <div className="evaluation-stats"><div><small>Accuracy</small><strong>{(data.test.accuracy*100).toFixed(2)}%</strong></div><div><small>Balanced accuracy</small><strong>{percentage(data.test.balanced_accuracy)}</strong></div><div><small>Macro F1</small><strong>{data.test.macro_f1.toFixed(3)}</strong></div><div><small>Matched test records</small><strong>{data.test_count.toLocaleString()}</strong></div></div>
 <div className="table-scroll"><table><thead><tr><th>Class</th><th>Precision</th><th>Recall</th><th>Test support</th></tr></thead><tbody>{['Pathogenic','Likely pathogenic','Uncertain significance','Likely benign','Benign'].map(label=>{const row=data.test.per_class[label];return <tr key={label}><td>{label}</td><td>{percentage(row.precision)}</td><td>{percentage(row.recall)}</td><td>{row.support.toLocaleString()}</td></tr>;})}</tbody></table></div>
 <details><summary>Evaluation scope and limitations</summary><p>{data.scope}</p><p>{data.upstream_overlap}</p><p>{data.calibration_count.toLocaleString()} separate calibration records. Test labels were not used to select the model. {data.test_policy || "Internal validation does not establish clinical accuracy."}</p></details>
 </>}
 </section>;
}
