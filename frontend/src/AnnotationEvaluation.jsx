import React from 'react';
import useCurrentResource from './useCurrentResource.js';
import {percentage} from './PredictionResult.jsx';
import ConfusionPerformance from './ConfusionPerformance.jsx';
const metricPercent=value=>value==null?'Not calculated':`${(value*100).toFixed(2)}%`;
function CorrectedComparison({data}){
 const candidate=data.corrected_candidate;
 if(!candidate)return null;
 const newer=candidate.validation;
 const partitions=candidate.partition_evaluation?.scores;
 return <article>
 <h3>New corrected model — validation comparison</h3>
 <div className="evaluation-stats" aria-label="Latest candidate validation metrics"><div><small>Overall accuracy · validation</small><strong>{metricPercent(newer.accuracy)}</strong></div><div><small>Balanced accuracy · validation</small><strong>{metricPercent(newer.balanced_accuracy)}</strong></div><div><small>Macro F1 · validation</small><strong>{newer.macro_f1.toFixed(3)}</strong></div><div><small>Validation records</small><strong>{newer.n.toLocaleString()}</strong></div></div>
 <h3>Latest model — training, validation and testing</h3>
 <div className="table-scroll"><table aria-label="Latest candidate partition scores"><thead><tr><th>Partition</th><th>Overall accuracy</th><th>Balanced accuracy</th><th>Macro F1</th><th>Status</th></tr></thead><tbody>
 <tr><th scope="row">Training</th><td>{metricPercent(partitions?.training?.accuracy)}</td><td>{metricPercent(partitions?.training?.balanced_accuracy)}</td><td>{partitions?.training?.macro_f1?.toFixed(3)??'Not calculated'}</td><td>{partitions?`${partitions.training.n.toLocaleString()} fitting records · Seen during training.`:'Saved ensemble training scores are not available.'}</td></tr>
 <tr><th scope="row">Validation</th><td>{metricPercent(newer.accuracy)}</td><td>{metricPercent(newer.balanced_accuracy)}</td><td>{newer.macro_f1.toFixed(3)}</td><td>{newer.n.toLocaleString()} records · Used to select this candidate.</td></tr>
 <tr><th scope="row">Testing · historical internal</th><td>{partitions?metricPercent(partitions.internal_test.accuracy):'Not evaluated'}</td><td>{partitions?metricPercent(partitions.internal_test.balanced_accuracy):'Not evaluated'}</td><td>{partitions?.internal_test?.macro_f1?.toFixed(3)??'Not evaluated'}</td><td>{partitions?`${partitions.internal_test.n.toLocaleString()} records · Previously inspected historical partition, not an untouched final test.`:'This candidate has no saved test evaluation.'}</td></tr>
 </tbody></table></div>
 <p className="footnote">All rows refer to {candidate.model_id}. Previous-model scores are not substituted for missing results.</p>
 <p className="footnote">New: {candidate.model_id} · Previous: {data.model_id} · Same {newer.n.toLocaleString()} validation records</p>
 <p>Precision improved for the two “Likely” classes, but recall decreased: fewer false predictions, more actual cases missed. The previous model is preserved; the new candidate is not deployed.</p>
 <div className="table-scroll"><table><thead><tr><th>Class</th><th>Previous precision</th><th>New precision</th><th>Previous recall</th><th>New recall</th></tr></thead><tbody>{newer.class_order.map(label=>{
 const before=data.validation.classification_report[label], after=newer.classification_report[label];
 return <tr key={label}><th scope="row">{label}</th><td>{metricPercent(before.precision)}</td><td><strong>{metricPercent(after.precision)}</strong></td><td>{metricPercent(before.recall)}</td><td>{metricPercent(after.recall)}</td></tr>;
 })}</tbody></table></div>
 <p>Overall accuracy counts correct predictions across all five classes. Balanced accuracy averages recall equally across those classes. These results belong to the new candidate, not the active prediction model. {partitions?'The frozen ensemble was evaluated without retraining or tuning on test results.':'Training and testing scores have not yet been calculated.'}</p>
 <details><summary>New model confusion matrix and per-class results</summary><ConfusionPerformance data={{validation:newer}}/></details>
 <p>Changes: 57 features including exact DNA/protein substitutions and gene identity; corrected event parsing; stopping and model selection prioritise macro F1.</p>
 </article>;
}
function PartitionScores({data,advanced=false}){
 const scores=data.partition_evaluation?.scores;
 const rows=[
 ['training','Training / fitting',data.data_audit.fit_rows,'Seen during fitting; not generalization performance'],
 ['early_stopping','Early stopping',data.data_audit.stopping_rows,'Used to choose iteration counts'],
 ['validation','Validation',data.data_audit.validation_rows,'Used to select the ensemble'],
 ['internal_test','Internal test',null,'Historical test partition; previously inspected in this project']];
 return <div className="table-scroll"><table><thead><tr><th>Partition</th><th>Records</th><th>Overall accuracy</th>{advanced&&<><th>Balanced accuracy</th><th>Macro F1</th></>}<th>Meaning</th></tr></thead><tbody>{rows.map(([key,label,count,note])=>{
 const result=scores?.[key]||(key==='validation'?data.validation:null);
 return <tr key={key}><td>{label}</td><td>{(result?.n??count)?.toLocaleString()??'Not evaluated'}</td><td>{metricPercent(result?.accuracy)}</td>{advanced&&<><td>{metricPercent(result?.balanced_accuracy)}</td><td>{result?.macro_f1?.toFixed(3)??'Not calculated'}</td></>}<td>{note}</td></tr>;
 })}</tbody></table></div>;
}
export default function AnnotationEvaluation(){
 const [view,setView]=React.useState('latest');
 const latest=useCurrentResource('/api/model/latest-training-evaluation');
 const {data,error,refresh,checkedAt}=useCurrentResource('/api/model/annotation-evaluation');
 if(view==='latest')return <section className="panel model-charts" id="annotation-model-quality">
 <p className="eyebrow">MODEL COMPARISON · CLASS-BY-CLASS PERFORMANCE</p>
 <h2>Model predictions and mistakes</h2>
 <p>Compare the new candidate with the preserved previous ensemble below. Individual variant predictions still use the active annotation model. Historical internal testing is not an untouched later or external evaluation.</p>
 <button className="secondary" onClick={latest.refresh}>Refresh metrics</button>{' '}
 <button className="secondary" onClick={()=>setView('active')}>View active model test results</button>
 {latest.error?<p role="alert">{latest.error}</p>:latest.data?.status==='not_configured'?<p>No completed training report is available.</p>:!latest.data?<p role="status">Loading training results…</p>:<>
 <CorrectedComparison data={latest.data}/>
 <details><summary>Previous ensemble — preserved results and training details</summary>
 <p className="footnote">Previous run: {latest.data.model_id}</p>
 <ConfusionPerformance data={latest.data}/>
 {latest.data.partition_evaluation&&<p className="footnote">Scores calculated: {new Date(latest.data.partition_evaluation.evaluated_at).toLocaleString()} · Model unchanged during evaluation · No overlapping genomic loci across train/validation/test.</p>}
 <p>Untouched later-dataset evaluation remains incomplete because that cohort had too few Benign and Likely benign examples. These internal scores do not establish clinical performance.</p>
 <details open><summary>Data cleaning and feature preparation</summary><ol>
 <li>Filter to supported GRCh38 germline/unknown SNVs; consolidate duplicate alleles and exclude conflicting records and unsupported labels.</li>
 <li>Match biological annotations by allele. Preserve missing numeric values, add frequency transformations and missingness flags, and encode missing categories explicitly.</li>
 <li>Use 51 population, conservation, sequence and biological-context features. Exclude unaudited pretrained pathogenicity scores from this run.</li>
 </ol><p>The newer-label audit found 27,308 changed labels and 11,737 now-conflicting development records (overlapping counts). These were flagged for review, not automatically corrected or removed. Missing assay evidence remains missing; negative measured assay scores are preserved. The proposed quality-based dataset reduction has not been performed.</p></details>
 <details open><summary>Model training and evaluation</summary><ol>
 <li>Keep genomic loci together. Sample approximately 20% of existing training and validation partitions; reserve approximately 10% of sampled training loci for early stopping.</li>
 <li>Fit eight CatBoost candidates with depths 6/8, class-weight powers 0–1, and up to 800 iterations with early stopping. No synthetic examples were used.</li>
 <li>Compare 84 pairwise ensembles on validation; select by balanced accuracy, then macro F1. Selected ensemble combines two models with 2:1 voting weights.</li>
 <li>Calculate partition scores using that saved ensemble with its weights and decision rule fixed. Training scores use fitting records only; the early-stopping subset is shown separately.</li>
 </ol><p>Original nominal split: 65% training, 15% validation, 10% calibration, 10% test, with legacy-locus exclusions. This run does not establish independent test or clinical performance.</p></details>
 <h3>Overall training and validation accuracy</h3>
 <p>These percentages count correct predictions across all classes, not only VUS. Read them together with the class-level results above.</p>
 <PartitionScores data={latest.data}/>
 <details><summary>Additional summary metrics</summary><p>Balanced accuracy averages recall equally across classes. Macro F1 averages each class's precision–recall F1 score. Neither replaces the confusion matrix.</p><PartitionScores data={latest.data} advanced/></details>
 </details>
 </>}
 </section>;
 if(data?.status==='not_configured')return null;
 return <section className="panel model-charts" id="annotation-model-quality">
 <p className="eyebrow">ACTIVE PREDICTION MODEL · SAVED TEST RESULTS</p><h2>Active model test performance</h2>
 <button className="secondary" onClick={()=>setView('latest')}>View latest ensemble partition scores</button>{' '}
 <button className="secondary" onClick={refresh}>Refresh metrics</button>
 {error?<p role="alert">{error}</p>:!data?<p role="status">Loading annotation evaluation…</p>:<>
 <p>These metrics belong to the active annotation model. They refresh every minute and when you return to this page. Accuracy changes after a newly evaluated model is installed, not after each prediction.</p>
 <p className="footnote">Model: {data.model_id} · Checked {checkedAt?.toLocaleTimeString()}</p>
 <div className="evaluation-stats"><div><small>Accuracy</small><strong>{(data.test.accuracy*100).toFixed(2)}%</strong></div><div><small>Balanced accuracy</small><strong>{percentage(data.test.balanced_accuracy)}</strong></div><div><small>Macro F1</small><strong>{data.test.macro_f1.toFixed(3)}</strong></div><div><small>Matched test records</small><strong>{data.test_count.toLocaleString()}</strong></div></div>
 <div className="table-scroll"><table><thead><tr><th>Class</th><th>Precision</th><th>Recall</th><th>Test support</th></tr></thead><tbody>{['Pathogenic','Likely pathogenic','Uncertain significance','Likely benign','Benign'].map(label=>{const row=data.test.per_class[label];return <tr key={label}><td>{label}</td><td>{percentage(row.precision)}</td><td>{percentage(row.recall)}</td><td>{row.support.toLocaleString()}</td></tr>;})}</tbody></table></div>
 <details><summary>Evaluation scope and limitations</summary><p>{data.scope}</p><p>{data.upstream_overlap}</p><p>{data.calibration_count.toLocaleString()} separate calibration records. Test labels were not used to select the model. {data.test_policy || "Internal validation does not establish clinical accuracy."}</p></details>
 </>}
 </section>;
}
