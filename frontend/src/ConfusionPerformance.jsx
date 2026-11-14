import React from 'react';
import {summarizeConfusion} from './performanceMetrics.mjs';
const percent=v=>v==null?'N/A':`${(100*v).toFixed(1)}%`;
const shortLabel=label=>label==='Uncertain significance'?'VUS':label;
export default function ConfusionPerformance({data}){
 const [partition,setPartition]=React.useState('validation');
 const [normalized,setNormalized]=React.useState(false);
 const available=data.partition_evaluation?.scores||{validation:data.validation};
 const result=available[partition];
 if(!result?.confusion_matrix)return <p>Class-level results are not available for this partition.</p>;
 const classes=result.class_order;
 const summary=summarizeConfusion(result.confusion_matrix,classes);
 const total=summary.reduce((n,r)=>n+r.support,0);
 return <div>
 <h3>Where does the model make mistakes?</h3>
 <label>Evaluation dataset: <select value={partition} onChange={e=>setPartition(e.target.value)}>
 {Object.entries({validation:'Validation — used for selection',training:'Training — fitting records',early_stopping:'Early stopping',internal_test:'Historical internal test'}).filter(([key])=>available[key]).map(([key,label])=><option key={key} value={key}>{label}</option>)}
 </select></label>{' '}
 <label><input type="checkbox" checked={normalized} onChange={e=>setNormalized(e.target.checked)}/> Show percentages within each actual class</label>
 <p>Rows = actual recorded classes. Columns = model predictions. The diagonal contains correct predictions; every other cell is a mistake. Counts and row percentages describe the same records.</p>
 <div className="table-scroll"><table><caption>Confusion matrix — {total.toLocaleString()} records</caption><thead><tr><th scope="col">Actual ↓ / Predicted →</th>{classes.map(c=><th scope="col" key={c}>{shortLabel(c)}</th>)}<th scope="col">Actual records</th></tr></thead>
 <tbody>{classes.map((label,i)=><tr key={label}><th scope="row">{shortLabel(label)}</th>{result.confusion_matrix[i].map((count,j)=><td key={classes[j]} style={{backgroundColor:i===j?'rgba(26,180,164,0.18)':undefined,fontWeight:i===j?700:400}} title={`${label} predicted as ${classes[j]}: ${count} records`}>
 {normalized?percent(summary[i].support?count/summary[i].support:null):count.toLocaleString()}</td>)}<td>{summary[i].support.toLocaleString()}</td></tr>)}</tbody></table></div>
 <h3>Performance for every class</h3>
 <p>Recall: how many actual examples were found. Precision: how many predictions of that class were correct. F1 combines both. These are different from overall accuracy.</p>
 <div className="table-scroll"><table><thead><tr><th>Class</th><th>Actual records</th><th>Predicted records</th><th>Correct</th><th>Precision</th><th>Recall</th><th>F1</th></tr></thead><tbody>{summary.map(r=><tr key={r.label}><th scope="row">{shortLabel(r.label)}</th><td>{r.support.toLocaleString()} ({percent(r.actualShare)})</td><td>{r.predicted.toLocaleString()} ({percent(r.predictedShare)})</td><td>{r.correct.toLocaleString()}</td><td>{percent(r.precision)}</td><td>{percent(r.recall)}</td><td>{r.f1==null?'N/A':r.f1.toFixed(3)}</td></tr>)}</tbody></table></div>
 <p>Class imbalance describes the dataset. Comparing actual and predicted counts helps identify over-prediction and under-prediction; neither imbalance nor a single score proves the cause of an error.</p>
 </div>;
}
