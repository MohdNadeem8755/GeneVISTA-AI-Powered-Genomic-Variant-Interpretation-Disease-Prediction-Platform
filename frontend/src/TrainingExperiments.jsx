import React from 'react';
import useCurrentResource from './useCurrentResource.js';
import {percentage} from './PredictionResult.jsx';
export default function TrainingExperiments(){
 const {data,error,refresh,checkedAt}=useCurrentResource('/api/model/research-progress');
 if(!error&&!data?.experiments?.length)return null;
 return <details className="legacy-experiment"><summary>Latest training experiments: CatBoost + protein language models</summary>
 <button className="secondary" onClick={refresh}>Refresh training results</button>
 {error?<p role="alert">{error}</p>:<><p>{data.note}</p><p>Checked {checkedAt?.toLocaleTimeString()}</p>
 <div className="table-scroll"><table><thead><tr><th>Experiment</th><th>Status</th><th>Validation accuracy</th><th>Validation balanced accuracy</th><th>Training records</th></tr></thead>
 <tbody>{data.experiments.map(row=><tr key={row.name}><td>{row.name}</td><td>{row.status.replaceAll('_',' ')}</td><td>{row.metrics?percentage(row.metrics.accuracy):'Pending'}</td><td>{row.metrics?percentage(row.metrics.balanced_accuracy):'Pending'}</td><td>{row.train_rows?.toLocaleString()||'Preparing'}</td></tr>)}</tbody></table></div></>}
 </details>;
}
