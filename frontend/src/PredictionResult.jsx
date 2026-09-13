import React from 'react';
export function percentage(value){return value>0&&value<.001?'<0.1%':`${(value*100).toFixed(1)}%`;}
export function PredictionSummary({prediction}){
 if(!prediction)return <p role="status">Select a variant to calculate its predicted class.</p>;
 if(prediction.status!=='available')return <p className="prediction-note">{prediction.reason}</p>;
 const change=prediction.change?.nucleotide?.label||prediction.change?.protein?.label;
 return <div className="prediction-result"><div><p className="eyebrow">PREDICTED RESULT · RESEARCH MODEL</p><h3>{prediction.predicted_class==='Uncertain significance'?'VUS — Uncertain significance':prediction.predicted_class}</h3><span className="prediction-confidence">{percentage(prediction.top_probability)} model class estimate</span>{change&&<p className="variant-change">Variant change: <strong>{change}</strong></p>}</div><div className="prediction-guidance"><p>Class probability estimates describe the model output, not the chance a patient will develop a disease.</p>{prediction.warnings.map(text=><p className="prediction-warning" key={text}>{text}</p>)}<small>{prediction.next_step}</small><p className="footnote">Model cohort: {prediction.cohort_membership.replaceAll('_',' ')}. Predictions for training records are not independent validation.</p></div></div>;
}
