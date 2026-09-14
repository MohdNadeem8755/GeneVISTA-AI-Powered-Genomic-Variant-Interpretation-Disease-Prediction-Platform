import React, {useEffect, useState} from 'react';
import {fetchApi, errorMessage} from './api.js';

export default function Brca1Evidence({selected}) {
  const id = selected?.data?.variation_id;
  const isBrca1 = selected?.data?.latest_details?.gene_symbol === 'BRCA1';
  const [state, setState] = useState(null);
  useEffect(() => {
    if (!id || !isBrca1) { setState(null); return; }
    const controller = new AbortController(); setState(null);
    fetchApi(`/api/variants/${id}/brca1-evidence`, controller.signal)
      .then(data => { if (!controller.signal.aborted) setState({id, data}); })
      .catch(error => { if (!controller.signal.aborted) setState({id, error: errorMessage(error)}); });
    return () => controller.abort();
  }, [id, isBrca1]);
  if (!id || !isBrca1) return null;
  const current = state?.id === id ? state : null;
  const data = current?.data;
  return <section id="brca1-evidence" className="panel" aria-label="BRCA1 functional evidence">
    <p className="eyebrow">BRCA1 · MaveDB FUNCTIONAL EVIDENCE</p>
    <h2>Measured functional assay evidence</h2>
    {!current ? <p role="status">Matching the public BRCA1 score set…</p> : current.error ? <p role="alert">{current.error}</p> : data.status !== 'available' ? <p>{data.reason}</p> : <>
      <p><strong>{data.variant.coding_change}</strong> · {data.score_set}</p>
      <p>This saturation-genome-editing score is an auxiliary functional result. It is not a five-class prediction, diagnosis, or patient disease probability.</p>
      <div className="classification-summary">
        <div><span>Functional score</span><h3>{data.score.toFixed(3)}</h3><p>{data.functional_class}</p></div>
        <div><span>Research criterion candidate</span><h3>{data.criterion_candidate || 'No criterion'}</h3><p>Requires expert review before any use.</p></div>
      </div>
      {data.rna_score !== null && <p>RNA score: <strong>{data.rna_score.toFixed(3)}</strong> · same assay family, not an independent evidence item.</p>}
      <details><summary>Source and limitations</summary><p>{data.limitations}</p>{data.sources.map(source => <p key={source.url}><a href={source.url} target="_blank" rel="noreferrer">{source.label} ↗</a></p>)}</details>
    </>}
  </section>;
}
