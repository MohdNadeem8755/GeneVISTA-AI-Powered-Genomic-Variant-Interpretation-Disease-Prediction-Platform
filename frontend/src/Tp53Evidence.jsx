import React, {useEffect, useState} from 'react';
import {fetchApi, errorMessage} from './api.js';

const criterion = value => value && value !== 'No evidence' ? value.replaceAll('_', ' ') : 'No criterion';

export default function Tp53Evidence({selected}) {
  const id = selected?.data?.variation_id;
  const isTp53 = selected?.data?.latest_details?.gene_symbol === 'TP53';
  const [state, setState] = useState(null);
  useEffect(() => {
    if (!id || !isTp53) {setState(null); return;}
    const controller = new AbortController();
    setState(null);
    fetchApi(`/api/variants/${id}/tp53-evidence`, controller.signal)
      .then(data => {if (!controller.signal.aborted) setState({id, data});})
      .catch(error => {if (!controller.signal.aborted) setState({id, error: errorMessage(error)});});
    return () => controller.abort();
  }, [id, isTp53]);
  if (!id || !isTp53) return null;
  const current = state?.id === id ? state : null;
  const data = current?.data;
  return <section id="tp53-evidence" className="panel" aria-label="TP53 evidence review">
    <p className="eyebrow">TP53 · FUNCTIONAL EVIDENCE PILOT</p>
    <h2>Evidence behind this DNA change</h2>
    {!current ? <p role="status">Matching published TP53 evidence…</p> : current.error ? <p role="alert">{current.error}</p> : data.status !== 'available' ? <p>{data.reason}</p> : <>
      <p><strong>{data.variant.transcript}:{data.variant.coding_change}</strong> · {data.variant.protein_change}<br/>
        Condition under review: {data.condition.name} · {data.condition.inheritance}</p>
      <p>Published computational and laboratory findings help review this variant. They do not provide a final clinical classification or a patient’s chance of disease.</p>
      <div className="classification-summary">
        <div><span>Computational criterion candidate</span>
          <h3>{criterion(data.computational.candidate.code)}</h3>
          <p>{data.computational.candidate.reason}</p>
          <small>TP53 specification {data.rule_version} · Not applied</small>
        </div>
        <div><span>Published functional criterion</span>
          <h3>{data.functional.status === 'matched' ? criterion(data.functional.published_code) : 'No unambiguous assay match'}</h3>
          <p>{data.functional.criterion_status === 'withheld' ? 'Withheld for splice/dependency review' : 'Requires review against current functional rules'}</p>
          <small>Published preliminary code · Not applied</small>
        </div>
      </div>
      {data.functional.hold_reasons.map(reason => <p className="prediction-warning" key={reason}>{reason}</p>)}
      {data.conflicting_evidence_directions && <p className="prediction-warning">Computational and functional evidence point in different directions. Resolve this difference before classification.</p>}
      <div className="tp53-predictors">
        <span>Align-GVGD <strong>{data.computational.align_gvgd}</strong></span>
        <span>BayesDel <strong>{data.computational.bayesdel.toFixed(4)}</strong></span>
        <span>SpliceAI delta <strong>{data.computational.spliceai.toFixed(2)}</strong></span>
      </div>
      <p className="evidence-note">These are predictor scores, not disease probabilities. PP3/BP4 describe computational evidence; PS3/BS3 describe functional evidence.</p>
      {data.functional.status === 'matched' && <div className="table-scroll"><table>
        <caption>Published assay findings · {data.functional.match_basis}</caption>
        <thead><tr><th scope="col">Assay</th><th scope="col">Reported functional result</th></tr></thead>
        <tbody>{Object.entries(data.functional.assays).map(([name, result]) => <tr key={name}><th scope="row">{name}</th><td>{result === 'NA' ? 'No assay result' : result === 'LOF' ? 'Loss of function' : result === 'noLOF' ? 'No loss of function reported' : result}</td></tr>)}</tbody>
      </table></div>}
      <details><summary>Evidence sources and remaining review</summary>
        <p>Computational data: {data.computational.source_sheet}, row {data.computational.source_row}.
          {data.functional.rows.length > 0 && ` Functional data: Table S3, rows ${data.functional.rows.map(row => row.source_row).join(', ')}.`}</p>
        <p>Transcript coding sequences .5 and .6 were checked against NCBI and match. Published functional codes still need review against the updated v2.4 flowchart.</p>
        <ul>{data.review_required.map(item => <li key={item}>{item}</li>)}</ul>
        <p>{data.limitations}</p>
        <p>{data.sources.map(source => <React.Fragment key={source.url}><a href={source.url} target="_blank" rel="noreferrer">{source.label} ↗</a><br/></React.Fragment>)}</p>
        <p>Tables: Fortuno et al. (2025), CC BY 4.0. GeneVISTA adds allele matching and a computational-criterion preview.</p>
      </details>
    </>}
  </section>;
}
