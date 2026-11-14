// Confusion-matrix rows are actual labels; columns are predicted labels.
export function summarizeConfusion(matrix, classes) {
 if (!Array.isArray(matrix)||matrix.length!==classes.length||matrix.some(row=>row.length!==classes.length||row.some(v=>!Number.isFinite(v)||v<0))) throw new Error('Invalid confusion matrix');
 const total=matrix.flat().reduce((a,b)=>a+b,0);
 return classes.map((label,i)=>{
  const support=matrix[i].reduce((a,b)=>a+b,0);
  const predicted=matrix.reduce((sum,row)=>sum+row[i],0);
  const correct=matrix[i][i];
  const precision=predicted?correct/predicted:null;
  const recall=support?correct/support:null;
  return {label,support,predicted,correct,precision,recall,
   f1:support+predicted?2*correct/(support+predicted):null,
   actualShare:total?support/total:null,predictedShare:total?predicted/total:null};
 });
}
