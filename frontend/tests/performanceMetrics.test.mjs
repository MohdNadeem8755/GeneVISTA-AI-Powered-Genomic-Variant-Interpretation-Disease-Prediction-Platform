import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {summarizeConfusion} from '../src/performanceMetrics.mjs';
test('rows are actual; columns are predictions; minority errors stay visible',()=>{
 const r=summarizeConfusion([[80,20],[40,60]],['a','b']);
 assert.equal(r[0].support,100); assert.equal(r[0].predicted,120);
 assert.equal(r[0].recall,.8); assert.equal(r[0].precision,80/120);
});
test('absent actual/predicted classes are not called perfect',()=>{
 const r=summarizeConfusion([[0,0],[0,2]],['a','b']);
 assert.equal(r[0].precision,null); assert.equal(r[0].recall,null);
});
test('displayed metrics reproduce saved evaluation for every partition',()=>{
 const report=JSON.parse(readFileSync(new URL('../../reports/biological_retuning_fresh_20261031/partition_evaluation.json',import.meta.url)));
 for(const r of Object.values(report.scores)){
  const s=summarizeConfusion(r.confusion_matrix,r.class_order);
  assert.equal(s.reduce((n,c)=>n+c.support,0),r.n);
  assert.ok(Math.abs(s.reduce((n,c)=>n+c.correct,0)/r.n-r.accuracy)<1e-12);
  for(const c of s){
   const saved=r.classification_report[c.label];
   assert.ok(Math.abs(c.precision-saved.precision)<1e-12);
   assert.ok(Math.abs(c.recall-saved.recall)<1e-12);
  }
 }
});
