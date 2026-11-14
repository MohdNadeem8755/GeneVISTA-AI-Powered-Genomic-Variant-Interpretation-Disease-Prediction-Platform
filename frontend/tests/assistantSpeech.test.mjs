import test from 'node:test';
import assert from 'node:assert/strict';
import {spokenText,matchingVoices,speechParts} from '../src/assistantSpeech.mjs';
test('Hindi audio pronounces brand and isolated DNA letters without changing numbers',()=>{
  assert.equal(spokenText('GeneVISTA: A → T, 91.36%. DNA','hi-IN'),'जीन विस्टा: ए  से  टी, 91.36%. डी एन ए');
  assert.equal(spokenText('ATCG BRCA1','hi-IN'),'ATCG BRCA1');
});
test('English article is preserved and brand is not uppercase acronym',()=>{
  assert.equal(spokenText('A variant in GeneVISTA. A to T.','en-IN'),'A variant in Gene Vista. base A, to base T.');
});
test('Prefer Indian voice; do not replace Hindi with English or remote voice',()=>{
  const voices=[{lang:'en-US',localService:true},{lang:'hi-IN',localService:false},{lang:'en-IN',localService:true}];
  assert.equal(matchingVoices(voices,'en-IN')[0].lang,'en-IN');
  assert.deepEqual(matchingVoices(voices,'hi-IN'),[]);
});
test('Sentence splitting keeps decimal scores intact',()=>{
  const parts=speechParts('Accuracy is 91.36%. Next sentence.','en-IN');
  assert.ok(parts.some(p=>p.includes('91.36%')));
  assert.equal(parts.length,2);
});
