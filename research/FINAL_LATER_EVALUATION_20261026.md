# Frozen later-snapshot evaluation — insufficient support

The Step 4 candidate and evaluation protocol were frozen before downloading or
reading the later labels. Candidate SHA-256:
`ccc87eb46a54792254d096cb43fbecb3f17c61702a5dc21f019fdb3a7ee3f3c7`.

Source: official ClinVar `variant_summary.txt.gz`, last modified October 26, 2026.
The downloaded file passed the official MD5 check:
`c79b46f3c02ef544608e0d9a83cac46e`.

The pre-registered cohort required exact GRCh38 germline/unknown SNVs with one of
the five exact labels, no ambiguous label/identity, a VariationID and allele absent
from the October 5 snapshot, an allele absent from the original September cohort,
an exact dbNSFP match, and at least 30 examples in every class.

Before annotation or prediction, the frozen novel cohort contained:

| Class | Eligible records |
| --- | ---: |
| Benign | 2 |
| Likely benign | 11 |
| Likely pathogenic | 465 |
| Pathogenic | 1,515 |
| Uncertain significance | 3,258 |

Benign and Likely benign failed the pre-registered minimum. Therefore:

- no predictions were generated;
- no accuracy, balanced accuracy, macro-F1, thresholds or confusion matrix were computed;
- the cohort definition and support threshold were not changed after labels were read;
- the final test was not used for model selection or production replacement;
- the active model remains unchanged.

Computing five-class balanced accuracy from two Benign records would be unstable
and could falsely suggest success or failure. The valid next action is to wait for
a later pre-registered snapshot with at least 30 genuinely novel records in every
class. Repeatedly scoring this October 26 cohort would contaminate it as a future
final test.

Artifacts:

- `reports/final_later_20261026/frozen_protocol.json`
- `reports/final_later_20261026/prepare_report.json`
- `reports/final_later_20261026/insufficient_support.json`
- `reports/final_later_20261026/cohort.sqlite`

The source snapshot is retained with its checksum for auditability. No annotation
scan was run because support failure was established before model inputs were
needed.
