"""Read-only progress and validation metrics for the current research experiments."""
import json
import math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
RUNS=(('CatBoost baseline','catboost_baseline_v6_pilot'),
      ('CatBoost + protein language models','catboost_protein_v6_pilot'),
      ('Full-data protein-model candidate','catboost_protein_v6_full'))


def progress(root=ROOT):
    experiments=[]
    for label,folder in RUNS:
        directory=Path(root)/'reports'/folder
        if not (directory/'started.json').is_file():continue
        entry={'name':label,'status':'training','evaluation_split':'development validation'}
        try:
            report=json.loads((directory/'report.json').read_text(encoding='utf-8'))
            entry.update(status=report['status'],train_rows=report.get('train_rows'),validation_rows=report.get('validation_rows'))
            selected=report.get('selected')
            if selected:
                metrics={k:selected[k] for k in ('accuracy','balanced_accuracy','macro_f1')}
                if not all(isinstance(v,(int,float)) and math.isfinite(v) and 0<=v<=1 for v in metrics.values()):raise ValueError('Invalid metrics')
                entry['metrics']=metrics
        except (FileNotFoundError,json.JSONDecodeError):pass
        except (KeyError,ValueError,TypeError):entry={'name':label,'status':'invalid_report','evaluation_split':'development validation'}
        experiments.append(entry)
    return {'experiments':experiments,'production_model_changed':False,
            'note':'Development validation is used for model selection. These results do not replace the active model test evaluation. Protein language-model scores are cached ESM1b/MutFormer outputs, not chat-generated labels.'}
