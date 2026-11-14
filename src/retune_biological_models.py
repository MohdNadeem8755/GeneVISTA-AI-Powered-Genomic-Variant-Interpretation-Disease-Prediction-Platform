"""Fresh development search with locus-separated stopping and checkpointed fits."""
import argparse
import itertools
import json
from pathlib import Path
import joblib
import numpy as np
from catboost import CatBoostClassifier
from src.annotation_pipeline import bucket
from src.biological_features import load_development
from src.robust_variant_training import metrics, save_json, file_hash
from src.tune_biological_models import predict_bundle

ROOT = Path(__file__).resolve().parents[1]


def stopping_mask(groups):
    return np.asarray([bucket('step4-redo-stop|' + str(g), 10) == 0 for g in groups])


def passes_gate(candidate, baseline):
    if any(candidate[k] < baseline[k] for k in ('accuracy','balanced_accuracy','macro_f1')):
        return False
    for label in baseline['class_order']:
        for key in ('precision','recall'):
            if candidate['classification_report'][label][key] < baseline['classification_report'][label][key]:
                return False
    return all(candidate['classification_report'][label]['recall'] > baseline['classification_report'][label]['recall'] for label in ('Likely benign','Likely pathogenic'))


def run(output, resume=False, corrected=False, minority_all=False):
    if minority_all and not corrected:
        raise ValueError('Expanded run requires corrected features')
    out = Path(output)
    out.mkdir(parents=True, exist_ok=resume)
    configs = [dict(profile='all', depth=6, power=p) for p in (0., .25, .5, .75, 1.)]
    configs += [dict(profile='all', depth=8, power=.5),
                dict(profile='population_sequence', depth=6, power=.5),
                dict(profile='biological_context', depth=6, power=.5)]
    feature_profile='variant_identity_v2' if corrected else 'legacy'
    objective=['macro_f1','balanced_accuracy'] if corrected else ['balanced_accuracy','macro_f1']
    protocol = dict(configs=configs, modulus=5, max_iterations=800, patience=100,
                    minority_all=minority_all,
                    seed=20261030, learning_rate=.05, l2_leaf_reg=20, threads=16,
                    stopping='10% of train loci, unweighted macro F1',
                    selection=objective,feature_profile=feature_profile,
                    preprocessing_sha256=file_hash(ROOT/'src/biological_features.py'),
                    revised_features_sha256=file_hash(ROOT/'src/variant_feature_revision.py') if corrected else None,
                    validation='Previously inspected development data; not independent evaluation',
                    pretrained_scores=False, final_test_rows_read=0,
                    feature_report_sha256=file_hash(ROOT/'data/processed/biological_features_20261029/report.json'),
                    script_sha256=file_hash(Path(__file__)))
    if resume:
        if json.loads((out/'protocol.json').read_text()) != protocol:
            raise ValueError('Protocol changed: use a fresh directory')
        if (out/'selection.json').exists():
            print('Run already complete', flush=True)
            return
    else:
        save_json(out/'protocol.json', protocol)
    print('Loading larger development sample', flush=True)
    cohort=ROOT/'data/processed/annotation_v3/full_cohort.sqlite'
    cache=ROOT/'data/processed/biological_features_20261029/features.sqlite'
    x,y,cats,groups=load_development(cohort,cache,'train',5,return_groups=True,feature_profile=feature_profile,minority_all=minority_all)
    xv,yv,_,vgroups=load_development(cohort,cache,'validation',5,return_groups=True,feature_profile=feature_profile)
    if set(groups) & set(vgroups):
        raise ValueError('Train/validation locus overlap')
    stop=stopping_mask(groups)
    classes=np.unique(y)
    for subset in (y[stop], y[~stop], yv):
        if set(subset) != set(classes):
            raise ValueError('Missing class in development partition')
    for frame in (x,xv):
        for c in cats:
            frame[c]=frame[c].fillna('__MISSING__').astype(str)
    prefixes=('clingen_','inheritance_','brca1_','tp53_','assay_')
    profiles={'all':list(x),
              'population_sequence':[c for c in x if not c.startswith(prefixes)],
              'biological_context':[c for c in x if c.startswith(prefixes)]}
    save_json(out/'data_audit.json',dict(train_rows=len(y), fit_rows=int((~stop).sum()),
        stopping_rows=int(stop.sum()), validation_rows=len(yv), train_validation_locus_overlap=0,
        fit_support=y[~stop].value_counts().to_dict(), stop_support=y[stop].value_counts().to_dict(),
        validation_support=yv.value_counts().to_dict(), feature_profiles=profiles))
    print(f'Fit {(~stop).sum()}, stop {stop.sum()}, validation {len(yv)}',flush=True)
    trained=[]; probs=[]; results=[]
    for i,cfg in enumerate(configs):
        path=out/f'candidate_{i}.joblib'
        cols=profiles[cfg['profile']]; cat=[c for c in cats if c in cols]
        if resume and path.exists():
            member=joblib.load(path)
        else:
            print(f'Training {i}: {cfg}',flush=True)
            counts=y[~stop].value_counts()
            fit_y=y[~stop]
            weights=(len(fit_y)/(len(classes)*fit_y.map(counts).to_numpy()))**cfg['power']
            model=CatBoostClassifier(iterations=800, depth=cfg['depth'], learning_rate=.05,
                l2_leaf_reg=20, thread_count=16, random_seed=20261030,
                loss_function='MultiClass', eval_metric='TotalF1:average=Macro;use_weights=false',
                allow_writing_files=False)
            model.fit(x.loc[~stop,cols],y[~stop],cat_features=cat,sample_weight=weights,
                eval_set=(x.loc[stop,cols],y[stop]),early_stopping_rounds=100,verbose=100)
            member=dict(classifier=model,columns=cols,categorical=cat)
            joblib.dump(member,path,compress=3)
        p=predict_bundle(dict(models=[member],classes=classes,weights=[1]),xv)
        report=dict(kind='single',members=[i],weights=[1],config=cfg,
                    iterations=member['classifier'].tree_count_,
                    **metrics(yv,classes[p.argmax(axis=1)],classes))
        trained.append(member); probs.append(p); results.append(report)
        save_json(out/'development.json',results)
        print({k:report[k] for k in ('accuracy','balanced_accuracy','macro_f1','iterations')},flush=True)
    # Fixed grid of pairwise votes, including unweighted and feature specialists.
    ensembles=[]
    for a,b in itertools.combinations(range(len(trained)),2):
        for weights in ((1,1),(2,1),(1,2)):
            p=np.average([probs[a],probs[b]],axis=0,weights=weights)
            ensembles.append(dict(kind='ensemble',members=[a,b],weights=list(weights),
                **metrics(yv,classes[p.argmax(axis=1)],classes)))
    save_json(out/'ensembles.json',ensembles)
    selected=max(results+ensembles,key=lambda r:tuple(r[key] for key in objective))
    if minority_all:
        baseline_dir=ROOT/'reports/biological_corrected_20261109'
        baseline=json.loads((baseline_dir/'selection.json').read_text())
        assert file_hash(baseline_dir/'selected.joblib')==baseline['model_sha256']
        baseline_bundle=joblib.load(baseline_dir/'selected.joblib')
        bp=predict_bundle(baseline_bundle,xv)
        reproduced=metrics(yv,baseline_bundle['classes'][bp.argmax(axis=1)],baseline_bundle['classes'])
        for key in ('accuracy','balanced_accuracy','macro_f1'):
            assert abs(reproduced[key]-baseline['selected'][key])<1e-10
        eligible=[r for r in results+ensembles if passes_gate(r,baseline['selected'])]
        save_json(out/'acceptance.json',dict(accepted=bool(eligible),eligible_count=len(eligible),
            baseline_sha256=baseline['model_sha256'],baseline=baseline['selected'],
            rule='No decrease in overall accuracy, balanced accuracy, macro F1 or any class precision/recall; both Likely recalls strictly improve.',
            production_model_changed=False,test_rows_read=0,
            limitation='Repeatedly inspected validation; this gate does not prove absence of overfitting or independent generalization.'))
        if not eligible:
            save_json(out/'rejected_best.json',selected)
            print('COMPLETED: no candidate passed all no-decrease requirements; old model preserved.',flush=True)
            return
        selected=max(eligible,key=lambda r:tuple(r[key] for key in objective))
    best_accuracy=max(results+ensembles,key=lambda r:r['accuracy'])
    bundle=dict(models=[trained[i] for i in selected['members']],weights=selected['weights'],
        classes=classes,metadata=dict(status='development_candidate',test_evaluated=False,selected=selected,feature_profile=feature_profile))
    joblib.dump(bundle,out/'selected.joblib',compress=3)
    reloaded=joblib.load(out/'selected.joblib')
    np.testing.assert_allclose(predict_bundle(bundle,xv.iloc[:100]),predict_bundle(reloaded,xv.iloc[:100]))
    save_json(out/'selection.json',dict(selected=selected,highest_accuracy=best_accuracy,
        production_model_changed=False,test_evaluated=False,
        model_sha256=file_hash(out/'selected.joblib'),
        caution='Reused development selection metrics; not final-test accuracy.'))
    print('COMPLETED', {k:selected[k] for k in ('accuracy','balanced_accuracy','macro_f1')},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--minority-all',action='store_true')
    parser.add_argument('--corrected',action='store_true',help='New fit only: retain variant identity, fix event parsing, and align stopping/selection on macro F1')
    args=parser.parse_args()
    run(args.output,args.resume,args.corrected,args.minority_all)
