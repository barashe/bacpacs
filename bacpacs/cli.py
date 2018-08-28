import os
import pandas as pd
from sklearn.svm import LinearSVC
import bacpacs
import util


def run_cli(args):
    wd = args.working_directory
    bp = bacpacs.Bacpacs(wd)
    if args.mode == 'merge':
        bp.merge_genome_files()
    elif args.mode == 'reduce':
        pass
    elif args.mode == 'create_pfs':
        pass
    elif args.mode == 'genomes_to_pfs':
        pass
    elif args.mode == 'extract':
        pass


def init(args):
    wd = args.working_directory
    if args.pre_trained:
        bp, _ = util.load_trained_model(wd)
        bp.trained_clf_ = os.path.abspath(os.path.join(wd, 'trained_model/linearsvc_full.json'))
    else:
        bp = bacpacs.Bacpacs(wd)
    bp.to_json(os.path.join(wd, 'bp.json'))
    print 'bacpacs initialized in {}'.format(wd)


def merge(args):
    wd = args.working_directory
    bp = util.read_json(os.path.join(wd, 'bp.json'), wd)
    bp.merge_genome_files(args.genome_input_dir, args.output)
    bp.to_json(os.path.join(wd, 'bp.json'))


def reduce(args):
    wd = args.working_directory
    bp = util.read_json(os.path.join(wd, 'bp.json'), wd)
    bp.reduce(args.long_ratio, args.input, args.output)
    bp.to_json(os.path.join(wd, 'bp.json'))


def create_pfs(args):
    wd = args.working_directory
    bp = util.read_json(os.path.join(wd, 'bp.json'), wd)
    bp.create_pfs(args.memory, args.n_jobs, args.cdhit, args.input, args.output)
    bp.to_json(os.path.join(wd, 'bp.json'))


def genomes_vs_pfs(args):
    wd = args.working_directory
    bp = util.read_json(os.path.join(wd, 'bp.json'), wd)
    bp.genomes_vs_pfs(args.genome_input_dir, args.feats_type, args.memory, args.n_jobs, args.cdhit, args.pf_path,
                      args.output)
    bp.to_json(os.path.join(wd, 'bp.json'))


def extract_feats(args):
    wd = args.working_directory
    bp = util.read_json(os.path.join(wd, 'bp.json'), wd)
    X = bp.extract_features(args.feats_type, args.clusters_dir)
    if args.output is None:
        feats_path = os.path.join(wd, '{}_feats.csv'.format(args.feats_type))
    else:
        feats_path = args.output
    X.to_csv(feats_path)
    setattr(bp, '{}_feats_'.format(args.feats_type), feats_path)
    print '{} feats stored in {}'.format('Prediction' if args.feats_type == 'pred' else 'Training', feats_path)
    bp.to_json(os.path.join(wd, 'bp.json'))


def train(args):
    wd = args.working_directory
    bp = util.read_json(os.path.join(wd, 'bp.json'), wd)
    if args.labels_path is None:
        raise ValueError('Please provide a valid labels path in --labels_path')
    if args.feats_path is None:
        if not hasattr(bp, 'train_feats_'):
            raise ValueError('Please provide a valid features path in --feats_path')
        feats_path = bp.train_feats_
    else:
        feats_path = args.feats_path
    if args.clf is None:
        clf = LinearSVC('l1', dual=False, class_weight='balanced')
    else:
        clf = util.json_to_clf(args.clf)
    X_train = pd.read_csv(feats_path, converters={'genome_id': str}).set_index('genome_id')
    y_train = util.read_labels(args.labels_path, X_train)
    clf.fit(X_train, y_train)
    if args.output is None:
        clf_path = os.path.join(wd, 'trained_clf.json')
    else:
        clf_path = args.output
    util.clf_to_json(clf, clf_path)
    bp.trained_clf_ = os.path.abspath(clf_path)
    bp.to_json(os.path.join(wd, 'bp.json'))
    print 'Trained classifier is stored at {}'.format(clf_path)


def predict(args):
    wd = args.working_directory
    bp = util.read_json(os.path.join(wd, 'bp.json'), wd)
    if args.feats_path is None:
        if not hasattr(bp, 'pred_feats_'):
            raise ValueError('Please provide a valid features path in --feats_path')
        feats_path = bp.pred_feats_
    else:
        feats_path = args.feats_path
    if args.clf is None:
        if not hasattr(bp, 'trained_clf_'):
            raise ValueError('Please provide a valid trained classifier path in --clf or run "train"')
        clf_path = bp.trained_clf_
    else:
        clf_path = args.clf
    clf = util.json_to_clf(clf_path, LinearSVC)
    X_pred = pd.read_csv(feats_path, converters={'genome_id': str}).set_index('genome_id')
    y_pred = pd.Series(clf.predict(X_pred), X_pred.index)
    if args.output is None:
        pred_csv_path = os.path.join(wd, 'predictions.csv')
    else:
        pred_csv_path = args.output
    y_pred.to_csv(pred_csv_path)
    print 'Predictions stored at {}'.format(pred_csv_path)


modes = {'init': init,
         'merge': merge,
         'reduce': reduce,
         'create_pfs': create_pfs,
         'genomes_vs_pfs': genomes_vs_pfs,
         'extract_feats': extract_feats,
         'train': train,
         'predict': predict}

