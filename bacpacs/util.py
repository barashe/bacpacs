import subprocess
import glob
import os
import joblib
import urllib
import warnings
import tarfile
import json
import sklearn
import bacpacs
import numpy as np
import pandas as pd


def cdhit(input_path, output_path, m, t, cdhit_path):
    """A wrapper for CD-HIT. Runs CD-HIT greedy mode.

    Parameters
    ----------
    input_path : basestring
        Path to input faa file.
    output_path : basestring
        Output path.
    m : int
        Memory limit (in MB) for CD-HIT, 0 for unlimited.
    t : int
        Number of threads. 0 to use all CPUs.
    cdhit_path : basestring
        Path to CD-HIT.

    """
    if cdhit_path is None:
        cdhit_path = 'cd-hit'
    else:
        cdhit_path = os.path.join(cdhit_path, 'cd-hit')
    cmd_line = '{} -i {} -o {} -c 0.4 -n 2 -M {} -T {}'.format(cdhit_path, input_path, output_path, m, t)
    print cmd_line
    subprocess.check_call(cmd_line.split())
    os.unlink(output_path + '.clstr')


def cdhit_2d(org_path, clusters_path, output_dir, m, t, cdhit_path):
    """A wrapper for CD-HIT-2D

    Parameters
    ----------
    org_path : basestring
        Path to org's faa file.
    clusters_path : basestring
        Path to protein clusters file.
    output_dir : basestring
        Path to output directory. Resulting file will have a .clstr extension.
    m : int
        Memory limit (in MB) for CD-HIT-2D, 0 for unlimited.
    t : int
        Number of threads. 0 to use all CPUs.
    cdhit_path : basestring
        Path to CD-HIT.

    """
    if cdhit_path is None:
        cdhit_path = 'cd-hit-2d'
    else:
        cdhit_path = os.path.join(cdhit_path, 'cd-hit-2d')
    org_name = get_file_name(org_path)
    output_path = os.path.join(output_dir, org_name)
    cmd_line = '{} -i {} -i2 {} -o {} -c 0.4 -n 2 -d 0 -M {} -T {} -g 1'
    cmd_line = cmd_line.format(cdhit_path, clusters_path, org_path, output_path, m, t)
    subprocess.check_call(cmd_line.split())
    # We are only interested in output_path.clstr, which is automatically created by CD-HIT-2D
    os.unlink(output_path)


def orgs_to_vecs(feat_list, clusters_dir):
    """Generated a feature pandas.Dataframe, with genome id as index, and protein families features as columns.

    Parameters
    ----------
    feat_list : list
        Feature names, for resulting feature dataframe columns.
    clusters_dir : basestring
        Directory in which the relevant CD-HIT-2D clusters are found.
    output_path : basestring
        Path to dump resulting features.

    """

    cluster_files = glob.glob(os.path.join(clusters_dir, '*.clstr'))
    genome_ids = pd.Series([get_file_name(cluster_file) for cluster_file in cluster_files], name='genome_id')
    ids_and_clusters = pd.Series(cluster_files, genome_ids)
    features = pd.DataFrame(0, index=genome_ids, columns=feat_list, dtype=np.bool)
    for genome_id, feat_vec in features.iterrows():
        genome_cluster_file = ids_and_clusters[genome_id]
        tmp_feat = ''
        with open(genome_cluster_file) as f:
            for line in f:
                if line[0] == '0':
                    tmp_feat = line.split()[2][1:-3]
                elif line[0] == '1':
                    feat_vec.loc[tmp_feat] = True
    return features


def get_file_name(path):
    return os.path.splitext(os.path.basename(path))[0]


def read_labels(csv_path, X=None):
    """Generates a pandas.Series of pathogenicity labels, given a csv file. the file should include two
    columns; genome id (first), and pathogenicity label (second).

    Parameters
    ----------
    csv_path : basestring
        Path to csv file. The file should include two columns; genome id (first), and pathogenicity label (second).
    X : pandas.Dataframe, optional
        The corresponding feature matrix, with genome ids as index. This is used for ordering purposes. If X is
        given, the labels returned are sorted according to it.

    Returns
    -------
    y : pandas.Series
        Pathogenicity labels, sorted according to X, if given.

    """
    if isinstance(csv_path, basestring):
        y = pd.read_csv(csv_path, header=None, dtype={0: np.object}).set_index(0)[1]
        y.index.name = 'genome_id'
    else:
        raise ValueError("'csv_path' is expected to be a string.")
    if X is not None:
        if not isinstance(X, pd.DataFrame):
            raise ValueError("'X' is expected to be a pandas.Dataframe")
        y = y.loc[X.index]
    if y.isnull().any():
        print 'Warning: Labels include null values'
    return y


def load_trained_model(output_dir):
    """Downloads and returns the Bacpacs and sklearn.svm.LinearSVC used to train
    the official bacpacs model.

    Parameters
    ----------
    output_dir : basestring
        Output directory in which bacpacs will cache files, and store resulting features.

    Returns
    -------
    bp : Bacpacs
        Bacpacs object to generate prediction features.
    svc : sklearn.svm.LinearSVC
        Trained estimator to predict new organisms.

    """

    github_path = 'https://github.com/barashe/bacpacs/raw/develop/trained/{}'
    # github_path = 'https://github.com/barashe/bacpacs/raw/master/trained/{}'
    file_names = ['full_bacpacs.json', 'linearsvc_full.json', 'protein_families']
    local_dir = os.path.join(output_dir, 'trained_model')
    if os.path.isdir(output_dir):
        warnings.warn('Directory {} already exists'.format(output_dir))
    else:
        os.mkdir(output_dir)
    if os.path.isdir(local_dir):
        warnings.warn('Directory {} already exists'.format(local_dir))
    else:
        os.mkdir(local_dir)
    print 'Retrieving files from github'
    for file_name in file_names:
        local_path = os.path.join(local_dir, file_name)
        if not os.path.isfile(local_path):
            print 'Downloading {}'.format(file_name)
            urllib.urlretrieve(github_path.format(file_name), local_path)
        else:
            print '{} exists. Skipping.'.format(file_name)
    bp = read_json(os.path.join(local_dir, 'full_bacpacs.json'), output_dir)
    bp.pf_path_ = os.path.abspath(os.path.join(local_dir, 'protein_families'))
    svc = json_to_clf(os.path.join(local_dir, 'linearsvc_full.json'))
    return bp, svc


def download_toy_data(target_directory):
    """Downloads toy data set from github, and stores it in target_directory/toy.

    Parameters
    ----------
    target_directory : basestring
        Directory to download the toy set directory to

    """
    local_tar_path = os.path.join(target_directory, 'toy.tar.gz')
    if os.path.isdir(os.path.join(target_directory, 'toy')):
        raise ValueError('{} already exists'.format(os.path.join(target_directory, 'toy')))
    urllib.urlretrieve('https://github.com/barashe/bacpacs/raw/master/toy.tar.gz', local_tar_path)
    with tarfile.open(local_tar_path) as f:
        f.extractall()
    os.remove(local_tar_path)
    print 'Toy data stored in {}'.format(os.path.join(target_directory, 'toy'))


def clf_to_json(clf, path):
    """Dumps a Scikit-Learn classifier to a JSON file

    Parameters
    ----------
    clf : sklearn classifier
        Classifier to dump. Classifier must have a self.get_params() method.
    path : basestring
        Path for resulting JSON file.

    """
    attrs = dict((attr, getattr(clf, attr).tolist()) for attr in dir(clf) if attr.endswith('_') and '__' not in attr)
    data = dict()
    data['params'] = clf.get_params()
    data['attr'] = attrs
    json.dump(data, open(path, 'wb'))


def json_to_clf(path, clf_class=None):
    """Loads a Scikit-Learn classifier from a JSON file.

    Parameters
    ----------
    path : basestring
        Path of JSON file
    clf_class :
        Classifier class to load, matching the stored classifier class.
    Returns
    -------
    clf : sklearn classifier
    """
    if clf_class is None:
        clf_class = sklearn.svm.LinearSVC
    data = json.load(open(path))
    clf = clf_class(**data['params'])
    for attr_name, attr in data['attr'].items():
        setattr(clf, attr_name, np.array(attr))
    return clf


def read_pickle(path, output_dir):
    """Load pickled Bacpacs object.

    Parameters
    ----------
    path : basestring
        File path where the pickled object will be loaded.
    output_dir : basestring
        Output directory in which bacpacs will cache files, and store resulting features.

    Returns
    -------
    bacpacs : a Bacpacs object

    """
    bp = joblib.load(path)
    if os.path.isdir(output_dir):
        warnings.warn('Directory {} already exists'.format(output_dir))
    else:
        os.mkdir(output_dir)
    bp._output_dir = output_dir
    return bp


def read_json(path, output_dir):
    """Load JSON Bacpacs object.

        Parameters
        ----------
        path : basestring
            File path where the JSON string will be loaded.
        output_dir : basestring
            Output directory in which bacpacs will cache files, and store resulting features.

        Returns
        -------
        bacpacs : a Bacpacs object

        """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bp = bacpacs.Bacpacs(output_dir)
    data = json.load(open(path))
    for attr_name, attr in data.items():
        if not isinstance(attr, basestring):
            attr = np.array(attr)
        setattr(bp, attr_name, attr)
    return bp
