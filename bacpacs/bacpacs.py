import glob
import warnings
import joblib
import pandas as pd
from os import mkdir
from os.path import join, isdir
from Bio import SeqIO
from util import cdhit, cdhit_2d, orgs_to_vecs


class Bacpacs(object):
    def __init__(self, output_dir):
        """Returns a bacpacs object which processes features from raw faa files.

        Parameters
        ----------
        output_dir : basestring
            Output directory in which bacpacs will cache files, and store resulting features.
        """
        self._output_dir = output_dir
        if isdir(output_dir):
            warnings.warn('Directory {} already exists'.format(output_dir))
        else:
            mkdir(output_dir)

    def merge_genome_files(self, genomes_dir, output_path=None):
        """Merges the raw training faa files.

        Parameters
        ----------
        genomes_dir : basestring
            Path to a directory containing genome .faa files.
        output_path : basestring, optional
            Output path. If None, saves 'reduced.faa' in self.output_dir.

        """
        if output_path is None:
            output_path = join(self._output_dir, 'merged.faa')
        raw_seqs_paths = glob.glob(join(genomes_dir, '*.faa'))
        unique = set()
        with open(output_path, 'wb') as out_file:
            for genome_path in raw_seqs_paths:
                for rec in SeqIO.parse(genome_path, 'fasta'):
                    if rec.id not in unique:
                        SeqIO.write(rec, out_file, 'fasta')
                        unique.add(rec.id)
        print 'Saving merged proteins as {}'.format(output_path)
        self.merged_path_ = output_path

    def reduce(self, long_percent=10, merged_path=None, output_path=None):
        """Samples the longest 10% proteins from the merged fasta file

        Parameters
        ----------
        long_percent : float, optional
            Determines the percentage of long proteins to be used for creating protein families
        merged_path : basestring, optional
            Path to merged faa file. If None, the path used by self.merge_genomes_files is used.
        output_path : basestring, optional
            Output path. If None, saves 'reduced.faa' in self.output_dir.

        """
        if output_path is None:
            output_path = join(self._output_dir, 'reduced.faa')
        if merged_path is None:
            if not hasattr(self, 'merged_path_'):
                raise ValueError('No merged fasta file')
            merged_path = self.merged_path_
        lens_and_ids = sorted([(len(rec), rec.id) for rec in SeqIO.parse(merged_path, 'fasta')], reverse=True)
        ids = [id for (length, id) in lens_and_ids]
        del lens_and_ids
        ids = ids[: len(ids) // long_percent]
        rec_index = SeqIO.index(merged_path, 'fasta')
        with open(output_path, 'wb') as out_file:
            for id in ids:
                SeqIO.write(rec_index[id], out_file, 'fasta')
        print 'Saving reduced proteins as {}'.format(output_path)
        self.reduced_path_ = output_path

    def cluster(self, memory=800, n_jobs=1, cdhit_path=None, reduced_path=None, output_path=None):
        """

        Parameters
        ----------
        memory : int, optional
            Memory limit (in MB) for CD-HIT, 0 for unlimited.
        n_jobs : int, optional
            Number of threads for CD-HIT. 0 to use all CPUs.
        cdhit_path : basestring, optional
            Path to CD-HIT
        reduced_path : basestring, optional
            Path to reduced faa file. If None, the path used by self.reduced_path_ is used (created in self.reduce()).
        output_path : basestring, optional
            Output path. If None, 'protein_families' is saved in self.output_dir).

        """
        if output_path is None:
            output_path = join(self._output_dir, 'protein_families')
        if reduced_path is None:
            if not hasattr(self, 'reduced_path_'):
                raise ValueError('No input fasta file')
            reduced_path = self.reduced_path_
        print 'Clustering genomes.'
        cdhit(reduced_path, output_path, memory, n_jobs, cdhit_path)
        print 'Clustering finished successfully. Protein families dumped in {}'.format(output_path)
        self.pf_path_ = output_path

    def extract_features(self, genomes_dir, feats_type='pred', memory=800, n_jobs=1, cdhit_path=None,
                         pf_path=None, output_clusters_dir=None, output_features_path=None):
        """Creates feature vectors for training/predicting genomes. Runs CD-HIT-2D for every genome, against the
        previously created protein families.

        Parameters
        ----------
        genomes_dir : basestring
            Path to a directory containing genome .faa files.
        feats_type : {'train', 'pred'}, optional
            Indication whether genomes are used for training, or for prediction.
        memory : int, optional
            Memory limit (in MB) for CD-HIT-2D, 0 for unlimited.
        n_jobs : int, optional
            Number of threads for CD-HIT-2D. 0 to use all CPUs.
        cdhit_path : basestring, optional
            Path to CD-HIT
        pf_path : basestring, optional
            Path to protein families file. If None, the path used by self.pf_path_ is used (created in self.cluster()).
        output_clusters_dir : basestring, optional
            Output dir for resulting genome clusters. If None, creates 'train_clusters'/'pred_clusters'
            (depending on 'feats_type') in self.output_dir.
        output_features_path: basestring, optional
            Output path for the resulting features. If None, 'pred_feats.pkl'/'train_feats.pkl' (depending on
            'feats_type') is saved in self.output_dir.

        """
        if feats_type == 'pred':
            training = False
        elif feats_type == 'train':
            training = True
        else:
            raise ValueError("'{}' is not a valid feats_type. Use 'pred' or 'train'")
        if output_clusters_dir is None:
            dir_name = 'train_clusters' if training else 'pred_clusters'
            output_clusters_dir = join(self._output_dir, dir_name)
        try:
            mkdir(output_clusters_dir)
        except Exception:
            raise IOError('Output dir {} already exists.'.format(output_clusters_dir))
        if output_features_path is None:
            features_file = 'train_feats.pkl' if training else 'pred_feats.pkl'
            output_features_path = join(self._output_dir, features_file)
        if pf_path is None:
            if not hasattr(self, 'pf_path_'):
                raise ValueError("Please specify a valid 'pf_path' or run self.cluster()")
            pf_path = self.pf_path_
        if not hasattr(self, 'feat_list_'):
            feat_list = [rec.id for rec in SeqIO.parse(pf_path, 'fasta')]
            self.feat_list_ = feat_list
        feat_list = self.feat_list_
        print 'Running genomes against protein families representatives'
        for genome in glob.glob(join(genomes_dir, '*.faa')):
            cdhit_2d(genome, pf_path, output_clusters_dir, memory, n_jobs, cdhit_path)
        print 'Organizing features'
        orgs_to_vecs(feat_list, output_clusters_dir, output_features_path)
        print 'Done. Features dumped to {}'.format(output_features_path)
        if training:
            self.train_feats_path_ = output_features_path
        else:
            self.pred_feats_path_ = output_features_path

    def get_features(self, feats_type='pred', feats_path=None, feat_order=None):
        """Get features matrix X (pandas.Dataframe) for training/prediction. If 'labels' is not None,
        a series of pathogenicity labels,y (pandas.Series) is returned as well.

        Parameters
        ----------
        feats_type : {'train', 'pred'}, optional
            Indication whether genomes are used for training, or for prediction.
        feats_path : basestring, optional
            Path to training/prediction (defined in 'feats_type') features. If None, the path used by
            self.train_feats/self.pred_feats is used (created in self.extract_features()).
        feat_order : list, optional
            Order of features for the resulting feature matrix. If not included, the feature list created in the
            training feature extraction process is used. If this list does not exist (e.g. when loading a model),
            the train set feature order is recorded, and used for the prediction set. If this doesn't exist, order is
            not validated and not guaranteed to be consistent.

        Returns
        -------
        X or (X, y) : pandas.Dataframe or  (pd.Dataframe, pd.Series)
            If labels is None, only features (X) is returned. Otherwise, both features and labels are return in a tuple,
            indexed by genome ids.
        """
        if feats_type == 'train':
            train = True
        elif feats_type == 'pred':
            train = False
        else:
            raise ValueError("'feats_type should be 'train' or 'pred'")
        if feats_path is None:
            if train:
                feats_path = self.train_feats_path_
            else:
                feats_path = self.pred_feats_path_
        if feat_order is None:
            if hasattr(self, 'feat_list_'):
                feat_order = self.feat_list_
            else:
                warnings.warn("Features might not have consistent order. Run self.order_feats() and "
                              "self.extract_features().")
        X = pd.read_pickle(feats_path)
        if feat_order is not None:
            X = X.loc[:, feat_order]
        return X

    def to_pickle(self, path):
        """Dumps the Bacpacs object to a pickle file

        Parameters
        ----------
        path : basestring
            File path where the pickled object will be stored.
        Returns
        -------

        """
        joblib.dump(self, path)


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
    if isdir(output_dir):
        warnings.warn('Directory {} already exists'.format(output_dir))
    else:
        mkdir(output_dir)
    bp._output_dir = output_dir
    return bp


def get_labels(csv_path, X=None):
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
        y = pd.read_csv(csv_path, index_col=0, dtype=np.bool)
    else:
        raise ValueError("'csv_path' is expected to be a string.")
    if X is not None:
        if not isinstance(X, pd.DataFrame):
            raise ValueError("'X' is expected to be a pandas.Dataframe")
        y = y.loc[X.index]
    if y.isnull().any():
        warnings.warn('Labels include some Null values')
    return y


def get_trained_model(output_dir):
    """Returns the Bacpacs and sklearn.svm.LinearSVC used to train the official bacpacs model.

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
    bp = read_pickle(r'..\trained\full_bacpacs.pkl', output_dir)
    svc = joblib.load(r'..\trained\linearsvc_full.pkl')
    return bp, svc
