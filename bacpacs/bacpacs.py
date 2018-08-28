import sys
import glob
import warnings
import argparse
import joblib
import os
import json
from Bio import SeqIO

import util
import cli


class Bacpacs(object):
    def __init__(self, output_dir):
        """Returns a bacpacs object which processes features from raw faa files.

        Parameters
        ----------
        output_dir : basestring
            Output directory in which bacpacs will cache files, and store resulting features.
        """
        self._output_dir = output_dir
        if os.path.isdir(output_dir):
            warnings.warn('Directory {} already exists'.format(output_dir))
        else:
            os.mkdir(output_dir)

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
            output_path = os.path.join(self._output_dir, 'merged.faa')
        raw_seqs_paths = glob.glob(os.path.join(genomes_dir, '*.faa'))
        unique = set()
        with open(output_path, 'wb') as out_file:
            for genome_path in raw_seqs_paths:
                for rec in SeqIO.parse(genome_path, 'fasta'):
                    if rec.id not in unique:
                        SeqIO.write(rec, out_file, 'fasta')
                        unique.add(rec.id)
        print 'Saving merged proteins as {}'.format(output_path)
        self.merged_path_ = os.path.abspath(output_path)

    def reduce(self, long_percent=10, merged_path=None, output_path=None):
        """Selects the longest 10% proteins from the merged fasta file

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
            output_path = os.path.join(self._output_dir, 'reduced.faa')
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
        self.reduced_path_ = os.path.abspath(output_path)

    def create_pfs(self, memory=800, n_jobs=1, cdhit_path=None, reduced_path=None, output_path=None):
        """Runs CD-HIT to cluster the merged and reduced fasta file to protein families.

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
            output_path = os.path.join(self._output_dir, 'protein_families')
        if reduced_path is None:
            if not hasattr(self, 'reduced_path_'):
                raise ValueError('No input fasta file')
            reduced_path = self.reduced_path_
        print 'Clustering genomes.'
        util.cdhit(reduced_path, output_path, memory, n_jobs, cdhit_path)
        print 'Clustering finished successfully. Protein families dumped in {}'.format(output_path)
        self.pf_path_ = os.path.abspath(output_path)

    def genomes_vs_pfs(self, genomes_dir, feats_type='pred', memory=800, n_jobs=1, cdhit_path=None,
                       pf_path=None, output_clusters_dir=None):
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

        """
        if feats_type == 'pred':
            training = False
        elif feats_type == 'train':
            training = True
        else:
            raise ValueError("'{}' is not a valid feats_type. Use 'pred' or 'train'")
        if output_clusters_dir is None:
            dir_name = 'train_clusters' if training else 'pred_clusters'
            output_clusters_dir = os.path.join(self._output_dir, dir_name)
        if not os.path.isdir(output_clusters_dir):
            os.mkdir(output_clusters_dir)
        if pf_path is None:
            if not hasattr(self, 'pf_path_'):
                raise ValueError("Please specify a valid 'pf_path' or run self.cluster()")
            pf_path = self.pf_path_
        if not hasattr(self, 'feat_list_'):
            self.order_features(pf_path)
        print 'Running genomes against protein families representatives'
        for genome in glob.glob(os.path.join(genomes_dir, '*.faa')):
            print 'Processing {}'.format(genome)
            util.cdhit_2d(genome, pf_path, output_clusters_dir, memory, n_jobs, cdhit_path)
        print 'Genome cluster files are stored in {}'.format(output_clusters_dir)
        if training:
            self.train_clusters_dir_ = os.path.abspath(output_clusters_dir)
        else:
            self.pred_clusters_dir_ = os.path.abspath(output_clusters_dir)

    def extract_features(self, feats_type='pred', clusters_dir=None):
        """Get features matrix X (pandas.Dataframe) for training/prediction. If 'labels' is not None,
        a series of pathogenicity labels,y (pandas.Series) is returned as well.

        Parameters
        ----------
        feats_type : {'train', 'pred'}, optional
            Indication whether genomes are used for training, or for prediction.
        clusters_dir : basestring, optional
            Path to training/prediction (defined in 'feats_type') clusters. If None, the path used by
            self.train_clusters_dir_/self.pred_clusters_dir_ is used (created in self.genomes_vs_pfs()).

        Returns
        -------
        X or (X, y) : pandas.Dataframe or  (pd.Dataframe, pd.Series)
            If labels is None, only features (X) is returned. Otherwise, both features and labels are return in a tuple,
            indexed by genome ids.
        """

        if feats_type == 'train':
            training = True
        elif feats_type == 'pred':
            training = False
        else:
            raise ValueError("'feats_type should be 'train' or 'pred'")
        if clusters_dir is None:
            if training:
                if hasattr(self, 'train_clusters_dir_'):
                    clusters_dir = self.train_clusters_dir_
                else:
                    raise ValueError('No clusters specified. Specify clusters_dir or run self.extract_features()')
            else:
                if hasattr(self, 'pred_clusters_dir_'):
                    clusters_dir = self.pred_clusters_dir_
                else:
                    raise ValueError('No clusters specified. Specify clusters_dir or run self.extract_features()')
        if not hasattr(self, 'feat_list_'):
            raise ValueError("Features are not correctly defined. Run self.extract_features() first or "
                             "self.order_features()")
        X = util.orgs_to_vecs(self.feat_list_, clusters_dir)
        return X

    def order_features(self, pf_path=None):
        """This method uses a protein families file (pf_path) to make sure features are created in consistent order.
        It should only be used when running self.extract_features() without first running self.genomes_vs_pfs().

        Parameters
        ----------
        pf_path : basestring, optional
            Path to protein families file. If None, the path used by self.pf_path_ is used (created in self.cluster()).

        """
        if pf_path is None:
            if hasattr(self, 'pf_path_'):
                pf_path = self.pf_path_
            else:
                raise ValueError('Protein families path unknown. Please specify pf_path')
        self.feat_list_ = [rec.id for rec in SeqIO.parse(pf_path, 'fasta')]

    def to_pickle(self, path):
        """Dumps the Bacpacs object to a pickle file

        Parameters
        ----------
        path : basestring
            File path where the pickled object will be stored.

        """
        joblib.dump(self, path)

    def to_json(self, path):
        """Dumps the Bacpacs object to a JSON file

        Parameters
        ----------
        path : basestring
            File path where the JSON file will be stored.

        """
        data = dict()
        for attr_name in dir(self):
            if attr_name.endswith('_') and '__' not in attr_name:
                attr = getattr(self, attr_name)
                if isinstance(attr, (basestring, list)):
                    data[attr_name] = attr
                else:
                    data[attr_name] = attr.tolist()
        json.dump(data, open(path, 'wb'))


modes_string = 'bacpacs operating mode. init: Initiates a bacpacs working directory. Will create a "bp.json" file, ' \
               'which stores previous operations history. merge: Merges the raw training faa files. Reduce: Selects' \
               ' the longest 10 precent proteins from the merged fasta file. create_pfs: Runs CD-HIT to cluster the' \
               ' merged and reduced fasta file to protein families. genomes_vs_pf: Creates feature vectors for ' \
               'training/predicting genomes. Runs CD-HIT-2D for every genome, against the previously created' \
               ' protein families. extract_feats: Get features matrix X (pandas.DataFrame) for training/prediction. ' \
               'train: Trains a sklearn.svm.LinearSVC model on the extracted feats. predict: Using either the trained' \
               ' classifier trained in "train", or a classifier from a JSON file ' \
               '(created by bacpacs.util.clf_to_json) a prediction is made and stored in a csv file.'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-m', '--mode', choices=cli.modes.keys(), required=True, help=modes_string)
    s = 'Working directory in which bacpacs will cache files, and store resulting features.'
    required.add_argument('-w', '--working_directory', required=True, help=s)
    optional.add_argument('-i', '--input', help='Input file path')
    s = 'Working directory in which bacpacs will cache files, and store resulting features.'
    optional.add_argument('--genome_input_dir', help=s)
    s = 'Path to protein families file, created in "create_pfs". Applies to "pf_vs_genomes". If not specified, the' \
        ' last created pf_file is used.'
    optional.add_argument('--pf_path', help=s)
    s = 'Output file path. Applies to all modes except "init". If not specified, default paths in the working ' \
        'directory are used and printed to the screen.'
    optional.add_argument('-o', '--output', help=s)
    s = 'Indication whether genomes are used for training, or for prediction. Applies to "genomes_vs_pfs" and ' \
        '"extract_feats"'
    optional.add_argument('-t', '--feats_type', choices=['pred', 'train'], default='pred', help=s)
    optional.add_argument('-r', '--long_ratio', help='Ratio of long proteins to use in "reduce".', type=int, default=10)
    s = 'Path to training/prediction (defined in --feats_type) clusters. Applies to "extract_feats". If not ' \
        'specified, the directory used by "genomes_vs_pfs" to store clusters is used.'
    optional.add_argument('-c', '--clusters_dir', help=s)
    s = 'Path to training/prediction csv file. Applies to "train" and "predict". If not specified, the path used to' \
        ' store features in "extract_feats" is used.'
    optional.add_argument('-f', '--feats_path', help=s)
    optional.add_argument('-l', '--labels_path', help='Path to labels csv file. Applies to "train".')
    s = 'Path to scikit-learn classifier, stored in JSON format, using bacpacs.util.clf_to_json. If not supplied, a ' \
        'new sklearn.svm.LinearSVC is used.'
    optional.add_argument('--clf', help=s)
    s = 'Path to CD-HIT. Only required if CD-HIT not in environmental path. Applies to "create_pfs" and' \
        ' "genomes_vf_pfs".'
    optional.add_argument('--cdhit', help=s)
    s = 'Number of threads for CD-HIT-2D. 0 to use all CPUs. Applies to "create_pfs" and "genomes_vf_pfs".'
    optional.add_argument('--n_jobs', help=s, type=int, default=1)
    s = 'Memory limit (in MB) for CD-HIT, 0 for unlimited. Applies to "create_pfs" and "genomes_vf_pfs".'
    optional.add_argument('--memory', help=s, type=int, default=800)
    optional.add_argument('--pre_trained', action='store_true',
                          help='Using the pre-trained bacpacs model. Applies to init.')
    args = parser.parse_args()
    sys.exit(cli.modes[args.mode](args))
