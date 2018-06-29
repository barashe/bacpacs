import glob
import warnings
import joblib
import os
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
        self.reduced_path_ = output_path

    def create_pfs(self, memory=800, n_jobs=1, cdhit_path=None, reduced_path=None, output_path=None):
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
            output_path = os.path.join(self._output_dir, 'protein_families')
        if reduced_path is None:
            if not hasattr(self, 'reduced_path_'):
                raise ValueError('No input fasta file')
            reduced_path = self.reduced_path_
        print 'Clustering genomes.'
        cdhit(reduced_path, output_path, memory, n_jobs, cdhit_path)
        print 'Clustering finished successfully. Protein families dumped in {}'.format(output_path)
        self.pf_path_ = output_path

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
            cdhit_2d(genome, pf_path, output_clusters_dir, memory, n_jobs, cdhit_path)
        print 'Genome cluster files are stored in {}'.format(output_clusters_dir)
        if training:
            self.train_clusters_dir_ = output_clusters_dir
        else:
            self.pred_clusters_dir_ = output_clusters_dir

    def extract_features(self, feats_type='pred', clusters_dir=None):
        """Get features matrix X (pandas.Dataframe) for training/prediction. If 'labels' is not None,
        a series of pathogenicity labels,y (pandas.Series) is returned as well.

        Parameters
        ----------
        feats_type : {'train', 'pred'}, optional
            Indication whether genomes are used for training, or for prediction.
        clusters_dir : basestring, optional
            Path to training/prediction (defined in 'feats_type'). If None, the path used by
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
        X = orgs_to_vecs(self.feat_list_, clusters_dir)
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
    if os.path.isdir(output_dir):
        warnings.warn('Directory {} already exists'.format(output_dir))
    else:
        os.mkdir(output_dir)
    bp._output_dir = output_dir
    return bp

