import glob
from os import mkdir
from os.path import join
from Bio import SeqIO

from util import cdhit, cdhit_2d, orgs_to_vecs


class Bacpacs(object):
    def __init__(self, output_dir):
        self._output_dir = output_dir

    def merge_genome_files(self, genomes_dir, output_path=None):
        """Merges the raw training faa files.

        Parameters
        ----------
        genomes_dir
        output_path

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
            Output path. If None, saves 'reduced.faa' in the initial output dir (created in self.merge()).

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

    def cluster(self, memory=800, n_jobs=1, cdhit_path='.', reduced_path=None, output_path=None):
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

    def extract_features(self, genomes_dir, feats_type='pred', memory=800, n_jobs=1, cdhit_path='.',
                         pf_path=None, output_clusters_dir=None, output_features_path=None):
        """Creates feature vectors for training/predicting genomes. Runs CD-HIT-2D for every genome, against the
        previously created protein families.

        Parameters
        ----------
        genomes_dir : basestring
            Path to a directory containing genome .faa files.
        feats_type : ['train', 'pred']
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
        if training:
            feat_list = [rec.id for rec in SeqIO.parse(pf_path, 'fasta')]
            self.feat_list_ = feat_list
        else:
            if not hasattr(self, 'features_'):
                raise ValueError("Training features were not yet created. Run with 'feats_type='train''")
        feat_list = self.feat_list_
        print 'Running genomes against protein families representatives'
        for genome in glob.glob(join(genomes_dir, '*.faa')):
            cdhit_2d(genome, pf_path, output_clusters_dir, memory, n_jobs, cdhit_path)
        print 'Organizing features'
        orgs_to_vecs(feat_list, output_clusters_dir, output_features_path)
        print 'Done. Features dumped to {}'.format(output_features_path)

    def train(self):
        pass

    def predict(self):
        pass
