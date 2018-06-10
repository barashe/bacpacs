import glob

from os.path import join
from Bio import SeqIO

from cdhit import run_cdhit

class Bacpacs(object):
    def __init__(self, output_dir):
        self._output_dir = output_dir

    def merge_genome_files(self, genomes_dir, output_path=None):
        """Merges all the raw training faa files.
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
        if output_path is None:
            output_path = join(self._output_dir, 'protein_families')
        if reduced_path is None:
            if not hasattr(self, 'reduced_path_'):
                raise ValueError('No input fasta file')
            reduced_path = self.reduced_path_
        print 'Clustering genomes.'
        run_cdhit(reduced_path, output_path, memory, n_jobs, cdhit_path)
        print 'Clustering finished successfully. Protein families dump in {}'.format(output_path)
        self.pf_path_ = output_path
