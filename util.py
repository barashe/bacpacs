import subprocess
import glob
import numpy as np
import pandas as pd

from os import unlink
from os.path import basename, splitext, join


def get_file_name(path):
    return splitext(basename(path))[0]


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
    cmd_line = '{}cd-hit -i {} -o {} -c 0.4 -n 2 -M {} -T {}'.format(cdhit_path + '/', input_path, output_path, m, t)
    subprocess.check_call(cmd_line.split())
    unlink(output_path + '.clstr')


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
    org_name = get_file_name(org_path)
    output_path = join(output_dir, org_name)
    cmd_line = '{}cd-hit-2d -i {} -i2 {} -o {} -c 0.4 -n 2 -d 0 -M {} -T {} -g 1'
    cmd_line = cmd_line.format(cdhit_path + '/', clusters_path, org_path, output_path, m, t)
    subprocess.check_call(cmd_line.split())
    # We are only interested in output_path.clstr, which is automatically created by CD-HIT-2D
    unlink(output_path)


def orgs_to_vecs(feat_list, clusters_dir, output_path):
    # returns a train/test set in the following format
    # [[org_ids][org_feats][org_pathogenicities]]
    cluster_files = glob.glob(join(clusters_dir, '*.clstr'))
    genome_ids = [get_file_name(cluster_file) for cluster_file in cluster_files]
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
    features.to_pickle(output_path)
