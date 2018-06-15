import subprocess

from os import unlink
from os.path import basename, splitext, join

def run_cdhit(input_path, output_path, m, t, cdhit_path):
    """A wrapper for CD-HIT. Runs CD-HIT greedy mode.

    Parameters
    ----------
    input_path : basestring
        Path to input faa file.
    output_path : basestring
        Output path.
    m : int
        Memory limit (in MB) for the CD-HIT, 0 for unlimited.
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
        Memory limit (in MB) for the CD-HIT, 0 for unlimited.
    t : int
        Number of threads. 0 to use all CPUs.
    cdhit_path : basestring
        Path to CD-HIT.

    """
    org_name = splitext(basename(org_path))[0]
    output_path = join(output_dir, org_name)
    cmd_line = '{}cd-hit-2d -i {} -i2 {} -o {} -c 0.4 -n 2 -d 0 -M {} -T {} -g 1'
    cmd_line = cmd_line.format(cdhit_path + '/', clusters_path, org_path, output_path, m, t)
    subprocess.check_call(cmd_line.split())
    # We are only interested in output_path.clstr, which is automatically created by CD-HIT-2D
    unlink(output_path)
