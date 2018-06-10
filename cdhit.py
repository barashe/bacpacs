import subprocess

from os import unlink

def run_cdhit(input_path, output_path, m, t, cdhit_path):
    """The actual CD-HIT wrapper. Runs CD-HIT accurate mode.
    :param m: Memory allocation for CD-HIT
    :param t: Number of CPUs to use for CD-HIT
    """
    cmd_line = '{}cd-hit -i {} -o {} -c 0.4 -n 2 -M {} -T {}'.format(cdhit_path + '/', input_path, output_path, m, t)
    subprocess.check_call(cmd_line.split())
    unlink(output_path + '.clstr')