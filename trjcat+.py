#! /usr/bin/env python

"""Functionize it and make it act like a gromacs tool"""

import os
import argparse
import subprocess

# import numpy as np
# from MDAnalysis import Universe

from utils import main_load, backup

def main(args):
    main_load()

    output = args.optf
    if output is None:
        outputfile = '{0:s}.output.xtc'.format(args.tprf)
    else:
        outputfile = output
    backup(outputfile)

    # do calculation
    trjcat_plus(args.xtcf, args.tprf, outputfile)

def trjcat_plus(xtcf, tprf, output):
    xtcf = sorted(xtcf)
    usable_xtcf = []
    for f in xtcf:
        print "%" * 79
        print "PROCESSING {0}".format(f)
        print "%" * 79
        returncode = subprocess.call(
            "gmxcheck -f {0}".format(f),
            shell=True)
        if returncode != 0:
            print "%" * 79
            print "CRASHED: {0}, trying to fix".format(f)
            print "%" * 79
            bkf = f + '.bk.xtc'
            os.rename(f, bkf)
            subprocess.call(
                "echo 'System' | trjconv -f {0} -s {1} -o {2}".format(bkf, tprf, f), 
                shell=True)
        usable_xtcf.append(f)
    subprocess.call(
        "trjcat -f {0} -o {1}".format(' '.join(xtcf), output),
        shell=True)

def parse_cmd():
    usage='used to calculate sequence spacing'
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('-f', '--xtcf', type=str, dest='xtcf', default=None, required=True, nargs="+",
                        help='Trajectory: xtc')
    parser.add_argument('-s', '--tprf', type=str, dest='tprf', default=None, required=True,
                        help='Structure: gro')
    parser.add_argument('-o', '--optf', type=str, dest='optf', default=None, 
                        help='Output file name (OPT.)')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_cmd()
    main(args)
