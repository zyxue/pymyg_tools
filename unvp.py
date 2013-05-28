#! /usr/bin/env python

"""Functionize it and make it act like a gromacs tool"""

import os
import sys
import time
import numpy as np
from MDAnalysis import Universe
import argparse

from utils import main_load, write_header, write_footer, backup_old_output

def main(ARGS):
    if ARGS.debug:
        main_load()

    # check the validity of output file name, do backup
    output = ARGS.optf
    if output is None:
        outputfile = '{0:s}.output.xvg'.format(ARGS.grof)
    else:
        outputfile = output
    backup_old_output(outputfile)

    # Do some logging at the beginning
    outputf = open(outputfile, 'w')
    beginning_time = write_header(outputf)

    # do calculation
    result = count_interactions(
        ARGS.grof, 
        ARGS.xtcf, 
        ARGS.btime,
        ARGS.cutoff * 10,                          # convert to angstrom from nm
        ARGS.debug)

    # write results to the outputfile
    outputf.write('# {0:>10s}{1:>8s}\n'.format('time', 'num'))
    for r in result:
        outputf.write(r)

    # Do some logging at the end
    write_footer(outputf, beginning_time)

    outputf.close()


def count_interactions(grof, xtcf, btime, cutoff, debug):
    u = Universe(grof, xtcf)
    un_query = ('(resname PRO and (name CB or name CG or name CD)) or'
                '(resname VAL and (name CG1 or name CG2)) or'
                '(resname GLY and name CA) or'
                '(resname ALA and name CB)')
    vp_query = ('name OW')
    # MDAnalysis will convert the unit of length to angstrom, though in Gromacs the unit is nm
    un_atoms = u.selectAtoms(un_query)
    for ts in u.trajectory:
        if ts.time >= btime:
            numcount = 0
            tropo_vp_atoms = u.selectAtoms(
                '({0}) and around 8 ({1})'.format(vp_query, un_query))
            # different from when calculating unun, there is no overlap atom
            # between un_atoms & tropo_vp_atoms
            for ai in un_atoms:
                for aj in tropo_vp_atoms:
                    d = np.linalg.norm(ai.pos - aj.pos)
                    if d <= cutoff:
                        numcount += 1
            yield '{0:10.0f}{1:8d}\n'.format(ts.time,  numcount)
        # per 100 frames, num of frames changes with the size of xtc file, for debugging
        if debug and ts.frame % 2 == 0: 
            print "time: {0:10.0f}; step: {1:10d}; frame: {2:10d}".format(ts.time, ts.step, ts.frame)


def parse_cmd():
    parser = argparse.ArgumentParser(usage='used to calculate unvp')
    parser.add_argument('-f', '--xtcf', type=str, dest='xtcf', default=None,
                        help='Trajectory: xtc')
    parser.add_argument('-s', '--grof', type=str, dest='grof', default=None,
                        help='Structure: gro')
    parser.add_argument('-o', '--optf', type=str, dest='optf', default=None,
                        help='Output file name (OPT.)')
    parser.add_argument('-b', '--btime', type=int, dest='btime', default=0,
                        help='beginning time in ps (confirmed) xtc file records time in unit of ps')
    parser.add_argument('-c', '--cutoff', type=float, dest='cutoff', default=9999,
                        help='cutoff in nm (will be converted to angstrom in the code so as to work with MDAnalysis), default it 9999 (nm)')
    parser.add_argument('--debug', dest='debug', default=False, action='store_true',
                        help='for debugging, verbose info will be printted to the screen')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_cmd()
    main(args)
