#! /usr/bin/env python

"""Functionize it and make it act like a gromacs tool"""

# import os
# import sys
# import time
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
    result = calc_rg(
        ARGS.grof, 
        ARGS.xtcf, 
        ARGS.btime,
        ARGS.debug)

    # write results to the outputfile
    outputf.write('# {0:>10s}{1:>8s}\n'.format('time', 'rg(A)'))
    for r in result:
        outputf.write(r)

    # Do some logging at the end
    write_footer(outputf, beginning_time)

    outputf.close()


def calc_rg(grof, xtcf, btime, debug):
    u = Universe(grof, xtcf)
    query = 'name CA'
    # MDAnalysis will convert the unit of length to angstrom, though in Gromacs the unit is nm
    atoms = u.selectAtoms(query)
    natoms = atoms.numberOfAtoms()
    for ts in u.trajectory:
        if ts.time >= btime:
            com = atoms.centerOfMass()                                # center of mass
            _sum = sum((sum(i**2 for i in (a.pos - com)) for a in atoms))
            rg = np.sqrt(_sum / natoms)
            yield '{0:10.0f}{1:15.6f}\n'.format(ts.time, rg)
        # per 100 frames, num of frames changes with the size of xtc file, for debugging
        if debug and ts.frame % 2 == 0: 
            print "time: {0:10.0f}; step: {1:10d}; frame: {2:10d}".format(ts.time, ts.step, ts.frame)

def parse_cmd():
    usage = 'myrg only calculates the radius of gyration of CA'
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('-f', '--xtcf', type=str, dest='xtcf', default=None,
                        help='Trajectory: xtc')
    parser.add_argument('-s', '--grof', type=str, dest='grof', default=None,
                        help='Structure: gro')
    parser.add_argument('-o', '--optf', type=str, dest='optf', default=None,
                        help='Output file name (OPT.)')
    parser.add_argument('-b', '--btime', type=int, dest='btime', default=0,
                        help='beginning time in ps (confirmed) xtc file records time in unit of ps')
    parser.add_argument('--debug', dest='debug', default=False, action='store_true',
                        help='for debugging, verbose info will be printted to the screen')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_cmd()
    main(args)
