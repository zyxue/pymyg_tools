#! /usr/bin/env python

"""Functionize it and make it act like a gromacs tool"""

import time
import os
import sys
import numpy as np
import argparse

import tables
from MDAnalysis import Universe

def timeit(method):
    def timed(*args, **kw):
        sys.stdout.write('{0}\n'.format("#" * 79))
        sys.stdout.write("{0:^79s}\n".format("SUCCESSFULLY LOAD main()! by Zhuyi Xue zhuyi.xue@utoronto.ca"))
        sys.stdout.write('{0}\n'.format("#" * 79))

        ts = time.time()
        sys.stdout.write('# main calculations BEGIN at: {0:s}\n'.format(time.ctime()))
        sys.stdout.write('# {0:s}\n\n'.format(' '.join(sys.argv[:])))
        sys.stdout.write('{0:s}\n\n'.format('#' * 50))

        res = method(*args, **kw)

        te = time.time()
        delta_time = te - ts
        sys.stdout.write('\n{0:s}\n\n'.format('#' * 50))
        sys.stdout.write('# main calculations END at: {0:s}\n'.format(time.ctime()))
        sys.stdout.write('# time consumed in TOTOAL: {0}\n'.format(
                time.strftime('%H:%M:%S', time.gmtime(delta_time))))
        return res
    return timed

@timeit
def main(args):
    A = args

    if not A.h5: 
        raise IOError('h5 output file not specified')

    # *10: convert to angstrom from nm
    result = count_interactions(A)
    path = os.path.join('/', os.path.dirname(A.xtcf))
    tb_name = os.path.join(path, 'unun_map')

    h5 = tables.openFile(A.h5, mode='a')
    if h5.__contains__(tb_name):
        _ = h5.getNode(tb_name)
        _.remove()
    h5.createArray(where=path, name='unun_map', object=result)

    h5.close()

def count_interactions(A):
    univ = Universe(A.grof)
    pro_atoms = univ.selectAtoms('protein and not resname ACE and not resname NH2')
    pl = pro_atoms.residues.numberOfResidues()
    # +1: for missing resname ACE, such that it's easier to proceed in the next
    # step
    unun_maps = []
    u = Universe(A.grof, A.xtcf)
    query = ('(resname PRO and (name CB or name CG or name CD)) or'
             '(resname VAL and (name CG1 or name CG2)) or'
             '(resname GLY and name CA) or'
             '(resname ALA and name CB)')
    # MDAnalysis will convert the unit of length to angstrom, though in Gromacs the unit is nm
    atoms = u.selectAtoms(query)
    cutoff = A.cutoff * 10
    for ts in u.trajectory:
        if ts.time >= A.btime:
            map_ = np.zeros((pl+1, pl+1))                   # map for a single frame
            for i, ai in enumerate(atoms):
                for j, aj in enumerate(atoms):
                    # to avoid counting the same pair twices,
                    # the 2 resid cannot be neigbors
                    if i < j and abs(ai.resid - aj.resid) >= 2: 
                        d = np.linalg.norm(ai.pos - aj.pos)
                        if d <= cutoff:
                            # -1: resid in MDAnalysis starts from 1
                            map_[ai.resid-1][aj.resid-1] += 1
            unun_maps.append(map_)
        # per 100 frames, num of frames changes with the size of xtc file, for
        # verbosity
        if ts.frame % 2 == 0: 
            sys.stdout.write("\r[38;5;226mtime: {0:10.0f}; step: {1:10d}; frame: {2:10d}".format(ts.time, ts.step, ts.frame))
            sys.stdout.flush()
    sys.stdout.write("\n")                                  # to get back to the normal color
    return np.array(unun_maps).mean(axis=0)

def get_args():
    parser = argparse.ArgumentParser(usage='used to calculate unun, neighbour residues are excluded')
    parser.add_argument('-f', '--xtcf', type=str, dest='xtcf', default=None,
                        help='Trajectory: xtc')
    parser.add_argument('-s', '--grof', type=str, dest='grof', default=None,
                        help='Structure: gro')
    parser.add_argument('--h5', help='written to hdf5 file directory')
    parser.add_argument('-b', '--btime', type=int, dest='btime', default=0,
                        help='beginning time in ps (confirmed) xtc file records time in unit of ps')
    parser.add_argument('-c', '--cutoff', default=0.55, type=float,
                        help=('cutoff in nm (will be converted to angstrom in the code '
                              'so as to work with MDAnalysis)'))
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)
