#! /usr/bin/env python

"""Functionize it and make it act like a gromacs tool"""

import os
import sys
import numpy as np
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level='DEBUG')

import tables
from MDAnalysis import Universe

import utils

def main(cmd_args):
    args = get_args(cmd_args)

    utils.main_load()
    
    output = args.optf
    if output is None:
        # it's a log since the results are written to the h5 file directly
        outputfile = '{0:s}.output.log'.format(args.grof)
    else:
        outputfile = output

    utils.backup(outputfile)

    outputf = open(outputfile, 'w')
    beginning_time = utils.write_header(outputf)

    A = args
    if not os.path.exists(A.h5): 
        raise IOError('{0} does not exist'.format(A.h5))

    # *10: convert to angstrom from nm
    result = count_interactions(A)
    path = os.path.join('/', os.path.dirname(A.xtcf))
    tb_name = os.path.join(path, 'unun_map')

    h5 = tables.openFile(A.h5, mode='a')
    if h5.__contains__(tb_name):
        logger.info('found {0} already in {0}, replacing with new calculated values'.format(tb_name, A.h5))
        _ = h5.getNode(tb_name)
        _.remove()
    h5.createArray(where=path, name='unun_map', object=result)
    h5.close()

    utils.write_footer(outputf, beginning_time)
    outputf.close()

def count_interactions(A):
    logger.debug('loading {0}'.format(A.grof))
    univ = Universe(A.grof)
    logger.debug('loaded {0}'.format(A.grof))

    pro_atoms = univ.selectAtoms('protein and not resname ACE and not resname NH2')
    pl = pro_atoms.residues.numberOfResidues()
    # +1: for missing resname ACE, such that it's easier to proceed in the next
    # step

    logger.debug('loading {0}, {1}'.format(A.grof, A.xtcf))
    u = Universe(A.grof, A.xtcf)
    logger.debug('loaded {0}, {1}'.format(A.grof, A.xtcf))

    # Just for reference to the content of query when then code was first
    # written and used
    # query = ('(resname PRO and (name CB or name CG or name CD)) or'
    #          '(resname VAL and (name CG1 or name CG2)) or'
    #          '(resname GLY and name CA) or'
    #          '(resname ALA and name CB)')

    query = A.query
    atoms = u.selectAtoms(query)
    logger.info('Number of atoms selected: {0}'.format(atoms.numberOfAtoms()))

    # MDAnalysis will convert the unit of length to angstrom, though in Gromacs
    # the unit is nm
    cutoff = A.cutoff * 10
    nres_away = A.nres_away
    btime = A.btime
    etime = A.etime
    nframe = 0
    unun_map = None
    for ts in u.trajectory:
        if btime > ts.time:
            continue
        if etime > 0 and etime < ts.time:
            break

        nframe += 1
        map_ = np.zeros((pl+1, pl+1))                   # map for a single frame
        for i, ai in enumerate(atoms):
            ai_resid = ai.resid
            for j, aj in enumerate(atoms):
                aj_resid = aj.resid
                # to avoid counting the same pair twices,
                # the 2 resid cannot be neigbors
                if i < j and aj_resid - ai_resid >= nres_away:
                    d = np.linalg.norm(ai.pos - aj.pos)
                    if d <= cutoff:
                        # -1: resid in MDAnalysis starts from 1
                        map_[ai_resid-1][aj_resid-1] += 1
        if unun_map is None:
            unun_map = map_
        else:
            unun_map = unun_map + map_
        utils.print_progress(ts)
    sys.stdout.write("\n")
    return unun_map / float(nframe)

def get_args(cmd_args):
    parser = utils.get_basic_parser(usage='used to calculate unun, neighbour residues are excluded')

    parser.add_argument('--query', required=True, 
                        help='query for selecting atoms used for calculating non-polar interactions')
    parser.add_argument('--h5', required=True, 
                        help='written to hdf5 file directory')
    parser.add_argument('-c', '--cutoff', default=0.55, type=float,
                        help='cutoff in nm (will be converted to angstrom in the code so as to work with MDAnalysis)')
    parser.add_argument('--nres-away', default=2, type=int,
                        help='the atoms should be nres away to avoid neigbor interactions)')

    args = parser.parse_args(cmd_args)
    return args

if __name__ == "__main__":
    main(sys.argv[1:])
