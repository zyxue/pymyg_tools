#! /usr/bin/env python

"""
This is trying to do the same thing as g_mindist_excl1 for unun, but without
duplication count, the drawback is that python code is much much slower!
"""

import numpy as np
from MDAnalysis import Universe

import utils
# from utils import main_load, write_header, write_footer, backup

def main(args):
    utils.main_load()
    outputfile =  args.optf if args.optf else '{0:s}.unun.xvg'.format(args.grof)
    utils.backup(outputfile)
    outputf = open(outputfile, 'w')
    beginning_time = utils.write_header(outputf)

    result = count_interactions(args.grof, args.xtcf, args.btime, args.etime, args.cutoff)

    # write headers
    outputf.write('# {0:>10s}{1:>8s}\n'.format('time', 'num'))
    # write results to the outputfile
    for r in result:
        outputf.write(r)

    # Do some logging at the end
    utils.write_footer(outputf, beginning_time)
    outputf.close()

def count_interactions(grof, xtcf, btime, etime, cutoff):
    cutoff = cutoff * 10 # * 10: convert from nm to angstrom to work with MDAnalysis
    u = Universe(grof, xtcf)
    query = ('(resname PRO and (name CB or name CG or name CD)) or'
             '(resname VAL and (name CG1 or name CG2)) or'
             '(resname GLY and name CA) or'
             '(resname ALA and name CB)')
    # MDAnalysis will convert the unit of length to angstrom, though in Gromacs the unit is nm
    atoms = u.selectAtoms(query)
    for ts in u.trajectory:
        if btime > ts.time:
            continue
        if etime > 0 and etime < ts.time:
            break

        numcount = 0
        for i, ai in enumerate(atoms):
            for j, aj in enumerate(atoms):
                # to avoid counting the same pair twices,
                # the 2 resid cannot be neigbors
                if i < j and abs(ai.resid - aj.resid) >= 2: 
                    d = np.linalg.norm(ai.pos - aj.pos)
                    if d <= cutoff:
                        numcount += 1
        yield '{0:10.0f}{1:8d}\n'.format(ts.time,  numcount)
        utils.print_progress(ts)

def get_args():
    p = utils.get_basic_parser(usage='used to calculate unun, neighbour residues are excluded')
    p.add_argument('-c', '--cutoff', type=float, required=True, 
                   help='cutoff in nm (will be converted to angstrom in the code so as to work with MDAnalysis)')
    p.add_argument('--debug', dest='debug', default=False, action='store_true',
                   help='for debugging, verbose info will be printted to the screen')
    args = p.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)
