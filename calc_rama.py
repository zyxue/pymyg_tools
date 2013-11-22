#! /usr/bin/env python

"""
This is a make up for g_rama from gromacs since the former cannot identify rama
from sr[1-3] in CHARMM force field. Look at the gro/pdb file of sr[1-3] for
details about why. Because the N-ter NH2 group is part of the center residue in
CHARMM force field.

WARNING: this script only works for sr[1-3] for when it was written on
2013-07-24
"""

from MDAnalysis import Universe

import utils as U

def main(args):
    U.main_load()
    outputfile =  args.optf if args.optf else '{0:s}.rama.xvg'.format(args.grof)
    U.backup(outputfile)
    outputf = open(outputfile, 'w')
    beginning_time = U.write_header(outputf)

    result = calc_rama(args.grof, args.xtcf, args.btime, args.etime)

    # write headers
    outputf.write('# {0:>10s}{1:>8s}\n'.format('phi', 'psi', 'resname-resid'))
    # write results to the outputfile
    for r in result:
        outputf.write(r)

    # Do some logging at the end
    U.write_footer(outputf, beginning_time)
    outputf.close()

def calc_rama(grof, xtcf, btime, etime):
    u = Universe(grof, xtcf)

    resname_query = 'resname GLY or resname VAL or resname PRO'
    atoms = u.selectAtoms(resname_query)
    resname = atoms.resnames()[0] # [0] because .resnames() returns a list of one element
    resid = atoms.resids()[0] # [0] because .resnames() returns a list of one element

    phi_query = ('(resname ACE and name C) or '
                 '(resname GLY or resname VAL or resname PRO and '
                 '(name N or name CA or name C))')

    psi_query = ('(resname GLY or resname VAL or resname PRO and (name N or name CA or name C or name NT)) or '
                 '(resname NH2 and name N)')

    # MDAnalysis will convert the unit of length to angstrom, though in Gromacs the unit is nm
    phi = u.selectAtoms(phi_query)
    psi = u.selectAtoms(psi_query)

    for _ in phi.atoms:
        print _

    for _ in psi.atoms:
        print _


    for ts in u.trajectory:
        if btime > ts.time:
            continue
        if etime > 0 and etime < ts.time:
            break

        yield '{0:.3f}  {1:.3f}  {2}-{3}\n'.format(
            phi.dihedral(), psi.dihedral(), resname, resid)
        U.print_progress(ts)

def get_args():
    p = U.get_basic_parser(usage='used to calculate unun, neighbour residues are excluded')
    args = p.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)
