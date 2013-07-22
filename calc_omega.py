#!/usr/bin/env python

import numpy as np
from MDAnalysis import Universe
from MDAnalysis.core.AtomGroup import AtomGroup

import utils

def main(args):
    utils.main_load()

    # check the validity of output file name, do backup
    output = args.optf
    if output is None:
        outputfile = '{0:s}.output.xvg'.format(args.grof)
    else:
        outputfile = output

    utils.backup(outputfile)

    # Do some logging at the beginning
    outputf = open(outputfile, 'w')
    beginning_time = utils.write_header(outputf)

    # do calculation
    result = calc_dihedral(args.grof, args.xtcf, args.btime, args.etime)

    # write results to the outputfile
    for r in result:
        outputf.write(r)

    # Do some logging at the end
    utils.write_footer(outputf, beginning_time)

    outputf.close()

def verify_atom_order(ca1, c, n, ca2):
    assert ca1.resid == c.resid
    assert n.resid == ca2.resid
    assert ca2.resid - ca1.resid == 1

def select_dihedrals(univer):
    CA = univer.selectAtoms('name CA and not resname ACE and not resname NH2')
    C  = univer.selectAtoms('name C  and not resname ACE and not resname NH2')
    N  = univer.selectAtoms('name N  and not resname ACE and not resname NH2')

    ##########ILLUSTRATION: 2 DIHEDRAL ANGLE FORMD BY 3 RESIDUES############
    # dihedral: // or \\
    #             O
    #             "
    #   Ca   N    C    Ca   O-H
    #  /  \ //\  / \\ /  \ /
    # N    C   Ca    N    C
    #      "              "
    #      O              O

    tets = []
    for ca1, c, n, ca2 in zip(CA[:-1], C[:-1], N[1:], CA[1:]):
        verify_atom_order(ca1, c, n, ca2)
        tet = [ca1, c, n, ca2]
        tets.append(AtomGroup(tet))

    # TEST CASE: angle around 180
    # ca1 = univer.selectAtoms("resid 18 and name CA")[0]
    # c   = univer.selectAtoms("resid 18 and name C")[0]
    # n   = univer.selectAtoms("resid 19 and name N")[0]
    # ca2 = univer.selectAtoms("resid 19 and name CA")[0]
    # tets = [AtomGroup([ca1, c, n, ca2])]
    return tets

def calc_dih(tets, t=0):
    dihs = []
    for tet in tets:
        dih = np.float64(np.nan_to_num(tet.dihedral()))
        # for some reason numpy.float32, the type of which dihedral
        # angle is returned in cannot be formatted by {0:f} either
        # (python 2.7.2 + numpy 1.6.1) or (python 2.7.2 + numpy 1.6.1)
        # works, use convert to float

        # equivalent to delta = abs(abs(dih) - 180) - abs(dih - 0)
        delta = abs(abs(dih) - 180) - abs(dih)
        if delta <= 0:                                   # closer to 180: trans
            dih = 0
        else:
            dih = 1                                       # closer to 0  : cis
        dihs.append(dih)
    return ' {0:<8.0f}{1}\n'.format(t, ' '.join('{0:<4d}'.format(d) for d in dihs))

def calc_dihedral(grof, xtcf, btime, etime):
    # xtcf=None, so if only gro file is parsed, it still works
    univer = Universe(grof, xtcf)
    
    tets = select_dihedrals(univer)     # there should be a better name for tet

    # Write headers, hdr: header
    hdrs = []
    for k, tet in enumerate(tets):
        # ca1 + ca2 + (NO. peptide-bond)
        hdr = (utils.swap_aa_name(tet[0].resname) +
               utils.swap_aa_name(tet[-1].resname) + 
               "{0:02d}".format(k+1))
        # hdrs.append("{0:8s}".format(hdr))
        hdrs.append("{0:4s}".format(hdr))
    yield '#{0:8s}{1}\n'.format('t(ps)', ' '.join(hdrs))

    if not xtcf:
        yield calc_dih(tets)
    else:
        for ts in univer.trajectory:
            if btime > ts.time:
                continue
            if etime > 0 and etime < ts.time:
                break

            res = calc_dih(tets, ts.time)
            yield res
            utils.print_progress(ts)

if __name__ == "__main__":
    parser = utils.get_basic_parser()
    args = parser.parse_args()
    main(args)
