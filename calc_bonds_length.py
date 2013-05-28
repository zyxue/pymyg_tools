#!/usr/bin/env python

import argparse

import numpy as np
from MDAnalysis import Universe

from utils import main_load, write_header, write_footer, backup_old_output

def main(args):
    if args.debug:
        main_load()

    # check the validity of output file name, do backup
    output = args.optf
    if output is None:
        outputfile = '{0:s}.output.xvg'.format(args.grof)
    else:
        outputfile = output
    backup_old_output(outputfile)

    # Do some logging at the beginning
    outputf = open(outputfile, 'w')
    beginning_time = write_header(outputf)

    # do calculation
    result = calc_bond_length(
        args.grof,
        args.xtcf,
        args.btime,
        args.etime,
        args.debug)

    # write results to the outputfile
    for r in result:
        outputf.write(r)

    # Do some logging at the end
    write_footer(outputf, beginning_time)

    outputf.close()

# grof = '/scratch/p/pomes/zyxue/mono_su_as/w300/sq1w/sq1w00/sq1w00_pro.gro'
# xtcf = '/scratch/p/pomes/zyxue/mono_su_as/w300/sq1w/sq1w00/sq1w00_pro.xtc'

def calc_bond_length(grof, xtcf, btime, etime, debug):
    # thebonds contains all the bonds that I am interested
    thebonds = { #atom names should be UNIQUE within each residue for this script
        'BACKBONE_INTRA': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ],       # backbone, intramolecular interactions
        # PB: peptide bond, which is the only intermolecular bonds that I am interested
        'PB':  [('C', 'N'),],

        'GLY': [('CA', 'HA1'),],
        'PRO': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'N' )],
        'VAL': [('CA', 'CB'), ('CB', 'CG1'), ('CB', 'CG2')],
        
        'MeO': [('C', 'OA'), ('C', 'H'), ('OA', 'HO')],
        'SOL': [('OW', 'HW1')],
        }

    aas = ['GLY', 'PRO', 'VAL']                                      # rl: residue list
    solvents = ['MeO', 'SOL']

    # initialize ibonds
    ibonds = {}                    # interested bonds, not very legible to human
    for k in thebonds:
        ibonds[k] = {}
        if k in aas:
            for kk in thebonds[k] + thebonds['BACKBONE_INTRA']:
                ibonds[k][tuple(sorted(kk))] = []
        elif k in solvents:
            for kk in thebonds[k]:
                ibonds[k][tuple(sorted(kk))] = []
            

    ibonds['PB'] = {}
    ibonds['PB'][('C', 'N')] = []

    # data structure would be (to do)
    # ibonds = {
    #     'PRO': {
    #         (a1, b1):[(atom_object1, atom_object2), (atom_object3, atom_object4), ... , ],
    #         (a2, b2):[(atom_object1, atom_object2), (atom_object3, atom_object4), ... , ],
    #         ...
    #         },
    #     'VAL': {
    #         (a1, b1):[(atom_object1, atom_object2), (atom_object3, atom_object4), ... , ],
    #         (a2, b2):[(atom_object1, atom_object2), (atom_object3, atom_object4), ... , ],
    #         ...
    #         },
    #     'GLY': {
    #         (a1, b1):[(atom_object1, atom_object2), (atom_object3, atom_object4), ... , ],
    #         (a2, b2):[(atom_object1, atom_object2), (atom_object3, atom_object4), ... , ],
    #         ...
    #         },
    #     }

    univer = Universe(grof, xtcf)

    atom_selection = "not resname ACE and not resname NH2"            # get rid of the ends
    # atom_selection = "resname MeO and resid 3000"
    atoms = univer.selectAtoms(atom_selection)

    # initialize ibonds data structure
    # a bondname is composed of readable plain text
    # a bond is composed of Atom object
    for ki, ai in enumerate(atoms):
        for kj, aj in enumerate(atoms):
            if ki < kj:
                if ai.resid == aj.resid: # collecting intramolecular bonds associated with real atom objects
                    resname= ai.resname                               # will also equal aj.resname
                    bondname = tuple(sorted([ai.name, aj.name]))
                    if bondname in ibonds[resname]:
                        bond = [ai, aj]
                        ibonds[resname][bondname].append(bond)
                elif ai.resid - aj.resid == -1: # collecting itermolecular bonds: i.e. peptide bond
                    bondname = tuple([ai.name, aj.name])
                    if bondname == ('C', 'N'):
                        bond = [ai, aj]
                        ibonds['PB'][bondname].append(bond)

################################################################################
# VERIFICATION STATUS: ibonds initiation verified for
# sq1w00_md.gro & sq1m00_md.gro
# 2012-04-25
#     for i in ibonds:
#         for j in ibonds[i]:
#             print i, j, len(ibonds[i][j])
    
#     from pprint import pprint as pp
#     pp(ibonds)

# VAL ('CB', 'CG2') 14
# VAL ('C', 'CA') 14
# VAL ('CA', 'N') 14
# VAL ('CB', 'CG1') 14
# VAL ('C', 'O') 14
# VAL ('CA', 'CB') 14
# PRO ('CD', 'CG') 7
# PRO ('C', 'CA') 7
# PRO ('CA', 'N') 7
# PRO ('CA', 'CB') 7
# PRO ('C', 'O') 7
# PRO ('CD', 'N') 7
# PRO ('CB', 'CG') 7
# SOL ('HW1', 'OW') 0
# PB ('C', 'N') 34
# GLY ('CA', 'N') 14
# GLY ('C', 'O') 14
# GLY ('CA', 'HA1') 14
# GLY ('C', 'CA') 14
# MeO ('HO', 'OA') 0
# MeO ('C', 'OA') 0
# MeO ('C', 'H') 0

#     import sys
#     sys.exit()

################################################################################

    # Just for Printing the Header
    sorted_resname = sorted(ibonds.keys()) # sort to keep the value in the right order
    partial_header = []
    for resname in sorted_resname:
        resname_header = []                         # the header specific to residue
        for bondname in sorted(ibonds[resname].keys()):
            # bn: since bondname has been used in previous codes
            bn = '{0}|{1}'.format(resname[0], '-'.join(bondname))
            resname_header.append('{0:9s}'.format(bn))
        partial_header.extend(resname_header)
    yield '#{0:8s}{1}\n'.format('t(ps)', ''.join(partial_header))

    # import sys
    # sys.exit()

    # Production Calculation
    # use < when for formatting values to align with headers, and the width will
    # be 1 col narrower than that in the corresponding header
    for ts in univer.trajectory:
        # for debugging only
        if debug and ts.frame % 2 == 0:
            print "time: {0:10.0f}; step: {1:10d}; frame: {2:10d}".format(ts.time, ts.step, ts.frame)

        if etime > ts.time >= btime:
            partial_yield = []
            for resname in sorted_resname:
                resname_yield = []
                for bondname in sorted(ibonds[resname].keys()):
                    bonds = ibonds[resname][bondname]
                    ds = []
                    for bond in bonds:
                        r = bond[0].pos - bond[1].pos # end-to-end vector from atom positions
                        d = np.linalg.norm(r)  # distance
                        ds.append(d)
                    resname_yield.append('{0:<8.3f}'.format(np.average(ds))) #, np.std(ds))
                partial_yield.extend(resname_yield)
            # a space in order to align with # in the header
            yield ' {0:<8.0f}{1}\n'.format(ts.time, ' '.join(partial_yield))

def parse_cmd():
    usage = 'calc_bond_length.py analyze the bond lengths of interested bonds'
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('-f', '--xtcf', type=str, dest='xtcf', default=None,
                        help='Trajectory: xtc')
    parser.add_argument('-s', '--grof', type=str, dest='grof', default=None,
                        help='Structure: gro')
    parser.add_argument('-o', '--optf', type=str, dest='optf', default=None,
                        help='Output file name (OPT.)')
    parser.add_argument('-b', '--btime', type=int, dest='btime', default=0,
                        help='beginning time in ps (confirmed) xtc file records time in unit of ps')
    parser.add_argument('-e', '--etime', type=int, dest='etime', default=0,
                        help='ending time in ps (confirmed) xtc file records time in unit of ps')
    parser.add_argument('--debug', dest='debug', default=False, action='store_true',
                        help='for debugging, verbose info will be printted to the screen')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_cmd()
    main(args)
