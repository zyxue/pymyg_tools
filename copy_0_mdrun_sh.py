#! /usr/bin/env python

"""Functionize it and make it act like a gromacs tool"""

import os
import sys
import time
import numpy as np
from MDAnalysis import Universe
from optparse import OptionParser

def main(OPTIONS):
    print "#" * 79
    print "{0:^79s}".format("SUCCESSFULLY LOAD main()! by Zhuyi Xue zhuyi.xue@utoronto.ca")
    print "#" * 79

    # check the validity of output file name, do backup
    output = OPTIONS.optf
    if output is None:
        outpufile = '{0:s}_output.xvg'.format(gro)
    else:
        outputfile = output
    backup_old_output(outputfile)

    if OPTIONS.atom_sel is None:
        raise ValueError("atom_selection must be specified, check --atom_selection option!")

    # Do some logging at the beginning
    outputfile = open(output, 'w')
    beginning_time = time.time()
    outputfile.write('# main calculations BEGIN at: {0:s}\n\n'.format(time.ctime()))
    outputfile.write('# {0:s}\n\n'.format(' '.join(sys.argv[:])))

    # do calculation
    ijdist_dict = sequence_spacing(OPTIONS.pf, 
                                   OPTIONS.grof, 
                                   OPTIONS.xtcf, 
                                   OPTIONS.peptide_length, 
                                   OPTIONS.atom_sel,
                                   outputfile)

    # Do some logging at the end
    ending_time = time.time()
    delta_time = ending_time - beginning_time
    outputfile.write('# main calculations END at: {0:s}\n'.format(time.ctime()))
    outputfile.write('# time consumed in TOTOAL: {0:.5f} (seconds)\n'.format(delta_time))
    outputfile.write('{0:s}\n\n'.format('#' * 50))

    # write results to the outputfile
    outputfile.write('{0:8s}{1:20s}{2:20s}{3:10s}\n'.format('# i-j', 'average', 'std', 'num'))
    for k in sorted(ijdist_dict.keys()):
        data = np.array(ijdist_dict[k])
        mean = data.mean()                      # mean of ijdist
        std = data.std()                        # standard deviation of ijdist
        num = len(data)                         # num of data in that ijdist
        outputfile.write('{0:8d}{1:20.8f}{2:20.8f}{3:10d}\n'.format(k, mean, std, num))

    outputfile.close()

def backup_old_output(output):
    if os.path.exists(output):
        count = 1
        rn_to = '#' + output + '.{0}#'.format(count)                 # rename to
        while os.path.exists(rn_to):
            count += 1
            rn_to = '#' + output + '.{0}#'.format(count)                 # rename to
        os.rename(output, rn_to)

def sequence_spacing(pf, grof, xtcf, peptide_length, atom_sel, output=None):
    u = Universe(grof, xtcf)
    # this selection part should be better customized
    # here, only have been backbone atoms are used, u.selectAtoms doesn't
    # include Hydrogen atoms
    # REMMEMBER: OPTIONS verification should be done in main ONLY!
    residues = [u.selectAtoms(atom_sel.format(i)) for i in range(2, peptide_length)]
    ijdist_dict = {}
    for ts in u.trajectory:
        for i, resi in enumerate(residues):
            for j, resj in enumerate(residues):
                if i < j:
                    resi_pos = resi.centerOfGeometry()                # residue i position
                    resj_pos = resj.centerOfGeometry()                # residue j position
                    ijdist = np.linalg.norm(resi_pos - resj_pos)
                    dij = j - i                                       # distance between i and j
                    if dij not in ijdist_dict.keys():
                        ijdist_dict[dij] = [dij]
                    else:
                        ijdist_dict[dij].append(ijdist)
        if ts.step % 2000000 == 0:                             # 2000ps
            print "time step: {0:d}".format(ts.step)
    return ijdist_dict

def parse_cmd():
    parser = OptionParser(usage='used to calculate sequence spacing"')
    parser.add_option('--pf', type='str', dest='pf', default=None,
                      help='python code needs it due to not so good design')
    parser.add_option('-f', '--xtcf', type='str', dest='xtcf', default=None,
                      help='Trajectory: xtc')
    parser.add_option('-s', '--grof', type='str', dest='grof', default=None,
                      help='Structure: gro')
    parser.add_option('-o', '--optf', type='str', dest='optf', default=None,
                      help='Output file name (OPT.)')
    parser.add_option('-l', '--peptide_length', type='int', dest='peptide_length', default=None,
                      help='specify the peptide_length')
    parser.add_option('--atom-selection', type='str', dest='atom_sel', default=None,
                      help="atom selection, e.g. \'resid {0} and not type H': {0} will be substituted with a number, the atom selection option only applies to each resid! I know, it's WIRED and UGLY, let me know if you have a better way. Thanks! 2011-12-12")
    OPTIONS, args = parser.parse_args()
    return OPTIONS

if __name__ == "__main__":
    options = parse_cmd()
    main(options)
