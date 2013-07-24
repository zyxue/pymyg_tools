#! /usr/bin/env python

"""Reference: Andreas Vitalis, Xiaoling Wang and Rohi V.Pappu 2008 JMB"""
import numpy as np
from MDAnalysis import Universe

import utils
 # import main_load, write_header, write_footer, backup, print_progress

def main(args):
    utils.main_load()
    outputfile =  args.optf if args.optf else '{0:s}.sespacing.xvg'.format(args.grof)
    utils.backup(outputfile)
    outputf = open(outputfile, 'w')
    beginning_time = utils.write_header(outputf)

    # This line will be used when there is a better code design
    # if ARGS.atom_sel is None:
    #     raise ValueError("atom_selection must be specified, check --atom_selection option!")

    # do calculation
    ijdist_dict = sequence_spacing(args.grof, args.xtcf, args.btime, args.etime,
                                   args.peptide_length, args.atom_sel)

    # cannot yield from sequence_spacing function because the result cannot be
    # calculated until all frames have been looped through

    # write headers
    outputf.write('# {0:8s}{1:20s}{2:20s}{3:10s}\n'.format('i-j', 'average', 'std', 'num_of_data_points'))
    # write results to the outputfile
    for k in sorted(ijdist_dict.keys()):
        data = np.array(ijdist_dict[k])
        mean = data.mean()                      # mean of ijdist
        std = data.std()                        # standard deviation of ijdist
        num = len(data)                         # num of data in that ijdist
        outputf.write('{0:8d}{1:20.8f}{2:20.8f}{3:10d}\n'.format(k, mean, std, num))

    # Do some logging at the end
    utils.write_footer(outputf, beginning_time)
    outputf.close()

def sequence_spacing(grof, xtcf, btime, etime, peptide_length, atom_sel):
    u = Universe(grof, xtcf)
    # this selection part should be better customized
    # here, only have been backbone atoms are used, u.selectAtoms doesn't
    # include Hydrogen atoms
    # REMMEMBER: ARGS verification should be done in main ONLY!
    # range works like this:

    # in MDAnalysis, resid starts from 1, in sequence_spacing.py, we don't count
    # the C- and N- termini, so it's from 2 to peptide_len+2
    residues = [u.selectAtoms(atom_sel.format(i)) for i in range(2, peptide_length + 2)]
    ijdist_dict = {}
    for ts in u.trajectory:
        # btime, etime defaults to 0, if etime is 0, loop till the end of the
        # trajectory
        if btime > ts.time:
            continue
        if etime > 0 and etime < ts.time:
            break

        # the good stuff
        for i, resi in enumerate(residues):
            for j, resj in enumerate(residues):
                # to remove duplicate since resi & resj are within the same peptide 
                if i < j:
                    dij = abs(i - j)
                    d_atomi_atomj = []
                    # loop through every atom in both residues
                    for atomi in resi:
                        for atomj in resj:
                            d_atomi_atomj.append(
                                np.linalg.norm(atomi.pos - atomj.pos))
                # add the result to the dictionary
                    ij_dist = np.average(d_atomi_atomj)   # distance between i and j
                    if dij not in ijdist_dict.keys():
                        ijdist_dict[dij] = [ij_dist]
                    else:
                        ijdist_dict[dij].append(ij_dist)
        utils.print_progress(ts)

    return ijdist_dict

def parse_cmd():
    usage='used to calculate sequence spacing'
    p = utils.get_basic_parser(usage=usage)
    p.add_argument('--pl', type=int, dest='peptide_length', default=None, required=True,
                        help='specify the peptide_length')
    p.add_argument('--atom-selection', type=str, dest='atom_sel', default='resid {0} and not type H',
                        help=("atom selection, the default is 'resid {0} and not type H': {0}"
                              "will be substituted with a number, the atom selection option only"
                              "applies to each resid! I know, it's WIRED and UGLY, let me know"
                              "if you have a better way. Thanks! 2011-12-12"))
    args = p.parse_args()
    return args

if __name__ == "__main__":
    args = parse_cmd()
    main(args)
