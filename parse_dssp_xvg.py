#! /usr/bin/env python

import os
import re
import sys

"""needs to be substantially improved"""

def main(infile):
    ss_map = {'B-Sheet'  : 'E',                             # secondary structure map
              'A-Helix'  : 'H',
              'Turn'     : 'T',
              'Coil'     : 'C',
              '3-Helix'  : 'G',
              '5-Helix'  : 'I',
              'Bend'     : 'S',
              'B-Bridge' : 'B',
              'Structure': 'X'}

    regex = re.compile(r'@ s(\d) legend "(.*)"')
    # temp = re.compile(r'@ s\d legend ')
    

    # map ss to column, also providing info about available secondary structure
    col_map = {} 

    # dssp_count data structure
    # dssp_count = {
    #     "E": [[t1, count1], [t2, count2], [t3, count3] ... ]
    #     "H": [[t2, count2], [t2, count2], [t3, count3] ... ]
    #     ...
    #     }
    dssp_count = {}
    for ss in ss_map.values():                              # initialization
        dssp_count[ss] = []

    with open(infile) as inf:
        for line in inf:
            match = regex.search(line)
            if match:
                col_map[ss_map[match.group(2)]] = int(match.group(1)) + 1
            elif line[0] != '@' and line[0] != '#':
                convert_line_to_ss_count(line, dssp_count, col_map)

    outputfilepattern = infile.replace(".xvg", "_{0}.xvg")

    write_to_file(dssp_count, outputfilepattern)


def write_to_file(dssp_count, outputfilepattern):
    for ss in dssp_count:
        dirname = os.path.dirname(outputfilepattern)
        outputdir = os.path.join(dirname, 'r_dssp_{0}'.format(ss))
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)
    
        outputfile = os.path.join(outputdir, os.path.basename(outputfilepattern).format(ss))

        with open(outputfile, 'w') as opf:
            opf.writelines('{0:10s}{1:10s}\n'.format(entry[0], entry[1]) for entry in dssp_count[ss])

def convert_line_to_ss_count(line, dssp_count, col_map):
    sl = line.split()
    for ss in dssp_count:
        if ss in col_map:
            dssp_count[ss].append([sl[0], sl[col_map[ss]]])
        else:
            dssp_count[ss].append([sl[0], '0'])

if __name__ == "__main__":
    infiles = sys.argv[1:]
    for infile in infiles:
        print infile
        main(infile)
    # parser = my_basic_parser()
    # parser.add_argument('--fpp', dest='fpp', required=True, # fpp: file path pattern
    #                     help=("specify the file path pattern, e.g. r_dssp/{s}{c}{n}_dssp.xvg,"
    #                           "{s}{c}{n} are the variables"))
    # args = parser.parse_args()
    # for s in args.SEQS:
    #     for c in args.CDTS:
    #         for n in args.NUMS:
    #             infile = args.fpp.format(s=s, c=c, n=n)
    #             if not os.path.exists(infile):
    #                 print "{0} doesn't exist, YOU KONW THIS, RIGHT? ".format(infile)
    #             else:
    #                 main(infile)
    #                 print "{0} is done".format(infile)
