#! /usr/bin/env python

import os
import sys

"""needs to be substantially improved"""

def main(infile):
    opfs = {}
    with open(infile) as inf:
        for line in inf:
            if line[0] == '@' or line[0] == '#':
                continue
            if not line.strip():
                continue

            sl = line.split()
            resname = sl[2][:3]
            res_dir = 'r_rama_{0}'.format(resname)
            if not os.path.exists(res_dir):
                os.mkdir(res_dir)
            if resname not in opfs:
                opf_name = os.path.join(
                    res_dir, infile.replace('.xvg', '_{0}.xvg'.format(resname)))
                opfs[resname] = open(opf_name, 'w')
            opfs[resname].write(line)

if __name__ == "__main__":
    infiles = sys.argv[1:]
    for infile in infiles:
        print infile
        main(infile)
