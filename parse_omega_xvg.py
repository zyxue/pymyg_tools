#! /usr/bin/env python


import sys
import os
import re

import numpy as np

PROP1 = 'omega_x_pro_percent'
PROP2 = 'omega_x_y_percent'     # y != pro

def main(xvgf):
    XP_RE = '[ACDEFGHIKLMNPQRSTVWY]P[0-9]+'
    omega_list, data = process_xvg(xvgf)
    x_pro_list = [_ for _ in omega_list if re.match(XP_RE, _)]
    x_y_list = [_ for _ in omega_list if not re.match(XP_RE, _)]
    print 'x_pro_list: {0}'.format(x_pro_list)
    print 'x_y_list: {0}'.format(x_y_list)

    if data.shape == (0,):
        raise ValueError("no data in {0}, please check.".format(xvgf))

    time = data[:, 0]
    # count number of frames from time.shape. Surely it can be obtained from
    # any element of omega_list, as well.
    nframes = time.shape[0] 
    omegas = data[:, 1:].transpose()

    # trans: 0,  cis: 1
    cis_x_pro, cis_x_y = [], []
    for om_name, om  in zip(omega_list, omegas):    # om: omega
        match = re.search(XP_RE, om_name)  # X-Pro
        if match:
            cis_x_pro.append(om.sum())
        else:
            cis_x_y.append(om.sum())

    cis_x_y_percent = sum(cis_x_y) / float(len(x_y_list) * nframes)
    if x_pro_list:
        cis_x_pro_percent = sum(cis_x_pro) / float(len(x_pro_list) * nframes)
    else:
        cis_x_pro_percent = 0

    write_to_file(xvgf, cis_x_pro_percent, cis_x_y_percent)

def process_xvg(xvgf):
    data = []
    with open(xvgf) as inf:
        for line in inf:
            if line.startswith('#t(ps)'):
                omega_list = line.split()[1:]
            elif not line.strip() or line.startswith('#') or line.startswith('@'):
                pass
            else:
                data.append([float(i) for i in line.split()])
    data = np.array(data)
    return omega_list, data

def write_to_file(xvgf, cis_x_pro, cis_x_y):
    dirname = os.path.dirname(xvgf)
    dir_x_pro = os.path.join(dirname, 'r_{0}'.format(PROP1))
    dir_x_y = os.path.join(dirname, 'r_{0}'.format(PROP2))

    for _ in [dir_x_pro, dir_x_y]:
        if not os.path.exists(_):
            os.mkdir(_)
        
    basename = os.path.basename(xvgf)
    opf_x_pro = os.path.join(dir_x_pro, basename.replace('omega', PROP1))
    opf_x_y = os.path.join(dir_x_y, basename.replace('omega', PROP2))
    with open(opf_x_pro, 'w') as opf:
        opf.writelines('# {0:<20s}{1:<20s}{2:<20s}\n'.format(
                'replica_id', 'trans_x_pro', 'cis_x_pro'))
                
        # use it as a title
        opf.writelines('  {0:<20s}{1:<20.5f}{2:<20.5f}\n'.format(
                xvgf, 1-cis_x_pro, cis_x_pro))

    with open(opf_x_y, 'w') as opf:
        opf.writelines('# {0:<20s}{1:<20s}{2:<20s}\n'.format(
                'replica_id', 'trans_x_y', 'cis_x_y'))
                
        opf.writelines('  {0:<20s}{1:<20.5f}{2:<20.5f}\n'.format(
                xvgf, 1-cis_x_y, cis_x_y))

if __name__ == "__main__":
    infiles = sys.argv[1:]
    for infile in infiles:
        main(infile)
