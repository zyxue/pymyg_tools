#! /usr/bin/env python

import os
import sys
import time
import argparse

def main_load():
    print "#" * 79
    print "{0:^79s}".format("SUCCESSFULLY LOAD main()! by Zhuyi Xue (zhuyi.xue@mail.utoronto.ca)")
    print "#" * 79

def write_header(outputfile):
    beginning_time = time.time()
    outputfile.write('# main calculations BEGIN at: {0:s}\n\n'.format(time.ctime()))
    outputfile.write('# {0:s}\n\n'.format(' '.join(sys.argv[:])))
    outputfile.write('{0:s}\n\n'.format('#' * 50))
    return beginning_time                                             # for calculating the time used

def write_footer(outputfile, beginning_time):
    ending_time = time.time()
    delta_time = ending_time - beginning_time
    outputfile.write('\n{0:s}\n\n'.format('#' * 50))
    outputfile.write('# main calculations END at: {0:s}\n'.format(time.ctime()))
    outputfile.write('# time consumed in TOTOAL: {0}\n'.format(
            time.strftime('%H:%M:%S', time.gmtime(delta_time))))

def backup(output):
    if os.path.exists(output):
        dirname = os.path.dirname(output)
        basename = os.path.basename(output)
        count = 1
        newbasename =  '#' + basename + '.{0}#'.format(count)
        rn_to = os.path.join(dirname, newbasename)                 # rename to
        while os.path.exists(rn_to):
            count += 1
            newbasename =  '#' + basename + '.{0}#'.format(count)
            rn_to = os.path.join(dirname, newbasename)                 # rename to
        os.rename(output, rn_to)

def print_progress(ts):
    if ts.frame % 2 == 0:
        sys.stdout.write("\r[38;5;226mtime: {0:10.0f}; step: {1:10d}; frame: {2:10d}".format(
                ts.time, ts.step, ts.frame))
        sys.stdout.flush()

def debug_print(ts, debug_flag):
    # for debugging only, verbose print
    if debug_flag and ts.frame % 2 == 0:
        print 'time: {0:10.0f}; step: {1:10d}; frame: {2:10d}'.format(ts.time, ts.step, ts.frame)

def get_basic_parser(usage=None):
    p = basic_parser = argparse.ArgumentParser(usage=usage)
    p.add_argument('-f', '--xtcf', type=str, default=None,
                   help='Trajectory: xtc')
    p.add_argument('-s', '--grof', type=str, default=None,
                   help='Structure: gro')
    p.add_argument('-b', '--btime', type=float, default=0,
                   help='beginning time in ps (confirmed) xtc file records time in unit of ps')
    p.add_argument('-e', '--etime', type=float, default=0,
                   help='ending time in ps (confirmed) xtc file records time in unit of ps')
    p.add_argument('-o', '--optf', type=str, default=None,
                   help='Output file name (OPT.)')
    return basic_parser

def swap_aa_name(abbr):
    map3 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
        "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
        "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
        "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        }
    map1 = {
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP",
        "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY",
        "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
        "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",
        "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
        }
    if len(abbr) == 3:
        return map3[abbr]
    elif len(abbr) == 1:
        return map1[abbr]
        
