#!/bin/env python

#######################################################################
# Copyright (C) 2023 Vinh Tran
#
#  This file is part of FAS. This script is used to convert FAS output
#  from TSV (tab-delimited) to JSON formar. Only proteins IDs and their
#  (pairwise) FAS scores will be stored in the JSON output.
#
#  FAS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FAS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with FAS.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################


import os
import sys
import argparse
import json
import glob
from pkg_resources import get_distribution


def read_file(file):
    """ Read a file and return list of lines"""
    if os.path.exists(file):
        with open(file, 'r') as f:
            lines = f.read().splitlines()
            f.close()
            return(lines)
    else:
        sys.exit('%s not found' % file)


def tsv_to_json(inFile, outName, outPath):
    """ Convert .tsv output to .json output
    """
    json_dict = {}
    tsv_file = read_file(inFile)
    for line in tsv_file:
        if not line.startswith('Seed'):
            tmp = line.split()
            json_dict[f'{tmp[0]}_{tmp[1]}'] = [f'{tmp[2].split("/")[0]}', f'{tmp[2].split("/")[1]}']

    final_out = f'{os.path.abspath(outPath)}/{outName}.json'
    if not os.path.exists(final_out):
        with open(final_out, "w") as f:
            json.dump(json_dict, f)
    else:
        sys.exit(f'Output file {final_out} exists!')


def main():
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-i", "--input", default='.', type=str, required=True,
                          help="Input file in TSV format")
    optional.add_argument("-n", "--outName", default=None, type=str, required=False,
                          help="Name of the output json")
    optional.add_argument("-o", "--outPath", default='.', type=str, required=False,
                          help="Path to output directory")
    args = parser.parse_args()

    if not args.outName:
        outName = args.input.split('/')[-1]
    else:
        outName = args.outName
    tsv_to_json(args.input, outName, args.outPath)
    print(f'DONE! Final output saved in {os.path.abspath(args.outPath)}/{outName}.json')


if __name__ == '__main__':
    main()
