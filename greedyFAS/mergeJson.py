#!/bin/env python

#######################################################################
# Copyright (C) 2023 Vinh Tran
#
# This file is part of FAS.
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


def merge_JsonFiles(json_files, outName, outPath):
    out_dict = {}
    # Load each JSON file
    for json_file in json_files:
        with open(json_file, "r") as f:
            dd = json.load(f)
            for id in list(dd):
                tmp = id.split('_')
                if not len(tmp) == 2:
                    sys.exit(f'{json_file} seems not to be a valid FAS output!')
                if f'{tmp[0]}_{tmp[1]}' in out_dict or f'{tmp[1]}_{tmp[0]}' in out_dict:
                    del dd[id]
            out_dict.update(dd)
    # Dump into a single JSON file
    final_out = f'{os.path.abspath(outPath)}/{outName}.json'
    if not os.path.exists(final_out):
        with open(final_out, "w") as f:
            json.dump(out_dict, f)
    else:
        sys.exit(f'Output file {final_out} exists!')


def main():
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.' +
                                     ' This script is used to merge multiple calculated FAS json output files!' +
                                     ' Note: if one protein pair has different scores in different files, only' +
                                     ' the scores in the first file will be kept. If you want to update the FAS scores,' +
                                     ' please use fas.updateJson function!',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-i", "--input", default='.', type=str, required=True,
                          help="Path to input folder contatining json files")
    required.add_argument("-n", "--outName", default=None, type=str, required=True,
                          help="Name of the output json")
    optional.add_argument("-o", "--outPath", default='.', type=str, required=False,
                          help="Path to output directory")
    args = parser.parse_args()

    path = f'{os.path.abspath(args.input)}/*.json'
    files = glob.glob(path)
    merge_JsonFiles(files, args.outName, args.outPath)
    print(f'DONE! Final output saved in {os.path.abspath(args.outPath)}/{args.outName}.json')


if __name__ == '__main__':
    main()
