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
from greedyFAS.annoFAS import annoModules
from pkg_resources import get_distribution


def update_JsonFile(js_old, js_new, out_file, force):
    if os.path.exists(f'{js_old}.updated'):
        if not force:
            sys.exit(f'ERROR: Output file {js_old}.updated exists!')
    # Load each JSON file
    with open(js_new, 'r') as f1:
        fas_new = json.load(f1)
        print(f'Pairs in new file: {len(fas_new)}')
    with open(js_old, 'r') as f2:
        fas_old = json.load(f2)
        print(f'Pairs in old file: {len(fas_old)}')
        c = 0
        for id in list(fas_old):
            tmp = id.split('_')
            if not len(tmp) == 2:
                sys.exit(f'{js_old} seems not to be a valid FAS output!')
            if f'{tmp[0]}_{tmp[1]}' in fas_new or f'{tmp[1]}_{tmp[0]}' in fas_new:
                c += 1
                del fas_old[id]
        fas_old.update(fas_new)
        print(f'Updated pairs: {c}')
        print(f'Total pairs in final file: {len(fas_old)}')
    # Dump into a single JSON file
    if not out_file:
        out_file = f'{js_old}.updated'
    with open(out_file, "w") as f:
        json.dump(fas_old, f)
    return(out_file)


def main():
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.' +
                                     ' This script is used to update an old FAS json file by a new one!',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("--old", default='.', type=str, required=True,
                          help="Old FAS json file")
    required.add_argument("--new", default=None, type=str, required=True,
                          help="New FAS json file")
    optional.add_argument("--out", default=None, type=str,
                          help="Output file")
    optional.add_argument('--force', help='Force override output files', action='store_true')
    args = parser.parse_args()

    annoModules.checkFileExist(args.new)
    annoModules.checkFileExist(args.old)

    out_file = update_JsonFile(args.old, args.new, args.out, args.force)
    print(f'DONE! Final output saved in {out_file}')


if __name__ == '__main__':
    main()
