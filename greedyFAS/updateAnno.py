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
from pkg_resources import get_distribution
from greedyFAS.annoFAS import annoModules


def update_annoFile(js_old, js_new, force):
    if os.path.exists(f'{js_old}.updated'):
        if not force:
            sys.exit(f'ERROR: Output file {js_old}.updated exists!')
    # Load each JSON file
    with open(js_new, 'r') as f1:
        anno_new = json.load(f1)
    with open(js_old, 'r') as f2:
        anno_old = json.load(f2)
        for i in anno_old:
            anno_old[i].update(anno_new[i])
    # Dump into a single JSON file
    with open(f'{js_old}.updated', "w") as f:
        json.dump(anno_old, f)


def main():
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.' +
                                     ' This script is used to update an old annotation json file by a new one!',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("--old", default='.', type=str, required=True,
                          help="Old annotation json file")
    required.add_argument("--new", default=None, type=str, required=True,
                          help="New annotation json file")
    optional.add_argument('--force', help='Force override output files', action='store_true')
    args = parser.parse_args()

    annoModules.checkFileExist(args.new)
    annoModules.checkFileExist(args.old)

    update_annoFile(args.old, args.new, args.force)
    print(f'DONE! Final output saved in {os.path.abspath(args.old)}.updated')


if __name__ == '__main__':
    main()
