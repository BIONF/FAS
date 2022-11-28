#!/bin/env python

#######################################################################
# Copyright (C) 2022 Julian Dosch
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


import argparse
from pkg_resources import get_distribution
from greedyFAS.mainFAS.fasInput import read_json


def get_anno_version(path):
    anno = read_json(path)
    if 'version' in anno:
        for tool in anno['version']:
            print(tool)
            for x in anno['version'][tool]:
                print(x + ': ' + str(anno['version'][tool][x]))
            print('-----------------------------------')
    else:
        print('This seems to be a pre 1.15 annotation file. Specific version data is only for annotations with FAS version '
              '1.15 and higher.')


def main():
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", default=None, type=str, required=True,
                          help="path to annotation json")
    args = parser.parse_args()
    get_anno_version(args.input)


if __name__ == '__main__':
    main()
