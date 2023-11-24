#!/bin/env python

#######################################################################
# Copyright (C) 2019 Julian Dosch
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
from greedyFAS.annoFAS.annoModules import mergeNestedDic
from greedyFAS.mainFAS.fasInput import read_json
from greedyFAS.annoFAS.annoModules import save2json
from pkg_resources import get_distribution


def merge_anno(pathlist, outpath, name):
    feature = []
    main = []
    for path in pathlist:
        input = read_json(path)
        feature.append(input['feature'])
        main.append({'clan': input['clan'], 'count': input['count']})
    mergedfeature = mergeNestedDic(feature)
    mergedmain = mergeNestedDic(main)
    mergedmain['feature'] = mergedfeature
    save2json(mergedmain, name, outpath)


def main():
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.' +
                                     ' This script is used to merge all annotation json files into one file!' +
                                     ' Note: if one protein has different annotations in different files, only' +
                                     ' the first annotation will be kept. If you want to update the annotations,' +
                                     ' please use fas.updateAnno function!',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-i", "--input", nargs='*', default=None, type=str, required=True,
                          help="path to input jsons, seperated by a space")
    required.add_argument("-o", "--outPath", default='.', type=str, required=True,
                          help="path to output directory.")
    optional.add_argument("-n", "--outName", default='merged', type=str, required=False,
                          help="Name of the output json.")
    args = parser.parse_args()
    merge_anno(args.input, args.outPath, args.outName)


if __name__ == '__main__':
    main()
