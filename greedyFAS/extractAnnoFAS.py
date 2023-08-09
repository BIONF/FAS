#!/bin/env python

#######################################################################
# Copyright (C) 2020 Julian Dosch
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
import sys
import json
from tqdm import tqdm
from time import sleep
from greedyFAS.mainFAS.fasInput import read_json
from pkg_resources import get_distribution

def get_options():
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    parser.add_argument('--version', action='version', version=str(version))
    required.add_argument("-x", "--extract_ids", default=None, type=str, required=True,
                          help="path to file, which contains a list of ids that should be extracted")
    required.add_argument("-i", "--input", default=None, type=str, required=True,
                          help="path to input json")
    required.add_argument("-o", "--output", default=None, type=str, required=True,
                          help="path to output json")
    arguments = parser.parse_args()
    return arguments


def get_ids(path):
    ids = []
    with open(path, 'r') as infile:
        for line in infile.readlines():
            pid = line.rstrip('\n').lstrip('>')
            ids.append(pid)
    return ids


def extract_architectures(ids, annotation):
    clan = annotation['clan']
    interprokeys = {}
    if 'interproID' in annotation:
        interprokeys = annotation['interproID']
    phmm = {}
    if 'length' in annotation:
        phmm = annotation['length']
    feature = {}
    count = {}
    missing = []
    progress = tqdm(total=len(ids), file=sys.stdout)
    for pid in ids:
        if pid in annotation['feature']:
            feature[pid] = annotation['feature'][pid]
            for tool in feature[pid]:
                if not tool == 'length':
                    for fid in feature[pid][tool]:
                        if fid in count:
                            count[fid] += len(feature[pid][tool][fid]['instance'])
                        else:
                            count[fid] = len(feature[pid][tool][fid]['instance'])
        else:
            missing.append(pid)
        progress.update(1)
    progress.refresh()
    progress.close()
    sleep(1.0)
    new_dict = {'feature': feature, 'clan': clan, 'count': count, 'length': phmm, 'interproID': interprokeys}
    if len(missing) > 0:
        print('\nWARNING: The following architectures could not be extracted:')
        for pid in missing:
            print(pid)
    return new_dict


def savejson(dict2save, outpath):
    jsonOut = json.dumps(dict2save, ensure_ascii=False)
    f = open(outpath, 'w')
    f.write(jsonOut)
    f.close()


def main():
    options = get_options()
    ids = get_ids(options.extract_ids)
    annotation = read_json(options.input)
    new_dict = extract_architectures(ids, annotation)
    savejson(new_dict, options.output)


if __name__ == '__main__':
    main()
