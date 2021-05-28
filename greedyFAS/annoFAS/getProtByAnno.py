# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2021 Vinh Tran
#
#  This file is part of FAS.
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

import sys
import os
import argparse
from pathlib import Path
import json
import re

def checkCompleteAnno(featureList, jsonFile):
    with open(jsonFile) as jf:
        dt = json.load(jf)
        # get list of proteins that contain selected features
        protList = []
        for prot in list(dt['feature'].keys()):
            for item in dt['feature'][prot]:
                if isinstance(dt['feature'][prot][item], dict):
                    if any(re.search(feat+';', ';'.join(list(dt['feature'][prot][item].keys()))+';', re.I)
                        for feat in featureList):
                        protList.append(prot)
        # get all features for found proteins
        out = {}
        for prot in protList:
            out[prot] = []
            for item in dt['feature'][prot]:
                if isinstance(dt['feature'][prot][item], dict):
                    features = list(dt['feature'][prot][item].keys())
                    if len(features) > 0:
                        out[prot].append(';'.join(features))
        return(out)

def main():
    version = '0.0.1'
    parser = argparse.ArgumentParser(description='You are running getProtByAnno version ' + str(version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-a', '--annoFile', help='Input annotation file in json format', action='store', default='',
                          required=True)
    required.add_argument('-f', '--features', help='List of features of interest, separated by comma. E.g. pfam_ig,pfam_TMA7', action='store', default='',
                          required=True)
    optional.add_argument('-i', '--idOnly', help='Get only protein IDs', action='store_true')

    args = parser.parse_args()
    featureList = str(args.features).split(",")
    if len(featureList) == 0:
        sys.exit('No feature given! Please specify feature of interest using --features option!')
    jsonFile = os.path.abspath(args.annoFile)
    if not os.path.exists(jsonFile):
        sys.exit('%s not found!' % jsonFile)
    idOnly = args.idOnly

    out = checkCompleteAnno(featureList, jsonFile)
    if not idOnly:
        for prot in out:
            print('%s\t%s' % (prot, ';'.join(out[prot])))
    else:
        print('\n'.join(out))

if __name__ == '__main__':
    main()
