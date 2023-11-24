# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2022 Vinh Tran
#
#  This file is part of FAS.
#  Update annotation file with InterPro Acc, phmm length and tool versions
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
import greedyFAS.annoFAS.annoModules as annoModules
from pkg_resources import get_distribution


def main():
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.' +
                                     ' This script os used to add interpro IDs and tool versions to an existing annotation file',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-a', '--annoFile', help='Input annotation file in json format', action='store', default='',
                          required=True)
    optional.add_argument('-o', '--outPath', help='Output directory', action='store', default='')
    optional.add_argument('--silent', help='Turn off terminal output', action='store_true')
    optional.add_argument('--annoToolFile', help='Path to files contains annotation tool names',
                          action='store', default='')

    args = parser.parse_args()

    annoFile = args.annoFile
    annoModules.checkFileExist(annoFile)
    annoFile = os.path.abspath(annoFile)
    outPath = os.path.abspath(args.outPath)
    Path(outPath).mkdir(parents=True, exist_ok=True)
    silent = args.silent
    annoToolFile = args.annoToolFile
    annoModules.checkFileExist(annoToolFile)

    try:
        annoModules.updateAnnoFile(annoFile)
        annoModules.printMsg(silent, '%s updated successfully!' % annoFile)
    except:
        print('ERROR: Cannot update %s! You can check with fas.checkAnno or contact us!' % annoFile)

if __name__ == '__main__':
    main()
