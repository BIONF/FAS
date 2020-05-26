#!/bin/env python

#######################################################################
# Copyright (C) 2019 Vinh Tran
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

import os
import sys
from sys import platform
from pathlib import Path
import re
import subprocess
import argparse
import readline
import glob
import inspect
import greedyFAS
from os.path import expanduser

home = expanduser("~")

def main():
    version = "1.1.0"
    parser = argparse.ArgumentParser(description="You are running annoFAS version " + str(version) + ".")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--fasta', help='Input sequence(s) in fasta format', action='store', default='',
                          required=True)
    required.add_argument('-o', '--outpath', help='Output directory', action='store', default='', required=True)
    required.add_argument('-n', '--name', help='Name of output annotation folder', action='store', default='',
                          required=True)
    required.add_argument('-t', '--toolPath', help='Path to annotation tools', action='store', default='',
                          required=True)
    optional.add_argument('-e', '--extract', help='Path to save the extracted annotation for input sequence',
                          action='store', default='')
    optional.add_argument('-r', '--redo', help='Re-annotation the sequence with flps|coils|seg|pfam|signalp|smart|tmhmm. '
                                         'Only one selection allowed!', action='store', default=0)
    optional.add_argument('-f', '--force', help='Force override annotations', action='store_true')
    optional.add_argument('-c', '--cores', help='number of cores', action='store', default='')
    optional.add_argument('--getToolPath', help='Get path to annotation tools', action='store_true')

    args = parser.parse_args()

    options = {
        'fasta': args.fasta,
        'path': args.path,
        'name': args.name,
        'extract': args.extract,
        'redo': args.redo,
        'force': args.force,
        'cores': args.cores,
    }
    # if args.getToolPath:
    #     perl_script = get_path() + '/annoFAS.pl'
    #     get_toolPath(perl_script)
    # else:
    #     call_annoFAS_perl(options)


if __name__ == '__main__':
    main()
