#!/bin/env python

#######################################################################
# Copyright (C) 2024 Vinh Tran
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
from pathlib import Path
import sys
import argparse
from pkg_resources import get_distribution


def split_domainFile(domain_file, outdir, idlist):
    # get list of input protein ids (if given):
    ids = ''
    if os.path.exists(idlist):
        with open(idlist, 'r') as idf:
            ids = idf.read().splitlines()
    # parse domain file
    domain_dict = {}
    header = ''
    if os.path.exists(domain_file):
        with open(domain_file, 'r') as df:
            for l in df.read().splitlines():
                if l.startswith('#'):
                    header = l
                else:
                    id = l.split('#')[0]
                    if (len(ids) > 0 and id in ids) or len(ids) == 0:
                        if not id in domain_dict:
                            domain_dict[id] = [l]
                        else:
                            domain_dict[id].append(l)
    else:
        sys.exit(f'ERROR: {domain_file} not found!')
    # write single domain files to outdir
    if len(domain_dict) > 0:
        if not os.path.exists(os.path.abspath(outdir)):
            Path(os.path.abspath(outdir)).mkdir(parents = True, exist_ok = True)
        for id in domain_dict:
            with open(f'{outdir}/{id}.domains', 'w') as o:
                o.write('%s\n%s' % (header, '\n'.join(domain_dict[id])))


def main():
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.' +
                                     ' This script is used to split a domain file into individual file for each protein',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-i", "--input", default='.', type=str, required=True,
                          help="Input domain file")
    required.add_argument("-o", "--outdir", default='.', type=str, required=True,
                          help="Output directory")
    optional.add_argument("--idlist", default='', type=str, required=False,
                          help="File containing list of protein IDs that need to be splitted")
    args = parser.parse_args()


    domain_dict = split_domainFile(args.input, args.outdir, args.idlist)
    print(f'DONE! Output files saved in {os.path.abspath(args.outdir)}')


if __name__ == '__main__':
    main()
