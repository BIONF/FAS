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


from greedyFAS.mainFAS import fasInput
import argparse
import os


def get_options():
    version = '1.5.3'
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.',
                                     epilog="This script allows you to create domain input files for phyloprofile "
                                            "without doing a FAS calculation.")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    parser.add_argument('--version', action='version', version=str(version))
    required.add_argument("-j", "--json", default=None, type=str, required=True,
                          help="path to protein architecture json file")
    required.add_argument("-p", "--protein_id", default=None, nargs='*', type=str,
                          help="Choose the protein ids for which you want to create a domain output (multiple ids "
                               "possible, divide with a space)")
    required.add_argument("-o", "--outfile", default=None, type=str,
                          help="path and name of outfile")
    optional.add_argument("-d", "--featuretypes", default=None, type=str,
                        help="inputfile that contains the tools/databases used to predict features. Please look at the "
                             "FAS wiki pages for templates of the the featuretypes input file")
    optional.add_argument("-n", "--groupname", default='Group', type=str,
                          help="Name of the protein group in the domain file")
    optional.add_argument("--toolPath", dest="toolPath", default=None, type=str,
                            help="Path to Annotion tool directory created with setupFAS")
    args = parser.parse_args()
    return args


def write_domain_file(path, idlist, outpath, tools, groupname):
    proteome = fasInput.read_json(path)["feature"]
    with open(outpath, 'w') as out:
        for pid in idlist:
            for tool in tools:
                for feature in proteome[pid][tool]:
                    for instance in proteome[pid][tool][feature]["instance"]:
                        out.write(groupname + "#" + pid + "\t" + pid + "\t" + str(proteome[pid]["length"]) +
                                  "\t" + feature + "\t" + str(instance[0]) + "\t" + str(instance[1]) + "\tNA\tN\n")


def main():
    args = get_options()
    toolpath = args.toolPath
    if toolpath is None:
        pathconfigfile = os.path.realpath(__file__).replace('domainFAS.py', 'pathconfig.txt')
        with open(pathconfigfile) as f:
            toolpath = f.readline().strip()
    if args.featuretypes is not None:
        option_dict = fasInput.featuretypes(args.featuretypes, {})
    else:
        option_dict = fasInput.featuretypes(toolpath + '/' + 'annoTools.txt', {})
    write_domain_file(args.json, args.protein_id, args.outfile, option_dict["input_linearized"] +
                      option_dict["input_normal"], args.groupname)


if __name__ == '__main__':
    main()
