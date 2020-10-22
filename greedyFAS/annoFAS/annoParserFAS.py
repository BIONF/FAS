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
import json


def read_infile_single(ignore_lines, feature_columns, tool_name, path):
    infile = open(path, 'r')
    for i in range(ignore_lines):
        infile.readline()
    line = infile.readline()
    feature_dict = {}
    while line:
        cells = line.rstrip('\n').split('\t')
        p_id, length = cells[feature_columns[0]], int(cells[feature_columns[1]])
        start, stop = int(cells[feature_columns[3]]), int(cells[feature_columns[4]])
        f_id = tool_name + '_' + cells[feature_columns[2]]
        if p_id not in feature_dict:
            feature_dict[p_id] = {'length': length, tool_name: {}}
        if f_id not in feature_dict[p_id][tool_name]:
            feature_dict[p_id][tool_name][f_id] = {'evalue': 'NA', 'instance': []}
        feature_dict[p_id][tool_name][f_id]['instance'].append([start, stop, 'NA'])
        line = infile.readline()
    return feature_dict


def read_infile_multiple(ignore_lines, feature_columns, tool_names, tool_column, path):
    infile = open(path, 'r')
    for i in range(ignore_lines):
        infile.readline()
    line = infile.readline()
    feature_dict = {}
    while line:
        cells = line.rstrip('\n').split('\t')
        p_id, length = cells[feature_columns[0]], int(cells[feature_columns[1]])
        start, stop = int(cells[feature_columns[3]]), int(cells[feature_columns[4]])
        tool_name = cells[tool_column].lower()
        f_id = tool_name + '_' + cells[feature_columns[2]]
        if tool_name not in tool_names:
            raise Exception('tsv file contains the tool: "' + tool_name + '", which was not given in '
                                                                          '[-t] [--tool_names]')
        if p_id not in feature_dict:
            feature_dict[p_id] = {'length': length}
            for tool in tool_names:
                feature_dict[p_id][tool.lower()] = {}
        if f_id not in feature_dict[p_id][tool_name]:
            feature_dict[p_id][tool_name][f_id] = {'evalue': 'NA', 'instance': []}
        feature_dict[p_id][tool_name][f_id]['instance'].append([start, stop, 'NA'])
        line = infile.readline()
    return feature_dict


def count_features(feature_dict):
    count = {}
    for p_id in feature_dict:
        for tool in feature_dict[p_id]:
            if not tool == 'length':
                for f_id in feature_dict[p_id][tool]:
                    if f_id not in count:
                        count[f_id] = 0
                    count[f_id] += len(feature_dict[p_id][tool][f_id]['instance'])
    return count


def write_json(dict2save, path):
    jsonOut = json.dumps(dict2save, ensure_ascii=False)
    f = open(path, 'w')
    f.write(jsonOut)
    f.close()


def create_json(args):
    if len(args.tool_names) > 1 and args.tool_column:
        tool_names = []
        for i in args.tool_names:
            tool_names.append(str(i).lower())
        feature_dict = read_infile_multiple(args.ignore_lines, args.feature_columns, tool_names,
                                            args.tool_column, args.input)
    elif len(args.tool_names) > 1:
        raise Exception("Multiple tool names given. Please select the column containing tool names with "
                        "[-c][--tool_column]")
    else:
        feature_dict = read_infile_single(args.ignore_lines, args.feature_columns, args.tool_names[0].lower(),
                                          args.input)
    count = count_features(feature_dict)
    out_dict = {'feature': feature_dict, 'clan': {}, 'count': count}
    write_json(out_dict, args.output)


def main():
    parser = argparse.ArgumentParser(description="Parses .tsv formatted files into the json format used by FAS")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-i", "--input", type=str, default=None, required=True,
                          help="Path to input file, needs to be in tsv format")
    required.add_argument("-o", "--output", type=str, default=None, required=True,
                          help="Path to output file, should be a .json file")
    optional.add_argument("--ignore_lines", type=int, help="skip the first n lines (for example headers)", default=0)
    required.add_argument("-t", "--tool_names", type=str, nargs='+', required=True,
                          help="name of the annotation tool(s) in the tsv file, if multiple are given, the argument "
                               "--tool_name must point to the column containing the tool name and the tool names must "
                               "be identical to the ones in the tsv file")
    required.add_argument("-f", "--feature_columns", type=int, nargs=5, required=True,
                          help="needs 5 integer values that point to the columns that contain: "
                               "(1) the protein id, (2) the protein length, (3) the feature id, "
                               "(4) the start position of the feature, (5) the end position of the feature; the column "
                               "indices start at 0 so the first column has index 0, the second index 1, etc.")
    optional.add_argument("-c", "--tool_column", type=int, default=None,
                          help="The index (integer) of column that contains the annotation tool name, only necessary "
                               "if the tsv file contains multiple annotation tools like the tsv output of InterPro")
    args = parser.parse_args()
    create_json(args)
