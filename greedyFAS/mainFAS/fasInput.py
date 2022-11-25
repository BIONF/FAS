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

import json


def featuretypes(path):
    """
    Input function,
    reads the tools/featuretypes input-file and stores information in option

    :param path: path to input-file
    :param option: dictionary that contains the main option variables of FAS
    :return: option
    """
    input_linearized = []
    input_normal = []
    ifile = open(path, "r")
    lines = ifile.readlines()
    mode = "NULL"
    for line in lines:
        tmp = line.rstrip("\n").lower()
        if tmp == "#linearized":
            mode = "lin"
        elif tmp == "#normal":
            mode = "nor"
        elif tmp == "#checked":
            mode = "ignore"
        elif mode == "NULL":
            raise Exception(path + " is not a valid input file")
        elif mode == "lin" and len(tmp) > 0:
            input_linearized.append(tmp)
        elif mode == "nor" and len(tmp) > 0:
            input_normal.append(tmp)
    ifile.close()
    return input_linearized, input_normal


def constraints_in(path):
    """
    Input function,
    reads the constraints file

    :param path: path to input-file
    :return: constraints
    """
    constraints = {}
    cfile = open(path, "r")
    lines = cfile.readlines()
    i = 1
    if lines[0][0] == "#":
        while lines[i][0] != "#":
            split = (lines[i].rstrip("\n")).split("\t")
            if split[1] != "N":
                constraints[split[0]] = float(split[1])
            i += 1
    else:
        raise Exception(path + " might be in the wrong format. Please see the sample file in config directory.")
    if lines[i][0] == "#":
        i += 1
        while i < len(lines):
            split = (lines[i].rstrip("\n")).split("\t")
            constraints[split[0]] = float(split[1])
            i += 1
    cfile.close()
    return constraints


def check_version(version, proteome, v_warning):
    if 'version' in proteome:
        if version and not proteome['version'] == version:
            v_warning = True
        elif not version:
            version = proteome['version']
    else:
        if version and not version == 'NA':
            v_warning = True
        elif not version:
            version = 'NA'
    return version, v_warning


def read_pairwise(path):
    with open(path, 'r') as infile:
        pairwise = []
        line = infile.readline()
        while line:
            pairwise.append(line.rstrip('\n').split('\t'))
            line = infile.readline()
    return pairwise


def read_json(path):
    with open(path, 'r') as infile:
        in_dict = json.loads(infile.read())
    return in_dict
