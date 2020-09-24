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

import xml.etree.ElementTree as ElTre
import os
import json


def xmlreader(path, mode, tool, assess, proteome, protein_lengths, clan_dict, option):
    """
    Input function,
    read input-files for seed, query and reference

    :param path: path to input file
    :param mode: 0 (seed), 1 (query), 2 (reference)
    :param tool: tool/database name (pfam, seg, ...)
    :param assess: e-value assess
    :param proteome: proteome var to write in [dict]
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param clan_dict: dictionary that maps features to clans
    :param option: dictionary that contains the main option variables of FAS
    :return: proteome, protein_lengths, clan_dict
    """
    start = 0
    if os.path.exists(path):
        try:
            xmltree = ElTre.parse(path)
        except ElTre.ParseError:
            raise Exception(path + ' is not a legit xml file')
        root = xmltree.getroot()
        for protein in root:
            p_id = protein.attrib["id"]
            plength = protein.attrib["length"]
            # set up of protein IDs, differentiate between proteins from different files
            if mode == 1:
                if(not option["query_id"]) or p_id in option["query_id"]:
                    protein_lengths["query_" + str(p_id)] = float(plength)
            elif mode == 0:
                if (not option["seed_id"]) or p_id in option["seed_id"]:
                    protein_lengths["seed_" + str(p_id)] = float(plength)
            elif mode == 2:
                protein_lengths[p_id] = float(plength)
            # set up of datastructure to store annotations
            if mode == 2 or (mode == 1 and ((not option["query_id"]) or p_id in option["query_id"])) or (
               mode == 0 and ((not option["seed_id"]) or p_id in option["seed_id"])):
                if not (p_id in proteome):
                    proteome[p_id] = {}

                for feature in protein:
                    if len(feature) > 0:
                        ftype = tool + "_" + feature.attrib["type"]
                        feat_eval = 'NULL'

                        # evalue check: family/ftype based
                        if 'evalue' in feature.attrib and float(feature.attrib["evalue"]) > option["efilter"]:
                            # skip current ftype and continue with the next one
                            continue
                        else:
                            # keep feature type bases evalue
                            if 'evalue' in feature.attrib:
                                feat_eval = float(feature.attrib["evalue"])

                            if 'clan' in feature.attrib:
                                fclan = feature.attrib["clan"]
                            else:
                                fclan = "NO_CLAN"
                            proteome[p_id][ftype] = []
                            proteome[p_id][ftype].append(assess)
                            proteome[p_id][ftype].append(feat_eval)

                            i = 0
                            # counting appended instances
                            inst_count = 0
                            for instance in feature:
                                inst_eval = 'NULL'
                                # XMLcase 1: feature instance contains evalue information (XML field inst_eval)
                                if 'inst_eval' in instance.attrib:
                                    # print tool + " instance evalue: "+ str(instance.attrib)
                                    inst_eval = float(instance.attrib["inst_eval"])
                                    start = int(instance.attrib["start"])
                                    end = int(instance.attrib["end"])

                                    if inst_eval <= option["inst_efilter"]:
                                        proteome[p_id][ftype].append((inst_eval, start, end))
                                        inst_count += 1

                                # XMLcase 2: feature instance contains NO evalue information (XML field inst_eval)
                                else:
                                    # NO instance based evalue information --> no instances can be rejected:
                                    # set inst_count = 1
                                    inst_count = 1
                                    if len(instance.attrib) == 2:
                                        start = int(instance.attrib["start"])
                                        end = int(instance.attrib["end"])
                                        proteome[p_id][ftype].append((inst_eval, start, end))

                                    else:
                                        if i == 0:
                                            start = int(instance.attrib["start"])
                                            i = 1
                                        else:
                                            end = int(instance.attrib["end"])
                                            proteome[p_id][ftype].append((inst_eval, start, end))
                                            i = 0
                            # any instance appended?
                            if inst_count < 1:
                                # delete feature type
                                proteome[p_id].pop(ftype)
                            if ftype not in clan_dict:
                                clan_dict[ftype] = fclan
    else:
        raise Exception(path + " does not exist")
    return proteome, protein_lengths, clan_dict


def featuretypes(path, option):
    """
    Input function,
    reads the tools/featuretypes input-file and stores information in option

    :param path: path to input-file
    :param option: dictionary that contains the main option variables of FAS
    :return: option
    """
    option["input_linearized"] = []
    option["input_normal"] = []
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
            option["input_linearized"].append(tmp)
        elif mode == "nor" and len(tmp) > 0:
            option["input_normal"].append(tmp)
    ifile.close()
    return option


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
            split = (lines[i].rstrip("\n")).split(" ")
            if split[1] != "N":
                constraints[split[0]] = float(split[1])
            i += 1
    else:
        raise Exception(path + " might be in the wrong format. Please see the sample file in config directory.")
    if lines[i][0] == "#":
        i += 1
        while i < len(lines):
            split = (lines[i].rstrip("\n")).split(" ")
            constraints[split[0]] = float(split[1])
            i += 1
    cfile.close()
    return constraints


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
