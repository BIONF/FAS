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


def write_tsv_out(outpath, bidirectional, results):
    out = open(outpath + ".tsv", "w")
    outdict = {}
    out.write("Seed\tQuery\tScore(Forward/Reverse)\tMS(Forward/Reverse)\tPS(Forward/Reverse)\tCS(Forward/Reverse)"
              "\tLS(Forward/Reverse)\tMethod\n")
    for result in results[0]:
        outdict[result[0], result[1]] = (result[2], ("NA", "NA", "NA", "NA", "NA"), result[3])
    if bidirectional:
        for result in results[1]:
            outdict[result[1], result[0]] = (outdict[result[1], result[0]][0], result[2], outdict[result[1],
                                                                                                  result[0]][2])
    for pair in outdict:
        out.write(pair[0] + "\t" + pair[1] + "\t" + f"{outdict[pair][0][0]:.4}" + "/" +
                  f"{outdict[pair][1][0]:.4}" + "\t" + f"{outdict[pair][0][1]:.4}" + "/" +
                  f"{outdict[pair][1][1]:.4}" + "\t" + f"{outdict[pair][0][2]:.4}" + "/" +
                  f"{outdict[pair][1][2]:.4}" + "\t" + f"{outdict[pair][0][3]:.4}" + "/" +
                  f"{outdict[pair][1][3]:.4}" + "\t" + f"{outdict[pair][0][4]:.4}" + "/" +
                  f"{outdict[pair][1][4]:.4}" + "\t" + outdict[pair][2] + "\n")
    out.close()


def write_domain_out(seed_proteome, query_proteome, seed, query, weights, scale, seedpath, querypath, out, option):
    tools = option["input_linearized"] + option["input_normal"]
    uni_weight = None
    if option['reverse']:
        groupname = query
    else:
        groupname = seed
    if option["MS_uni"] and len(seedpath) > 0:
        uni_weight = round(1.0 / len(seedpath), 4)
    for tool in tools:
        for feature in seed_proteome[seed][tool]:
            if feature in seedpath:
                if option["MS_uni"]:
                    weight = uni_weight
                else:
                    weight = round(weights[feature] * scale, 4)
                for instance in seed_proteome[seed][tool][feature]["instance"]:
                    if (instance[0], instance[1]) in seedpath[feature]:
                        inpath = "Y"
                    else:
                        inpath = "N"
                    if option["reverse"]:
                        out.write(groupname + "#" + seed + "\t" + seed + "\t" + str(seed_proteome[seed]["length"]) +
                                  "\t" + feature + "\t" + str(instance[0]) + "\t" + str(instance[1]) + "\t" +
                                  str(weight) + "\t" + inpath + "\n")
                    else:
                        out.write(groupname + "#" + query + "\t" + seed + "\t" + str(seed_proteome[seed]["length"]) +
                                  "\t" + feature + "\t" + str(instance[0]) + "\t" + str(instance[1]) + "\t" +
                                  str(weight) + "\t" + inpath + "\n")
            else:
                for instance in seed_proteome[seed][tool][feature]["instance"]:
                    if option["reverse"]:
                        out.write(groupname + "#" + seed + "\t" + seed + "\t" + str(seed_proteome[seed]["length"]) +
                                  "\t" + feature + "\t" + str(instance[0]) + "\t" + str(instance[1]) + "\tNA\tN\n")
                    else:
                        out.write(groupname + "#" + query + "\t" + seed + "\t" + str(seed_proteome[seed]["length"]) +
                                  "\t" + feature + "\t" + str(instance[0]) + "\t" + str(instance[1]) + "\tNA\tN\n")
    if not seed == query:
        for tool in tools:
            for feature in query_proteome[query][tool]:
                if feature in seedpath:
                    if option["MS_uni"]:
                        weight = uni_weight
                    else:
                        weight = round(weights[feature] * scale, 4)
                else:
                    weight = "NA"
                if feature in querypath:
                    for instance in query_proteome[query][tool][feature]["instance"]:
                        if (instance[0], instance[1]) in querypath[feature]:
                            inpath = "Y"
                        else:
                            inpath = "N"
                        if option["reverse"]:
                            out.write(groupname + "#" + seed + "\t" + query + "\t" +
                                      str(query_proteome[query]["length"]) + "\t" + feature + "\t" + str(instance[0]) +
                                      "\t" + str(instance[1]) + "\t" + str(weight) + "\t" + inpath + "\n")
                        else:
                            out.write(groupname + "#" + query + "\t" + query + "\t" +
                                      str(query_proteome[query]["length"]) + "\t" + feature + "\t" + str(instance[0]) +
                                      "\t" + str(instance[1]) + "\t" + str(weight) + "\t" + inpath + "\n")
                else:
                    for instance in query_proteome[query][tool][feature]["instance"]:
                        if option["reverse"]:
                            out.write(groupname + "#" + seed + "\t" + query + "\t" +
                                      str(query_proteome[query]["length"]) + "\t" + feature + "\t" + str(instance[0]) +
                                      "\t" + str(instance[1]) + "\t" + str(weight) + "\tN\n")
                        else:
                            out.write(groupname + "#" + query + "\t" + query + "\t" +
                                      str(query_proteome[query]["length"]) + "\t" + feature + "\t" + str(instance[0]) +
                                      "\t" + str(instance[1]) + "\t" + str(weight) + "\tN\n")


def phyloprofile_out(outpath, bidirectional, mapping_file, results):
    with open(mapping_file) as infile:
        map = {}
        for line in infile.readlines():
            cells = line.rstrip("\n").split("\t")
            map[cells[0]] = cells[1]
    out = open(outpath + ".phyloprofile", "w")
    out.write("geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n")
    outdict = {}

    for result in results[0]:
        outdict[result[0], result[1]] = (result[2][0], 0.0)
    if bidirectional:
        for result in results[1]:
            outdict[result[1], result[0]] = (outdict[result[1], result[0]][0], result[2][0])
    for pair in outdict:
        try:
            out.write(pair[0] + "\t" + map[pair[1]] + "\t" + pair[1] + "\t" + str(outdict[pair][0]) + "\t" +
                      str(outdict[pair][1]) + "\n")
        except KeyError:
            raise Exception(pair[1] + " not in mapping file")
    out.close()
