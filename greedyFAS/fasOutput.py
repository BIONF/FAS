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


def bidirectionout(outpath):
    """Output function for bidirectional mode: This function summarizes the scores of the both scoring directions into a
    table (csv format).

    :param outpath: path to the out directory
    :return: no returns
    """
    outdict = {}
    forwardtree = ElTre.parse(outpath + ".xml")
    reversetree = ElTre.parse(outpath + "_reverse.xml")
    forwardroot = forwardtree.getroot()
    reverseroot = reversetree.getroot()
    for query in forwardroot:
        query_id = query.attrib["id"]
        for seed in query:
            seed_id, forward_score, seed_mode = seed.attrib["id"], seed.attrib["score"], seed.attrib["mode"]
            forward_s_path = {}
            forward_q_path = {}
            for path in seed:
                if path.tag == "template_path":
                    for feature in path:
                        forward_s_path[feature.attrib["type"]] = []
                        for instance in feature:
                            forward_s_path[feature.attrib["type"]].append((instance.attrib["start"],
                                                                           instance.attrib["end"]))
                if path.tag == "query_path":
                    for feature in path:
                        forward_q_path[feature.attrib["type"]] = []
                        for instance in feature:
                            forward_q_path[feature.attrib["type"]].append((instance.attrib["start"],
                                                                           instance.attrib["end"]))
            for node in reverseroot:
                if node.attrib["id"] == seed_id:
                    for child in node:
                        if child.attrib["id"] == query_id:
                            reverse_score, query_mode = child.attrib["score"], child.attrib["mode"]
                            reverse_s_path = {}
                            reverse_q_path = {}
                            for path in child:
                                if path.tag == "query_path":
                                    for feature in path:
                                        reverse_s_path[feature.attrib["type"]] = []
                                        for instance in feature:
                                            reverse_s_path[feature.attrib["type"]].append((instance.attrib["start"],
                                                                                           instance.attrib["end"]))
                                if path.tag == "template_path":
                                    for feature in path:
                                        reverse_q_path[feature.attrib["type"]] = []
                                        for instance in feature:
                                            reverse_q_path[feature.attrib["type"]].append((instance.attrib["start"],
                                                                                           instance.attrib["end"]))
                            consistence = "yes"
                            for feature in forward_q_path:
                                if feature in reverse_q_path:
                                    for instance in forward_q_path[feature]:
                                        if instance not in reverse_q_path[feature]:
                                            consistence = "no"
                                else:
                                    consistence = "no"
                            for feature in forward_s_path:
                                if feature in reverse_s_path:
                                    for instance in forward_s_path[feature]:
                                        if instance not in reverse_s_path[feature]:
                                            consistence = "no"
                                else:
                                    consistence = "no"
                            outdict[(seed_id, query_id)] = (forward_score, reverse_score, consistence,
                                                            seed_mode + "/" + query_mode)
    out = open(outpath + "_table.csv", "w")
    out.write("seedID,queryID,forward,reverse,Path_consistency,linearization_mode\n")
    for pair in outdict:
        out.write(pair[0] + "," + pair[1] + "," + outdict[pair][0] + "," + outdict[pair][1] + "," +
                  outdict[pair][2] + "," + outdict[pair][3] + "\n")
    out.close()


def domain_out(outpath, bidirectional, extendedout):
    arc = {}
    if extendedout:
        arctree = ElTre.parse(outpath + "_architecture.xml")
        arcroot = arctree.getroot()
        for architecture in arcroot:
            pid = architecture.attrib["id"]
            arc[pid] = {}
            for feature in architecture[0]:
                type = feature.attrib["type"]
                arc[pid][type] = []
                for inst in feature:
                    arc[pid][type].append((inst.attrib["start"], inst.attrib["end"]))
    # outdict = {}
    # groupname = outpath.split("/")[-1]
    d0_out = open(outpath + "_forward.domains", "w")
    forwardtree = ElTre.parse(outpath + ".xml")
    forwardroot = forwardtree.getroot()

    for query in forwardroot:
        query_id, query_length = query.attrib["id"], query.attrib["length"]
        for seed in query:
            seed_id, forward_score, seed_length = seed.attrib["id"], seed.attrib["score"], seed.attrib["length"]
            # outdict[(seed_id, query_id)] = (forward_score, "0.0")
            if extendedout:
                # weights = {}
                forward_s_path = {}
                forward_q_path = {}
                for path in seed:
                    if path.tag == "template_path":
                        for feature in path:
                            # if noref:
                            #     weights[feature.attrib["type"]] = str(1.0 / len(path))
                            # else:
                            #     weights[feature.attrib["type"]] = feature.attrib["corrected_weight"]
                            forward_s_path[feature.attrib["type"]] = []
                            for instance in feature:
                                forward_s_path[feature.attrib["type"]].append((instance.attrib["start"],
                                                                               instance.attrib["end"]))
                        for feature in arc[seed_id]:
                            if feature in forward_s_path:
                                for inst in arc[seed_id][feature]:
                                    if inst in forward_s_path[feature]:
                                        d0_out.write(seed_id + "#" + query_id + "\t" + seed_id + "\t" + seed_length +
                                                     "\t" + feature + "\t" + inst[0] + "\t" + inst[1] + "\t" +
                                                     "NA\tY\n")  # weights[feature] + "\tY\n")
                                    else:
                                        d0_out.write(seed_id + "#" + query_id + "\t" + seed_id + "\t" + seed_length +
                                                     "\t" + feature + "\t" + inst[0] + "\t" + inst[1] + "\t" +
                                                     "NA\tN\n")  # weights[feature] + "\tN\n")
                            else:
                                for inst in arc[seed_id][feature]:
                                    d0_out.write(seed_id + "#" + query_id + "\t" + seed_id + "\t" + seed_length +
                                                 "\t" + feature + "\t" + inst[0] + "\t" + inst[1] + "\tNA\tN\n")

                    if path.tag == "query_path":
                        for feature in path:
                            forward_q_path[feature.attrib["type"]] = []
                            for instance in feature:
                                forward_q_path[feature.attrib["type"]].append((instance.attrib["start"],
                                                                               instance.attrib["end"]))
                        for feature in arc[query_id]:
                            if feature in forward_q_path and feature in forward_s_path:
                                for inst in arc[query_id][feature]:
                                    if inst in forward_q_path[feature]:
                                        d0_out.write(seed_id + "#" + query_id + "\t" + query_id + "\t" + query_length
                                                     + "\t" + feature + "\t" + inst[0] + "\t" + inst[1] + "\t"
                                                     + "NA\tY\n")  # weights[feature] + "\tY\n")
                                    else:
                                        d0_out.write(seed_id + "#" + query_id + "\t" + query_id + "\t" + query_length
                                                     + "\t" + feature + "\t" + inst[0] + "\t" + inst[1] + "\t"
                                                     + "NA\tN\n")  # weights[feature] + "\tN\n")
                            elif feature in forward_q_path:
                                for inst in arc[query_id][feature]:
                                    if inst in forward_q_path[feature]:
                                        d0_out.write(seed_id + "#" + query_id + "\t" + query_id + "\t" + query_length
                                                     + "\t" + feature + "\t" + inst[0] + "\t" + inst[1] + "\tNA\tY\n")
                                    else:
                                        d0_out.write(seed_id + "#" + query_id + "\t" + query_id + "\t" + query_length
                                                     + "\t" + feature + "\t" + inst[0] + "\t" + inst[1] + "\tNA\tN\n")
                            else:
                                for inst in arc[query_id][feature]:
                                    d0_out.write(seed_id + "#" + query_id + "\t" + query_id + "\t" + query_length +
                                                 "\t" + feature + "\t" + inst[0] + "\t" + inst[1] + "\tNA\tN\n")
    if extendedout:
        d0_out.close()
    if bidirectional:
        reversetree = ElTre.parse(outpath + "_reverse.xml")
        reverseroot = reversetree.getroot()
        d1_out = open(outpath + "_reverse.domains", "w")
        for seed in reverseroot:
            seed_id, seed_length = seed.attrib["id"], seed.attrib["length"]
            for query in seed:
                query_id, reverse_score, query_length = query.attrib["id"], query.attrib["score"], \
                                                        query.attrib["length"]
                # outdict[(seed_id, query_id)] = (outdict[(seed_id, query_id)][0], reverse_score)
                if extendedout:
                    # weights = {}
                    reverse_q_path = {}
                    reverse_s_path = {}
                    for path in query:
                        if path.tag == "template_path":
                            for feature in path:
                                # weights[feature.attrib["type"]] = feature.attrib["corrected_weight"]
                                reverse_q_path[feature.attrib["type"]] = []
                                for instance in feature:
                                    reverse_q_path[feature.attrib["type"]].append((instance.attrib["start"],
                                                                                   instance.attrib["end"]))
                            for feature in arc[query_id]:
                                if feature in reverse_q_path:
                                    for inst in arc[query_id][feature]:
                                        if inst in reverse_q_path[feature]:
                                            d1_out.write(seed_id + "#" + query_id + "\t" + query_id + "\t" +
                                                         query_length + "\t" + feature + "\t" + inst[0] + "\t" +
                                                         inst[1] + "\t" + "NA\tY\n")  # weights[feature] + "\tY\n")
                                        else:
                                            d1_out.write(seed_id + "#" + query_id + "\t" + query_id + "\t" +
                                                         query_length + "\t" + feature + "\t" + inst[0] + "\t" +
                                                         inst[1] + "\t" + "NA\tY\n")  # weights[feature] + "\tN\n")
                                else:
                                    for inst in arc[query_id][feature]:
                                        d1_out.write(seed_id + "#" + query_id + "\t" + query_id + "\t" + query_length
                                                     + "\t" + feature + "\t" + inst[0] + "\t" + inst[1] + "\tNA\tN\n")

                        if path.tag == "query_path":
                            for feature in path:
                                reverse_s_path[feature.attrib["type"]] = []
                                for instance in feature:
                                    reverse_s_path[feature.attrib["type"]].append((instance.attrib["start"],
                                                                                   instance.attrib["end"]))
                            for feature in arc[seed_id]:
                                if feature in reverse_s_path and feature in reverse_q_path:
                                    for inst in arc[seed_id][feature]:
                                        if inst in reverse_s_path[feature]:
                                            d1_out.write(seed_id + "#" + query_id + "\t" + seed_id + "\t" +
                                                         seed_length + "\t" + feature + "\t" + inst[0] + "\t" +
                                                         inst[1] + "\t" + "NA\tY\n")  # weights[feature] + "\tY\n")
                                        else:
                                            d1_out.write(seed_id + "#" + query_id + "\t" + seed_id + "\t" +
                                                         seed_length + "\t" + feature + "\t" + inst[0] + "\t" +
                                                         inst[1] + "\t" + "NA\tY\n")  # weights[feature] + "\tN\n")
                                elif feature in reverse_s_path:
                                    for inst in arc[seed_id][feature]:
                                        if inst in reverse_s_path[feature]:
                                            d1_out.write(seed_id + "#" + query_id + "\t" + seed_id + "\t" +
                                                         seed_length + "\t" + feature + "\t" + inst[0] + "\t" +
                                                         inst[1] + "\tNA\tY\n")
                                        else:
                                            d1_out.write(seed_id + "#" + query_id + "\t" + seed_id + "\t" +
                                                         seed_length + "\t" + feature + "\t" + inst[0] + "\t" +
                                                         inst[1] + "\tNA\tN\n")
                                else:
                                    for inst in arc[seed_id][feature]:
                                        d1_out.write(seed_id + "#" + query_id + "\t" + seed_id + "\t" + seed_length +
                                                     "\t" + feature + "\t" + inst[0] + "\t" + inst[1] + "\tNA\tN\n")
        if extendedout:
            d1_out.close()


def phyloprofile_out(outpath, bidirectional, mapping_file):
    with open(mapping_file) as infile:
        map = {}
        for line in infile.readlines():
            cells = line.rstrip("\n").split("\t")
            map[cells[0]] = cells[1]
    outdict = {}
    groupname = outpath.split("/")[-1]
    forwardtree = ElTre.parse(outpath + ".xml")
    forwardroot = forwardtree.getroot()
    for query in forwardroot:
        query_id, query_length = query.attrib["id"], query.attrib["length"]
        for seed in query:
            seed_id, forward_score, seed_length = seed.attrib["id"], seed.attrib["score"], seed.attrib["length"]
            outdict[(seed_id, query_id)] = (forward_score, "0.0")
    if bidirectional:
        reversetree = ElTre.parse(outpath + "_reverse.xml")
        reverseroot = reversetree.getroot()
        for seed in reverseroot:
            seed_id, seed_length = seed.attrib["id"], seed.attrib["length"]
            for query in seed:
                query_id, reverse_score, query_length = query.attrib["id"], query.attrib["score"], \
                                                        query.attrib["length"]
                outdict[(seed_id, query_id)] = (outdict[(seed_id, query_id)][0], reverse_score)

    out = open(outpath + ".phyloprofile", "w")
    out.write("geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n")
    for pair in outdict:
        try:
            out.write(groupname + "\t" + map[pair[1]] + "\t" + pair[1] + "\t" + outdict[pair][0] + "\t" +
                      outdict[pair][1] + "\n")
        except KeyError:
            raise Exception(pair[1] + " not in mapping file")
    out.close()
