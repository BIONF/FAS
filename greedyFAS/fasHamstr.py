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


import multiprocessing
from sys import version_info
import xml.etree.ElementTree as ElTre
import argparse
if version_info.major == 3:
    from greedyFAS import greedyFAS
    from greedyFAS.processAnnotation import get_data
    from greedyFAS.fasInput import xmlreader
elif version_info.major == 2:
    import greedyFAS
    from processAnnotation import get_data
    from fasInput import xmlreader


def main():
    args = get_options()
    joblist = read_extended_fa(args.extendended_fa)
    jobdict = create_pairs(joblist)
    features = [["pfam", "smart"], ["flps", "coils", "seg", "signalp", "tmhmm"]]
    manage_jobpool(jobdict, args.seed_path, args.weight_dir, args.seed_spec, args.tmp_dir, args.cores, features,
                   args.bidirectional)
    write_phyloprofile(jobdict, args.tmp_dir, args.out_dir, args.bidirectional, args.groupname, args.seed_spec)


def read_extended_fa(path):
    joblist = []
    with open(path, 'r') as infile:
        for line in infile.readlines():
            if line[0] == '>':
                joblist.append(line.lstrip('>').rstrip('\n').split('|'))
    return joblist


def create_pairs(joblist):
    jobdict = {}
    for entry in joblist:
        spec = entry[1]
        prot_id = entry[2]
        if spec in jobdict:
            jobdict[spec].append(prot_id)
        else:
            jobdict[spec] = [prot_id]
    return jobdict


def manage_jobpool(jobdict, seed_path, weight_dir, seed_spec, tmp_path, cores, features, bidirectional):
    tmp_data = get_data(weight_dir + "/" + seed_spec + "/" + seed_spec + ".json")
    seed_weight = tmp_data["domain_count"]
    data = []
    for spec in jobdict:
        data.append([spec, jobdict[spec],
                    {"weight_const": False, "version": '1.0.0', "seed_id": None, "query_id": None,
                     "priority_mode": True, "priority_threshold": 50, "max_cardinality": 5000, "efilter": 0.001,
                     "cores": 1, "inst_efilter": 0.01, "e_output": True, "feature_info": None,
                     "bidirectional": bidirectional,
                     "max_overlap": 0, "classicMS": False, "timelimit": 7200, "ref_2": None, "phyloprofile": None,
                     "score_weights": (0.7, 0.0, 0.3), "output": 0, "max_overlap_percentage": 0.0, "domain": False,
                     "pairwise": None, "weight_correction": "loge", "outpath": tmp_path + "/" + spec,
                     "input_linearized": features[0], "input_normal": features[1]}, seed_path, seed_weight, weight_dir])
    jobpool = multiprocessing.Pool(processes=cores)
    pool_results = jobpool.map_async(run_fas, data)
    jobpool.close()
    jobpool.join()
    pool_results.wait()


def run_fas(data, bidirectional):
    tmp_data = get_data(data[5] + "/" + data[0] + "/" + data[0] + ".json")
    query_proteome = {}
    protein_lengths = {}
    clan_dict = tmp_data["clan_dict"]
    seed_proteome = {}
    for i in data[1]:
        query_proteome[i] = tmp_data["proteome"][i]
        protein_lengths["query_" + i] = tmp_data["protein_lengths"][i]
    for ftype in data[2]["input_linearized"]:
        seed_proteome, protein_lengths, clan_dict = xmlreader(data[3] + "/" + ftype + ".xml", 0, ftype, True,
                                                              seed_proteome, protein_lengths, clan_dict, data[2])
    for ftype in data[2]["input_normal"]:
        seed_proteome, protein_lengths, clan_dict = xmlreader(data[3] + "/" + ftype + ".xml", 0, ftype, False,
                                                              seed_proteome, protein_lengths, clan_dict, data[2])
    greedyFAS.fc_main([], 0, tmp_data["domain_count"], seed_proteome, query_proteome, protein_lengths, clan_dict,
                      data[3])
    if bidirectional:
        tmp = {}
        for protein in protein_lengths:
            if protein[0:4] == "seed":
                tmp["query_" + protein[5:]] = protein_lengths[protein]
            else:
                tmp["seed_" + protein[6:]] = protein_lengths[protein]
        data[2]["e_output"] = False
        data[2]["outpath"] += "_reverse"
        greedyFAS.fc_main([], 0, data[4], query_proteome, seed_proteome, tmp, clan_dict, data[3])


def write_phyloprofile(jobdict, tmp_path, out_path, bidirectional, groupname, seed_spec):
    out = open(out_path + "/" + groupname + ".phyloprofile", "w")
    out.write("geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n")
    out_f = open(out_path + "/" + groupname + "_forward.domains", "w")
    out_b = None
    if bidirectional:
        out_b = open(out_path + "/" + groupname + "_reverse.domains", "w")
    for spec in jobdict:
        outdict = {}
        arc = {}
        ncbi = spec.split("@")[1]
        arctree = ElTre.parse(tmp_path + "/" + spec + "_architecture.xml")
        arcroot = arctree.getroot()
        for architecture in arcroot:
            pid = architecture.attrib["id"]
            arc[pid] = {}
            for feature in architecture[0]:
                f_type = feature.attrib["type"]
                arc[pid][f_type] = []
                for inst in feature:
                    arc[pid][f_type].append((inst.attrib["start"], inst.attrib["end"]))
        forwardtree = ElTre.parse(tmp_path + "/" + spec + ".xml")
        forwardroot = forwardtree.getroot()
        for query in forwardroot:
            query_id, query_length = query.attrib["id"], query.attrib["length"]
            for seed in query:
                seed_id, forward_score, seed_length = seed.attrib["id"], seed.attrib["score"], seed.attrib["length"]
                outdict[(seed_id, query_id)] = (forward_score, "0.0")
                weights = {}
                forward_s_path = {}
                forward_q_path = {}
                for path in seed:
                    if path.tag == "template_path":
                        for feature in path:
                            weights[feature.attrib["type"]] = feature.attrib["corrected_weight"]
                            forward_s_path[feature.attrib["type"]] = []
                            for instance in feature:
                                forward_s_path[feature.attrib["type"]].append((instance.attrib["start"],
                                                                               instance.attrib["end"]))
                        for feature in arc[seed_id]:
                            if feature in forward_s_path:
                                for inst in arc[seed_id][feature]:
                                    out_f.write(
                                        groupname + "#" + groupname + "|" + spec + "|" + query_id + "|" +
                                        spec.split("|")[2] + "\t" + groupname + "|" + seed_spec + "|" + seed_id + "\t" +
                                        seed_length + "\t" + feature + "\t" + inst[0] + "\t" + inst[1] + "\t" +
                                        weights[feature])
                                    if inst in forward_s_path[feature]:
                                        out_f.write("\tY\n")
                                    else:
                                        out_f.write("\tN\n")
                            else:
                                for inst in arc[seed_id][feature]:
                                    out_f.write(groupname + "#" + groupname + "|" + spec + "|" + query_id + "|" +
                                                spec.split("|")[2] + "\t" + groupname + "|" + seed_spec + "|" +
                                                seed_id + "\t" + seed_length + "\t" + feature + "\t" + inst[0] + "\t" +
                                                inst[1] + "\tNA\tN\n")

                    if path.tag == "query_path":
                        for feature in path:
                            forward_q_path[feature.attrib["type"]] = []
                            for instance in feature:
                                forward_q_path[feature.attrib["type"]].append((instance.attrib["start"],
                                                                               instance.attrib["end"]))
                        for feature in arc[query_id]:
                            if feature in forward_q_path and feature in forward_s_path:
                                for inst in arc[query_id][feature]:
                                    out_f.write(
                                        groupname + "#" + groupname + "|" + spec + "|" + query_id + "|" +
                                        spec.split("|")[2] + "\t" + groupname + "|" + spec + "|" + query_id + "|" +
                                        spec.split("|")[2] + "\t" + query_length + "\t" + feature + "\t" + inst[0] +
                                        "\t" + inst[1] + "\t" + weights[feature])
                                    if inst in forward_q_path[feature]:
                                        out_f.write("\tY\n")
                                    else:
                                        out_f.write("\tN\n")
                            elif feature in forward_q_path:
                                for inst in arc[query_id][feature]:
                                    out_f.write(
                                        groupname + "#" + groupname + "|" + spec + "|" + query_id + "|" +
                                        spec.split("|")[2] + "\t" + groupname + "|" + spec + "|" + query_id + "|" +
                                        spec.split("|")[2] + "\t" + query_length + "\t" + feature + "\t" + inst[0] +
                                        "\t" + inst[1])
                                    if inst in forward_q_path[feature]:
                                        out_f.write("\tNA\tY\n")
                                    else:
                                        out_f.write("\tNA\tN\n")
                            else:
                                for inst in arc[query_id][feature]:
                                    out_f.write(groupname + "#" + groupname + "|" + spec + "|" + query_id + "|" +
                                                spec.split("|")[2] + "\t" + groupname + "|" + spec + "|" + query_id +
                                                "\t" + query_length + "\t" + feature + "\t" + inst[0] + "\t" + inst[1] +
                                                "\tNA\tN\n")
        if bidirectional:
            reversetree = ElTre.parse(tmp_path + "/" + spec + "_reverse.xml")
            reverseroot = reversetree.getroot()
            for seed in reverseroot:
                seed_id, seed_length = seed.attrib["id"], seed.attrib["length"]
                for query in seed:
                    query_id, reverse_score, query_length = query.attrib["id"], query.attrib["score"], \
                                                            query.attrib["length"]
                    outdict[(seed_id, query_id)] = (outdict[(seed_id, query_id)][0], reverse_score)
                    weights = {}
                    reverse_q_path = {}
                    reverse_s_path = {}
                    for path in query:
                        if path.tag == "template_path":
                            for feature in path:
                                weights[feature.attrib["type"]] = feature.attrib["corrected_weight"]
                                reverse_q_path[feature.attrib["type"]] = []
                                for instance in feature:
                                    reverse_q_path[feature.attrib["type"]].append((instance.attrib["start"],
                                                                                   instance.attrib["end"]))
                            for feature in arc[query_id]:
                                if feature in reverse_q_path:
                                    for inst in arc[query_id][feature]:
                                        out_b.write(
                                            groupname + "#" + groupname + "|" + spec + "|" + query_id + "|" +
                                            spec.split("|")[2] + "\t" + groupname + "|" + spec + "|" + query_id + "|" +
                                            spec.split("|")[2] + "\t" + query_length + "\t" + feature + "\t" + inst[0] +
                                            "\t" + inst[1] + "\t" + weights[feature])
                                        if inst in reverse_q_path[feature]:
                                            out_b.write("\tY\n")
                                        else:
                                            out_b.write("\tN\n")
                                else:
                                    for inst in arc[query_id][feature]:
                                        out_b.write(groupname + "#" + groupname + "|" + spec + "|" + query_id + "|" +
                                                    spec.split("|")[2] + "\t" + groupname + "|" + spec + "|" +
                                                    query_id + "|" + spec.split("|")[2] + "\t" + query_length + "\t" +
                                                    feature + "\t" + inst[0] + "\t" + inst[1] + "\tNA\tN\n")

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
                                            out_b.write(
                                                groupname + "#" + groupname + "|" + spec + "|" + query_id + "|" +
                                                spec.split("|")[2] + "\t" + groupname + "|" + seed_spec + "|" +
                                                seed_id + "\t" + seed_length + "\t" + feature + "\t" + inst[0] + "\t" +
                                                inst[1] + "\t" + weights[feature])
                                            if inst in reverse_s_path[feature]:
                                                out_f.write("\tY\n")
                                            else:
                                                out_f.write("\tN\n")
                                elif feature in reverse_s_path:
                                    for inst in arc[seed_id][feature]:
                                        out_b.write(
                                            groupname + "#" + groupname + "|" + spec + "|" + query_id + "|" +
                                            spec.split("|")[2] + "\t" + groupname + "|" + seed_spec + "|" + seed_id +
                                            "\t" + query_length + "\t" + feature + "\t" + inst[0] + "\t" + inst[1])
                                        if inst in reverse_s_path[feature]:
                                            out_f.write("\tNA\tY\n")
                                        else:
                                            out_f.write("\tNA\tN\n")
                                else:
                                    for inst in arc[seed_id][feature]:
                                        out_b.write(groupname + "#" + groupname + "|" + spec + "|" + query_id + "|" +
                                                    spec.split("|")[2] + "\t" + groupname + "|" + seed_spec + "|" +
                                                    seed_id + "\t" + query_length + "\t" + feature + "\t" + inst[0] +
                                                    "\t" + inst[1] + "\tNA\tN\n")
        for pair in outdict:
            out.write(groupname + "\tncbi" + ncbi + "\t" + groupname + "|" + spec + "|" + pair[1] + "|" +
                      spec.split("|")[2] + "\t" + outdict[pair][0] + "\t" + outdict[pair][1] + "\n")
    out.close()


def get_options():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-i", "--extended_fa", default=None, type=str, required=True,
                          help="path to extended.fa file")
    required.add_argument("-n", "--groupname", default=None, type=str, required=True,
                          help="groupname")
    required.add_argument("-w", "--weight_dir", default=None, type=str, required=True,
                          help="path to weight_dir of Hamstr")
    required.add_argument("-t", "--tmp_dir", default=None, type=str,
                          help="Path to working directory")
    required.add_argument("-o", "--out_dir", default=None, type=str,
                          help="path to out_dir")
    required.add_argument("-s", "--seed_path", default=None, type=str,
                          help="path to seed annotation")
    required.add_argument("-a", "--seed_spec", default=None, type=str,
                          help="name of the seed species (in weight dir)")
    optional.add_argument("--bidirectional", action="store_true",
                          help="calculate both scoring directions")
    optional.add_argument('--cores', action='store', type=int, default=1,
                          help='number of cores')
    arguments = parser.parse_args()
    return arguments


if __name__ == '__main__':
    main()
