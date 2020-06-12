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
import xml.etree.ElementTree as ElTre
import argparse
import os
from greedyFAS import greedyFAS
from greedyFAS.fasInput import read_json
from greedyFAS.fasWeighting import w_weight_correction


def main():
    args = get_options()
    joblist = read_extended_fa(args.extended_fa)
    jobdict, namedict = create_jobdict(joblist)
    features = [["pfam", "smart"], ["flps", "coils2", "seg", "signalp", "tmhmm"]]
    manage_jobpool(jobdict, args.seed_name, args.weight_dir, args.seed_spec, args.tmp_dir, args.cores, features,
                   args.bidirectional)
    write_phyloprofile(jobdict, args.tmp_dir, args.out_dir, args.bidirectional, args.groupname, args.seed_spec,
                       namedict)


def read_extended_fa(path):
    joblist = []
    with open(path, 'r') as infile:
        for line in infile.readlines():
            if line[0] == '>':
                joblist.append(line.lstrip('>').rstrip('\n').split('|'))
    return joblist


def create_jobdict(joblist):
    jobdict = {}
    namedict = {}
    for entry in joblist:
        spec = entry[1]
        prot_id = entry[2]
        namedict[prot_id] = entry
        if spec in jobdict:
            jobdict[spec].append(prot_id)
        else:
            jobdict[spec] = [prot_id]
    return jobdict, namedict


def manage_jobpool(jobdict, seed_name, weight_dir, seed_spec, tmp_path, cores, features, bidirectional):
    try:
        tmp_data = read_json(weight_dir + "/" + seed_spec + "/" + seed_spec + ".json")
    except FileNotFoundError:
        raise Exception('Taxon: "' + seed_spec + '" is missing in the weight_dir')
    seed_weight = w_weight_correction("loge", tmp_data["count"])
    seed_proteome = tmp_data["feature"]
    if seed_name not in seed_proteome:
        raise Exception('The protein: "' + seed_name + '" is missing in taxon: "' + seed_spec + '". The annotations ' +
                        'in weight_dir should contain all proteins from the genome_dir.')
    clan_dict = tmp_data["clan"]
    data = []
    for spec in jobdict:
        data.append([spec, jobdict[spec],
                    {"weight_const": False, "version": '1.0.0', "seed_id": [seed_name], "query_id": None,
                     "priority_mode": True, "priority_threshold": 50, "max_cardinality": 5000, "eFeature": 0.001,
                     "cores": 1, "eInstance": 0.01, "e_output": True, "feature_info": None,
                     "bidirectional": bidirectional,
                     "max_overlap": 0, "classicMS": False, "timelimit": 7200, "ref_2": None, "phyloprofile": None,
                     "score_weights": (0.7, 0.0, 0.3), "output": 0, "max_overlap_percentage": 0.0, "domain": False,
                     "pairwise": None, "weight_correction": "loge", "outpath": tmp_path + "/" + spec,
                     "input_linearized": features[0], "input_normal": features[1], "MS_uni": 0,
                     "ref_proteome": [spec + '.json']},
                     seed_proteome, seed_weight, weight_dir, clan_dict])
    jobpool = multiprocessing.Pool(processes=cores)
    jobpool.map(run_fas, data)
    jobpool.close()
    jobpool.join()


def run_fas(data):
    try:
        tmp_data = read_json(data[5] + "/" + data[0] + "/" + data[0] + ".json")
    except FileNotFoundError:
        raise Exception('Taxon: "' + data[0] + '" is missing in the weight_dir')
    query_proteome = {}
    clan_dict = data[6]
    clan_dict.update(tmp_data["clan"])
    seed_proteome = data[3]
    weight = w_weight_correction("loge", tmp_data["count"])
    for i in data[1]:
        try:
            query_proteome[i] = tmp_data["feature"][i]
        except KeyError:
            raise Exception('The protein: "' + i + '" is missing in taxon: "' + data[0] + '". The annotations ' +
                            'in weight_dir should contain all proteins from the genome_dir.')

    greedyFAS.fc_main(weight, seed_proteome, query_proteome, clan_dict, data[2])
    if data[2]["bidirectional"]:
        data[2]["e_output"] = False
        data[2]["outpath"] += "_reverse"
        data[2]["query_id"] = data[2]["seed_id"]
        data[2]["seed_id"] = None
        greedyFAS.fc_main(data[4], query_proteome, seed_proteome, clan_dict, data[2])


def write_phyloprofile(jobdict, tmp_path, out_path, bidirectional, groupname, seed_spec, namedict):
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
                                        groupname + "#" + spec + "|" + query_id + "|" + namedict[query_id][3] + "\t" +
                                        seed_spec + "|" + seed_id + "\t" + seed_length + "\t" + feature + "\t" +
                                        inst[0] + "\t" + inst[1] + "\t" + weights[feature])
                                    if inst in forward_s_path[feature]:
                                        out_f.write("\tY\n")
                                    else:
                                        out_f.write("\tN\n")
                            else:
                                for inst in arc[seed_id][feature]:
                                    out_f.write(groupname + "#" + spec + "|" + query_id + "|" + namedict[query_id][3] +
                                                "\t" + seed_spec + "|" + seed_id + "\t" + seed_length + "\t" +
                                                feature + "\t" + inst[0] + "\t" + inst[1] + "\tNA\tN\n")

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
                                        groupname + "#" + spec + "|" + query_id + "|" + namedict[query_id][3] + "\t" +
                                        spec + "|" + query_id + "|" + namedict[query_id][3] + "\t" + query_length +
                                        "\t" + feature + "\t" + inst[0] + "t" + inst[1] + "\t" + weights[feature])
                                    if inst in forward_q_path[feature]:
                                        out_f.write("\tY\n")
                                    else:
                                        out_f.write("\tN\n")
                            elif feature in forward_q_path:
                                for inst in arc[query_id][feature]:
                                    out_f.write(
                                        groupname + "#" + spec + "|" + query_id + "|" + namedict[query_id][3] + "\t" +
                                        spec + "|" + query_id + "|" + namedict[query_id][3] + "\t" + query_length +
                                        "\t" + feature + "\t" + inst[0] + "\t" + inst[1])
                                    if inst in forward_q_path[feature]:
                                        out_f.write("\tNA\tY\n")
                                    else:
                                        out_f.write("\tNA\tN\n")
                            else:
                                for inst in arc[query_id][feature]:
                                    out_f.write(groupname + "#" + groupname + "|" + spec + "|" + query_id + "|" +
                                                namedict[query_id][3] + "\t" + spec + "|" + query_id + "|" +
                                                namedict[query_id][3] + "\t" + query_length + "\t" + feature + "\t" +
                                                inst[0] + "\t" + inst[1] + "\tNA\tN\n")
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
                                            groupname + "#" + spec + "|" + query_id + "|" + namedict[query_id][3] +
                                            "\t" + spec + "|" + query_id + "|" + namedict[query_id][3] + "\t" +
                                            query_length + "\t" + feature + "\t" + inst[0] + "\t" + inst[1] + "\t" +
                                            weights[feature])
                                        if inst in reverse_q_path[feature]:
                                            out_b.write("\tY\n")
                                        else:
                                            out_b.write("\tN\n")
                                else:
                                    for inst in arc[query_id][feature]:
                                        out_b.write(groupname + "#" + spec + "|" + query_id + "|" +
                                                    namedict[query_id][3] + "\t" + spec + "|" + query_id + "|" +
                                                    namedict[query_id][3] + "\t" + query_length + "\t" +
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
                                                groupname + "#" + spec + "|" + query_id + "|" +
                                                namedict[query_id][3] + "\t" + seed_spec + "|" +
                                                seed_id + "\t" + seed_length + "\t" + feature + "\t" + inst[0] + "\t" +
                                                inst[1] + "\t" + weights[feature])
                                            if inst in reverse_s_path[feature]:
                                                out_b.write("\tY\n")
                                            else:
                                                out_b.write("\tN\n")
                                elif feature in reverse_s_path:
                                    for inst in arc[seed_id][feature]:
                                        out_b.write(
                                            groupname + "#" + spec + "|" + query_id + "|" + namedict[query_id][3] +
                                            "\t" + seed_spec + "|" + seed_id + "\t" + query_length + "\t" + feature +
                                            "\t" + inst[0] + "\t" + inst[1])
                                        if inst in reverse_s_path[feature]:
                                            out_b.write("\tNA\tY\n")
                                        else:
                                            out_b.write("\tNA\tN\n")
                                else:
                                    for inst in arc[seed_id][feature]:
                                        out_b.write(groupname + "#" + spec + "|" + query_id + "|" +
                                                    namedict[query_id][3] + "\t" + seed_spec + "|" + seed_id + "\t" +
                                                    query_length + "\t" + feature + "\t" + inst[0] + "\t" + inst[1] +
                                                    "\tNA\tN\n")
            os.remove(tmp_path + "/" + spec + "_reverse.xml")
        os.remove(tmp_path + "/" + spec + "_architecture.xml")
        os.remove(tmp_path + "/" + spec + ".xml")
        for pair in outdict:
            out.write(groupname + "\tncbi" + ncbi + "\t" + spec + "|" + pair[1] + "|" + namedict[pair[1]][3] + "\t" +
                      outdict[pair][0] + "\t" + outdict[pair][1] + "\n")
    out.close()


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--extended_fa", default=None, type=str, required=True,
                        help="path to extended.fa file")
    parser.add_argument("-n", "--groupname", default=None, type=str, required=True,
                        help="name of the ortholog group")
    parser.add_argument("-w", "--weight_dir", default=None, type=str, required=True,
                        help="path to weight_dir of Hamstr")
    parser.add_argument("-t", "--tmp_dir", default=None, type=str,
                        help="Path to working directory")
    parser.add_argument("-o", "--out_dir", default=None, type=str,
                        help="path to out directory")
    parser.add_argument("-s", "--seed_name", default=None, type=str,
                        help="name of the seed protein as it appears in the species fasta and species json")
    parser.add_argument("-a", "--seed_spec", default=None, type=str,
                        help="name of the seed species in genome dir and weight dir")
    parser.add_argument("--bidirectional", action="store_true",
                        help="calculate both scoring directions")
    parser.add_argument('--cores', action='store', type=int, default=1,
                        help='number of cores used for parallel calculation')
    arguments = parser.parse_args()
    return arguments


if __name__ == '__main__':
    main()
