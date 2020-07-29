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
import argparse
import os
from greedyFAS.mainFAS import greedyFAS
from greedyFAS.mainFAS.fasInput import read_json
from greedyFAS.mainFAS.fasWeighting import w_weight_correction
from pathlib import Path
from tqdm import tqdm


def main():
    args = get_options()
    joblist = read_extended_fa(args.extended_fa)
    jobdict, namedict = create_jobdict(joblist)
    features = [["pfam", "smart"], ["flps", "coils2", "seg", "signalp", "tmhmm"]]
    if not args.groupname:
        groupname = args.extended_fa.split('/')[-1].split('.')[0]
    else:
        groupname = args.groupname
    print('calculating FAS scores for ' + groupname + '...')
    Path(args.tmp_dir + '/' + groupname).mkdir(parents=True, exist_ok=True)
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    results = manage_jobpool(jobdict, args.seed_name, args.weight_dir, args.seed_spec, args.tmp_dir + '/' + groupname,
                   args.cores, features, args.bidirectional)
    print('writing phyloprofile output...')
    write_phyloprofile(results, args.out_dir, groupname, namedict)
    join_domain_out(jobdict, args.tmp_dir + "/" + groupname, args.out_dir, args.bidirectional, groupname,
                    args.seed_spec, namedict)
    print('hamstrFAS finished!')
    if args.out_dir[0] == '/':
        out_dir = args.out_dir
    else:
        out_dir = os.getcwd() + '/' + args.out_dir
    if args.bidirectional:
        print('Output files: ' + groupname + '.phyloprofile, ' + groupname + '_forward.domains, ' +
              groupname + '_reverse.domains in ' + out_dir)
    else:
        print('Output files: ' + groupname + '.phyloprofile, ' + groupname + '_forward.domains in' +
              out_dir)


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
        prot_id = '|'.join(entry[2:-1])
        namedict[prot_id] = entry[-1]
        if spec in jobdict:
            jobdict[spec].append(prot_id)
        else:
            jobdict[spec] = [prot_id]
    return jobdict, namedict


def manage_jobpool(jobdict, seed_name, weight_dir, seed_spec, tmp_path, cores, features, bidirectional):
    try:
        tmp_data = read_json(weight_dir + "/" + seed_spec + ".json")
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
                     "priority_mode": True, "priority_threshold": 30, "max_cardinality": 500, "eFeature": 0.001,
                     "cores": 1, "eInstance": 0.01, "e_output": True, "feature_info": None,
                     "bidirectional": bidirectional, "raw": False, "silent": False, "reverse": False,
                     "max_overlap": 0, "classicMS": False, "timelimit": 0, "ref_2": None, "phyloprofile": None,
                     "score_weights": (0.7, 0.0, 0.3), "output": 0, "max_overlap_percentage": 0.0, "domain": True,
                     "pairwise": None, "weight_correction": "loge", "outpath": tmp_path + "/" + spec,
                     "input_linearized": features[0], "input_normal": features[1], "MS_uni": 0,
                     "ref_proteome": [spec + '.json'], "progress": False},
                     seed_proteome, seed_weight, weight_dir, clan_dict])
    jobpool = multiprocessing.Pool(processes=cores)
    results = []
    for _ in tqdm(jobpool.imap_unordered(run_fas, data), total=len(jobdict)):
        results.append(_)
    jobpool.close()
    jobpool.join()
    return results


def run_fas(data):
    try:
        tmp_data = read_json(data[5] + "/" + data[0] + ".json")
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

    f_results = greedyFAS.fc_main(weight, seed_proteome, query_proteome, clan_dict, data[2])
    outdata = {}
    for result in f_results:
        outdata[result[1]] = (result[2][0], 0.0)
    if data[2]["bidirectional"]:
        data[2]["reverse"] = True
        data[2]["query_id"] = data[2]["seed_id"]
        data[2]["seed_id"] = None
        r_results = greedyFAS.fc_main(data[4], query_proteome, seed_proteome, clan_dict, data[2])
        for result in r_results:
            outdata[result[0]] = (outdata[result[0]][0], result[2][0])
    return outdata, data[0]


def join_domain_out(jobdict, tmp_path, out_path, bidirectional, groupname, seed_spec, namedict):
    out_f = open(out_path + "/" + groupname + "_forward.domains", "w")
    out_r = None
    if bidirectional:
        out_r = open(out_path + "/" + groupname + "_reverse.domains", "w")
    for spec in jobdict:
        with open(tmp_path + "/" + spec + "_forward.domains", "r") as infile:
            for line in infile.readlines():
                cells = line.split("\t")
                q_id = cells[0].split("#")[1]
                if not cells[1] == q_id:
                    p_id = seed_spec + "|" + cells[1]
                else:
                    p_id = spec + "|" + cells[1] + "|" + namedict[q_id]
                out_f.write(groupname + "#" + groupname + "|" + spec + "|" + q_id + "|" + namedict[q_id] + "\t" +
                            groupname + "|" + p_id + "\t" + "\t".join(cells[2:]))
        os.remove(tmp_path + "/" + spec + "_forward.domains")
        if bidirectional:
            with open(tmp_path + "/" + spec + "_reverse.domains", "r") as infile:
                for line in infile.readlines():
                    cells = line.split("\t")
                    q_id = cells[0].split("#")[1]
                    if not cells[1] == q_id:
                        p_id = seed_spec + "|" + cells[1]
                    else:
                        p_id = spec + "|" + cells[1] + "|" + namedict[q_id]
                    out_r.write(groupname + "#" + groupname + "|" + spec + "|" + q_id + "|" + namedict[q_id] + "\t" +
                                groupname + "|" + p_id + "\t" + "\t".join(cells[2:]))
            os.remove(tmp_path + "/" + spec + "_reverse.domains")
    out_f.close()
    if bidirectional:
        out_r.close()


def write_phyloprofile(results, out_path, groupname, namedict):
    out = open(out_path + "/" + groupname + ".phyloprofile", "w+")
    out.write("geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n")
    for result in results:
        spec = result[1]
        ncbi = spec.split("@")[1]
        for q_id in result[0]:
            out.write(groupname + "\tncbi" + ncbi + "\t" + groupname + "|" + spec + "|" + q_id + "|" + namedict[q_id] +
                      "\t" + str(result[0][q_id][0]) + "\t" + str(result[0][q_id][1]) + "\n")
    out.close()


def get_options():
    version = '1.2.6'
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-i", "--extended_fa", default=None, type=str, required=True,
                          help="path to extended.fa file")
    optional.add_argument("-n", "--groupname", default=None, type=str,
                          help="name of the ortholog group")
    required.add_argument("-w", "--weight_dir", default=None, type=str, required=True,
                          help="path to weight_dir of Hamstr")
    optional.add_argument("-t", "--tmp_dir", type=str, default='tmp/',
                          help="Path to working directory")
    optional.add_argument("-o", "--out_dir", type=str, default='out/',
                          help="path to out directory")
    required.add_argument("-s", "--seed_name", default=None, type=str, required=True,
                          help="name of the seed protein as it appears in the species fasta and species json")
    required.add_argument("-a", "--seed_spec", default=None, type=str, required=True,
                          help="name of the seed species in genome dir and weight dir")
    optional.add_argument("--bidirectional", action="store_true",
                          help="calculate both scoring directions")
    optional.add_argument('--cores', action='store', type=int, default=1,
                          help='number of cores used for parallel calculation')
    arguments = parser.parse_args()
    return arguments


if __name__ == '__main__':
    main()
