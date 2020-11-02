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
import shutil


def main():
    args = get_options()
    if not args.outname:
        outname = args.extended_fa.split('/')[-1].split('.')[0]
    else:
        outname = args.outname
    if not args.out_dir:
        out_dir = '/'.join(os.path.abspath(args.extended_fa).split('/')[0:-1]) + '/'
    else:
        out_dir = args.out_dir + '/'
    if not args.tmp_dir:
        tmp_dir = out_dir
    else:
        tmp_dir = args.tmp_dir + '/'
    print(out_dir)
    joblist = read_extended_fa(args.extended_fa, args.groupnames)
    jobdict, groupdict, seedspec = create_jobdict(joblist)
    if args.no_lin:
        features = [[], ["pfam", "smart", "flps", "coils2", "seg", "signalp", "tmhmm"]]
    else:
        features = [["pfam", "smart"], ["flps", "coils2", "seg", "signalp", "tmhmm"]]
    print('calculating FAS scores for ' + outname + '...')
    Path(tmp_dir + '/' + outname).mkdir(parents=True, exist_ok=True)
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    results = manage_jobpool(jobdict, groupdict, seedspec, args.weight_dir, tmp_dir + '/' + outname, args.cores,
                             features, args.bidirectional)
    print('writing phyloprofile output...')
    write_phyloprofile(results, out_dir, outname, groupdict)
    join_domain_out(jobdict, tmp_dir + "/" + outname, out_dir, args.bidirectional, outname,
                    seedspec, groupdict)
    shutil.rmtree(tmp_dir + "/" + outname, ignore_errors=True)
    print('fdogFAS finished!')
    if args.bidirectional:
        print('fdogFAS files: ' + outname + '.phyloprofile, ' + outname + '_forward.domains, ' +
              outname + '_reverse.domains in ' + out_dir)
    else:
        print('fdogFAS files: ' + outname + '.phyloprofile, ' + outname + '_forward.domains in' +
              out_dir)


def read_extended_fa(path, grouplist):
    joblist = {}
    with open(path, 'r') as infile:
        for line in infile.readlines():
            if line[0] == '>':
                cells = line.lstrip('>').rstrip('\n').split('|')
                if grouplist:
                    if cells[0] in grouplist and cells[0] not in joblist:
                        joblist[cells[0]] = []
                elif cells[0] not in joblist:
                    joblist[cells[0]] = []
                if cells[0] in joblist:
                    joblist[cells[0]].append(cells[1:])
    return joblist


def create_jobdict(joblist):# check jobdict generation
    jobdict = {}
    groupdict = {}
    seedspec = None
    for entry in joblist:
        seed = '|'.join(joblist[entry][0][1:-1])
        if seedspec and not seedspec == joblist[entry][0][0]:
            raise Exception(
                'There seem to be multiple seed species in the extended.fa but fdogFAS only supports a single one')
        elif not seedspec:
            seedspec = joblist[entry][0][0]
        if seed not in groupdict:
            groupdict[seed] = {entry: {}}
        else:
            groupdict[seed][entry] = {}
        for query in joblist[entry]:
            prot_id = '|'.join(query[1:-1])
            spec = query[0]
            groupdict[seed][entry][prot_id] = query[-1]
            if spec not in jobdict:
                jobdict[spec] = [(seed, prot_id)]
            elif (seed, prot_id) not in jobdict[spec]:
                jobdict[spec].append((seed, prot_id))
    return jobdict, groupdict, seedspec


def manage_jobpool(jobdict, seed_names, seed_spec, weight_dir, tmp_path, cores, features, bidirectional):
    try:
        tmp_data = read_json(weight_dir + "/" + seed_spec + ".json")
    except FileNotFoundError:
        raise Exception('Taxon: "' + seed_spec + '" is missing in the weight_dir')
    seed_weight = w_weight_correction("loge", tmp_data["count"])
    seed_proteome = tmp_data["feature"]
    for seed_name in seed_names:
        if seed_name not in seed_proteome:
            raise Exception('The protein: "' + seed_name + '" is missing in taxon: "' + seed_spec +
                            '". The annotations in weight_dir should contain all proteins from the genome_dir.')
    clan_dict = tmp_data["clan"]
    data = []
    for spec in jobdict:
        data.append([spec,
                    {"weight_const": False, "seed_id": None, "query_id": None,
                     "priority_mode": True, "priority_threshold": 30, "max_cardinality": 500, "eFeature": 0.001,
                     "cores": 1, "eInstance": 0.01, "e_output": True, "feature_info": None,
                     "bidirectional": bidirectional, "raw": False, "silent": False, "reverse": False,
                     "max_overlap": 0, "classicMS": False, "timelimit": 0, "ref_2": None, "phyloprofile": None,
                     "score_weights": (0.7, 0.0, 0.3), "output": 0, "max_overlap_percentage": 0.0, "domain": True,
                     "pairwise": jobdict[spec], "weight_correction": "loge", "outpath": tmp_path + "/" + spec,
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
        tmp_data = read_json(data[4] + "/" + data[0] + ".json")
    except FileNotFoundError:
        raise Exception('Taxon: "' + data[0] + '" is missing in the weight_dir')
    query_proteome = {}
    clan_dict = data[5]
    clan_dict.update(tmp_data["clan"])
    seed_proteome = data[2]
    weight = w_weight_correction("loge", tmp_data["count"])
    for i in data[1]['pairwise']:
        try:
            query_proteome[i[1]] = tmp_data["feature"][i[1]]
        except KeyError:
            raise Exception('The protein: "' + i[1] + '" is missing in taxon: "' + data[0] + '". The annotations ' +
                            'in weight_dir should contain all proteins from the genome_dir.')

    f_results = greedyFAS.fc_main(weight, seed_proteome, query_proteome, clan_dict, data[1])
    outdata = {}
    for result in f_results:
        outdata[result[0], result[1]] = (result[2][0], 0.0)
    if data[1]["bidirectional"]:
        data[1]["reverse"] = True
        pairtmp = []
        for pair in data[1]['pairwise']:
            pairtmp.append((pair[1], pair[0]))
        data[1]['pairwise'] = pairtmp
        r_results = greedyFAS.fc_main(data[3], query_proteome, seed_proteome, clan_dict, data[1])
        for result in r_results:
            outdata[result[1], result[0]] = (outdata[result[1], result[0]][0], result[2][0])
    return outdata, data[0]


def join_domain_out(jobdict, tmp_path, out_path, bidirectional, outname, seed_spec, groupdict):
    out_f = open(out_path + outname + "_forward.domains", "w")
    out_r = None
    if bidirectional:
        out_r = open(out_path + outname + "_reverse.domains", "w")
    for spec in jobdict:
        with open(tmp_path + "/" + spec + "_forward.domains", "r") as infile:
            for line in infile.readlines():
                cells = line.split("\t")
                s_id, q_id = cells[0].split("#")
                for seed in groupdict[s_id]:
                    if q_id in groupdict[s_id][seed]:
                        if not cells[1] == q_id:
                            p_id = seed_spec + "|" + cells[1]
                        else:
                            p_id = spec + "|" + cells[1] + "|" + groupdict[s_id][seed][q_id]
                        out_f.write(seed + "#" + seed + "|" + spec + "|" + q_id + "|" + groupdict[s_id][seed][q_id] +
                                    "\t" + seed + "|" + p_id + "\t" + "\t".join(cells[2:]))
        os.remove(tmp_path + "/" + spec + "_forward.domains")
        if bidirectional:
            with open(tmp_path + "/" + spec + "_reverse.domains", "r") as infile:
                for line in infile.readlines():
                    cells = line.split("\t")
                    s_id, q_id = cells[0].split("#")
                    for seed in groupdict[s_id]:
                        if q_id in groupdict[s_id][seed]:
                            if not cells[1] == q_id:
                                p_id = seed_spec + "|" + cells[1]
                            else:
                                p_id = spec + "|" + cells[1] + "|" + groupdict[s_id][seed][q_id]
                            out_r.write(seed + "#" + seed + "|" + spec + "|" + q_id + "|" + groupdict[s_id][seed][q_id]
                                        + "\t" + seed + "|" + p_id + "\t" + "\t".join(cells[2:]))
            os.remove(tmp_path + "/" + spec + "_reverse.domains")
    out_f.close()
    if bidirectional:
        out_r.close()


def write_phyloprofile(results, out_path, outname, groupdict):
    out = open(out_path + outname + ".phyloprofile", "w+")
    out.write("geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n")
    for result in results:
        spec = result[1]
        ncbi = spec.split("@")[1]
        for pair in result[0]:
            for seed in groupdict[pair[0]]:
                if pair[1] in groupdict[pair[0]][seed]:
                    out.write(seed + "\tncbi" + ncbi + "\t" + seed + "|" + spec + "|" + pair[1] + "|" +
                          groupdict[pair[0]][seed][pair[1]] + "\t" + str(result[0][pair][0]) + "\t" +
                          str(result[0][pair][1]) + "\n")
    out.close()


def get_options():
    version = '1.4.9'
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    parser.add_argument('--version', action='version', version=str(version))
    required.add_argument("-i", "--extended_fa", default=None, type=str, required=True,
                          help="path to extended.fa file")
    required.add_argument("-w", "--weight_dir", default=None, type=str, required=True,
                          help="path to weight_dir of Hamstr")
    optional.add_argument("-n", "--outname", default=None, type=str,
                          help="name of the ortholog group")
    optional.add_argument("-t", "--tmp_dir", type=str, default=None,
                          help="Path to working directory (temporary files are stored here)")
    optional.add_argument("-o", "--out_dir", type=str, default=None,
                          help="path to out directory")
    optional.add_argument("-s", "--groupnames", default=None, type=str, nargs='*',
                          help="specify which groups in the extended.fa will be calculated")
    optional.add_argument("--bidirectional", action="store_false",
                          help="deactivate bidirectional scoring")
    optional.add_argument("--cores", action="store", type=int, default=1,
                          help="number of cores used for parallel calculation")
    optional.add_argument("--no_lin", action='store_true',
                          help="deactivate linearization for pfam/smart")
    arguments = parser.parse_args()
    return arguments


if __name__ == '__main__':
    main()
