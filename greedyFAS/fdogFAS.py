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
import fnmatch
from pathlib import Path
import shutil
from tqdm import tqdm
from greedyFAS.mainFAS import greedyFAS
from greedyFAS.mainFAS.fasInput import read_json, featuretypes
from greedyFAS.mainFAS.fasWeighting import w_weight_correction
from greedyFAS.annoFAS.checkAnno import doAnnoForMissing
from pkg_resources import get_distribution

def main():
    args = get_options()
    if not args.outname:
        outname = args.extended_fa.split('/')[-1].replace('.extended.fa', '')  #.split('.')[0]
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
    joblist, fasta = read_extended_fa(args.extended_fa, args.groupnames, args.redo_anno)
    if len(joblist) == 0:
        raise Exception(args.extended_fa + " is empty!")
    jobdict, groupdict, seedspec = create_jobdict(joblist)
    pathconfigfile = os.path.realpath(__file__).replace('fdogFAS.py', 'pathconfig.txt')
    with open(pathconfigfile) as f:
        toolpath = f.readline().strip()
    if args.featuretypes:
        features = featuretypes(args.featuretypes)
    else:
        features = featuretypes(toolpath + '/' + 'annoTools.txt')
    if args.no_lin:
        features = ([], features[1] + features[0])
    print('calculating FAS scores for ' + outname + '...')
    Path(tmp_dir + '/' + outname).mkdir(parents=True, exist_ok=True)
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    results = manage_jobpool(jobdict, groupdict, seedspec, args.weight_dir, tmp_dir + '/' + outname, args.cores,
                             features, args.bidirectional, fasta, args.featuretypes)
    print('writing phyloprofile output...')
    write_phyloprofile(results, out_dir, outname, groupdict)
    join_domain_out(jobdict, tmp_dir + "/" + outname, out_dir, args.bidirectional, outname,
                    seedspec, groupdict)
    if args.redo_anno:
        err_taxa = []
        files = os.listdir(tmp_dir + "/" + outname)
        pattern = "*.json"
        for entry in files:
            if fnmatch.fnmatch(entry, pattern):
                err_taxa.append(entry)
        if err_taxa:
            print('\033[93mWARNING: There were missing annotations in the following taxa:')
            for i in err_taxa:
                print(i)
    shutil.rmtree(tmp_dir + "/" + outname, ignore_errors=True)
    print('\033[0mfdogFAS finished!')


def read_extended_fa(path, grouplist, reanno):
    joblist = {}
    seq = ''
    fasta = {}
    cells = []
    with open(path, 'r') as infile:
        for line in infile.readlines():
            if line[0] == '>':
                if seq and reanno:
                    if cells == []:
                        raise Exception('Please check the fasta file for errors')
                    elif cells[1] not in fasta:
                        fasta[cells[1]] = {}
                    fasta[cells[1]][cells[2]] = seq
                seq = ''
                cells = line.lstrip('>').rstrip('\n').split('|')
                if grouplist:
                    if cells[0] in grouplist and cells[0] not in joblist:
                        joblist[cells[0]] = []
                elif cells[0] not in joblist:
                    joblist[cells[0]] = []
                if cells[0] in joblist:
                    joblist[cells[0]].append(cells[1:])
            else:
                seq = seq + line
        if seq and reanno:
            if cells[1] not in fasta:
                fasta[cells[1]] = {}
            fasta[cells[1]][cells[2]] = seq
    return joblist, fasta


def create_jobdict(joblist):  # check jobdict generation
    jobdict = {}
    groupdict = {}
    seedspec = None
    for entry in joblist:
        seed = joblist[entry][0][1]
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
            prot_id = query[1]
            spec = query[0]
            groupdict[seed][entry][prot_id] = query[-1]
            if spec not in jobdict:
                jobdict[spec] = [(seed, prot_id)]
            elif (seed, prot_id) not in jobdict[spec]:
                jobdict[spec].append((seed, prot_id))
    return jobdict, groupdict, seedspec


def manage_jobpool(jobdict, seed_names, seed_spec, weight_dir, tmp_path, cores, features, bidirectional, fasta, annoToolFile):
    missing = []
    for spec in jobdict:
        if not os.path.exists(weight_dir + "/" + spec + ".json"):
            missing.append(spec)
    try:
        tmp_data = read_json(weight_dir + "/" + seed_spec + ".json")
    except FileNotFoundError:
        missing.append(seed_spec)
    if missing:
        raise Exception('The following taxa are missing in the weight_dir:\n' + '\n'.join(missing))
    seed_weight = w_weight_correction("loge", tmp_data["count"])
    seed_proteome = tmp_data["feature"]
    missing = []
    for seed_name in seed_names:
        if seed_name not in seed_proteome:
            if fasta:
                missing.append(seed_name)
            else:
                raise Exception('The protein: "' + seed_name + '" is missing in taxon: "' + seed_spec +
                                '". The annotations in weight_dir should contain all proteins from the genome_dir.')
    if missing:
        missinseq = []
        for i in missing:
            missinseq.append('>' + i + '\n' + fasta[seed_spec][i])
        doAnnoForMissing(seed_spec, missinseq, weight_dir + "/" + seed_spec + ".json", tmp_path + "/", cores, True,
                         annoToolFile)
        tmp_data = read_json(tmp_path + "/" + seed_spec + ".json")
        seed_weight = w_weight_correction("loge", tmp_data["count"])
        seed_proteome = tmp_data["feature"]
    clan_dict = tmp_data["clan"]
    interprokeys = {}
    phmm = {}
    if 'interproID' in tmp_data:
        interprokeys.update(tmp_data['interproID'])
    if 'length' in tmp_data:
        phmm.update(tmp_data['length'])
    data = []
    for spec in jobdict:
        data.append([spec,
                    {"weight_const": False, "seed_id": None, "query_id": None, "empty_as_1": False,
                     "priority_mode": True, "priority_threshold": 30, "max_cardinality": 500, "eFeature": 0.001,
                     "cores": 1, "eInstance": 0.01, "e_output": True, "feature_info": None,
                     "bidirectional": bidirectional, "raw": False, "silent": False, "reverse": False,
                     "max_overlap": 0, "classicMS": False, "timelimit": 0, "ref_2": None, "phyloprofile": None,
                     "score_weights": (0.7, 0.0, 0.3), "output": 0, "max_overlap_percentage": 0.0, "domain": True,
                     "pairwise": jobdict[spec], "weight_correction": "loge", "outpath": tmp_path + "/" + spec,
                     "input_linearized": features[0], "input_normal": features[1], "MS_uni": 0,
                     "ref_proteome": [spec + '.json'], "progress": False, "old_json": False, "paths_limit": 0},
                     seed_proteome, seed_weight, weight_dir, clan_dict, fasta, tmp_path, interprokeys, phmm])
    jobpool = multiprocessing.Pool(processes=cores)
    results = []
    for _ in tqdm(jobpool.imap_unordered(run_fas, data), total=len(jobdict)):
        results.append(_)
    jobpool.close()
    jobpool.join()
    return results


def run_fas(data):
    tmp = True
    try:
        tmp_data = read_json(data[7] + "/" + data[0] + ".json")
    except FileNotFoundError:
        try:
            tmp_data = read_json(data[4] + "/" + data[0] + ".json")
            tmp = False
        except FileNotFoundError:
            raise Exception('Taxon: "' + data[0] + '" is missing in the weight_dir')
    query_proteome = {}
    missing = []
    for i in data[1]['pairwise']:
        try:
            query_proteome[i[1]] = tmp_data["feature"][i[1]]
        except KeyError:
            if data[6]:
                missing.append(i[1])
            else:
                raise Exception('The protein: "' + i[1] + '" is missing in taxon: "' + data[0] + '". The annotations ' +
                            'in weight_dir should contain all proteins from the genome_dir.')
    if missing:
        missinseq = []
        for i in missing:
            missinseq.append('>' + i + '\n' + data[6][data[0]][i])
        if tmp:
            weight_dir = data[7]
        else:
            weight_dir = data[4]
        doAnnoForMissing(data[0], missinseq, weight_dir + "/" + data[0] + ".json", data[7] + "/", 1, True, "None")
        tmp_data = read_json(data[7] + "/" + data[0] + ".json")
        for i in missing:
            query_proteome[i] = tmp_data["feature"][i]
    interprokeys = data[8]
    phmm = data[9]
    if 'interproID' in tmp_data:
        interprokeys.update(tmp_data['interproID'])
    if 'length' in tmp_data:
        phmm.update(tmp_data['length'])
    clan_dict = data[5]
    clan_dict.update(tmp_data["clan"])
    seed_proteome = data[2]
    weight = w_weight_correction("loge", tmp_data["count"])
    f_results = greedyFAS.fc_main(weight, seed_proteome, query_proteome, clan_dict, data[1], interprokeys, phmm)
    outdata = {}
    for result in f_results:
        outdata[result[0], result[1]] = (result[2][0], 0.0)
    if data[1]["bidirectional"]:
        data[1]["reverse"] = True
        pairtmp = []
        for pair in data[1]['pairwise']:
            pairtmp.append((pair[1], pair[0]))
        data[1]['pairwise'] = pairtmp
        r_results = greedyFAS.fc_main(data[3], query_proteome, seed_proteome, clan_dict, data[1], interprokeys, phmm)
        for result in r_results:
            outdata[result[1], result[0]] = (outdata[result[1], result[0]][0], result[2][0])
    return outdata, data[0]


def join_domain_out(jobdict, tmp_path, out_path, bidirectional, outname, seed_spec, groupdict):
    out_f = open(out_path + outname + "_forward.domains", "w")
    out_f.write('# pairID\torthoID\tseqLen\tfeature\tfStart\tfEnd\tfWeight\tfPath\tinterProID\te-value\tbitScore'
                + '\tpStart\tpEnd\tpLen\n')
    out_r = None
    if bidirectional:
        out_r = open(out_path + outname + "_reverse.domains", "w")
        out_r.write('# pairID\torthoID\tseqLen\tfeature\tfStart\tfEnd\tfWeight\tfPath\tinterProID\te-value\tbitScore'
                    + '\tpStart\tpEnd\tpLen\n')
    for spec in jobdict:
        with open(tmp_path + "/" + spec + "_forward.domains", "r") as infile:
            for line in infile.readlines():
                if not line[0] == '#':
                    cells = line.split("\t")
                    s_id, q_id = cells[0].split("#")
                    for seed in groupdict[s_id]:
                        if q_id in groupdict[s_id][seed]:
                            if not cells[1] == q_id:
                                p_id = seed_spec + "|" + cells[1]
                            else:
                                p_id = spec + "|" + cells[1] + "|" + groupdict[s_id][seed][q_id]
                            out_f.write(seed + "#" + seed + "|" + spec + "|" + q_id + "|" + groupdict[s_id][seed][q_id]
                                        + "\t" + seed + "|" + p_id + "\t" + "\t".join(cells[2:]))
        os.remove(tmp_path + "/" + spec + "_forward.domains")
        if bidirectional:
            with open(tmp_path + "/" + spec + "_reverse.domains", "r") as infile:
                for line in infile.readlines():
                    if not line[0] == '#':
                        cells = line.split("\t")
                        s_id, q_id = cells[0].split("#")
                        for seed in groupdict[s_id]:
                            if q_id in groupdict[s_id][seed]:
                                if not cells[1] == q_id:
                                    p_id = seed_spec + "|" + cells[1]
                                else:
                                    p_id = spec + "|" + cells[1] + "|" + groupdict[s_id][seed][q_id]
                                out_r.write(seed + "#" + seed + "|" + spec + "|" + q_id + "|"
                                            + groupdict[s_id][seed][q_id] + "\t" + seed + "|" + p_id + "\t"
                                            + "\t".join(cells[2:]))
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
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    parser.add_argument('--version', action='version', version=str(version))
    required.add_argument("-i", "--extended_fa", default=None, type=str, required=True,
                          help="path to extended.fa file")
    required.add_argument("-w", "--weight_dir", default=None, type=str, required=True,
                          help="path to weight_dir of fdog")
    optional.add_argument("-r", "--redo_anno", action="store_true",
                          help="enable automatic annotation of missing proteins")
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
                          help="deactivate linearization (for pfam/smart)")
    optional.add_argument("-d", "--featuretypes", type=str, default=None,
                          help="inputfile that contains the tools/databases used for comparison. Please look at the "
                               "FAS wiki pages for templates of the the featuretypes input file")
    arguments = parser.parse_args()
    return arguments


if __name__ == '__main__':
    main()
