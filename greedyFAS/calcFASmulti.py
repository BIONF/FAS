#!/bin/env python

#######################################################################
#  Copyright (C) 2023 Vinh Tran
#
#  This file is part of FAS. It is used to calculate FAS scores for
#  input that contains pairwise proteins from different pair of taxa
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


import os
import sys
import argparse
import multiprocessing as mp
from sys import argv
from greedyFAS.annoFAS import annoFAS
from greedyFAS.annoFAS import annoModules
from greedyFAS.mainFAS import fasInput, fasOutput, greedyFAS
from greedyFAS import calcFAS
from pkg_resources import get_distribution
import multiprocessing as mp
import subprocess
from tqdm import tqdm
import json
import shutil
import glob
from datetime import date
# import linecache
import random


def get_options():
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    general = parser.add_argument_group('general arguments')
    inargs = parser.add_argument_group('input arguments')
    outargs = parser.add_argument_group('output arguments')
    annotation = parser.add_argument_group('annotation arguments')
    weighting = parser.add_argument_group('weighting arguments')
    thresholds = parser.add_argument_group('threshold arguments')
    obscure = parser.add_argument_group('obscure arguments')
    other = parser.add_argument_group('other arguments')
    parser.add_argument('--version', action='version', version=str(version))
    required.add_argument("--input", default=None, type=str, required=True,
                          help="Input file (see https://github.com/BIONF/FAS/wiki/File-Formats)")
    required.add_argument("-a", "--annotation_dir", default=None, type=str, required=True,
                          help='working directory, all annotations are be stored here')
    required.add_argument("-o", "--out_dir", default=None, type=str, required=True,
                          help="output directory, all outputfiles will be stored here")
    general.add_argument("--bidirectional", action="store_true",
                         help="calculate both scoring directions")
    general.add_argument('--cpus', help='number of cores', type=int, action='store', default=0)
    annotation.add_argument("-f", "--eFeature", default="0.001", type=float,
                            help="eValue cutoff for PFAM/SMART domain, applied during annotation but also during "
                                 "calculation")
    annotation.add_argument("-i", "--eInstance", default="0.01", type=float,
                            help="eValue cutoff for PFAM/SMART instances, applied during annotation but also during "
                                 "calculation")
    annotation.add_argument('--eFlps', help='eValue cutoff for fLPS', action='store', default=0.0000001, type=float)
    annotation.add_argument('--org', help='Organism of input sequence(s) for SignalP search',
                            choices=['euk', 'gram+', 'gram-'], action='store', default='euk', type=str)
    annotation.add_argument("--toolPath", dest="toolPath", default='', type=str,
                            help="Path to Annotion tool directory created with fas.setup")
    weighting.add_argument("-g", "--weight_correction", default="loge", type=str,
                           help="Function applied to the frequency of feature types during weighting, options are "
                                "linear(no function), loge(natural logarithm[Default]), log10(base-10 logarithm), "
                                "root4(4th root) and root8(8th root).")
    weighting.add_argument("-x", "--weight_constraints", default=None, type=str,
                           help="Apply weight constraints via constraints file, by default there are no constraints. "
                                "Please look at the FAS wiki pages for templates for the constraints file")
    inargs.add_argument("--oldJson", default='', type=str,
                        help="Input old FAS output in JSON format. Only new pairs will be considered!")
    inargs.add_argument("-d", "--featuretypes", default='', type=str,
                        help="inputfile that contains the tools/databases used to predict features. Please look at the "
                             "FAS wiki pages for templates of the the featuretypes input file")
    inargs.add_argument("--extra_annotation", default=None, nargs='*', type=str,
                        help="give naming conventions for extra annotation files, eg. if the file name for the seed is "
                             "someseed.json than the extra annotation files should be named someseed_"
                             "[EXTRA_ANNOTATION].json, these extra files need to exist for seed, query and the "
                             "references (if given)")
    inargs.add_argument("--pairLimit", default=0, type=int,
                        help="Limit the number of pairs to be considered")
    outargs.add_argument("--raw", dest="raw", action="store_true",
                         help="print FAS scores to terminal")
    outargs.add_argument("--tsv", dest="tsv", action="store_true",
                         help="deactivates creation of the tsv output, either use together with --raw or "
                              "--phyloprofile")
    outargs.add_argument("--json", dest="json", action="store_true",
                         help="create json output")
    outargs.add_argument("--mergeJson", action="store_true",
                        help="merge multiple json outputs")
    outargs.add_argument("--outName", default=None, type=str,
                        help="name for merged output file, if none is given the name will be created by the date")
    outargs.add_argument("--phyloprofile", dest="phyloprofile", default=None, type=str,
                         help="activate phyloprofile output, needs mapping file for all query proteins")
    outargs.add_argument("--domain", dest="domain", action="store_false",
                         help="deactivate .domains output")
    outargs.add_argument("--no_config", dest="config", action="store_false",
                         help="deactivate *_config.yml output")
    outargs.add_argument('--force', help='Force override output files', action='store_true')
    outargs.add_argument('--silentOff', help='Turn on STDOUT (not recommended for parallel run)', action='store_true')
    outargs.add_argument('--keep', help='Keep split input files', action='store_true')
    thresholds.add_argument("-c", "--max_overlap", dest="max_overlap", default=0, type=int,
                            help="maximum size overlap allowed, default is 0 amino acids")
    thresholds.add_argument("--max_overlap_percentage", dest="max_overlap_percentage", default=0.4, type=float,
                            help="defines how much percent of a feature the overlap is allowed to cover, default "
                                 "is 0.4 (40%%)")
    thresholds.add_argument("-t", "--priority_threshold", type=int, default=30,
                            help="Change to define the feature number threshold for activating priority mode in the "
                                 "path evaluation. default=30")
    thresholds.add_argument("-m", "--max_cardinality", default=500, type=int,
                            help="Change to define the threshold for the maximal cardinality (number) of feature paths "
                                 "in a graph. If max. cardinality is exceeded the priority mode will be used to for "
                                 "the path evaluation. default=500")
    obscure.add_argument("--priority_mode", action='store_false',
                         help="deactivates the greedy strategy priority mode for larger architectures, NOT RECOMMENDED")
    obscure.add_argument("--timelimit", default=3600, type=int,
                         help="Sets a soft time limit in seconds for the calculation between a pair of proteins,"
                              "This limit does not necessarily represent the actual runtime. It only stops the "
                              "exhaustive strategy through the architecture and activate priority mode. "
                              "This option is only relevant if the option [--priority_mode] is set")
    obscure.add_argument("-w", "--score_weights", nargs=3, default=[0.7, 0.0, 0.3], type=float,
                         help="Defines how the three scores MS, CS and PS are weighted, takes three float arguments, "
                              "sum should be 1.0, the default is 0.7, 0.0, 0.3")
    obscure.add_argument("--empty_as_1", action='store_true',
                         help="If both proteins have no features, score 1.0 instead of the default 0.0")
    other.add_argument("--check_anno", action='store_true',
                         help="Checking for existing annotations")
    args = parser.parse_args()
    return args

def check_anno(in_file, annotation_dir):
    # get list of taxa and their protein IDs
    input_dict = {}
    with open(in_file, 'r') as fr:
        for line in fr:
            tmp = line.strip().split()
            if len(tmp) == 4:
                if not tmp[1] in input_dict:
                    input_dict[tmp[1]] = []
                input_dict[tmp[1]].append(tmp[0])
                if not tmp[3] in input_dict:
                    input_dict[tmp[3]] = []
                input_dict[tmp[3]].append(tmp[2])
    # check missing anno
    missing_dict = {}
    for taxon in input_dict:
        with open(f'{annotation_dir}/{taxon}.json') as f:
            anno_dict = json.load(f)
            for id in input_dict[taxon]:
                if not id in anno_dict['feature']:
                    if not taxon in missing_dict:
                        missing_dict[taxon] = []
                    missing_dict[taxon].append(id)
    return(missing_dict)


def create_jobs(in_file, args, annotation_dir, out_dir, toolpath):
    tax_dict = {}
    os.makedirs(f'{out_dir}/split_inputs', exist_ok=True)
    with open(in_file, 'r') as fr:
        lines = fr.read().splitlines()
        if args.pairLimit > 0:
            if args.silentOff:
                print(f"NOTE: only random {args.pairLimit} pairs will be calculated!")
            if args.pairLimit < len(lines):
                lines = random.sample(lines, args.pairLimit)
        for line in lines:
            tmp = line.strip().split()
            if len(tmp) == 4:
                taxa_pair = f'{tmp[1]}#{tmp[3]}'
                if taxa_pair not in tax_dict:
                    fp = open(f'{out_dir}/split_inputs/{taxa_pair}.txt', 'w')
                    fp.write(f'{tmp[0]}\t{tmp[2]}\n')
                    fp.close()
                    tax_dict[taxa_pair] = [f'{tmp[0]}\t{tmp[2]}']
                else:
                    fp = open(f'{out_dir}/split_inputs/{taxa_pair}.txt', 'a')
                    fp.write(f'{tmp[0]}\t{tmp[2]}\n')
                    fp.close()
                    tax_dict[taxa_pair].append(f'{tmp[0]}\t{tmp[2]}')
    jobs = []
    for tax_pair in tax_dict:
        args_dict_new = {}
        args_dict = vars(args)
        for item in args_dict:
            args_dict_new[item] = args_dict[item]
        args_m = argparse.Namespace(**args_dict_new)
        tmp = tax_pair.split('#')
        args_m.seed = f'{tmp[0]}.json'
        args_m.query = f'{tmp[1]}.json'
        args_m.pairwise = f'{out_dir}/split_inputs/{tax_pair}.txt'
        args_m.out_name = tax_pair
        args_m.query_id = None
        args_m.seed_id = None
        args_m.ref_proteome = None
        args_m.ref_2 = None
        if args.silentOff:
            args_m.silent = False
        else:
            args_m.silent = True
        # check existing output files
        out_file = f'{args_m.out_dir}/{args_m.out_name}.tsv'
        if args.json:
            out_file = f'{args_m.out_dir}/{args_m.out_name}.json'
        if os.path.exists(out_file):
            if not args.force:
                print(f'Output file {os.path.abspath(out_file)} exists! Use --force if you want to overwrite it!')
                continue
        jobs.append([args_m, toolpath])
    return(jobs)


def main():
    args = get_options()
    toolpath = args.toolPath
    if toolpath == '':
        pathconfigfile = os.path.realpath(__file__).replace('calcFASmulti.py', 'pathconfig.txt')
        with open(pathconfigfile) as f:
            toolpath = f.readline().strip()
    else:
        toolpath = os.path.abspath(args.toolPath)
    if not os.path.isdir(args.annotation_dir):
        sys.exit(f'ERROR: {args.annotation_dir} not found!')
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    if args.force:
        print('WARNING: --force option is used. Old output files will be overwritten!')

    if args.check_anno:
        print('==> checking annotations...')
        missing_dict = check_anno(args.input, args.annotation_dir)
        if len(missing_dict) > 0:
            sys.exit(f'ERROR: Annotations for the following proteins are missing!\n{missing_dict}')

    print('==> preparing jobs...')
    jobs = create_jobs(args.input, args, args.annotation_dir, args.out_dir, toolpath)
    if args.cpus == 0:
        cpus = mp.cpu_count()-1
    else:
        if args.cpus > mp.cpu_count():
            cpus = mp.cpu_count()-1
        else:
            cpus = args.cpus

    print('==> calculating FAS scores...')
    pool = mp.Pool(cpus)
    for _ in tqdm(pool.imap_unordered(calcFAS.fas, jobs), total=len(jobs)):
        pass
    pool.close()

    if args.mergeJson:
        out_name = args.outName
        out_dir = os.path.abspath(args.out_dir)
        if not args.outName:
            today = date.today()
            out_name = f"fas_{today.strftime('%y%m%d')}"
        print(f'==> merge outputs into {out_name}...')
        if os.path.exists(f'{out_dir}/{out_name}.json'):
            if args.force:
                os.remove(f'{out_dir}/{out_name}.json')
        os.makedirs(f'{out_dir}/tmp', exist_ok = True)
        for json_file in glob.glob(os.path.join(out_dir, '*.json')):
            shutil.move(json_file, f"{out_dir}/tmp/{json_file.split('/')[-1]}")
        if args.oldJson:
            if os.path.exists(os.path.abspath(args.oldJson)):
                os.symlink(os.path.abspath(args.oldJson), f'{out_dir}/tmp/oldJson.json')
        cmd = f'fas.mergeJson -i {out_dir}/tmp/ -n {out_name} -o {out_dir}'
        try:
            mergedOut = subprocess.run([cmd], shell=True, capture_output=True, check=True)
        except:
            sys.exit(f'Error running\n{cmd}')

    if not args.keep:
        shutil.rmtree(f'{args.out_dir}/split_inputs')
        if args.mergeJson:
            shutil.rmtree(f'{args.out_dir}/tmp')
    print('==> DONE!')

if __name__ == '__main__':
    main()
