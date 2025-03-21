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


import os
import argparse
import multiprocessing as mp
from sys import argv
from greedyFAS.annoFAS import annoFAS
from greedyFAS.annoFAS import annoModules
from greedyFAS.mainFAS import fasInput, fasOutput, fadMain
from greedyFAS.calcFAS import anno
from importlib.metadata import version, PackageNotFoundError

def get_options():
    fas_version = version("greedyFAS")
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(fas_version) + '.',
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
    parser.add_argument('--version', action='version', version=str(version))
    required.add_argument("-s", "--seed", default=None, type=str, required=True,
                          help="path to seed protein fasta file")
    required.add_argument("-q", "--query", default=None, type=str, required=True,
                          help="path to query protein fasta file")
    required.add_argument("-a", "--annotation_dir", default=None, type=str, required=True,
                          help='working directory, all annotations are be stored here')
    required.add_argument("-o", "--out_dir", default=None, type=str, required=True,
                          help="output directory, all outputfiles will be stored here")
    general.add_argument("--bidirectional", action="store_true",
                         help="calculate both scoring directions")
    general.add_argument('--cpus', help='number of cores', type=int, action='store', default=0)
    annotation.add_argument('--force', help='Force override annotations', action='store_true')
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
    weighting.add_argument("-r", "--weights", default=None, type=str,
                           help="Path to file containing weighting. [json] format")
    inargs.add_argument("--query_id", default=None, nargs='*', type=str,
                        help="Choose specific proteins (ids divided by spaces) from the query input for calculation, "
                             "by default this is off (all proteins are used for calculation)")
    inargs.add_argument("--seed_id", default=None, nargs='*', type=str,
                        help="Choose specific proteins (ids divided by spaces) from the seed input for calculation, "
                             "by default this is off (all proteins are used for calculation)")
    inargs.add_argument("--pairwise", dest="pairwise", default=None, type=str,
                        help="deactivate all against all comparison, needs a pairing file with the ids that should be"
                             " compared (one pair per line tab seperated), please look at the FAS wiki pages for "
                             "templates")
    inargs.add_argument("-d", "--featuretypes", default=None, type=str,
                        help="inputfile that contains the tools/databases used to predict features. Please look at the "
                             "FAS wiki pages for templates of the the featuretypes input file")
    inargs.add_argument("--extra_annotation", default=None, nargs='*', type=str,
                        help="give naming conventions for extra annotation files, eg. if the file name for the seed is "
                             "someseed.json than the extra annotation files should be named someseed_"
                             "[EXTRA_ANNOTATION].json, these extra files need to exist for seed, query and the "
                             "references (if given)")
    outargs.add_argument("-n", "--out_name", default=None, type=str,
                         help="name for outputfiles, if none is given the name will be created from the seed and "
                              "query names")
    outargs.add_argument("--raw", dest="raw", action="store_true",
                         help="print FAS scores to terminal")
    outargs.add_argument("--tsv", dest="tsv", action="store_true",
                         help="deactivates creation of the tsv output, either use together with --raw or "
                              "--phyloprofile")
    outargs.add_argument("--phyloprofile", dest="phyloprofile", default=None, type=str,
                         help="activate phyloprofile output, needs mapping file for all query proteins")
    outargs.add_argument("--domain", dest="domain", action="store_false",
                         help="deactivate .domains output")
    outargs.add_argument("--no_config", dest="config", action="store_false",
                         help="deactivate *_config.yml output")
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
    args = parser.parse_args()
    return args


def fad(args, toolpath):
    option_dict = {"seed_id": args.seed_id, "query_id": args.query_id, "priority_mode": args.priority_mode,
                   "priority_threshold": args.priority_threshold, "max_cardinality": args.max_cardinality,
                   "cores": int(args.cpus), "raw": args.raw, "bidirectional": args.bidirectional,
                   "max_overlap": args.max_overlap, "timelimit": args.timelimit, "phyloprofile": args.phyloprofile,
                   "score_weights": [], "tsv": args.tsv, "max_overlap_percentage": 0.0, "domain": args.domain,
                   "pairwise": None, "eInstance": args.eInstance, "eFeature": args.eFeature, "progress": True,
                   "weight_path": args.weights}

    seedname = '.'.join(args.seed.split('/')[-1].split('.')[:-1])
    option_dict["p_path"] = [args.annotation_dir + '/' + seedname + '.json']
    queryname = '.'.join(args.query.split('/')[-1].split('.')[:-1])
    option_dict["s_path"] = [args.annotation_dir + '/' + queryname + '.json']
    if args.raw:
        option_dict["progress"] = False
    if option_dict["cores"] == 0:
        option_dict["cores"] = mp.cpu_count()-1
    if args.extra_annotation:
        for i in args.extra_annotation:
            option_dict["p_path"].append(args.annotation_dir + '/' + i + '.json')
            option_dict["s_path"].append(args.annotation_dir + '/' + i + '.json')
    try:
        test = 0.0
        for weight in args.score_weights:
            option_dict["score_weights"].append(weight)
            test += float(weight)
        if test != 1.0:
            raise Exception("The sum of all values in -w [--score_weights] should be 1.0")
        option_dict["score_weights"] = tuple(option_dict["score_weights"])
    except ValueError:
        raise Exception(str(args.score_weights) + " is not a valid input for -w [--score_weights]")
    if 0.0 <= float(args.max_overlap_percentage) <= 1.0:
        option_dict["max_overlap_percentage"] = float(args.max_overlap_percentage)
    else:
        raise Exception("[--max_overlap_percentage] should be between 0.0 and 1.0")
    if args.out_name:
        option_dict['outpath'] = args.out_dir.rstrip('/') + '/' + args.out_name
    else:
        option_dict['outpath'] = args.out_dir.rstrip('/') + '/' + seedname + '_' + queryname
    if args.featuretypes:
        option_dict['input_linearized'], option_dict['input_normal'] = fasInput.featuretypes(args.featuretypes)
    else:
        option_dict['input_linearized'], option_dict['input_normal'] = fasInput.featuretypes(toolpath
                                                                                             + '/' + 'annoTools.txt')
    if args.pairwise:
        option_dict["pairwise"] = fasInput.read_pairwise(args.pairwise)
    else:
        option_dict["pairwise"] = None
    print('Calculating FAD score...')
    fadMain.fc_start(option_dict)
    if args.config:
        fasOutput.write_metadata(option_dict['outpath'] + '_config.yml', args, ' '.join(argv),
                                 str(version("greedyFAS")))
    print('done!')


def main():
    args = get_options()
    toolpath = args.toolPath
    if toolpath == '':
        pathconfigfile = os.path.realpath(__file__).replace('calcFAD.py', 'pathconfig.txt')
        with open(pathconfigfile) as f:
            toolpath = f.readline().strip()
    else:
        toolpath = os.path.abspath(args.toolPath)
    if not os.path.isdir(args.annotation_dir):
        os.mkdir(args.annotation_dir)
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    annojobs = [args.seed]
    if args.query not in annojobs:
        annojobs.append(args.query)
    anno(annojobs, args, toolpath)
    fad(args, toolpath)


if __name__ == '__main__':
    main()
