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
import sys
import argparse
import multiprocessing as mp
from sys import argv
from greedyFAS.annoFAS import annoFAS
from greedyFAS.annoFAS import annoModules
from greedyFAS.mainFAS import fasInput, fasOutput, greedyFAS
from greedyFAS.mainFAS.fasInput import read_json
from pkg_resources import get_distribution


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
    annotation.add_argument('--forceAnno', help='Force override annotations', action='store_true')
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
    weighting.add_argument("-r", "--ref_proteome", default=None, type=str,
                           help="Path to a reference proteome which can be used for the weighting of features, "
                                "by default there is no reference proteome used, the weighting will be uniform")
    weighting.add_argument("--ref_2", dest="ref_2", default=None, type=str,
                           help="Give a second reference for bidirectional mode, does not do anything if bidirectional "
                                "mode is not active or if no main reference was given")
    weighting.add_argument("-g", "--weight_correction", default="loge", type=str,
                           help="Function applied to the frequency of feature types during weighting, options are "
                                "linear(no function), loge(natural logarithm[Default]), log10(base-10 logarithm), "
                                "root4(4th root) and root8(8th root).")
    weighting.add_argument("-x", "--weight_constraints", default=None, type=str,
                           help="Apply weight constraints via constraints file, by default there are no constraints. "
                                "Please look at the FAS wiki pages for templates for the constraints file")
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
    outargs.add_argument("-n", "--out_name", default=None, type=str,
                         help="name for outputfiles, if none is given the name will be created from the seed and "
                              "query names")
    outargs.add_argument("--raw", dest="raw", action="store_true",
                         help="print FAS scores to terminal")
    outargs.add_argument("--tsv", dest="tsv", action="store_true",
                         help="deactivates creation of the tsv output, either use together with --raw or "
                              "--phyloprofile")
    outargs.add_argument("--json", dest="json", action="store_true",
                         help="create json output")
    outargs.add_argument("--phyloprofile", dest="phyloprofile", default=None, type=str,
                         help="activate phyloprofile output, needs mapping file for all query proteins")
    outargs.add_argument("--domain", dest="domain", action="store_false",
                         help="deactivate .domains output")
    outargs.add_argument("--no_config", dest="config", action="store_false",
                         help="deactivate *_config.yml output")
    outargs.add_argument('--force', help='Force override output files', action='store_true')
    outargs.add_argument('--silent', help='Turn off STDOUT', action='store_true')
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
    thresholds.add_argument("--paths_limit", default=0, type=int,
                            help="Specify number of maximum paths to be considered (10^n). If this threshold is exceeded, "
                                "the corresponding protein will be ignored. Default: 0 for no limit")
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
    args = parser.parse_args()
    return args


def anno(annojobs, args, toolpath):
    eflps = args.eFlps
    signalporg = args.org
    efeature = args.eFeature
    einstance = args.eInstance
    hmmcores = 1
    cpus = args.cpus
    if cpus == 0:
        cpus = mp.cpu_count()-1
    for annojob in annojobs:
        name = '.'.join(annojob.split('/')[-1].split('.')[:-1])
        outpath = os.path.abspath(args.annotation_dir)
        annotate = True
        seqfile = annojob
        if not args.silent:
            print('Check annotation for "' + name + '"...')
        if os.path.exists(os.path.abspath(args.annotation_dir + '/' + name + '.json')):
            if not args.silent:
                print('Annotation for "' + name + '" already exists.')
            if args.forceAnno:
                if not args.silent:
                    print('Overwriting...')
                annoModules.checkFileExist(seqfile)
            else:
                annotate = False
        else:
            annoModules.checkFileExist(seqfile)
        if annotate:
            annoFAS.runAnnoFas(
                [seqfile, outpath, toolpath, args.forceAnno, name, eflps, signalporg, efeature, einstance, hmmcores, '',
                 '', '', cpus, args.featuretypes])

def fas(opts):
    (args, toolpath) = opts
    option_dict = {
                   "weight_const": False, "seed_id": args.seed_id, "query_id": args.query_id,
                   "priority_mode": args.priority_mode, "priority_threshold": args.priority_threshold,
                   "max_cardinality": args.max_cardinality, "cores": int(args.cpus), "raw": args.raw,
                   "bidirectional": args.bidirectional, "max_overlap": args.max_overlap,
                   "timelimit": args.timelimit, "phyloprofile": args.phyloprofile, "score_weights": [],
                    "tsv": args.tsv, "json": args.json, "max_overlap_percentage": 0.0, "domain": args.domain, "pairwise": None,
                    "eInstance": args.eInstance, "eFeature": args.eFeature, "progress": True,
                    "empty_as_1": args.empty_as_1, "silent": args.silent, "paths_limit": 10**args.paths_limit
                   }
    seedname = '.'.join(args.seed.split('/')[-1].split('.')[:-1])
    option_dict["p_path"] = [args.annotation_dir + '/' + seedname + '.json']
    queryname = '.'.join(args.query.split('/')[-1].split('.')[:-1])
    option_dict["s_path"] = [args.annotation_dir + '/' + queryname + '.json']
    if args.ref_proteome:
        name = '.'.join(args.ref_proteome.split('/')[-1].split('.')[:-1])
        option_dict["ref_proteome"] = [args.annotation_dir + '/' + name + '.json']
    else:
        option_dict["ref_proteome"] = None
    if args.raw:
        option_dict["progress"] = False
    if args.silent:
        option_dict["progress"] = False
    if option_dict["cores"] == 0:
        option_dict["cores"] = mp.cpu_count()-1
    if args.ref_2:
        r2name = '.'.join(args.ref_2.split('/')[-1].split('.')[:-1])
        option_dict["ref_2"] = [args.annotation_dir + '/' + r2name + '.json']
    else:
        option_dict["ref_2"] = None
    if args.extra_annotation:
        for i in args.extra_annotation:
            option_dict["p_path"].append(args.annotation_dir + '/' + i + '.json')
            option_dict["s_path"].append(args.annotation_dir + '/' + i + '.json')
            if args.ref_proteome:
                option_dict["ref_proteome"].append(args.annotation_dir + '/' + i + '.json')
            if args.ref_2:
                option_dict["ref_2"].append(args.annotation_dir + '/' + i + '.json')
    if not option_dict["ref_proteome"]:
        option_dict["MS_uni"] = 1
        option_dict["weight_correction"] = None
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
    if not option_dict["ref_proteome"]:
        option_dict["MS_uni"] = 1
        option_dict["weight_correction"] = None
    else:
        option_dict["MS_uni"] = 0
        if args.weight_correction == "log10":
            option_dict["weight_correction"] = "log10"
        elif args.weight_correction == "loge":
            option_dict["weight_correction"] = "loge"
        elif args.weight_correction == "root4":
            option_dict["weight_correction"] = "root4"
        elif args.weight_correction == "root8":
            option_dict["weight_correction"] = "root8"
        elif args.weight_correction == "linear":
            option_dict["weight_correction"] = None
        else:
            raise Exception(str(args.weightcorrection) + "is not a valid input for -g [--weightcorrection]")
    if args.ref_proteome and args.weight_constraints:
        option_dict["weight_const"] = True
        option_dict["constraints"] = fasInput.constraints_in(args.weight_constraints)
    elif args.weight_constraints:
        raise Exception("[--weight_constraints] only works with a reference proteome")
    if 0.0 <= float(args.max_overlap_percentage) <= 1.0:
        option_dict["max_overlap_percentage"] = float(args.max_overlap_percentage)
    else:
        raise Exception("[--max_overlap_percentage] should be between 0.0 and 1.0")
    if args.out_name:
        option_dict['outpath'] = args.out_dir.rstrip('/') + '/' + args.out_name
        out_file = f'{option_dict["outpath"]}.tsv'
        if args.json:
            out_file = f'{option_dict["outpath"]}.json'
        if os.path.exists(out_file):
            if not args.force:
                sys.exit(f'Output file {os.path.abspath(out_file)} exists! Use --force if you want to overwrite it!')
            else:
                os.remove(f'{os.path.abspath(out_file)}')
    else:
        option_dict['outpath'] = args.out_dir.rstrip('/') + '/' + seedname + '_' + queryname
    if args.featuretypes:
        option_dict['input_linearized'], option_dict['input_normal'] = fasInput.featuretypes(args.featuretypes)
    else:
        option_dict['input_linearized'], option_dict['input_normal'] = fasInput.featuretypes(toolpath
                                                                                             + '/' + 'annoTools.txt')

    option_dict["old_json"] = False
    if args.oldJson:
        if os.path.exists(os.path.abspath(args.oldJson)):
            if not args.silent:
                print(f'### NOTE: existing output given ({os.path.abspath(args.oldJson)}). Only new pairs of proteins will be considered!')
            option_dict["old_json"] = os.path.abspath(args.oldJson)
        else:
            print(f'WARNING: {args.oldJson} not found!')

    if args.pairwise:
        option_dict["pairwise"] = fasInput.read_pairwise(args.pairwise)
        if args.oldJson:
            if os.path.exists(os.path.abspath(args.oldJson)):
                option_dict["pairwise"] = filter_oldJson(args.oldJson, fasInput.read_pairwise(args.pairwise))
        if len(option_dict['pairwise']) == 0:
            sys.exit(f'DONE! All pairwise FAS scores can be found in {args.oldJson}!')
    else:
        option_dict["pairwise"] = None

    if not args.silent:
        print('Calculating FAS score...')
    greedyFAS.fc_start(option_dict)
    if args.config:
        fasOutput.write_metadata(option_dict['outpath'] + '_config.yml', args, ' '.join(argv),
                                 str(get_distribution('greedyFAS').version))
    if not args.silent:
        print('done!')


def filter_oldJson(old_json, in_pairs):
    """ Function to check if a FAS scores for input protein pair already exists
    Return a list of new pairs only
    """
    old_results = read_json(old_json)
    new_pairs = []
    for pair in in_pairs:
        if not f'{pair[0]}_{pair[1]}' in old_results and not f'{pair[1]}_{pair[0]}' in old_results:
            new_pairs.append(pair)
    return(new_pairs)


def main():
    args = get_options()
    toolpath = args.toolPath
    if toolpath == '':
        pathconfigfile = os.path.realpath(__file__).replace('calcFAS.py', 'pathconfig.txt')
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
    if args.ref_proteome and args.ref_proteome not in annojobs:
        annojobs.append(args.ref_proteome)
    if args.ref_2 and args.ref_2 not in annojobs:
        annojobs.append(args.ref_2)
    anno(annojobs, args, toolpath)
    fas((args, toolpath))


if __name__ == '__main__':
    main()
