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
import logging
import inspect
from sys import version_info
import multiprocessing as mp
if version_info.major == 3:
    from greedyFAS import annoFAS
    from greedyFAS import greedyFAS
    from greedyFAS import annoModules
    from greedyFAS import fasInput
elif version_info.major == 2:
    import annoFAS
    import greedyFAS
    import annoModules
    import fasInput


def get_options():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    version = '1.1.0'
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.')
    parser.add_argument("-s", "--seed", default=None, type=str, required=True,
                          help="path to seed protein fasta file")
    parser.add_argument("-q", "--query", default=None, type=str, required=True,
                          help="path to query protein fasta file")
    parser.add_argument("--annotation_dir", default=None, type=str, required=True,
                          help='working directory, all annotations are be stored here')
    parser.add_argument("-o", "--out_dir", default=None, type=str, required=True,
                          help="output directory, all outputfiles will be stored here")
    parser.add_argument("-n", "--out_name", default=None, type=str,
                          help="name for outputfiles, if none is given the name will be created from the seed and "
                               "query names")
    parser.add_argument("-r", "--ref_proteome", default=None, type=str,
                          help="Path to a reference proteome which can be used for the weighting of features, "
                               "by default there is no reference proteome used, the weighting will be uniform")
    parser.add_argument('--force', help='Force override annotations', action='store_true')
    parser.add_argument("--query_id", default=None, nargs='*', type=str,
                          help="Choose specific proteins from the query input for calculation, by default this is off "
                               "(all proteins are used for calculation)")
    parser.add_argument("--seed_id", default=None, nargs='*', type=str,
                          help="Choose specific proteins from the seed input for calculation, by default this is off "
                               "(all proteins are used for calculation)")
    parser.add_argument("-w", "--score_weights", nargs=3, default=[0.7, 0.0, 0.3], type=float,
                          help="Defines how the three scores MS, CS and PS are weighted, takes three float arguments, "
                               "sum should be 1.0, the default is 0.7, 0.0, 0.3")
    parser.add_argument("-t", "--priority_threshold", type=int, default=30,
                          help="Change to define the feature number threshold for activating priority mode in the path "
                               "evaluation.")
    parser.add_argument("-m", "--max_cardinality", default=5000, type=int,
                          help="Change to define the threshold for the maximal cardinality of feature paths in a graph."
                               " If max. cardinality is exceeded the priority mode will be used to for the path "
                               "evaluation.")
    parser.add_argument("-f", "--eFeature", default="0.001", type=float,
                          help="eValue cutoff for PFAM/SMART domain")
    parser.add_argument("-i", "--eInstance", default="0.01", type=float,
                          help='eValue cutoff for PFAM/SMART instance')
    parser.add_argument("-g", "--weight_correction", default="loge", type=str,
                          help="Function applied to the frequency of feature types during weighting, options are "
                               "linear(no function), loge(natural logarithm[Default]), log10(base-10 logarithm), "
                               "root4(4th root) and root8(8th root).")
    parser.add_argument('--eFlps', help='eValue cutoff for fLPS', action='store', default=0.0000001, type=float)
    parser.add_argument('--org', help='Organism of input sequence(s) for SignalP search',
                          choices=['euk', 'gram+', 'gram-'], action='store', default='euk', type=str)
    parser.add_argument("-x", "--weight_constraints", default=None, type=str,
                          help="Apply weight constraints via constraints file, by default there are no constraints.")
    parser.add_argument("-e", "--no_arch", action="store_false",
                          help="deactivate creation of _architecture file")
    parser.add_argument("--bidirectional", action="store_true",
                          help="calculate both scoring directions (separate files), creates csv file with combined "
                               "scores")
    parser.add_argument("-c", "--max_overlap", dest="max_overlap", default=0, type=int,
                          help="maximum size overlap allowed, default is 0 amino acids")
    parser.add_argument("--max_overlap_percentage", dest="max_overlap_percentage", default=0.4, type=float,
                          help="defines how much percent of a feature the overlap is allowed to cover, default "
                               "is 0.4 (40%%)")
    parser.add_argument("--ref_2", dest="ref_2", default=None, type=str,
                          help="Give a second reference for bidirectional mode, does not do anything if bidirectional "
                               "mode is not active or if no main reference was given")
    parser.add_argument("--extra_annotation", default=None, nargs='*', type=str,
                          help="give naming conventions for extra annotation files, these files should be in the "
                               "corresponding directory in the annotation_dir")
    parser.add_argument('--cpus', help='number of cores', action='store', default=0)
    parser.add_argument("--pairwise", dest="pairwise", default=None, type=str,
                          help="deactivate all against all comparison, needs a pairing file with the ids that should be"
                               " compared (one pair per line tab seperated)")
    parser.add_argument("--toolPath", dest="toolPath", default=None, type=str, required=True,
                          help="Path to Annotion tool directory created with prepareFAS")
    parser.add_argument("-d", "--featuretypes", default=None, type=str,
                        help="inputfile that contains the tools/databases used to predict features")
    parser.add_argument("--phyloprofile", dest="phyloprofile", default=None, type=str,
                        help="activate phyloprofile output, needs mapping file for all query proteins, single seed "
                             "only, will run with more but output won't work without editing")
    parser.add_argument("--domain", dest="domain", action="store_true",
                        help="activate domain tabular output")
    args = parser.parse_args()
    return args


def anno(annojobs, args):
    toolpath = os.path.abspath(args.toolPath)
    eflps = args.eFlps
    signalporg = args.org
    efeature = args.eFeature
    einstance = args.eInstance
    hmmcores = 1
    cpus = args.cpus
    if cpus == 0:
        cpus = mp.cpu_count()-1

    for annojob in annojobs:
        name = ''.join(annojob.split('/')[-1].split('.')[:-1])
        outpath = os.path.abspath(args.annotation_dir + '/' + name)
        annotate = True
        seqfile = annojob
        if os.path.isdir(args.annotation_dir + '/' + name + '/' + name + '.json'):
            print('Annotation for "' + name + '" already exists.')
            if args.force:
                print('Overwriting...')
                annoModules.checkFileExist(seqfile)
            else:
                print('Skipping annotation...')
                annotate = False
        else:
            annoModules.checkFileExist(seqfile)
            print('Annotating "' + name + '"...')
        if annotate:
            annoFAS.runAnnoFas(
                [seqfile, outpath, toolpath, args.force, name, eflps, signalporg, efeature, einstance, hmmcores, '',
                 '', '', cpus])


def fas(args):
    version = "1.1.0"
    loglevel = "ERROR"
    option_dict = {
                   "weight_const": False, "version": version, "seed_id": args.seed_id, "query_id": args.query_id,
                   "priority_mode": True, "priority_threshold": args.priority_threshold,
                   "max_cardinality": args.max_cardinality, "cores": 1, "e_output": args.no_arch,
                    "bidirectional": args.bidirectional, "max_overlap": args.max_overlap, "classicMS": False,
                    "timelimit": 7200, "phyloprofile": args.phyloprofile, "score_weights": [], "output": 0,
                    "max_overlap_percentage": 0.0, "domain": args.domain, "pairwise": None
                   }
    name = ''
    r2name = ''
    seedname = ''.join(args.seed.split('/')[-1].split('.')[:-1])
    option_dict["p_path"] = [args.annotation_dir + '/' + seedname + '/' + seedname + '.json']
    queryname = ''.join(args.query.split('/')[-1].split('.')[:-1])
    option_dict["s_path"] = [args.annotation_dir + '/' + queryname + '/' + queryname + '.json']
    if args.ref_proteome:
        name = ''.join(args.ref_proteome.split('/')[-1].split('.')[:-1])
        option_dict["ref_proteome"] = [args.annotation_dir + '/' + name + '/' + name + '.json']
    else:
        option_dict["ref_proteome"] = None
    if args.ref_2:
        r2name = ''.join(args.ref_2.split('/')[-1].split('.')[:-1])
        option_dict["ref_2"] = [args.annotation_dir + '/' + r2name + '/' + r2name + '.json']
    else:
        option_dict["ref_2"] = None
    if args.extra_annotation:
        for i in args.extra_annotation:
            option_dict["p_path"].append(args.annotation_dir + '/' + seedname + '/' + i + '.json')
            option_dict["s_path"].append(args.annotation_dir + '/' + queryname + '/' + i + '.json')
            if args.ref_proteome:
                option_dict["ref_proteome"].append(args.annotation_dir + '/' + name + '/' + i + '.json')
            if args.ref_2:
                option_dict["ref_2"].append(args.annotation_dir + '/' + r2name + '/' + i + '.json')
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
    elif args.weight_constraints:
        raise Exception("[--weight_constraints] only works with a reference proteome")
    if 0.0 <= float(args.max_overlap_percentage) <= 1.0:
        option_dict["max_overlap_percentage"] = float(args.max_overlap_percentage)
    else:
        raise Exception("[--max_overlap_percentage] should be between 0.0 and 1.0")
    if args.out_name:
        option_dict['outpath'] = args.out_dir.rstrip('/') + '/' + args.out_name
    else:
        option_dict['outpath'] = args.out_dir.rstrip('/') + '/' + seedname + '_' + queryname
    if args.featuretypes is not None:
        option_dict = fasInput.featuretypes(args.featuretypes, option_dict)
    else:
        option_dict["input_linearized"] = ["pfam", "smart"]
        option_dict["input_normal"] = ["flps", "coils2", "seg", "signalp", "tmhmm"]

    if args.pairwise:
        option_dict["pairwise"] = greedyFAS.read_pairwise(args.pairwise)
    else:
        option_dict["pairwise"] = None

    logging.basicConfig(level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info(
        'greedyFAS.py started with options: priority_threshold=' + str(option_dict["priority_threshold"]) +
        ', log_level=' + str(loglevel))
    logging.info(
        'score_weights are set to: ' + str(option_dict["score_weights"][0]) + " " + str(option_dict["score_weights"][1])
        + " " + str(option_dict["score_weights"][2]))
    logging.info('ref_proteome is set to: ' + str(option_dict["ref_proteome"]))
    print('Calculating FAS score...')
    greedyFAS.fc_start(option_dict)
    print('done!')


def main():
    args = get_options()
    if not os.path.isdir(args.annotation_dir):
        os.mkdir(args.annotation_dir)
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    annojobs = [args.seed]
    if args.query not in annojobs:
        annojobs.append(nargs.query)
    if args.ref_proteome and args.ref_proteome not in annojobs:
        annojobs.append(args.ref_proteome)
    if args.ref_2 and args.ref_2 not in annojobs:
        annojobs.append(args.ref_2)
    anno(annojobs, args)
    fas(args)


if __name__ == '__main__':
    main()
