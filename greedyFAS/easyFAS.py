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
if version_info.major == 3:
    from greedyFAS import annoFAS
    from greedyFAS import greedyFAS
elif version_info.major == 2:
    import annoFAS
    import greedyFAS


def get_options():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-s", "--seed", default=None, type=str, required=True,
                          help="path to seed protein fasta file")
    required.add_argument("-q", "--query", default=None, type=str, required=True,
                          help="path to query protein fasta file")
    required.add_argument("--projectdir", default=None, type=str, required=True,
                          help='working directory, all annotations and the fas output will be stored here')
    optional.add_argument("-r", "--ref_proteome", default=None, type=str,
                          help="Path to a reference proteome which can be used for the weighting of features, "
                               "by default there is no reference proteome used, the weighting will be uniform")
    optional.add_argument('--force', help='Force override annotations', action='store_true')
    optional.add_argument("--query_id", default=None, nargs='*', type=str,
                          help="Choose a single protein from the query input for calculation, by default this is off "
                               "(all proteins are used for calculation)")
    optional.add_argument("--seed_id", default=None, nargs='*', type=str,
                          help="Choose a single protein from the seed input for calculation, by default this is off "
                               "(all proteins are used for calculation)")
    optional.add_argument("-w", "--score_weights", nargs=3, default=[0.7, 0.0, 0.3], type=float,
                          help="Defines how the three scores MS, CS and PS are weighted, takes three float arguments, "
                               "sum should be 1.0, the default is 0.7, 0.0, 0.3")
    optional.add_argument("-t", "--priority_threshold", type=int, default=50,
                          help="Change to define the feature number threshold for activating priority mode in the path "
                             "evaluation.")
    optional.add_argument("-m", "--max_cardinality", default=5000, type=int,
                          help="Change to define the threshold for the maximal cardinality of feature paths in a graph."
                               " If max. cardinality is exceeded the priority mode will be used to for the path "
                               "evaluation.")
    optional.add_argument("-f", "--efilter", default="0.001", type=float,
                          help="E-value filter for hmm based search methods (feature based/complete sequence).")
    optional.add_argument("-i", "--inst_efilter", default="0.01", type=float,
                          help="E-value filter for hmm based search methods (instance based/complete sequence).")
    optional.add_argument("-g", "--weightcorrection", default="loge", type=str,
                          help="Function applied to the frequency of feature types during weighting, options are "
                               "linear(no function), loge(natural logarithm[Default]), log10(base-10 logarithm), "
                               "root4(4th root) and root8(8th root).")
    optional.add_argument("-x", "--weight_constraints", default=None, type=str,
                          help="Apply weight constraints via constraints file, by default there are no constraints.")
    optional.add_argument("-e", "--no_arch", action="store_false",
                          help="deactivate creation of _architecture file")
    optional.add_argument("--feature_info", action="store_true",
                          help="create a file with information on the abundance of all seed and query features in the "
                               "reference")
    optional.add_argument("--bidirectional", action="store_true",
                          help="calculate both scoring directions (separate files), creates csv file with combined "
                               "scores")
    optional.add_argument("-c", "--max_overlap", dest="max_overlap", default=0, type=int,
                          help="maximum size overlap allowed, default is 0 amino acids")
    optional.add_argument("--max_overlap_percentage", dest="max_overlap_percentage", default=0.4, type=float,
                          help="defines how much percent of a feature the overlap is allowed to cover, default "
                               "is 0.4 (40%%)")
    optional.add_argument("--ref_2", dest="ref_2", default=None, type=str,
                          help="Give a second reference for bidirectional mode, does not do anything if bidirectional "
                               "mode is not active or if no main reference was given")
    optional.add_argument('--cores', help='number of cores', action='store', default='')
    optional.add_argument("--domain", dest="domain", action="store_true", help="activate domain tabular output")
    optional.add_argument("--pairwise", dest="pairwise", default=None, type=str,
                          help="deactivate all against all comparison, needs a pairing file with the ids that should be"
                             " compared (one pair per line tab seperated)")
    args = parser.parse_args()
    return args


def anno(annojobs, projectdir, force, cores):
    anno_options = {'prepare': 'n', 'path': projectdir.rstrip('/') + '/annotations/', 'annoPath': '', 'redo': 0,
                    'extract': '', 'force': force, 'cores': cores}
    for annojob in annojobs:
        name = ''.join(annojob.split('/')[-1].split('.')[:-1])
        annotate = True
        if os.path.isdir(projectdir + '/annotations/' + name):
            print('Annotation for "' + name + '" already exists.')
            if force:
                print('Overwriting...')
            else:
                print('Skipping annotation...')
                annotate = False
        else:
            print('Annotating "' + name + '"...')
        if annotate:
            anno_options['fasta'] = annojob
            anno_options['name'] = name
            annoFAS.easyfas_entry(anno_options)


def fas(args):
    expath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    version = "1.0.0"
    loglevel = "ERROR"
    option_dict = {
                   "weight_const": False, "version": version, "seed_id": args.seed_id, "query_id": args.query_id,
                   "priority_mode": True, "priority_threshold": args.priority_threshold,
                   "max_cardinality": args.max_cardinality, "efilter": args.efilter, "cores": 1,
                   "inst_efilter": args.inst_efilter, "e_output": args.no_arch, "feature_info": args.feature_info,
                   "bidirectional": args.bidirectional, "max_overlap": args.max_overlap, "classicMS": False,
                   "timelimit": 7200, "ref_2": args.ref_2, "phyloprofile": None, "score_weights": [], "output": 0,
                   "max_overlap_percentage": 0.0, "domain": args.domain, "pairwise": None
                   }

    seedname = ''.join(args.seed.split('/')[-1].split('.')[:-1])
    option_dict["p_path"] = args.projectdir.rstrip('/') + '/annotations/' + seedname
    queryname = ''.join(args.query.split('/')[-1].split('.')[:-1])
    option_dict["s_path"] = args.projectdir.rstrip('/') + '/annotations/' + queryname
    if args.ref_proteome:
        name = ''.join(args.ref_proteome.split('/')[-1].split('.')[:-1])
        option_dict["ref_proteome"] = args.projectdir.rstrip('/') + '/annotations/' + name
    else:
        option_dict["ref_proteome"] = None
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
        if args.weightcorrection == "log10":
            option_dict["weight_correction"] = "log10"
        elif args.weightcorrection == "loge":
            option_dict["weight_correction"] = "loge"
        elif args.weightcorrection == "root4":
            option_dict["weight_correction"] = "root4"
        elif args.weightcorrection == "root8":
            option_dict["weight_correction"] = "root8"
        elif args.weightcorrection == "linear":
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

    option_dict['outpath'] = args.projectdir.rstrip('/') + '/out/' + seedname + '_' + queryname
    option_dict["input_linearized"] = ["pfam", "smart"]
    option_dict["input_normal"] = ["flps", "coils", "seg", "signalp", "tmhmm"]

    if options.pairwise:
        option_dict["pairwise"] = greedyFAS.read_pairwise(options.pairwise)
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
    if not os.path.isdir(args.projectdir):
        os.mkdir(args.projectdir)
    if not os.path.isdir(args.projectdir + '/annotations/'):
        os.mkdir(args.projectdir + '/annotations/')
    if not os.path.isdir(args.projectdir + '/out/'):
        os.mkdir(args.projectdir + '/out/')
    annojobs = [args.seed, args.query]
    if args.ref_proteome:
        annojobs.append(args.ref_proteome)
    anno(annojobs, args.projectdir, args.force, args.cores)
    fas(args)


if __name__ == '__main__':
    main()
