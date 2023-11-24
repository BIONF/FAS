#!/bin/env python

#######################################################################
# Copyright (C) 2019 Julian Dosch
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


import argparse
import os
try:
    from graphviz import Digraph
except ModuleNotFoundError:
    pass
from greedyFAS.mainFAS import fasPathing
from greedyFAS.mainFAS import fasInput
from greedyFAS.mainFAS import greedyFAS
from pkg_resources import get_distribution

def get_options():
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.',
                                     epilog="This script allows you to assess the complexity (number of paths) of the "
                                            "feature architecture for linearization.")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    parser.add_argument('--version', action='version', version=str(version))
    required.add_argument("-i", "--input", default=None, type=str, required=True,
                          help="path to protein architecture json file")
    optional.add_argument("-p", "--id", default=None, nargs='*', type=str,
                          help="Choose specific proteins from the input for the complexity analysis")
    optional.add_argument("-d", "--featuretypes", default=None, type=str,
                          help="inputfile that contains the tools/databases used to predict features. Please look at "
                               "the FAS wiki pages for templates of the the featuretypes input file")
    optional.add_argument("-c", "--max_overlap", dest="max_overlap", default=0, type=int,
                          help="maximum size overlap allowed, default is 0 amino acids")
    optional.add_argument("--max_overlap_percentage", dest="max_overlap_percentage", default=0.4, type=float,
                          help="defines how much percent of a feature the overlap is allowed to cover, default "
                               "is 0.4 (40%%)")
    optional.add_argument("-f", "--eFeature", default="0.001", type=float,
                          help="eValue cutoff for PFAM/SMART domains")
    optional.add_argument("-e", "--eInstance", default="0.01", type=float,
                          help="eValue cutoff for PFAM/SMART instances")
    optional.add_argument("--show_graph", action='store_true',
                          help="visualize the feature graph, you may need to install graphviz for this feature")
    args = parser.parse_args()
    return args


def main():
    args = get_options()
    clan_dict = {}
    pathconfigfile = os.path.realpath(__file__).replace('complexityFAS.py', 'pathconfig.txt')
    option_dict = {}
    with open(pathconfigfile) as f:
        toolpath = f.readline().strip()
    if args.featuretypes is not None:
        option_dict['input_linearized'], option_dict['input_normal'] = fasInput.featuretypes(args.featuretypes)
    else:
        option_dict['input_linearized'], option_dict['input_normal'] = fasInput.featuretypes(toolpath + '/'
                                                                                             + 'annoTools.txt')
    option_dict["max_overlap"] = args.max_overlap
    if 0.0 <= float(args.max_overlap_percentage) <= 1.0:
        option_dict["max_overlap_percentage"] = float(args.max_overlap_percentage)
    else:
        raise Exception("[--max_overlap_percentage] should be between 0.0 and 1.0")
    option_dict["eFeature"] = args.eFeature
    option_dict["eInstance"] = args.eInstance
    proteome = fasInput.read_json(args.input)
    clan_dict.update(proteome["clan"])
    proteome = proteome["feature"]
    if args.id:
        proteins = args.id
    else:
        proteins = list(proteome.keys())
    print('Protein ID\t#Paths\tapproximate greedy complexity')
    for protein in proteins:
        (path_number, greedy_comp, graph) = calc_complex(protein, proteome, option_dict, clan_dict)
        print(f'{protein}\t{path_number}\t{greedy_comp}')
        if args.show_graph:
            show_graph(graph)


def calc_complex(protein, proteome, option, clan_dict):
    lin_set, features, a_f, clans, clan_dict = greedyFAS.su_lin_query_protein(
        protein, proteome, clan_dict, option)
    graph, path_number = fasPathing.pb_region_paths(fasPathing.pb_region_mapper(
        lin_set, features, option["max_overlap"], option["max_overlap_percentage"]))
    greedy_comp = approx_greedy_comp(graph)
    return((path_number, greedy_comp, graph))


def approx_greedy_comp(graph):
    node = 'START'
    complexity = 0
    while not node == 'END':
        if len(graph[node]) > 1:
            complexity += len(graph[node])
        node = graph[node][0]
    return complexity


def show_graph(graph):
    dot = Digraph(comment='feature graph')
    for feature in graph:
        dot.node(feature, feature)
    for feature in graph:
        for f2 in graph[feature]:
            dot.edge(feature, f2)
    dot.render('featuregraph.gv', view=True)


if __name__ == '__main__':
    main()
