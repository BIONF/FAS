#!/bin/env python

#######################################################################
# Copyright (C) 2023 Julian Dosch
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
import json
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from pkg_resources import get_distribution


def get_stats(annotation_path, tool1, tool2, threshold, svg):
    with open(annotation_path, 'r') as infile:
        in_dict = json.loads(infile.read())
    statistics = count_overlaps(in_dict['feature'], tool1, tool2, threshold)
    fig = make_subplots(rows=3, cols=2,
                        column_titles=[tool1, tool2],
                        row_titles=['Domains', 'Domain Instances (overlap)', 'Domain Instances (all)'],
                        specs=[[{'type': 'domain'}, {'type': 'domain'}],
                               [{'type': 'domain'}, {'type': 'domain'}],
                               [{'type': 'domain'}, {'type': 'domain'}]])
    fig.add_trace(
        go.Pie(labels=['Overlapping', 'Overlap Free'],
               values=[statistics[0], statistics[1]-statistics[0]]),
        row=1, col=1
    )
    fig.add_trace(
        go.Pie(labels=['Overlapping', 'Overlap Free'],
               values=[statistics[4], statistics[5]-statistics[4]]),
        row=1, col=2
    )
    fig.add_trace(
        go.Pie(labels=['Overlapping', 'Overlap Free'],
               values=[statistics[8], statistics[9]-statistics[8]]),
        row=2, col=1
    )
    fig.add_trace(
        go.Pie(labels=['Overlapping', 'Overlap Free'],
               values=[statistics[10], statistics[11]-statistics[10]]),
        row=2, col=2
    )
    fig.add_trace(
        go.Pie(labels=['Overlapping', 'Overlap Free'],
               values=[statistics[2], statistics[3]-statistics[2]]),
        row=3, col=1
    )
    fig.add_trace(
        go.Pie(labels=['Overlapping', 'Overlap Free'],
               values=[statistics[6], statistics[7]-statistics[6]]),
        row=3, col=2
    )
    fig.show()
    if svg:
        out = '.'.join(annotation_path.split('.')[0:-1])
        fig.write_image(file=out + "_Overlap_Statistics.svg", width=800, height=800)


def count_overlaps(proteome, t1, t2, threshold):
    t1_f_o = {}  # dictionary of features types with number of overlaps t1
    t1_f_a = {}  # set of all features in t2
    t1_o_total = 0  # total number of feature instances that have at least one overlap in t1
    t1_f_total = 0  # total number of feature instance in t1
    t1_rel_total = 0  # number of feature instances that have at least one overlap out of overlapping types
    t1_rel_a_total = 0  # number of feature instances out of overlapping types
    t2_f_o = {}
    t2_f_a = {}
    t2_f_total = 0
    t2_o_total = 0
    t2_rel_total = 0
    t2_rel_a_total = 0
    for protein in proteome:
        t1_f_a, t1_f_total, t1_f_o, t1_o_total = count_overlaps_sub(proteome[protein], t1_f_a, t1_f_total, t1_f_o,
                                                                    t1_o_total, t1, t2, threshold)
        t2_f_a, t2_f_total, t2_f_o, t2_o_total = count_overlaps_sub(proteome[protein], t2_f_a, t2_f_total, t2_f_o,
                                                                    t2_o_total, t2, t1, threshold)
    for feature in t1_f_o:
        t1_rel_total += t1_f_o[feature]
        t1_rel_a_total += t1_f_a[feature]
    for feature in t2_f_o:
        t2_rel_total += t2_f_o[feature]
        t2_rel_a_total += t2_f_a[feature]
    return (len(t1_f_o), len(t1_f_a), t1_o_total, t1_f_total, len(t2_f_o), len(t2_f_a), t2_o_total, t2_f_total,
            t1_rel_total, t1_rel_a_total, t2_rel_total, t2_rel_a_total)


def count_overlaps_sub(arc, f_a, f_total, f_o, o_total, t1, t2, threshold):
    for feature in arc[t1]:
        if feature not in f_a:
            f_a[feature] = 0
        for instance in arc[t1][feature]['instance']:
            f_a[feature] += 1
            f_total += 1
            start, stop = instance[0], instance[1]
            length = stop - start
            if count_overlaps_sub2(arc[t2], start, stop, length, threshold):
                if feature not in f_o:
                    f_o[feature] = 0
                f_o[feature] += 1
                o_total += 1
    return f_a, f_total, f_o, o_total


def count_overlaps_sub2(arc, start, stop, length, threshold):
    back = False
    for feature in arc:
        for instance in arc[feature]['instance']:
            length2 = instance[1] - instance[0]
            if instance[0] <= start <= instance[1]:
                o_length = instance[1] - start
                if (length2 * threshold) <= o_length >= (length * threshold):
                    back = True
            elif start <= instance[0] <= stop:
                o_length = stop - instance[0]
                if (length2 * threshold) <= o_length >= (length * threshold):
                    back = True
    return back


def main():
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-i", "--input", default=None, type=str, required=True,
                          help="path to input annotation file [.json]")
    optional.add_argument("-a", "--class1", default='pfam', type=str, required=False,
                          help="Evaluated feature class 1")
    optional.add_argument("-b", "--class2", default='smart', type=str, required=False,
                          help="Evaluated feature class 2")
    optional.add_argument("-t", "--threshold", default=0.8, type=float, required=False,
                          help="Threshold for overlap size in comparison to feature lengths. The overlap must at least "
                               "reach the given fraction [0.0 to 1.0] of both feature lengths to be considered.")
    optional.add_argument("-s", "--svg", action="store_true",
                          help="Save as .svg. Will be saved in same directory as input.")
    args = parser.parse_args()
    get_stats(args.input, args.class1, args.class2, args.threshold, args.svg)


if __name__ == '__main__':
    main()
