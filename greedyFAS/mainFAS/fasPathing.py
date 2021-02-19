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


def pb_region_paths(overlap_map):
    """Uses overlap information to build a directional, acyclic graph that contains all linear paths

    :param overlap_map: contains overlap information
    :return: graph
    """
    graph = {}
    reached_list = {}
    path_count = {'END': 1}
    for feature in overlap_map:
        reached = []
        links = []
        if not feature[0] == 'END':
            path_count[feature[0]] = 0
        for candidate in feature[1]:
            links.append(candidate)
            path_count[feature[0]] += path_count[candidate]
            reached.append(candidate)
            for i in reached_list[candidate]:
                if i in feature[1]:
                    feature[1].remove(i)
                if i not in reached:
                    reached.append(i)
        reached_list[feature[0]] = reached
        graph[feature[0]] = links
    return graph, path_count['START']


def pb_region_mapper(overlap_region, features, max_overlap, max_overlap_percentage):
    """ Finds all (relevant) overlaps in the feature architecture

    :param max_overlap_percentage: overlap_threshold in percentage
    :param max_overlap: overlap threshold in amino acids
    :param overlap_region: region that is looked at
    :param features: feature information [start/stop]
    :return: overlap_map
    """

    overlap_map = [("START", overlap_region + ['END'])]
    for i in range(0, len(overlap_region)):
        end = features[overlap_region[i]][3]
        length_i = features[overlap_region[i]][3] - features[overlap_region[i]][2]
        length_i = length_i * max_overlap_percentage
        overlap = []
        x = i + 1
        while x < len(overlap_region):
            length_x = features[overlap_region[x]][3] - features[overlap_region[x]][2]
            length_x = length_x * max_overlap_percentage
            overlap_size = end - features[overlap_region[x]][2] + 1
            if overlap_size <= max_overlap and overlap_size <= length_i and overlap_size <= length_x:
                overlap.append(overlap_region[x])
            x += 1
        overlap.append("END")
        overlap_map.append((overlap_region[i], overlap))
    overlap_map.append(("END", []))
    overlap_map.reverse()
    return overlap_map
