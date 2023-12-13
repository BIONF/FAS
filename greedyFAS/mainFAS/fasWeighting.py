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


import math


def w_weighting(protein, domain_count, proteome, option):
    """Calculates weights

    :param protein: protein id
    :param domain_count: feature counts of the reference
    :param proteome: seed or query proteome
    :return: weights, domain_count
    """
    weights = {}
    features = []
    scaling_factor = 0.0
    sum_of_features = 0.0
    for tool in option["input_linearized"]:
        for feature in proteome[protein][tool]:
            features.append(feature)
    for tool in option["input_normal"]:
        for feature in proteome[protein][tool]:
            features.append(feature)
    for feature in features:
        try:
            domain_count[feature]
        except KeyError:
            # print "feature not found in ref " + feature
            domain_count[feature] = 1
            sum_of_features += float(domain_count[feature])
        else:
            sum_of_features += float(domain_count[feature])
    for feature in features:
        weights[feature] = round(float(sum_of_features) / float(domain_count[feature]), 4)
        scaling_factor += round(float(sum_of_features) / float(domain_count[feature]), 4)
    for feature in features:
        weights[feature] = round(float(weights[feature]) / float(scaling_factor), 4)
    return weights, domain_count


def w_count_add_domains(protein, domain_count, proteome):
    features = []
    if domain_count == 'NA':
        return 0
    for feature in proteome[protein]:
        features.append(feature)
    for feature in features:
        try:
            domain_count[feature]
        except KeyError:
            # print "feature not found in ref " + feature
            domain_count[feature] = 1
    return domain_count


def w_weighting_constraints(protein, domain_count, proteome, option):
    """Calculates weights, while fulfilling the constraints

    :param protein: protein id
    :param domain_count: feature counts of the reference
    :param proteome: seed or query proteome
    :param option: dictionary that contains the main option variables of FAS
    :return: weights, domain_count
    """
    weights = {}
    tools = {}
    for tool in (option["input_linearized"] + option["input_normal"]):
        tools[tool] = []
    features = []
    single_constraints = []
    filled = 0.0
    for tool in tools:
        for feature in proteome[protein][tool]:
            features.append(feature)
    for feature in features:
        if feature in option["constraints"]:
            filled += option["constraints"][feature]
            weights[feature] = option["constraints"][feature]
            single_constraints.append(feature)
        elif feature.split('_')[0] in option["constraints"]:
            tools[feature.split('_')[0]].append(feature)
    for tool in tools:
        if len(tools[tool]) > 0:
            filled += option["constraints"][tool]
            sum_of_features = 0
            scaling_factor = 0.0
            for feature in tools[tool]:
                features.remove(feature)
                try:
                    domain_count[feature]
                except KeyError:
                    domain_count[feature] = 1
                    sum_of_features += float(domain_count[feature])
                else:
                    sum_of_features += float(domain_count[feature])
            for feature in tools[tool]:
                weights[feature] = round(float(sum_of_features) / float(domain_count[feature]), 4)
                scaling_factor += round(float(sum_of_features) / float(domain_count[feature]), 4)
            for feature in tools[tool]:
                weights[feature] = round(float(weights[feature]) / float(scaling_factor) * option["constraints"][tool],
                                         4)
    for feature in single_constraints:
        features.remove(feature)
    sum_of_features = 0.0
    for feature in features:
        try:
            domain_count[feature]
        except KeyError:
            # print "feature not found in ref " + feature
            domain_count[feature] = 1
            sum_of_features += float(domain_count[feature])
        else:
            sum_of_features += float(domain_count[feature])
    scaling_factor = 0.0
    for feature in features:
        weights[feature] = round(float(sum_of_features) / float(domain_count[feature]), 4)
        scaling_factor += round(float(sum_of_features) / float(domain_count[feature]), 4)
    for feature in features:
        weights[feature] = round(float(weights[feature]) / float(scaling_factor) * (1.0 - filled), 4)
    return weights, domain_count


def w_weight_correction(method, domain_count):
    """Rescales counts by applying one of the 4 functions (loge, log10, root4, root8)

    :param method: rescaling-function
    :param domain_count: feature counts of the reference
    :return: domain_count
    """
    if method == "loge":
        for feature in domain_count:
            domain_count[feature] = int(round(math.log(domain_count[feature]), 0) + 1)
    elif method == "log10":
        for feature in domain_count:
            domain_count[feature] = int(round(math.log10(domain_count[feature]), 0) + 1)
    elif method == "root4":
        for feature in domain_count:
            domain_count[feature] = int(round(math.pow(domain_count[feature], 0.25), 0))
    elif method == "root8":
        for feature in domain_count:
            domain_count[feature] = int(round(math.pow(domain_count[feature], 0.125), 0))
    return domain_count


def w_weight_const_rescale(path, weights, search_features, final, option):
    """Rescales weights to 1 after linearization

    :param path: Path to rescale
    :param weights: feature weights
    :param search_features: all features to be linearized in the seed protein
    :param final: final calculation? (1/0)
    :param option: dictionary that contains the main option variables of FAS
    :return: rescaled_weights
    """
    lindict = {}
    for ftype in option["input_linearized"]:
        lindict[ftype] = []
    if final:
        for ftype in option["input_normal"]:
            lindict[ftype] = []
    tmp = 0.0
    rescaled_weights = {}
    if final:
        for feature in path:
            db = feature.split('_')[0]
            if db == 'coils':
                db = 'coils2'
            if feature not in option["constraints"]:
                if feature not in lindict[db]:
                    lindict[db].append(feature)
        for ftype in option["input_linearized"]:
            if ftype in option["constraints"]:
                for feature in lindict[ftype]:
                    tmp += weights[feature]
                if option["constraints"][ftype] > tmp > 0.0:
                    scale = option["constraints"][ftype] / tmp
                    for feature in lindict[ftype]:
                        rescaled_weights[feature] = weights[feature] * scale
                tmp = 0.0
    else:
        for feature in path:
            if feature in search_features:
                if search_features[feature][0] not in option["constraints"]:
                    if feature not in lindict[search_features[feature][0].split('_')[0]]:
                        lindict[search_features[feature][0].split('_')[0]].append(feature)

        for ftype in option["input_linearized"]:
            if ftype in option["constraints"]:
                for feature in lindict[ftype]:
                    tmp += weights[search_features[feature][0]]
                if option["constraints"][ftype] > tmp > 0.0:
                    scale = option["constraints"][ftype] / tmp
                    for feature in lindict[ftype]:
                        rescaled_weights[search_features[feature][0]] = weights[search_features[feature][0]] * scale
                tmp = 0.0

    return rescaled_weights
