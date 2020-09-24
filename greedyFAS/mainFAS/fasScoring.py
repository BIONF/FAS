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


from greedyFAS.mainFAS.fasWeighting import w_weight_const_rescale


def sf_calc_score(path, protein, weights, search_features, query_features, seed_proteome, clan_dict, query_clans,
                  option):
    """Used to be main Scoring function, starts the functions for the scores and calculates the complete FAS score,
    only used during priority mode to make greedy decisions

    :param path: Path to score
    :param protein: protein id
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param search_features: all features to be linearized in the seed protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: (Pfam)-clans of the current query protein
    :param option: dictionary that contains the main option variables of FAS
    :return: score_ms, score_ps, score_cs, final_score
    """
    adjusted_weights = {}
    tmp_weight = {}
    if option["weight_const"]:
        adjusted_weights = w_weight_const_rescale(path, weights, search_features, False, option)
        tmp_weight = {}
        for i in adjusted_weights:
            tmp_weight[i] = weights[i]
            weights[i] = adjusted_weights[i]
    score_cs = round(sf_cs_score(path, clan_dict, query_clans, query_features), 4)
    tmp = sf_ms_score(path, protein, seed_proteome, query_features, option)
    score_ms = round(tmp[0], 4)
    score_ps = round(
        sf_ps_score(path, tmp[2], protein, query_features, seed_proteome, option), 4)
    final_score = (score_ms * option["score_weights"][0]) + (score_cs * option["score_weights"][1]) + (
        score_ps * option["score_weights"][2])
    if option["weight_const"]:
        for i in adjusted_weights:
            weights[i] = tmp_weight[i]
    return score_ms, score_ps, score_cs, final_score


def sf_entire_calc_score(path, query_path, weights, search_features, a_s_f, query_features, a_q_f, clan_dict, option):
    """Main Scoring function, starts the functions for the scores and calculates the complete FAS score

    :param path: Path to score
    :param query_path: Path to score (query)
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param search_features: all features to be linearized in the seed protein
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param option: dictionary that contains the main option variables of FAS
    :return: score_ms, score_ps, score_cs, final_score, path_weight, common_feature, score_ls
    """
    adjusted_weights = {}
    tmp_weight = {}
    if option["weight_const"]:
        adjusted_weights = w_weight_const_rescale(path, weights, search_features, False, option)
        tmp_weight = {}
        for i in adjusted_weights:
            tmp_weight[i] = weights[i]
            weights[i] = adjusted_weights[i]
    score_cs = round(sf_entire_cs_score(path, query_path, query_features, clan_dict, search_features), 4)
    tmp = sf_entire_ms_score(path, query_path, search_features, a_s_f, query_features, a_q_f, weights, option)
    score_ms = round(tmp[0], 4)
    path_weight = tmp[3]
    common_feature = tmp[4]
    ps_tmp = sf_entire_ps_score(path, tmp[2], query_path, search_features, a_s_f, query_features, a_q_f, weights,
                                option)
    score_ps = round(ps_tmp[0], 4)
    score_ls = round(ps_tmp[1], 4)
    final_score = (score_ms * option["score_weights"][0]) + (score_cs * option["score_weights"][1]) + (
        score_ps * option["score_weights"][2])
    if option["weight_const"]:
        for i in adjusted_weights:
            weights[i] = tmp_weight[i]
    return score_ms, score_ps, score_cs, float(final_score), path_weight, common_feature, score_ls


def sf_cs_score(path, clan_dict, query_clans, features):
    """Calculates clan score

    :param path: Path to score
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: Clans in the query
    :param features: features (seed or query) [dict]
    :return: score (PS)
    """

    counter = 0.0
    path_clans = {}
    score = 0.0

    # counting clans in path
    # path_clans: contains counts for clans in path
    for i in path:
        feature = features[i]

        if feature[0] in clan_dict:
            if clan_dict[feature[0]] in path_clans:
                path_clans[clan_dict[feature[0]]] += 1
            else:
                path_clans[clan_dict[feature[0]]] = 1
    for clan in path_clans:
        counter += 1.0
        if clan in query_clans:
            score += float(path_clans[clan] * query_clans[clan]) / float(
                max(path_clans[clan], query_clans[clan]) * max(path_clans[clan], query_clans[clan]))
    if counter == 0:
        score = 0.0
    else:
        score = score / counter
    return float(score)


def sf_entire_cs_score(path, query_path, query_features, clan_dict, search_features):
    """calculates clan score

    :param query_path:
    :param path: Path to score
    :param query_features: all features to be linearized in the query protein
    :param search_features: all features to be linearized in the seed protein
    :param clan_dict: dictionary that maps features to clans
    :return: score (PS)
    """

    counter = 0.0
    s_clans = {}  # clans from path
    q_clans = {}  # clans from query_path

    score = 0.0
    for i in path:
        if i in search_features:
            feature = search_features[i]
            if feature[0] in clan_dict:
                if clan_dict[feature[0]] in s_clans:
                    s_clans[clan_dict[feature[0]]] += 1
                else:
                    s_clans[clan_dict[feature[0]]] = 1
    for i in query_path:
        if i in query_features:
            feature = query_features[i]
            if feature[0] in clan_dict:
                if clan_dict[feature[0]] in q_clans:
                    q_clans[clan_dict[feature[0]]] += 1
                else:
                    q_clans[clan_dict[feature[0]]] = 1

    for clan in s_clans:
        counter += 1.0
        if clan in q_clans:
            score += float(s_clans[clan] * q_clans[clan]) / float(
                max(s_clans[clan], q_clans[clan]) * max(s_clans[clan], q_clans[clan]))
    if counter == 0:
        score = 0.0
    else:
        score = score / counter
    return float(score)


def sf_ms_score(path, protein, proteome, features, option):
    """Calculates multiplicity score, only used for priority mode now

    :param path : Path to score
    :param protein : protein id
    :param proteome : seed_proteome is a dictionary that contains the feature architecture of all seed proteins
    :param features : features (seed or query) [dict]
    :param option : specifies the behavior of the MS calculation
    :return: final_score (MS), search_domains, scale
    """
    domains = {}
    scale = 0
    scores = []
    final_score = 0.0
    tools = option["input_linearized"] + option["input_normal"]
    for i in path:
        feature = features[i]
        if feature[0] in domains:
            domains[feature[0]] += 1
        else:
            domains[feature[0]] = 1
    for feature in domains:
        scale += 1
        p_score = 0.0
        for tool in tools:
            if feature in proteome[protein][tool]:
                e_feature = False
                try:
                    if proteome[protein][tool][feature]["evalue"] <= option["eFeature"]:
                        e_feature = True
                except TypeError:
                    e_feature = True
                if e_feature:
                    s_length = 0
                    for instance in proteome[protein][tool][feature]["instance"]:
                        e_instance = False
                        try:
                            if instance[2] <= option["eInstance"]:
                                s_length += 1
                        except TypeError:
                            s_length += 1
                    p_score = min(float(domains[feature] * s_length) / float(domains[feature] * domains[feature]), 1.0)
        scores.append((feature, p_score))
    if scale > 0:
        scale = 1.0 / float(scale)
    for score in scores:
        final_score += score[1] * scale
    return float(final_score), domains, scale


def sf_entire_ms_score(path, query_path, search_features, a_s_f, query_features, a_q_f, weights, option):
    """Calculates multiplicity score

    :param path: Path to score
    :param query_path: Path to score (query)
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param search_features: all features to be linearized in the seed protein
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param option: dictionary that contains the main option variables of FAS
    :return: final_score (MS), search_domains, scale, final_weight, common_feature
    """
    search_domains = {}
    query_domains = {}
    scale = 0
    scores = []
    final_score = 0.0
    final_weight = 0.0
    main_features_s = {}
    main_features_q = []
    common_feature = False

    for i in path:
        if i in search_features:
            feature = search_features[i]
            main_features_s[feature[0]] = True
        else:
            feature = a_s_f[i]

        if feature[0] in search_domains:
            search_domains[feature[0]] += 1
        else:
            search_domains[feature[0]] = 1

    for i in query_path:
        if i in query_features:
            feature = query_features[i]
            main_features_q.append(feature[0])
        else:
            feature = a_q_f[i]
        if feature[0] in query_domains:
            query_domains[feature[0]] += 1
        else:
            query_domains[feature[0]] = 1

    # check for common features
    for y in main_features_q:
        try:
            if main_features_s[y]:
                common_feature = True
        except KeyError:
            pass
    for feature in search_domains:
        if option["MS_uni"] == 0:
            scale += weights[feature]
        else:
            scale += 1
        if feature in query_domains:
            s_length = query_domains[feature]
            p_score = min(float(search_domains[feature] * s_length) /
                          float(search_domains[feature] * search_domains[feature]), 1.0)
            scores.append((feature, p_score))
        else:
            scores.append((feature, 0.0))
    if scale > 0:
        scale = 1.0 / float(scale)
    for score in scores:
        if option["MS_uni"] == 0:
            final_score += score[1] * scale * weights[score[0]]
            final_weight += weights[score[0]]
        else:
            final_score += score[1] * scale
    return float(final_score), search_domains, scale, final_weight, common_feature


def sf_ps_score(path, scale, protein, features, seed_proteome, option):
    """Calculates positional score

    :param path: Path to score
    :param scale: contains scaling factor for each weight or number of feature types for uniform weighting
    :param protein: protein id
    :param features: features (seed or query) [dict]
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :return: final_score (PS)
    """
    count = {}
    final_score = 0.0
    scores = {}
    tools = option["input_linearized"] + option["input_normal"]
    best_match = 0.0
    for i in path:
        feature = features[i]
        for tool in tools:
            if feature[0] in seed_proteome[protein][tool]:
                e_feature = False
                try:
                    if seed_proteome[protein][tool][feature[0]]["evalue"] <= option["eFeature"]:
                        e_feature = True
                except TypeError:
                    e_eature = True
                if e_feature:
                    best_match = 0.0
                    if not feature[0] in scores:
                        scores[feature[0]] = 0.0
                        count[feature[0]] = 0
                    for instance in seed_proteome[protein][tool][feature[0]]["instance"]:
                        e_instance = False
                        try:
                            if instance[2] <= option["eInstance"]:
                                e_instance = True
                        except TypeError:
                            e_instance = True
                        if e_instance:
                            pos = (float(instance[1]) + float(instance[2])) / 2.0 / float(
                                seed_proteome[protein]["length"])
                            match = 1.0 - float(abs(feature[1]) - pos)
                            if best_match < match:
                                best_match = match
                scores[feature[0]] += best_match
                count[feature[0]] += 1

    for f_score in scores:
        final_score += scores[f_score] / count[f_score] * scale
    return float(final_score)


def sf_entire_ps_score(path, scale, query_path, search_features, a_s_f, query_features, a_q_f, weights, option):
    """Calculates positional score

    :param path: Path to score
    :param scale: contains scaling factor for each weight or number of feature types for uniform weighting
    :param query_path: Path to score (query)
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param search_features: all features to be linearized in the seed protein
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param option: dictionary that contains the main option variables of FAS
    :return: final_score (PS), final_ls_score
    """
    count = {}
    final_score = 0.0
    final_ls_score = 0.0
    scores = {}
    ls_scores = {}
    local_query_protein = {}
    ls = 0.0
    # get current features from query path
    for i in query_path:
        if i in query_features:
            feature = query_features[i]
        else:
            feature = a_q_f[i]
        if feature[0] in local_query_protein:
            length = feature[3] - feature[2] + 1.0
            local_query_protein[feature[0]].append((feature[1], length))
        else:
            local_query_protein[feature[0]] = [(feature[1], feature[3] - feature[2] + 1.0)]

    # compare features in path with features form query path
    for i in path:
        if i in search_features:
            feature = search_features[i]
        else:
            feature = a_s_f[i]
        if feature[0] in local_query_protein:
            best_match = 0.0
            if not feature[0] in scores:
                count[feature[0]] = 0
                scores[feature[0]] = 0.0
                ls_scores[feature[0]] = 0.0
            for instance in local_query_protein[feature[0]]:
                match = 1.0 - float(abs(feature[1] - instance[0]))
                if best_match < match:
                    best_match = match
                    if instance[1] <= (feature[3] - feature[2]):
                        ls = float(instance[1]) / float(feature[3] - feature[2] + 1.0)
                    else:
                        ls = float(feature[3] - feature[2] + 1.0) / float(instance[1])

            scores[feature[0]] += best_match
            ls_scores[feature[0]] += ls
            count[feature[0]] += 1

    for f_score in scores:
        if option["MS_uni"] == 0:
            final_score += scores[f_score] / count[f_score] * scale * weights[f_score]
            final_ls_score += ls_scores[f_score] / count[f_score] * scale * weights[f_score]
        else:
            final_score += scores[f_score] / count[f_score] * scale
            final_ls_score += ls_scores[f_score] / count[f_score] * scale

    return float(final_score), float(final_ls_score)
