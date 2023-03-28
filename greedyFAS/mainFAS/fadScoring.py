# done
def fad_entire_calc_score(path, query_path, weights, search_features, a_s_f, query_features, a_q_f, clan_dict, option):
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
    score_cs = round(fad_entire_cs_score(path, query_path, query_features, clan_dict, search_features), 4)
    tmp = fad_entire_ms_score(path, query_path, search_features, a_s_f, query_features, a_q_f, weights, option)
    score_ms = round(tmp[0], 4)
    scale = tmp[2]
    path_weight = tmp[3]
    common_feature = tmp[4]
    ps_tmp = fad_entire_ps_score(path, tmp[2], query_path, search_features, a_s_f, query_features, a_q_f, weights,
                                 option)
    score_ps = round(ps_tmp[0], 4)
    score_ls = round(ps_tmp[1], 4)
    final_score = ((score_ms * option["score_weights"][0]) + (score_cs * option["score_weights"][1])
                   + (score_ps * option["score_weights"][2]))
    return score_ms, score_ps, score_cs, float(final_score), path_weight, common_feature, score_ls, scale


# done
def fad_entire_cs_score(path, query_path, query_features, clan_dict, search_features):
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
    p_score = 0.0

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
            p_score += 1.0 - min(float(s_clans[clan] * q_clans[clan]) / float(s_clans[clan] * s_clans[clan]), 1.0)
    if counter == 0:
        score = 0.0
    else:
        score = 1.0 - p_score / counter
    return float(score)


# done
def fad_entire_ms_score(path, query_path, search_features, a_s_f, query_features, a_q_f, weights, option):
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
    final_score = 1.0
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
        if feature in weights['weights']:
            scale += weights['weights'][feature]
        else:
            scale += weights['default']
        if feature in query_domains:
            p_score = 1.0 - min(float(search_domains[feature] * query_domains[feature]) /
                                float(search_domains[feature] * search_domains[feature]), 1.0)
            scores.append((feature, p_score))
        else:
            scores.append((feature, 1.0))
    if scale > 1.0:
        scale = 1.0 / float(scale)
    else:
        scale = 1.0
    for score in scores:
        if score[0] in weights['weights']:
            weight = weights['weights'][score[0]]
        else:
            weight = weights['default']
        final_score -= score[1] * scale * weight
        final_weight += weight
    return float(final_score), search_domains, scale, final_weight, common_feature


#done
def fad_entire_ps_score(path, scale, query_path, search_features, a_s_f, query_features, a_q_f, weights, option):
    """Calculates positional score FAD variant

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
    final_score = 1.0
    final_ls_score = 1.0
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
            min_distance = 1.0
            if not feature[0] in scores:
                count[feature[0]] = 0
                scores[feature[0]] = 0.0
                ls_scores[feature[0]] = 0.0
            for instance in local_query_protein[feature[0]]:
                distance = float(abs(feature[1] - instance[0]))
                if min_distance > distance:
                    min_distance = distance
                    if instance[1] <= (feature[3] - feature[2]):
                        ls = 1 - float(instance[1]) / float(feature[3] - feature[2] + 1.0)
                    else:
                        ls = 1 - float(feature[3] - feature[2] + 1.0) / float(instance[1])

            scores[feature[0]] += min_distance
            ls_scores[feature[0]] += ls
            count[feature[0]] += 1
        else:
            scores[feature[0]] = 1.0
            ls_scores[feature[0]] = 1.0
            count[feature[0]] = 1
    for f_score in scores:
        if f_score in weights['weights']:
            weight = weights['weights'][f_score]
        else:
            weight = weights['default']
        final_score -= scores[f_score] / count[f_score] * scale * weight
        final_ls_score -= ls_scores[f_score] / count[f_score] * scale * weight
    return float(final_score), float(final_ls_score)
