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


from operator import itemgetter
import os
import multiprocessing
from functools import partial
from copy import deepcopy
from tqdm import tqdm
import itertools
import sys
from time import sleep
from greedyFAS.mainFAS.fasInput import read_json
from greedyFAS.mainFAS.fasInput import check_version
from greedyFAS.mainFAS.fasOutput import write_domain_out
from greedyFAS.mainFAS.fasOutput import write_tsv_out, write_json_out
from greedyFAS.mainFAS.fasOutput import phyloprofile_out
from greedyFAS.mainFAS.fasScoring import sf_calc_score
from greedyFAS.mainFAS.fasScoring import sf_entire_calc_score
from greedyFAS.mainFAS.fasWeighting import w_weight_const_rescale
from greedyFAS.mainFAS.fasWeighting import w_weight_correction
from greedyFAS.mainFAS.fasWeighting import w_weighting
from greedyFAS.mainFAS.fasWeighting import w_weighting_constraints
from greedyFAS.mainFAS.fasWeighting import w_count_add_domains
from greedyFAS.mainFAS.fasPathing import pb_region_mapper
from greedyFAS.mainFAS.fasPathing import pb_region_paths
from greedyFAS.annoFAS.annoModules import mergeNestedDic


# important vars #             ###  var looks ###
# naming
# seed_proteome = {}           #{("protein_id", {("domain_name", [("START", "STOP")])})}
# query_proteome = {}          #{("protein_id", {("domain_name", [("START", "STOP")])})}
# clan_dict = {}               #{("domain_name", "clan")}
# search_features = {}         #{("F_0", ("domain_name", "POSITION", "Start", "Stop" ))}
# a_s_f = {}                   # additional seed features [non linearized]
# query_features = {}          #{("F_0", ("domain_name", "POSITION", "Start", "Stop" ))}
# a_q_f = {}                   # additional query features [non linearized]
# query_protein = {}           #{("domain_name", ["POSITION_1", "POSITION_2"])}
# query_clans = {}             #{("clan", "INSTANCES")}
# weights = {}                 #{("domain_name", "weight")}
# domain_count = {}            #{("domain", "COUNT")}


# flow control functions


def fc_start(option):
    """Overhead function,
    this function manages the individual functions that read the input files and prepares the data for the main script.

    :param option: dictionary that contains the main option variables of FAS
    """
    clan_dict = {}
    domain_count = {}
    phmm = {}
    option["reverse"] = False
    interprokeys = {}
    v_warning = False
    version = None
    # MS_uni set to 0 when no weighting is conducted
    if option["MS_uni"] == 0:
        domain_count = {}
        for path in option["ref_proteome"]:
            proteome = read_json(path)
            domain_count.update(proteome["count"])
            version, v_warning = check_version(version, proteome, v_warning)
    proteome_list = []
    for path in option["p_path"]:
        proteome = read_json(path)
        proteome_list.append(proteome["feature"])
        clan_dict.update(proteome["clan"])
        if 'interproID' in proteome:
            interprokeys.update(proteome['interproID'])
        if 'length' in proteome:
            phmm.update(proteome['length'])
        version, v_warning = check_version(version, proteome, v_warning)
    seed_proteome = mergeNestedDic(proteome_list)
    proteome_list = []
    for path in option["s_path"]:
        proteome = read_json(path)
        proteome_list.append(proteome["feature"])
        clan_dict.update(proteome["clan"])
        if 'interproID' in proteome:
            interprokeys.update(proteome['interproID'])
        if 'length' in proteome:
            phmm.update(proteome['length'])
        version, v_warning = check_version(version, proteome, v_warning)
    query_proteome = mergeNestedDic(proteome_list)
    if v_warning:
        print('##########\nWARNING: There are version discrepancies in the annotation files. This means that there '
              'might be differences in tool/database versions of the annotation. You should check the annotation files '
              'with "fas.getAnnoVersion"!\n##########')
    if option["weight_correction"]:
        domain_count = w_weight_correction(option["weight_correction"], domain_count)
    for tool in option["input_linearized"]:
        if tool not in seed_proteome[list(seed_proteome)[0]]:
            raise Exception(tool + " is missing in the seed annotation")
        if tool not in query_proteome[list(query_proteome)[0]]:
            raise Exception(tool + " is missing in the query annotation")
    if option["seed_id"]:
        for protid in option["seed_id"]:
            if protid not in seed_proteome:
                raise Exception(protid + " is not in the seed annotation")
    if option["query_id"]:
        for protid in option["query_id"]:
            if protid not in query_proteome:
                raise Exception(protid + " is not in the query annotation")
    if option["bidirectional"]:
        if not option["silent"]:
            print("calculating forward scores...")
        f_results = fc_main(domain_count, seed_proteome, query_proteome, clan_dict, option, interprokeys, phmm)
        r_results = ()
        if not f_results[-1] == 'NA':
            if option["MS_uni"] == 0 and option["ref_2"]:
                domain_count_2 = {}
                for path in option["ref_2"]:
                    domain_count_2.update(read_json(path)["count"])
                if option["weight_correction"]:
                    domain_count_2 = w_weight_correction(option["weight_correction"], domain_count_2)
                option['ref_proteome'] = option['ref_2']
            else:
                domain_count_2 = domain_count
            id_tmp = option["seed_id"]
            option["reverse"] = True
            option["seed_id"] = option["query_id"]
            option["query_id"] = id_tmp
            if option["pairwise"]:
                pairtmp = []
                for pair in option["pairwise"]:
                    pairtmp.append([pair[1], pair[0]])
                option["pairwise"] = pairtmp
            if not option["silent"]:
                print("calculating backward scores...")
            r_results = fc_main(domain_count_2, query_proteome, seed_proteome, clan_dict, option, interprokeys, phmm)

        if option["phyloprofile"]:
            phyloprofile_out(option["outpath"], True, option["phyloprofile"], (f_results, r_results))
        if not option['tsv']:
            write_tsv_out(option["outpath"], True, (f_results, r_results))
        if option['json']:
            write_json_out(option["outpath"], True, (f_results, r_results))
    else:
        if not option["silent"]:
            print("calculating forward scores...")
        results = fc_main(domain_count, seed_proteome, query_proteome, clan_dict, option, interprokeys, phmm)
        if not option["tsv"]:
            write_tsv_out(option["outpath"], False, (results, None))
        if option["json"]:
            write_json_out(option["outpath"], False, (results, None))
        if option["phyloprofile"]:
            phyloprofile_out(option["outpath"], False, option["phyloprofile"], [results, None])


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


def fc_main(domain_count, seed_proteome, query_proteome, clan_dict, option, interprokeys, phmm):
    """Main function,
    manages linearization and scoring, creates the output
    Function calls: w_weighting_constraints(), w_weighting(), su_lin_query_protein(), pb_graphtraversal(),
                    pb_entire_graphtraversal_priority(),
                    sf_entire_calc_score(), pb_entire_main_nongreedy(), w_weight_const_rescale()

    :param domain_count: dictionary that contains the (adjusted) counts of each feature in the reference proteome
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param query_proteome: dictionary that contains the feature architecture of all query proteins
    :param clan_dict: dictionary that maps features to clans
    :param option: dictionary that contains the main option variables of FAS
    """

    progress = None
    domain_out = None
    results = []
    if option["domain"]:
        if option["reverse"]:
            domain_out = open(option["outpath"] + "_reverse.domains", "w")
        else:
            domain_out = open(option["outpath"] + "_forward.domains", "w")
        domain_out.write('# pairID\torthoID\tseqLen\tfeature\tfStart\tfEnd\tfWeight\tfPath\tinterProID\te-value'
                         + '\tbitScore\tpStart\tpEnd\tpLen\n')

    final_pairs = []
    if option['pairwise']:
        final_pairs = option['pairwise']
    else:
        if option["query_id"]:
            querylist = option["query_id"]
            for query in querylist:
                if query not in query_proteome:
                    raise Exception(query + ' is missing in annotation!')
        else:
            querylist = list(query_proteome.keys())
        if option["seed_id"]:
            seedlist = option["seed_id"]
            for seed in seedlist:
                if seed not in seed_proteome:
                    raise Exception(seed + ' is missing in annotation!')
        else:
            seedlist = list(seed_proteome.keys())
        # create pairwise list
        final_pairs = list(itertools.product(seedlist,querylist))
        if option['old_json']:
            final_pairs = filter_oldJson(option['old_json'], final_pairs)
    if len(final_pairs) > 0:
        if option['progress']:
             progress = tqdm(total=len(final_pairs), file=sys.stdout)

        for pair in final_pairs:
            query = pair[1]
            protein = pair[0]
            ### CHECK QUERY & PROTEIN IN JSON ################
            if query not in query_proteome:
                raise Exception(query + ' is missing in annotation!')
            if protein not in seed_proteome:
                raise Exception(protein + ' is missing in annotation!')
            tmp_query = fc_prep_query(query, domain_count, query_proteome, option, clan_dict)
            if not tmp_query == None:
                tmp_protein = fc_prep_query(protein, 'NA', seed_proteome, option, clan_dict)
                if not tmp_protein == None:
                    query_graph, all_query_paths, lin_query_set, query_features, a_q_f, query_clans, clan_dict = tmp_query[0:7]
                    go_priority, domain_count = tmp_query[7:9]

                    #### WRITE RESULTS ####################
                    results.append(fc_main_sub(protein, domain_count, seed_proteome, option, all_query_paths, query_features,
                                               go_priority, a_q_f, clan_dict, query_graph, query_proteome, query, query_clans,
                                               domain_out, interprokeys, phmm))
                    if option["progress"]:
                        progress.update(1)
                else:
                    results.append((protein, query, ('NA', 'NA', 'NA', 'NA', 'NA'), 'NA'))
            else:
                results.append((protein, query, ('NA', 'NA', 'NA', 'NA', 'NA'), 'NA'))

        if option["progress"]:
            progress.refresh()
            progress.close()
            sleep(1.0)
    else:
        print('==> No new protein pairs need to be calculated!')
    return results


def fc_prep_query(query, domain_count, query_proteome, option, clan_dict):
    go_priority = False
    if option["MS_uni"] == 0:
        domain_count = w_count_add_domains(query, domain_count, query_proteome)
    lin_query_set, query_features, a_q_f, query_clans, clan_dict = su_lin_query_protein(
        query, query_proteome, clan_dict, option)
    tmp_query_graph, path_number = pb_region_paths(pb_region_mapper(
        lin_query_set, query_features, option["max_overlap"], option["max_overlap_percentage"]))
    if option["paths_limit"] > 1 and option["paths_limit"] < path_number:
        return None
    # PRIORITY CHECK: checking for number of instances - assess complexity of the feature graph
    if (len(query_features) > option["priority_threshold"] or path_number > option["max_cardinality"]) and \
            option["priority_mode"]:
        go_priority = True
        all_query_paths = "PRIORITY"
    elif option["priority_mode"]:
        # creating all paths
        all_query_paths = pb_graphtraversal(tmp_query_graph, [], [], option)
    elif len(tmp_query_graph) == 0:
        all_query_paths = []
    else:
        all_query_paths = "NOPRIORITY"
    return tmp_query_graph, all_query_paths, lin_query_set, query_features, a_q_f, query_clans, clan_dict, \
        go_priority, domain_count


def fc_main_sub(protein, domain_count, seed_proteome, option, all_query_paths, query_features, go_priority, a_q_f,
                clan_dict, tmp_query_graph, query_proteome, query, query_clans, domain_out, interprokeys, phmm):
    go_priority_2 = False
    mode = 0
    pathcount = 0
    weights = None
    unsolved = True
    if option["MS_uni"] == 0:
        if option["weight_const"]:
            weights, domain_count = w_weighting_constraints(protein, domain_count, seed_proteome, option)
        else:
            weights, domain_count = w_weighting(protein, domain_count, seed_proteome, option)
    search_protein, search_features, a_s_f = su_search_protein(protein, seed_proteome, option)
    search_protein = tuple(search_protein)
    max_fixture = ([], (0.0, 0.0, 0.0, 0.0, 0, 0, False), [], protein)
    region = pb_region_mapper(list(search_protein), search_features, option["max_overlap"],
                              option["max_overlap_percentage"])
    search_graph, search_path_number = pb_region_paths(region)
    # check for available paths
    # query <--VS-- seed
    # case M2.1: empty(query)

    if int(len(all_query_paths)) == 0 or int(len(query_features) == 0):
        # case M2.1.1: empty(query)-empty(search)
        # should be the best fix independent from weight
        if int(len(search_features)) == 0:
            path = list(a_s_f.keys())
            query_architecture = list(a_q_f.keys())
            score_w = sf_entire_calc_score(path, query_architecture, weights, search_features, a_s_f,
                                           query_features, a_q_f, clan_dict, option)
            mode = 0
        else:
            # case M2.1.2: empty(query)-graph(search)
            tmp_path_score = pb_entire_main_nongreedy(protein, [], search_features, weights, query_features,
                                                      seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                                                      option, search_graph, search_path_number)
            path = tmp_path_score[0][0]
            score_w = tmp_path_score[0][1]
            query_architecture = tmp_path_score[0][2]
            mode = 0

        # set max_fixture according to
        max_fixture = (path, score_w, query_architecture, protein)
        unsolved = False

    # seed architecture == query architecture
    elif seed_proteome[protein] == query_proteome[query]:
        path = list(a_s_f.keys()) + list(search_features.keys())
        query_architecture = list(a_q_f.keys()) + list(query_features.keys())
        score_w = sf_entire_calc_score(path, query_architecture, weights, search_features, a_s_f,
                                       query_features, a_q_f, clan_dict, option)
        mode = 0
        max_fixture = (path, score_w, query_architecture, protein)
        unsolved = False

    # handle all paths
    # case M2.2: graph(query)
    elif not option["priority_mode"] and not int(len(query_features)) == 0:
        pool_results = None
        mode = 0
        stack, jobpaths = pb_create_jobs(tmp_query_graph, option)
        if len(stack) > 0:
            jobpool = multiprocessing.Pool(processes=option["cores"])
            func = partial(pb_graph_traversal_sub, protein, search_features, weights, query_features, seed_proteome,
                           a_s_f, a_q_f, clan_dict, query_clans, tmp_query_graph, option, search_graph,
                           search_path_number)
            try:
                pool_results = jobpool.map_async(func, stack).get(timeout=option["timelimit"])
            except multiprocessing.TimeoutError:
                go_priority_2 = True
                jobpool.terminate()
            if not go_priority_2:
                for pool_result in pool_results:
                    if pool_result[1][3] >= max_fixture[1][3]:
                        max_fixture = deepcopy(pool_result)
        if len(jobpaths) > 0 and not go_priority_2:
            # case M2.2.1: graph(query)-empty(search)
            if int(len(search_features)) == 0:
                for jobpath in jobpaths:
                    pathcount += 1
                    # special case: protein with no pfam or smart domains
                    # get score for a_s_f and query_path directly
                    path = list(a_s_f.keys())
                    query_path_ad = jobpath + list(a_q_f.keys())
                    score_w = sf_entire_calc_score(path, query_path_ad, weights, search_features, a_s_f,
                                                   query_features, a_q_f, clan_dict, option)

                    # check for max scoring fixture of path and query_path
                    if (score_w[4] >= max_fixture[1][4] and score_w[3] == max_fixture[1][3]) or \
                            score_w[3] >= max_fixture[1][3]:
                        max_fixture = (path, score_w, query_path_ad, protein)
            else:
                # case M2.2.2 graph(query)-graph(search)
                # regular traversal of graph based on search_protein
                jobpool = multiprocessing.Pool(processes=option["cores"])
                jobs = []
                jobcount = (int(len(jobpaths) / option["cores"]), len(jobpaths) % option["cores"])
                if jobcount[0] == 0:
                    for i in range(jobcount[1]):
                        jobs.append([jobpaths[-(i + 1)]])
                else:
                    counter = 0
                    for i in range(option["cores"]):
                        jobs.append(jobpaths[counter: (i + 1) * jobcount[0]])
                        counter = (i + 1) * jobcount[0]
                    for i in range(jobcount[1]):
                        jobs[i].append(jobpaths[-(i + 1)])
                func = partial(pb_calc_sub, protein, search_features, weights, query_features,
                               seed_proteome, a_s_f, a_q_f, clan_dict, query_clans, option, search_graph,
                               search_path_number)
                try:
                    pool_results = jobpool.map_async(func, jobs).get(timeout=option["timelimit"])
                except multiprocessing.TimeoutError:
                    go_priority_2 = True
                    jobpool.terminate()
                if not go_priority_2:
                    for pool_result in pool_results:
                        if pool_result[1][3] >= max_fixture[1][3]:
                            max_fixture = deepcopy(pool_result)
    if go_priority or go_priority_2 and unsolved:
        mode = 1
        all_query_paths = []
        priority_list = []
        for tool in option["input_linearized"]:
            for domain_type in query_proteome[query][tool]:
                priority_list.append(domain_type)
        priority_list.append("NONE")
        for domain_type in priority_list:
            all_query_paths += (
                pb_entire_graphtraversal_priority(tmp_query_graph, domain_type, protein, 1, search_features,
                                                  weights, query_features, seed_proteome, a_s_f, a_q_f,
                                                  clan_dict, query_clans, option))
            # max fixture of (search_path, score, query_path)
    if option["priority_mode"] or go_priority_2 and unsolved:
        for query_path in all_query_paths:
            pathcount += 1

            # case M2.2.1: graph(query)-empty(search)
            if int(len(search_features)) == 0:
                # special case: protein with no pfam or smart domains
                # get score for a_s_f and query_path directly
                path = list(a_s_f.keys())
                query_path_ad = query_path + list(a_q_f.keys())
                score_w = sf_entire_calc_score(path, query_path_ad, weights, search_features, a_s_f,
                                               query_features, a_q_f, clan_dict, option)
            else:
                # case M2.2.2 graph(query)-graph(search)
                # regular traversal of graph based on search_protein
                tmp_path_score = pb_entire_main_nongreedy(protein, query_path, search_features, weights,
                                                              query_features, seed_proteome, a_s_f, a_q_f, clan_dict,
                                                              query_clans, option, search_graph, search_path_number)

                path = tmp_path_score[0][0]
                score_w = tmp_path_score[0][1]
                query_path_ad = tmp_path_score[0][2]
                mode = tmp_path_score[1]

            # check for max scoring fixture of path and query_path
            if (score_w[4] >= max_fixture[1][4] and score_w[3] == max_fixture[1][3]) or \
                    score_w[3] >= max_fixture[1][3]:
                max_fixture = (path, score_w, query_path_ad, protein)

    score = max_fixture[1]
    try:
        if mode == 1:
            mode_out = "priority"
        else:
            mode_out = "exhaustive"
    except KeyError:
        mode_out = "default"
    if go_priority or go_priority_2:
        mode_out += "/priority"
    else:
        mode_out += "/exhaustive"

    best_template_path = []
    best_query_path = []

    for feature in max_fixture[0]:
        if feature in search_features:
            best_template_path.append((search_features[feature][0], search_features[feature][1],
                                       search_features[feature][2], search_features[feature][3]))
        else:
            best_template_path.append(
                (a_s_f[feature][0], a_s_f[feature][1], a_s_f[feature][2], a_s_f[feature][3]))

    for feature in max_fixture[2]:
        if feature in query_features:
            best_query_path.append((query_features[feature][0], query_features[feature][1],
                                    query_features[feature][2], query_features[feature][3]))
        else:
            best_query_path.append((a_q_f[feature][0], a_q_f[feature][1], a_q_f[feature][2], a_q_f[feature][3]))
    path_tmp = {}
    path_tmp_query = {}
    scale = 0
    if option["weight_const"] == 1:
        path_tmp2 = []
        for feature in best_template_path:
            if feature[0] not in path_tmp2:
                path_tmp2.append(feature[0])
        adjusted_weights = w_weight_const_rescale(path_tmp2, weights, search_features, True, option)
        for adj_feature in adjusted_weights:
            weights[adj_feature] = adjusted_weights[adj_feature]
    for feature in best_template_path:
        if feature[0] in path_tmp:
            path_tmp[feature[0]].append((feature[2], feature[3]))
        else:
            path_tmp[feature[0]] = [(feature[2], feature[3])]
            if option["MS_uni"] == 0:
                scale += weights[feature[0]]
    for feature in best_query_path:
        if feature[0] in path_tmp_query:
            path_tmp_query[feature[0]].append((feature[2], feature[3]))
        else:
            path_tmp_query[feature[0]] = [(feature[2], feature[3])]

    # unweighted case
    if option["MS_uni"] == 0:
        if scale > 0:
            scale = 1.0 / float(scale)
        else:
            scale = 1.0
    if option["domain"]:
        write_domain_out(seed_proteome, query_proteome, protein, query, weights, scale, path_tmp, path_tmp_query,
                             domain_out, option, interprokeys, phmm)
    if option["raw"]:
        print('#' + '\t' + protein + '\t' + query + '\t' + str(score[3]))
    return protein, query, (score[3], score[0], score[1], score[2], score[6]), mode_out


# Start Up Functions <su>


def su_lin_query_protein(protein_id, query_proteome, clan_dict, option):
    """Initializes variables for the current query protein (linearization version)

    :param option: dictionary that contains the main option variables of FAS
    :param protein_id: String that contains the identifier of the query protein
    :param query_proteome: dictionary that contains the feature architecture of all query proteins
    :param clan_dict: dictionary that maps features to clans
    :return:lin_query_protein, query_features, a_q_f(additional[not linearized] query features), query_clans, clan_dict
    """

    lin_query_protein = []
    query_clans = {}
    query_features = {}
    a_q_f = {}
    tmp = []
    i = 0

    for tool in option["input_linearized"]:
        for feature in query_proteome[protein_id][tool]:
            e_feature = False
            try:
                if query_proteome[protein_id][tool][feature]["evalue"] <= option["eFeature"]:
                    e_feature = True
            except TypeError:
                e_feature = True
            if e_feature:
                clan = None
                if feature in clan_dict:
                    clan = clan_dict[feature]
                tmp_clan = 0
                for instance in query_proteome[protein_id][tool][feature]["instance"]:
                    e_instance = False
                    try:
                        if instance[2] <= option["eInstance"]:
                            e_instance = True
                    except TypeError:
                        e_instance = True
                    if e_instance:
                        position = ((float(instance[0]) + float(instance[1])) / 2.0) / float(
                            query_proteome[protein_id]["length"])
                        position = round(position, 4)
                        key = "F_" + str(i)
                        query_features[key] = (feature, position, instance[0], instance[1])
                        tmp.append((key, instance[1]))
                        i += 1
                        tmp_clan += 1
                if clan in query_clans:
                    query_clans[clan] += tmp_clan
                elif tmp_clan > 0:
                    query_clans[clan] = tmp_clan
    for tool in option["input_normal"]:
        for feature in query_proteome[protein_id][tool]:
            e_feature = False
            try:
                if query_proteome[protein_id][tool][feature]["evalue"] <= option["eFeature"]:
                    e_feature = True
            except TypeError:
                e_feature = True
            if e_feature:
                clan = None
                if feature in clan_dict:
                    clan = clan_dict[feature]
                tmp_clan = 0
                for instance in query_proteome[protein_id][tool][feature]["instance"]:
                    e_instance = False
                    try:
                        if instance[2] <= option["eInstance"]:
                            e_instance = True
                    except TypeError:
                        e_instance = True
                    if e_instance:
                        position = ((float(instance[0]) + float(instance[1])) / 2.0) / float(
                            query_proteome[protein_id]["length"])
                        position = round(position, 4)
                        key = "O_" + str(i)
                        a_q_f[key] = (feature, position, instance[0], instance[1])
                        i += 1
                        tmp_clan += 1
                if clan in query_clans:
                    query_clans[clan] += tmp_clan
                elif tmp_clan > 0:
                    query_clans[clan] = tmp_clan

    # sort  instances
    tmp2 = sorted(tmp, key=itemgetter(1))
    for x in tmp2:
        lin_query_protein.append(x[0])

    return lin_query_protein, query_features, a_q_f, query_clans, clan_dict


def su_search_protein(protein_id, seed_proteome, option):
    """Initializes variables for the current seed protein

    :param option: dictionary that contains the main option variables of FAS
    :param protein_id: String that contains the identifier of the seed protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :return: search_protein, search_features, a_s_f (additional [not linearized] seed features)
    """
    search_features = {}
    search_protein = []
    a_s_f = {}
    tmp = []
    i = 0

    for tool in option["input_linearized"]:
        for feature in seed_proteome[protein_id][tool]:
            e_feature = False
            try:
                if seed_proteome[protein_id][tool][feature]["evalue"] <= option["eFeature"]:
                    e_feature = True
            except TypeError:
                e_feature = True
            if e_feature:
                for instance in seed_proteome[protein_id][tool][feature]["instance"]:
                    e_instance = False
                    try:
                        if instance[2] <= option["eInstance"]:
                            e_instance = True
                    except TypeError:
                        e_instance = True
                    if e_instance:
                        position = ((float(instance[0]) + float(instance[1])) / 2.0) / float(
                            seed_proteome[protein_id]["length"])
                        position = round(position, 4)
                        key = "F_" + str(i)
                        search_features[key] = (feature, position, instance[0], instance[1])
                        tmp.append((key, instance[1]))
                        i += 1
    for tool in option["input_normal"]:
        for feature in seed_proteome[protein_id][tool]:
            e_feature = False
            try:
                if seed_proteome[protein_id][tool][feature]["evalue"] <= option["eFeature"]:
                    e_feature = True
            except TypeError:
                e_feature = True
            if e_feature:
                for instance in seed_proteome[protein_id][tool][feature]["instance"]:
                    e_instance = False
                    try:
                        if instance[2] <= option["eInstance"]:
                            e_instance = True
                    except TypeError:
                        e_instance = True
                    if e_instance:
                        position = ((float(instance[0]) + float(instance[1])) / 2.0) / float(
                            seed_proteome[protein_id]["length"])
                        position = round(position, 4)
                        key = "O_" + str(i)
                        a_s_f[key] = (feature, position, instance[0], instance[1])
                        i += 1
    # sort  instances
    tmp2 = sorted(tmp, key=itemgetter(1))
    for x in tmp2:
        search_protein.append(x[0])
    return tuple(search_protein), search_features, a_s_f


# Path-building Functions # <pb>
# Used for the Pfam/Smart domains (default)
def pb_create_jobs(graph, option):
    """Does a breadth-first search to divide the feature architecture graph into several sub-problems so that it can be
    worked on by multiple cores. There are more jobs created than there are cores to counter the possibility of a core
    having nothing to do because of different running times for each job.

    :param graph: graph containing all paths through a protein architecture
    :param option: dictionary that contains the main option variables of FAS
    :return: queue (contains the starting points and sub-paths of the jobs), paths (empty unless b-f search reaches end
             node)
    """
    limit = 4 * option["cores"]
    queue = [("START", [])]
    paths = []
    while len(queue) < limit and queue:
        vertex, path = queue.pop(0)
        for next_vertex in graph[vertex]:
            if next_vertex == "END":
                paths.append(path)
            else:
                queue.append((next_vertex, path + [next_vertex]))
    return queue, paths


def pb_graph_traversal_sub(protein, search_features, weights, query_features, seed_proteome, a_s_f, a_q_f, clan_dict,
                           query_clans, query_graph, option, search_graph, search_path_number, stack):
    """Sub function for multi-processing: This function does a depth-first search in the feature graph. The starting
     point in the graph is given in the stack.

    :param search_path_number: number of paths in search protein
    :param search_graph: architecture graph of the search protein
    :param protein: String that contains the identifier of the seed protein
    :param search_features: all features to be linearized in the seed protein
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: (Pfam)-clans of the current query protein
    :param option: dictionary that contains the main option variables of FAS
    :param stack: contains the starting point in the graph and sub-path
    :param query_graph: graph of the query architecture
    :return: max_fixture
    """
    max_fixture = ([], (0.0, 0.0, 0.0, 0.0, 0, 0, False), [], protein)
    v_stack, p_stack = [stack[0]], [stack[1]]
    while v_stack:
        query_path, v_stack, p_stack = pb_graphtraversal(query_graph, v_stack, p_stack, option)
        # case M2.2.1: graph(query)-empty(search)
        if int(len(search_features)) == 0:
            # special case: protein with no pfam or smart domains
            # get score for a_s_f and query_path directly
            path = list(a_s_f.keys())
            query_path_ad = query_path + list(a_q_f.keys())
            score_w = sf_entire_calc_score(path, query_path_ad, weights, search_features, a_s_f,
                                           query_features, a_q_f, clan_dict, option)
        else:
            # case M2.2.2 graph(query)-graph(search)
            # regular traversal of graph based on search_protein
            tmp_path_score = pb_entire_main_nongreedy(protein, query_path, search_features, weights, query_features,
                                                      seed_proteome, a_s_f, a_q_f, clan_dict, query_clans, option,
                                                      search_graph, search_path_number)
            path = tmp_path_score[0][0]
            score_w = tmp_path_score[0][1]
            query_path_ad = tmp_path_score[0][2]

        # check for max scoring fixture of path and query_path
        if ((not max_fixture[1][5]) and (score_w[4] >= max_fixture[1][4]) and score_w[3] <= max_fixture[1][3]) or \
                score_w[3] >= max_fixture[1][3]:
            max_fixture = (path, score_w, query_path_ad, protein)
    return max_fixture


def pb_calc_sub(protein, search_features, weights, query_features, seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                option, search_graph, search_path_number, jobpaths):
    """Sub function for multiprocessing [2]: Goes through a list of (query) paths an evaluates them.

    :param protein: String that contains the identifier of the seed protein
    :param search_features: all features to be linearized in the seed protein
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: (Pfam)-clans of the current query protein
    :param option: dictionary that contains the main option variables of FAS
    :param jobpaths: a group of (query) paths to be evaluated
    :param search_path_number: number of paths in search protein
    :param search_graph: architecture graph of the search protein
    :return: max_fixture, tmpcheck - tmpstart (time taken for calculation)

    """
    max_fixture = ([], (0.0, 0.0, 0.0, 0.0, 0, 0, False), [], protein)
    for jobpath in jobpaths:
        tmp_path_score = pb_entire_main_nongreedy(protein, jobpath, search_features, weights, query_features,
                                                  seed_proteome, a_s_f, a_q_f, clan_dict, query_clans, option,
                                                  search_graph, search_path_number)
        path = tmp_path_score[0][0]
        score_w = tmp_path_score[0][1]
        query_path_ad = tmp_path_score[0][2]
        if ((not max_fixture[1][5]) and (score_w[4] >= max_fixture[1][4]) and score_w[3] <= max_fixture[1][3]) or \
                score_w[3] >= max_fixture[1][3]:
            max_fixture = (path, score_w, query_path_ad, protein)
    return max_fixture


def pb_entire_main_nongreedy(protein_id, query_path, search_features, weights, query_features, seed_proteome, a_s_f,
                             a_q_f, clan_dict, query_clans, option, search_graph, search_path_number):
    """Main Path-building function,
    creates graph with all paths(for seed/search protein), checks priority mode activation if necessary retrieves best
    path from graph traversal function, returns best path, score and mode
    Function calls: pb_region_mapper(), pb_region_paths(), pb_entire_priority_mode(),
                    pb_entire_graphtraversal()


    :param protein_id: String that contains the identifier of the seed protein
    :param query_path: currently evaluated query path
    :param search_features: all features to be linearized in the seed protein
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: (Pfam)-clans of the current query protein
    :param option: dictionary that contains the main option variables of FAS
    :param search_path_number: number of paths in search protein
    :param search_graph: architecture graph of the search protein
    :param search_path_number: number of paths in search protein
    :param search_graph: architecture graph of the search protein
    :return: path_score, mode[priority or extensive]

    """
    if (int(len(search_features)) >= option["priority_threshold"] or search_path_number > option['max_cardinality']) \
            and option["priority_mode"]:
        mode = 1
        path_score = pb_entire_priority_mode(protein_id, query_path, search_graph, search_features, weights,
                                             query_features, seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                                             option)
    else:
        mode = 0
        path_score = pb_entire_graphtraversal(search_graph, query_path, search_features, weights, query_features,
                                              a_s_f, a_q_f, clan_dict, option)
    return path_score, mode


def pb_entire_priority_mode(protein, query_path, search_graph, search_features, weights, query_features, seed_proteome,
                            a_s_f, a_q_f, clan_dict, query_clans, option):
    """Path-evaluation in priority mode

    :param protein: String that contains the identifier of the seed protein
    :param query_path: currently evaluated query path
    :param search_graph: graph containing all paths through the seed architecture
    :param search_features: all features to be linearized in the seed protein
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: (Pfam)-clans of the current query protein
    :param option: dictionary that contains the main option variables of FAS
    :return: best_path
    """

    # best_path (Path, (MS_score(float), PS_score(float), CS_score(float), final_score(float), path_weight(float),
    # common_feature(bool)), QPath)
    best_path = ("NULL", (0.0, 0.0, 0.0, 0.0, 0.0, False), "NULL")
    priority_list = []
    for tool in option["input_linearized"]:
        for domain_type in seed_proteome[protein][tool]:
            priority_list.append(domain_type)
    priority_list.append("NONE")
    for domain_type in priority_list:
        path_score_w = pb_entire_graphtraversal_priority(search_graph, domain_type, query_path, 0, search_features,
                                                         weights, query_features, seed_proteome, a_s_f, a_q_f,
                                                         clan_dict, query_clans, option)
        if path_score_w[1][5]:
            if path_score_w[1][3] >= best_path[1][3]:
                best_path = (path_score_w[0], path_score_w[1], path_score_w[2])
        else:
            # if so far no path with common features found AND the weight is higher
            if (not best_path[1][5]) and (path_score_w[1][4] >= best_path[1][4]):
                best_path = (path_score_w[0], path_score_w[1], path_score_w[2])
    return best_path


def pb_entire_graphtraversal(search_graph, query_path, search_features, weights, query_features, a_s_f, a_q_f,
                             clan_dict, option):
    """Traverses the feature architecture graph (seed) to score all paths

    :param search_graph: graph containing all paths through the seed architecture
    :param query_path: currently evaluated query path
    :param search_features: all features to be linearized in the seed protein
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param option: dictionary that contains the main option variables of FAS
    :return: best_path
    """

    v_stack = ["START"]
    p_stack = []
    best_path = ([], (0.0, 0.0, 0.0, 0.0, 0.0, False, 0.0), [])
    first = 1
    if len(query_path) > 0:
        while v_stack:
            vertex = v_stack.pop()
            if len(p_stack) == 0:
                path = []
            else:
                path = p_stack.pop()
            for next_vertex in search_graph[vertex]:
                if next_vertex == "END":
                    path_ad = path + list(a_s_f.keys())
                    query_path_ad = query_path + list(a_q_f.keys())
                    score_w = sf_entire_calc_score(path_ad, query_path_ad, weights, search_features, a_s_f,
                                                   query_features, a_q_f, clan_dict, option)
                    if (score_w[4] >= best_path[1][4] and score_w[3] == best_path[1][3]) or \
                       score_w[3] >= best_path[1][3]:
                        best_path = (path_ad, score_w, query_path_ad)
                else:
                    v_stack.append(next_vertex)
                    p_stack.append(path + [next_vertex])
    else:
        while v_stack and first:
            vertex = v_stack.pop()
            if len(p_stack) == 0:
                path = []
            else:
                path = p_stack.pop()
            for next_vertex in search_graph[vertex]:
                if next_vertex == "END":
                    path_ad = path + list(a_s_f.keys())
                    query_path_ad = query_path + list(a_q_f.keys())
                    score_w = sf_entire_calc_score(path_ad, query_path_ad, weights, search_features, a_s_f,
                                                   query_features, a_q_f, clan_dict, option)
                    if (score_w[4] >= best_path[1][4] and score_w[3] == best_path[1][3]) or \
                       score_w[3] >= best_path[1][3]:
                        best_path = (path_ad, score_w, query_path_ad)
                    first = 0
                else:
                    v_stack.append(next_vertex)
                    p_stack.append(path + [next_vertex])
    return best_path


def pb_graphtraversal(graph, v_stack, p_stack, option):
    """Traverses the feature architecture graph (query) to evaluate number of paths, returns all paths if # of paths is
     smaller than the max_cardinality

    :param graph: graph containing all paths through the (query) architecture
    :param option: dictionary that contains the main option variables of FAS
    :param v_stack: vertex stack
    :param p_stack: path stack
    :return: paths
    """
    paths = []
    if v_stack == [1]:
        v_stack = ["START"]
        p_stack = []
        mode = 1
    elif not v_stack:
        v_stack = ["START"]
        p_stack = []
        mode = 0
    else:
        mode = 1
    while v_stack:
        vertex = v_stack.pop()
        if len(p_stack) == 0:
            path = []
        else:
            path = p_stack.pop()
        for next_vertex in graph[vertex]:
            if next_vertex == "END":
                if mode == 0:
                    paths.append(path)
                    if len(paths) > option["max_cardinality"]:
                        return paths
                else:
                    return path, v_stack, p_stack
            else:
                v_stack.append(next_vertex)
                p_stack.append(path + [next_vertex])
    return paths


def pb_entire_graphtraversal_priority(search_graph, priority, query_path, mode, search_features, weights,
                                      query_features, seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                                      option):
    """Traverses the feature architecture graph in priority mode

    :param search_graph: graph containing all paths through the seed architecture
    :param priority: current priority feature
    :param query_path: graph containing all paths through the query architecture
    :param mode: query or seed linearization
    :param search_features: all features to be linearized in the seed protein
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: (Pfam)-clans of the current query protein
    :param option: dictionary that contains the main option variables of FAS
    :return: paths or best_path
    """

    v_stack = ["START"]
    p_stack = []
    protein = None
    paths = []
    best_path = None
    if mode == 0:
        best_path = ([], (0.0, 0.0, 0.0, 0.0, 0.0, False), [])
    else:
        protein = query_path
    while v_stack:
        vertex = v_stack.pop()
        if len(p_stack) == 0:
            path = []
        else:
            path = p_stack.pop()

        p_found = 0
        p_candidates = []
        for next_vertex in search_graph[vertex]:
            if next_vertex == "END":
                if mode == 1:
                    paths.append(path)
                    p_found = 2
                else:
                    path_ad = path + list(a_s_f.keys())
                    query_path_ad = query_path + list(a_q_f.keys())
                    score_w = sf_entire_calc_score(path_ad, query_path_ad, weights, search_features, a_s_f,
                                                   query_features, a_q_f, clan_dict, option)
                    if ((score_w[4] >= best_path[1][4] and score_w[3] == best_path[1][3])
                            or score_w[3] >= best_path[1][3]):
                        best_path = (path_ad, score_w, query_path_ad)
                    p_found = 2
            elif mode == 1:
                if query_features[next_vertex][0] == priority:
                    p_candidates.append(next_vertex)
                    p_found = 1
            elif mode == 0:
                p_candidates.append(next_vertex)
                p_found = 1
        if p_found == 1:
            if len(p_candidates) == 1:
                v_stack.append(p_candidates[0])
                p_stack.append(path + p_candidates)
            else:
                if mode == 1:
                    best_priority_bridger = ("NONE", (0.0, 0.0, 0.0, 0.0, 0.0, False))
                    for next_vertex in p_candidates:
                        score = sf_calc_score(path + [next_vertex], protein, weights, search_features, query_features,
                                              seed_proteome, clan_dict, query_clans, option)
                        if score[3] >= best_priority_bridger[1][3]:
                            best_priority_bridger = (next_vertex, score)
                    v_stack.append(best_priority_bridger[0])
                    p_stack.append(path + [best_priority_bridger[0]])
                elif mode == 0:
                    # greedy strategy: if feature type priority appears more than once
                    best_partial_path = ("NONE", (0.0, 0.0, 0.0, 0.0, 0.0, False))
                    for next_vertex in p_candidates:
                        path_ad = path + list(a_s_f.keys())
                        query_path_ad = query_path + list(a_q_f.keys())
                        score_w = sf_entire_calc_score(path_ad + [next_vertex], query_path_ad, weights,
                                                       search_features, a_s_f, query_features, a_q_f, clan_dict, option)
                        if (score_w[4] >= best_partial_path[1][4] and score_w[3] == best_partial_path[1][3]) \
                           or score_w[3] >= best_partial_path[1][3]:
                            best_partial_path = (next_vertex, score_w)

                    v_stack.append(best_partial_path[0])
                    p_stack.append(path + [best_partial_path[0]])

        elif p_found == 0:
            if mode == 1:
                best_priority_bridger = ("NONE", (0.0, 0.0, 0.0, 0.0, 0.0, False))
                for next_vertex in search_graph[vertex]:
                    score = sf_calc_score(path + [next_vertex], protein, weights, search_features, query_features,
                                          seed_proteome, clan_dict, query_clans, option)
                    if score[3] >= best_priority_bridger[1][3]:
                        best_priority_bridger = (next_vertex, score)
                v_stack.append(best_priority_bridger[0])
                p_stack.append(path + [best_priority_bridger[0]])

            elif mode == 0:
                # greedy strategy: if feature type (p priority) not found
                best_partial_path = ("NONE", (0.0, 0.0, 0.0, 0.0, 0.0, False))
                for next_vertex in search_graph[vertex]:
                    path_ad = path + list(a_s_f.keys())
                    query_path_ad = query_path + list(a_s_f.keys())
                    score_w = sf_entire_calc_score(path_ad + [next_vertex], query_path_ad, weights, search_features,
                                                   a_s_f, query_features, a_q_f, clan_dict, option)
                    if (score_w[4] >= best_partial_path[1][4] and score_w[3] == best_partial_path[1][3]) or \
                       score_w[3] >= best_partial_path[1][3]:
                        best_partial_path = (next_vertex, score_w)

                v_stack.append(best_partial_path[0])
                p_stack.append(path + [best_partial_path[0]])
    if mode == 1:
        return paths
    elif mode == 0:
        return best_path
