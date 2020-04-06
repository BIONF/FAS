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
import logging
import inspect
import os
import multiprocessing
import argparse
import time
from functools import partial
from copy import deepcopy
from greedyFAS.fasInput import xmlreader
from greedyFAS.fasInput import read_pairwise
from greedyFAS.fasInput import constraints_in
from greedyFAS.fasInput import featuretypes
from greedyFAS.fasOutput import bidirectionout
from greedyFAS.fasOutput import domain_out
from greedyFAS.fasOutput import phyloprofile_out
from greedyFAS.fasScoring import sf_calc_score
from greedyFAS.fasScoring import sf_entire_calc_score
from greedyFAS.fasWeighting import w_count
from greedyFAS.fasWeighting import w_count_ref
from greedyFAS.fasWeighting import w_weight_const_rescale
from greedyFAS.fasWeighting import w_weight_correction
from greedyFAS.fasWeighting import w_weighting
from greedyFAS.fasWeighting import w_weighting_constraints
from greedyFAS.fasWeighting import w_count_add_domains
from greedyFAS.fasPathing import pb_region_mapper
from greedyFAS.fasPathing import pb_region_paths


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
# protein_lengths = {}         #{("protein_id", "LENGTH")}
# domain_count = {}            #{("domain", "COUNT")}

# hidden options
# tab separated table as output file #
taciturn = 1

# Flow Control <fc>
# flow control functions


def fc_start(option):
    """Overhead function,
    this function manages the individual functions that read the input files and prepares the data for the main script.
    Function calls: xmlreader(), w_count_ref(), w_count(), fc_main()

    :param option: dictionary that contains the main option variables of FAS
    """
    clan_dict = {}
    seed_proteome = {}
    query_proteome = {}
    protein_lengths = {}
    domain_count = {}
    prot_count = 0
    ref_proteome = {}

    # MS_uni set to 0 when no weighting is conducted
    if option["MS_uni"] == 0:
        for ftype in option["input_linearized"]:
            ref_proteome, tmp, clan_dict = xmlreader(option["ref_proteome"] + "/" + ftype + ".xml", 2, ftype, True,
                                                     ref_proteome, protein_lengths, clan_dict, option)
        for ftype in option["input_normal"]:
            ref_proteome, tmp, clan_dict = xmlreader(option["ref_proteome"] + "/" + ftype + ".xml", 2, ftype, True,
                                                     ref_proteome, protein_lengths, clan_dict, option)

        prot_count, domain_count = w_count_ref(ref_proteome)
    for ftype in option["input_linearized"]:
        seed_proteome, protein_lengths, clan_dict = xmlreader(option["p_path"] + "/" + ftype + ".xml", 0, ftype, True,
                                                              seed_proteome, protein_lengths, clan_dict, option)
        query_proteome, protein_lengths, clan_dict = xmlreader(option["s_path"] + "/" + ftype + ".xml", 1, ftype, True,
                                                               query_proteome, protein_lengths, clan_dict, option)
    for ftype in option["input_normal"]:
        seed_proteome, protein_lengths, clan_dict = xmlreader(option["p_path"] + "/" + ftype + ".xml", 0, ftype, False,
                                                              seed_proteome, protein_lengths, clan_dict, option)
        query_proteome, protein_lengths, clan_dict = xmlreader(option["s_path"] + "/" + ftype + ".xml", 1, ftype, False,
                                                               query_proteome, protein_lengths, clan_dict, option)
    if option["MS_uni"] == 0:
        relevant_features, domain_count = w_count(prot_count, domain_count, seed_proteome, query_proteome)
    else:
        relevant_features = {}
    if option["weight_correction"]:
        domain_count = w_weight_correction(option["weight_correction"], domain_count)
    if option["bidirectional"]:
        fc_main(relevant_features, prot_count, domain_count, seed_proteome, query_proteome, protein_lengths, clan_dict,
                option)
        tmp = {}
        for protein in protein_lengths:
            if protein[0:4] == "seed":
                tmp["query_" + protein[5:]] = protein_lengths[protein]
            else:
                tmp["seed_" + protein[6:]] = protein_lengths[protein]
        if option["e_output"]:
            extmp = True
        else:
            extmp = False
        option["e_output"] = False
        org_outpath = option["outpath"]
        option["outpath"] += "_reverse"
        if option["MS_uni"] == 0 and option["ref_2"]:
            ref_proteome = {}
            for ftype in option["input_linearized"]:
                ref_proteome, tmp2, clan_dict = xmlreader(option["ref_2"] + "/" + ftype + ".xml", 2, ftype, True,
                                                          ref_proteome, protein_lengths, clan_dict, option)
            for ftype in option["input_normal"]:
                ref_proteome, tmp2, clan_dict = xmlreader(option["ref_2"] + "/" + ftype + ".xml", 2, ftype, True,
                                                          ref_proteome, protein_lengths, clan_dict, option)

            prot_count, domain_count_2 = w_count_ref(ref_proteome)
            if option["weight_correction"]:
                domain_count_2 = w_weight_correction(option["weight_correction"], domain_count_2)
        else:
            option["feature_info"] = False
            domain_count_2 = domain_count
        option['ref_proteome'] = option['ref_2']
        fc_main(relevant_features, prot_count, domain_count_2, query_proteome, seed_proteome, tmp, clan_dict, option)
        if option["domain"]:
            domain_out(org_outpath, True, extmp, option["MS_uni"])
        if option["phyloprofile"]:
            phyloprofile_out(org_outpath, True, option["phyloprofile"], extmp, option["MS_uni"])
        else:
            bidirectionout(org_outpath)
    else:
        fc_main(relevant_features, prot_count, domain_count, seed_proteome, query_proteome, protein_lengths, clan_dict,
                option)
        if option["domain"]:
            domain_out(option["outpath"], False, option["e_output"], option["MS_uni"])
        if option["phyloprofile"]:
            phyloprofile_out(option["outpath"], False, option["phyloprofile"], option["e_output"], option["MS_uni"])


def fc_main(relevant_features, prot_count, domain_count, seed_proteome, query_proteome, protein_lengths, clan_dict,
            option):
    """Main function,
    manages linearization and scoring, creates the output
    Function calls: w_weighting_constraints(), w_weighting(), su_lin_query_protein(), pb_graphtraversal(),
                    su_query_protein(), su_search_protein(), pb_entire_graphtraversal_priority(),
                    sf_entire_calc_score(), pb_entire_main_nongreedy(), w_weight_const_rescale()

    :param relevant_features: dictionary that contains the rarity of all features from seed and query
    :param prot_count: integer that stores the number of proteins in the reference
    :param domain_count: dictionary that contains the (adjusted) counts of each feature in the reference proteome
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param query_proteome: dictionary that contains the feature architecture of all query proteins
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param clan_dict: dictionary that maps features to clans
    :param option: dictionary that contains the main option variables of FAS
    """
    logging.info("fc_main")

    mode = {}
    weights = {}
    adjusted_weights = {}
    weight_tmp = {}
    out = None
    a_out = None
    if option["output"] == 0 or option["output"] == 2:
        if option["e_output"]:
            a_out = open(option["outpath"] + "_architecture.xml", "w+")
            a_out.write("<?xml version=\"1.0\"?>\n")
            a_out.write("<architectures FAS_version=\"" + str(option["version"]) + "\">\n")
        out = open(option["outpath"] + ".xml", "w+")
        out.write("<?xml version=\"1.0\"?>\n")
        settings_out = {"priority": "off", "weighting": "uniform", "constraints": "none", "overlap":
                        str(option["max_overlap"]) + "/" + str(option["max_overlap_percentage"]), "filters":
                        str(option["efilter"]) + "/" + str(option["inst_efilter"]), "lin": "", "norm": "", "ms":
                        "classic"}

        if not option["classicMS"]:
            settings_out["ms"] = "new"
        if option["timelimit"] == 0 or option["priority_mode"]:
            settings_out["time"] = "off"
        else:
            settings_out["time"] = str(option["timelimit"]) + "s"
        for tool in option["input_linearized"]:
            settings_out["lin"] = settings_out["lin"] + "/" + tool
        settings_out["lin"] = settings_out["lin"].lstrip("/")
        if len(settings_out["lin"]) == 0:
            settings_out["lin"] = "none"
        for tool in option["input_normal"]:
            settings_out["norm"] = settings_out["norm"] + "/" + tool
        settings_out["norm"] = settings_out["norm"].lstrip("/")
        if len(settings_out["norm"]) == 0:
            settings_out["norm"] = "none"
        if option["weight_const"] == 1:
            settings_out["constraints"] = option["constname"]
        if option["priority_mode"]:
            settings_out["priority"] = str(option["priority_threshold"]) + "/" + str(option["max_cardinality"])
        if option["MS_uni"] == 0 and option["weight_correction"]:
            settings_out["weighting"] = option["ref_proteome"].rstrip("/").split("/")[-1] + "/" + \
                                        option["weight_correction"]
        elif option["MS_uni"] == 0:
            settings_out["weighting"] = option["ref_proteome"].rstrip("/").split("/")[-1] + "/none"
        out.write("<out FAS_version=\"" + str(option["version"]) + "\" weighting=\"" + settings_out["weighting"] +
                  "\" constraints=\"" + settings_out["constraints"] + "\" MS=\"" + settings_out["ms"] + "\" priority=\""
                  + settings_out["priority"] + "\" overlap=\"" + settings_out["overlap"] + "\" efilters=\"" +
                  settings_out["filters"] + "\" timelimit=\"" + settings_out["time"] + "\" scoreweights=\"" +
                  str(option["score_weights"]) + "\" cores=\"" + str(option["cores"]) + "\" linearized=\"" +
                  settings_out["lin"] + "\" normal=\"" + settings_out["norm"] + "\">\n")
    if option['pairwise']:
        for pair in option['pairwise']:
            query = pair[1]
            protein = pair[0]
            if option["output"] == 0 or option["output"] == 2:
                out.write(
                    "\t<query id=\"" + query + "\" length=\"" + str(int(protein_lengths["query_" + query])) + "\">\n")
            elif taciturn == 1:
                out = open(option["outpath"], "w+")
            tmp_query = fc_prep_query(query, domain_count, query_proteome, option, clan_dict, protein_lengths)
            query_graph, all_query_paths, lin_query_set, query_features, a_q_f, query_clans, clan_dict = tmp_query[0:7]
            go_priority, domain_count = tmp_query[7:9]

#            query_protein, query_clans, clan_dict = su_query_protein(query, query_proteome, protein_lengths, clan_dict)
            ####
            fc_main_sub(protein, domain_count, seed_proteome, option, protein_lengths, all_query_paths, query_features,
                        out, go_priority, a_q_f, clan_dict, mode, query_graph, query_proteome, query, query_clans)
    else:
        for query in query_proteome:
            if option["output"] == 0 or option["output"] == 2:
                out.write("\t<query id=\"" + query + "\" length=\"" + str(int(protein_lengths["query_" + query])) +
                          "\">\n")
            elif taciturn == 1:
                out = open(option["outpath"], "w+")
            tmp_query = fc_prep_query(query, domain_count, query_proteome, option, clan_dict, protein_lengths)
            query_graph, all_query_paths, lin_query_set, query_features, a_q_f, query_clans, clan_dict = tmp_query[0:7]
            go_priority, domain_count = tmp_query[7:9]

#            query_protein, query_clans, clan_dict = su_query_protein(query, query_proteome, protein_lengths, clan_dict)
            for protein in seed_proteome:
                fc_main_sub(protein, domain_count, seed_proteome, option, protein_lengths, all_query_paths,
                            query_features, out, go_priority, a_q_f, clan_dict, mode, query_graph, query_proteome,
                            query, query_clans)
            if option["output"] == 0 or option["output"] == 2:
                out.write("\t</query>\n")
    if option["output"] == 0 or option["output"] == 2:
        out.write("</out>")
        out.close()
        if option["e_output"]:
            for protein in seed_proteome:
                a_out.write("\t<template id=\"" + protein + "\" length=\"" + str(
                    int(protein_lengths["seed_" + str(protein)])) + "\">\n")
                a_out.write("\t\t<architecture>\n")

                for feature in seed_proteome[protein]:
                    if option["MS_uni"] == 0:
                        a_out.write("\t\t\t<feature type=\"" + feature + "\" evalue=\"" + str(
                            seed_proteome[protein][feature][1]) + "\">\n")
                    else:
                        a_out.write("\t\t\t<feature type=\"" + feature + "\" evalue=\"" + str(
                            seed_proteome[protein][feature][1]) + "\">\n")
                    for instance in seed_proteome[protein][feature][2:]:
                        a_out.write("\t\t\t\t<instance inst_eval=\"" + str(instance[0]) + "\" start=\"" + str(
                            instance[1]) + "\" end=\"" + str(instance[2]) + "\"/>\n")
                    a_out.write("\t\t\t</feature>\n")
                a_out.write("\t\t</architecture>\n")
                a_out.write("\t</template>\n")
            for query in query_proteome:
                a_out.write(
                    "\t<query id=\"" + query + "\" length=\"" + str(int(protein_lengths["query_" + query])) + "\">\n")
                a_out.write("\t\t<architecture>\n")
                for feature in query_proteome[query]:
                    a_out.write("\t\t\t<feature type=\"" + feature + "\" evalue=\"" + str(
                        query_proteome[query][feature][1]) + "\">\n")
                    for instance in query_proteome[query][feature][2:]:
                        a_out.write("\t\t\t\t<instance inst_eval=\"" + str(instance[0]) + "\" start=\"" + str(
                            instance[1]) + "\" end=\"" + str(instance[2]) + "\"/>\n")

                    a_out.write("\t\t\t</feature>\n")
                a_out.write("\t\t</architecture>\n")
                a_out.write("\t</query>\n")
            a_out.write("</architectures>")
            a_out.close()
        if option["feature_info"]:
            f_out = open(option["outpath"] + "_features", "w+")
            f_out.write("#proteins: " + str(prot_count) + "\n")
            tmp = []
            for feature in relevant_features:
                tmp.append((feature, relevant_features[feature]))
            tmp = sorted(tmp, key=lambda x: x[1])
            for feature in tmp:
                f_out.write(feature[0] + "\t" + str(feature[1]) + "\n")
            f_out.close()


def fc_prep_query(query, domain_count, query_proteome, option, clan_dict, protein_lengths):
    go_priority = False
    if option["MS_uni"] == 0:
        domain_count = w_count_add_domains(query, domain_count, query_proteome)
    lin_query_set, query_features, a_q_f, query_clans, clan_dict = su_lin_query_protein(query, query_proteome,
                                                                                        protein_lengths,
                                                                                        clan_dict)
    tmp_query_graph, path_number = pb_region_paths(pb_region_mapper(lin_query_set, query_features,
                                                                    option["max_overlap"],
                                                                    option["max_overlap_percentage"]))
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


def fc_main_sub(protein, domain_count, seed_proteome, option, protein_lengths, all_query_paths, query_features, out,
                go_priority, a_q_f, clan_dict, mode, tmp_query_graph, query_proteome, query, query_clans):
    go_priority_2 = False
    calcstart = time.time()
    timeover = 0
    mode[protein] = 0
    pathcount = 0
    adjusted_weights = None
    weight_tmp = None
    weights = None
    if option["MS_uni"] == 0:
        if option["weight_const"]:
            weights, domain_count = w_weighting_constraints(protein, domain_count, seed_proteome, option)
        else:
            weights, domain_count = w_weighting(protein, domain_count, seed_proteome)
    search_protein, search_features, a_s_f = su_search_protein(protein, seed_proteome, protein_lengths)
    search_protein = tuple(search_protein)
    max_fixture = ([], (0.0, 0.0, 0.0, 0.0, 0, 0, False), [], protein)

    # check for available paths
    # query <--VS-- seed
    # case M2.1: empty(query)

    if int(len(all_query_paths)) == 0 or int(len(query_features) == 0):
        logging.warning("CASE M2.1: No paths (pfam or smart annotations) in query.")
        # case M2.1.1: empty(query)-empty(search)
        # should be the best fix independent from weight
        if int(len(search_features)) == 0:
            logging.warning("CASE M2.1.1: empty vs empty.")
            path = list(a_s_f.keys())
            query_architecture = list(a_q_f.keys())
            score_w = sf_entire_calc_score(path, query_architecture, weights, search_features, a_s_f,
                                           query_features, a_q_f, clan_dict, option)
            mode[protein] = 2
        else:
            # case M2.1.2: empty(query)-graph(search)
            logging.warning("CASE M2.1.2: empty vs graph.")
            tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, [], search_features, weights,
                                                      query_features, seed_proteome, a_s_f, a_q_f, clan_dict,
                                                      query_clans, protein_lengths, "OFF", option)
            path = tmp_path_score[0][0]
            score_w = tmp_path_score[0][1]
            query_architecture = tmp_path_score[0][2]
            mode[protein] = 2

        # set max_fixture according to
        max_fixture = (path, score_w, query_architecture, protein)

    # handle all paths
    # case M2.2: graph(query)
    elif not option["priority_mode"] and not int(len(query_features)) == 0:
        mode[protein] = 0
        stack, jobpaths = pb_create_jobs(tmp_query_graph, option)
        timelimit = 0.0
        if len(stack) > 0:
            jobpool = multiprocessing.Pool(processes=option["cores"])
            if option["timelimit"] > 0:
                timelimit = option["timelimit"] / len(stack)
                func = partial(pb_graph_traversal_sub, search_protein, protein, search_features, weights,
                               query_features, seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                               protein_lengths, timelimit, tmp_query_graph, option)
            else:
                func = partial(pb_graph_traversal_sub, search_protein, protein, search_features, weights,
                               query_features, seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                               protein_lengths, 0, tmp_query_graph, option)
            pool_results = jobpool.map_async(func, stack)
            jobpool.close()
            jobpool.join()
            pool_results.wait()
            for pool_result in pool_results.get():
                if pool_result[1] > timelimit and option["timelimit"] > 0:
                    timeover += 1
                    go_priority_2 = True
                if pool_result[0][1][3] >= max_fixture[1][3]:
                    max_fixture = deepcopy(pool_result[0])
        timecheck = time.time()
        if len(jobpaths) > 0:
            # case M2.2.1: graph(query)-empty(search)
            if int(len(search_features)) == 0:
                for jobpath in jobpaths:
                    logging.warning("CASE M2.2: graph.")
                    pathcount += 1
                    logging.warning("CASE M2.2.1: graph vs empty.")
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
                logging.warning("CASE M2.2.2: graph vs graph.")
                jobs = []
                jobcount = (int(len(jobpaths) / option["cores"]), len(jobpaths) % option["cores"])
                if jobcount[0] == 0:
                    timelimit = (option["timelimit"] - (timecheck - calcstart)) / float(jobcount[1])
                    for i in range(jobcount[1]):
                        jobs.append([jobpaths[-(i + 1)]])
                else:
                    timelimit = (option["timelimit"] - (timecheck - calcstart)) / float(option["cores"])
                    counter = 0
                    for i in range(option["cores"]):
                        jobs.append(jobpaths[counter: (i + 1) * jobcount[0]])
                        counter = (i + 1) * jobcount[0]
                    for i in range(jobcount[1]):
                        jobs[i].append(jobpaths[-(i + 1)])
                func = partial(pb_calc_sub, search_protein, protein, search_features, weights, query_features,
                               seed_proteome, a_s_f, a_q_f, clan_dict, query_clans, protein_lengths, timelimit,
                               option)
                pool_results = jobpool.map_async(func, jobs)
                jobpool.close()
                jobpool.join()
                pool_results.wait()
                for pool_result in pool_results.get():
                    if pool_result[1] > timelimit:
                        timeover += 1
                        go_priority_2 = True
                    if pool_result[0][1][3] >= max_fixture[1][3]:
                        max_fixture = deepcopy(pool_result[0])
    if go_priority or go_priority_2:
        mode[protein] = 1
        all_query_paths = []
        priority_list = []
        for domain_type in query_proteome[query]:
            if query_proteome[query][domain_type][0]:
                priority_list.append(domain_type)
        priority_list.append("NONE")
        for domain_type in priority_list:
            all_query_paths += (
                pb_entire_graphtraversal_priority(tmp_query_graph, domain_type, protein, 1, search_features,
                                                  weights, query_features, seed_proteome, a_s_f, a_q_f,
                                                  clan_dict, query_clans, protein_lengths, option))
            logging.info("domain_type: " + str(domain_type))
            logging.debug("all_query_paths: " + str(all_query_paths))
            # max fixture of (search_path, score, query_path)
    if option["priority_mode"] or go_priority_2:
        for query_path in all_query_paths:
            logging.warning("CASE M2.2: graph.")
            pathcount += 1

            # case M2.2.1: graph(query)-empty(search)
            if int(len(search_features)) == 0:
                logging.warning("CASE M2.2.1: graph vs empty.")
                # special case: protein with no pfam or smart domains
                # get score for a_s_f and query_path directly
                path = list(a_s_f.keys())
                query_path_ad = query_path + list(a_q_f.keys())
                score_w = sf_entire_calc_score(path, query_path_ad, weights, search_features, a_s_f,
                                               query_features, a_q_f, clan_dict, option)
            else:
                # case M2.2.2 graph(query)-graph(search)
                # regular traversal of graph based on search_protein
                logging.warning("CASE M2.2.2: graph vs graph.")
                if option["priority_mode"]:
                    tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, query_path,
                                                              search_features, weights, query_features,
                                                              seed_proteome, a_s_f, a_q_f, clan_dict,
                                                              query_clans, protein_lengths, "OFF", option)
                else:
                    tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, query_path,
                                                              search_features, weights, query_features,
                                                              seed_proteome, a_s_f, a_q_f, clan_dict,
                                                              query_clans, protein_lengths, "OVER", option)
                path = tmp_path_score[0][0]
                score_w = tmp_path_score[0][1]
                query_path_ad = tmp_path_score[0][2]
                mode[protein] = tmp_path_score[1]
                logging.debug("tmp_path_score " + str(tmp_path_score))  # search path, score, mode

            # check for max scoring fixture of path and query_path
            if (score_w[4] >= max_fixture[1][4] and score_w[3] == max_fixture[1][3]) or \
                    score_w[3] >= max_fixture[1][3]:
                max_fixture = (path, score_w, query_path_ad, protein)

        logging.info("Found: " + str(pathcount) + " path(s) for query.")
        logging.debug("Path fixture: " + str(max_fixture))

    timecheck = time.time()
    score = max_fixture[1]
    if option["output"] == 0 or option["output"] == 2:
        try:
            if mode[protein] == 2:
                mode_out = "greedy"
            elif mode[protein] == 1:
                mode_out = "priority"
            else:
                mode_out = "exhaustive"
        except KeyError:
            mode_out = "default"
        if go_priority or go_priority_2:
            mode_out += "/priority"
        else:
            mode_out += "/exhaustive"
        runtime = str(round(timecheck - calcstart, 4)) + "s"
        if timeover >= 1:
            runtime += "/exceeded_timelimit"
        out.write("\t\t<template id=\"" + protein + "\" score=\"" + str(score[3]) + "\" MS=\"" + str(
            score[0]) + "\" PS=\"" + str(score[1]) + "\" CS=\"" + str(score[2]) + "\" LS=\"" + str(
            score[6]) + "\" length=\"" + str(int(protein_lengths["seed_" + str(protein)])) + "\" mode=\"" +
                  mode_out + "\" calculationTime=\"" + runtime + "\" >\n")

    best_template_path = []
    best_query_path = []

    for feature in max_fixture[0]:
        if feature in search_features:
            logging.debug(str(search_features[feature][0]) + " " + str(search_features[feature][1]) + " " + str(
                search_features[feature][2]) + " " + str(search_features[feature][3]))
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
    if option["output"] == 0 or option["output"] == 2:
        path_tmp = {}
        path_tmp_query = {}
        scale = 0
        if option["weight_const"] == 1:
            path_tmp2 = []
            for feature in best_template_path:
                if feature[0] not in path_tmp2:
                    path_tmp2.append(feature[0])
            adjusted_weights = w_weight_const_rescale(path_tmp2, weights, search_features, True, option)
            weight_tmp = {}
            for adj_feature in adjusted_weights:
                weight_tmp[adj_feature] = weights[adj_feature]
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
        logging.debug("path_tmp: " + str(path_tmp))

        # unweighted case
        if option["MS_uni"] == 0:
            if scale > 0:
                scale = 1.0 / float(scale)
            else:
                scale = 1.0
        # print path
        out.write("\t\t\t<template_path>\n")
        for feature in path_tmp:
            if option["MS_uni"] == 0:
                out.write("\t\t\t\t<feature type=\"" + feature + "\" corrected_weight=\"" + str(
                    weights[feature] * scale) + "\">\n")
            else:
                out.write("\t\t\t\t<feature type=\"" + feature + "\">\n")
            for tmp_inst in path_tmp[feature]:
                out.write("\t\t\t\t\t<instance start=\"" + str(tmp_inst[0]) + "\" end=\"" + str(
                    tmp_inst[1]) + "\"/>\n")
            out.write("\t\t\t\t</feature>\n")
        out.write("\t\t\t</template_path>\n")
        out.write("\t\t\t<query_path>\n")
        for feature in path_tmp_query:
            out.write("\t\t\t\t<feature type=\"" + feature + "\">\n")
            for tmp_inst in path_tmp_query[feature]:
                out.write("\t\t\t\t\t<instance start=\"" + str(tmp_inst[0]) + "\" end=\"" + str(
                    tmp_inst[1]) + "\"/>\n")
            out.write("\t\t\t\t</feature>\n")
        out.write("\t\t\t</query_path>\n")
        out.write("\t\t</template>\n")
        if option["weight_const"] == 1:
            for adj_feature in adjusted_weights:
                weights[adj_feature] = weight_tmp[adj_feature]
    if option["output"] == 1 or option["output"] == 2:
        print(score[3])
        # hidden
        if taciturn == 1 and option["output"] != 2:
            out.write(protein + "\t" + str(score[3]) + "\n")

# Start Up Functions <su>


def su_query_protein(protein_id, query_proteome, protein_lengths, clan_dict):
    """Initializes variables for the current query protein

    :param protein_id: String that contains the identifier of the query protein
    :param query_proteome: dictionary that contains the feature architecture of all query proteins
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param clan_dict: dictionary that maps features to clans
    :return: query_protein, query_clans, clan_dict
    """
    query_protein = {}
    query_clans = {}

    for feature in query_proteome[protein_id]:
        query_protein[feature] = []
        clan = clan_dict[feature]
        if clan in query_clans:
            query_clans[clan] += len(query_proteome[protein_id][feature]) - 2
        else:
            query_clans[clan] = len(query_proteome[protein_id][feature]) - 2
        for instance in query_proteome[protein_id][feature][2:]:
            position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(
                protein_lengths["query_" + str(protein_id)])
            position = round(position, 4)
            query_protein[feature].append(position)
    return query_protein, query_clans, clan_dict


def su_lin_query_protein(protein_id, query_proteome, protein_lengths, clan_dict):
    """Initializes variables for the current query protein (linearization version)

    :param protein_id: String that contains the identifier of the query protein
    :param query_proteome: dictionary that contains the feature architecture of all query proteins
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param clan_dict: dictionary that maps features to clans
    :return:lin_query_protein, query_features, a_q_f(additional[not linearized] query features), query_clans, clan_dict
    """
    logging.info("su_lin_query_protein")

    lin_query_protein = []
    query_clans = {}
    query_features = {}
    a_q_f = {}
    tmp = []
    i = 0

    for feature in query_proteome[protein_id]:
        clan = clan_dict[feature]
        if clan in query_clans:
            query_clans[clan] += len(query_proteome[protein_id][feature]) - 2
        else:
            query_clans[clan] = len(query_proteome[protein_id][feature]) - 2

        if query_proteome[protein_id][feature][0]:
            for instance in query_proteome[protein_id][feature][2:]:
                position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(
                    protein_lengths["query_" + str(protein_id)])
                position = round(position, 4)
                key = "F_" + str(i)
                query_features[key] = (feature, position, instance[1], instance[2])
                tmp.append((key, instance[1]))
                i += 1
        else:
            for instance in query_proteome[protein_id][feature][2:]:
                position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(
                    protein_lengths["query_" + str(protein_id)])
                position = round(position, 4)
                key = "O_" + str(i)
                a_q_f[key] = (feature, position, instance[1], instance[2])
                i += 1

    # sort  instances
    tmp2 = sorted(tmp, key=itemgetter(1))
    for x in tmp2:
        lin_query_protein.append(x[0])

    return lin_query_protein, query_features, a_q_f, query_clans, clan_dict


def su_search_protein(protein_id, seed_proteome, protein_lengths):
    """Initializes variables for the current seed protein

    :param protein_id: String that contains the identifier of the seed protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :return: search_protein, search_features, a_s_f (additional [not linearized] seed features)
    """
    search_features = {}
    search_protein = []
    a_s_f = {}
    tmp = []
    i = 0
    for feature in seed_proteome[protein_id]:
        # True if feature is going to be linearized
        if seed_proteome[protein_id][feature][0]:
            for instance in seed_proteome[protein_id][feature][2:]:
                position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(
                    protein_lengths["seed_" + str(protein_id)])
                position = round(position, 4)
                key = "F_" + str(i)
                search_features[key] = (feature, position, instance[1], instance[2])
                tmp.append((key, instance[1]))
                i += 1
        else:
            for instance in seed_proteome[protein_id][feature][2:]:
                position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(
                    protein_lengths["seed_" + str(protein_id)])
                position = round(position, 4)
                key = "O_" + str(i)
                a_s_f[key] = (feature, position, instance[1], instance[2])
                i += 1

    # sort  instances
    tmp2 = sorted(tmp, key=itemgetter(1))
    for x in tmp2:
        search_protein.append(x[0])
    return tuple(search_protein), search_features, a_s_f


def su_set_path(jobname, path):
    """Checks default output path if only a name for job is given

    :param jobname: String contains name of output file
    :param path: String contains output path
    :return: jobname
    """
    if not os.path.exists(path + "/out/"):
        os.makedirs(path + "/out/")
    if os.path.exists(path + "/out/" + jobname + ".xml"):
        i = 1
        while os.path.exists(path + "/out/" + jobname + "_" + str(i) + ".xml"):
            i += 1
        jobname = jobname + "_" + str(i)
    return jobname


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


def pb_graph_traversal_sub(search_protein, protein, search_features, weights, query_features, seed_proteome, a_s_f,
                           a_q_f, clan_dict, query_clans, protein_lengths, timelimit, query_graph, option, stack):
    """Sub function for multi-processing: This function does a depth-first search in the feature graph. The starting
     point in the graph is given in the stack.

    :param search_protein: contains feature architecture of the current seed protein [dictionary]
    :param protein: String that contains the identifier of the seed protein
    :param search_features: all features to be linearized in the seed protein
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: (Pfam)-clans of the current query protein
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param timelimit: the amount of time this job is allowed to take before the calculation is stopped and priority
                      mode is used instead
    :param option: dictionary that contains the main option variables of FAS
    :param stack: contains the starting point in the graph and sub-path
    :param query_graph: graph of the query architecture
    :return: max_fixture, timecheck - tmp_calcstart (time taken for calculation)
    """
    tmp_calcstart = time.time()
    max_fixture = ([], (0.0, 0.0, 0.0, 0.0, 0, 0, False), [], protein)
    v_stack, p_stack = [stack[0]], [stack[1]]
    timecheck = time.time()
    while v_stack and (timecheck - tmp_calcstart <= timelimit or option["timelimit"] == 0):
        query_path, v_stack, p_stack = pb_graphtraversal(query_graph, v_stack, p_stack, option)
        # case M2.2.1: graph(query)-empty(search)
        timecheck = time.time()
        if int(len(search_features)) == 0:
            logging.warning("CASE M2.2.1: graph vs empty.")
            # special case: protein with no pfam or smart domains
            # get score for a_s_f and query_path directly
            path = list(a_s_f.keys())
            query_path_ad = query_path + list(a_q_f.keys())
            score_w = sf_entire_calc_score(path, query_path_ad, weights, search_features, a_s_f,
                                           query_features, a_q_f, clan_dict, option)
        else:
            # case M2.2.2 graph(query)-graph(search)
            # regular traversal of graph based on search_protein
            logging.warning("CASE M2.2.2: graph vs graph.")
            if option["timelimit"] >= 1:
                tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, query_path, search_features,
                                                          weights, query_features, seed_proteome, a_s_f, a_q_f,
                                                          clan_dict, query_clans, protein_lengths,
                                                          timelimit - (timecheck - tmp_calcstart), option)
            else:
                tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, query_path, search_features,
                                                          weights, query_features, seed_proteome, a_s_f, a_q_f,
                                                          clan_dict, query_clans, protein_lengths, "OFF", option)
            path = tmp_path_score[0][0]
            score_w = tmp_path_score[0][1]
            query_path_ad = tmp_path_score[0][2]
            logging.debug("tmp_path_score " + str(tmp_path_score))  # search path, score, mode

        # check for max scoring fixture of path and query_path
        if ((not max_fixture[1][5]) and (score_w[4] >= max_fixture[1][4]) and score_w[3] <= max_fixture[1][3]) or \
                score_w[3] >= max_fixture[1][3]:
            max_fixture = (path, score_w, query_path_ad, protein)
        timecheck = time.time()

    return max_fixture, timecheck - tmp_calcstart


def pb_calc_sub(search_protein, protein, search_features, weights, query_features, seed_proteome, a_s_f, a_q_f,
                clan_dict, query_clans, protein_lengths, timelimit, option, jobpaths):
    """Sub function for multiprocessing [2]: Goes through a list of (query) paths an evaluates them.

    :param search_protein: contains feature architecture of the current seed protein [dictionary]
    :param protein: String that contains the identifier of the seed protein
    :param search_features: all features to be linearized in the seed protein
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: (Pfam)-clans of the current query protein
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param timelimit: the amount of time this job is allowed to take before the calculation is stopped and priority
                      mode is used instead
    :param option: dictionary that contains the main option variables of FAS
    :param jobpaths: a group of (query) paths to be evaluated
    :return: max_fixture, tmpcheck - tmpstart (time taken for calculation)
    """
    max_fixture = ([], (0.0, 0.0, 0.0, 0.0, 0, 0, False), [], protein)
    tmpstart = time.time()
    for jobpath in jobpaths:
        tmpcheck = time.time()
        tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, jobpath, search_features, weights,
                                                  query_features, seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                                                  protein_lengths, timelimit - (tmpcheck - tmpstart), option)
        path = tmp_path_score[0][0]
        score_w = tmp_path_score[0][1]
        query_path_ad = tmp_path_score[0][2]
        logging.debug("tmp_path_score " + str(tmp_path_score))  # search path, score, mode
        if ((not max_fixture[1][5]) and (score_w[4] >= max_fixture[1][4]) and score_w[3] <= max_fixture[1][3]) or \
                score_w[3] >= max_fixture[1][3]:
            max_fixture = (path, score_w, query_path_ad, protein)
    tmpcheck = time.time()
    return max_fixture, tmpcheck - tmpstart


def pb_entire_main_nongreedy(search_protein, protein_id, query_path, search_features, weights, query_features,
                             seed_proteome, a_s_f, a_q_f, clan_dict, query_clans, protein_lengths, tmp_timelimit,
                             option):
    """Main Path-building function,
    creates graph with all paths(for seed/search protein), checks priority mode activation if necessary retrieves best
    path from graph traversal function, returns best path, score and mode
    Function calls: pb_region_mapper(), pb_region_paths(), pb_entire_priority_mode(),
                    pb_entire_graphtraversal()

    :param tmp_timelimit: timelimit for exhaustive mode
    :param search_protein: contains feature architecture of the current seed protein [dictionary]
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
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param option: dictionary that contains the main option variables of FAS
    :return: path_score, mode[priority or extensive]
    """
    logging.info("pb_entire_main_nongreedy")
    path_score = 0
    mode = 0
    priority_check = False

    region = pb_region_mapper(list(search_protein), search_features, option["max_overlap"],
                              option["max_overlap_percentage"])
    search_graph, path_number = pb_region_paths(region)
    logging.debug(region)
    logging.debug(search_graph)

    if (int(len(search_features)) >= int(option["priority_threshold"]) and option["priority_mode"]) \
            or tmp_timelimit == "OVER":
        mode = 1
        priority_check = True
        # checking best path for all domain types
        path_score = pb_entire_priority_mode(protein_id, query_path, search_graph, search_features, weights,
                                             query_features, seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                                             protein_lengths, option)
    elif option["priority_mode"]:
        # PRIORITY CHECK "2": for every protein in seed_proteome
        tmp_path_set_size = len(pb_graphtraversal(search_graph, [], [], option))
        logging.debug("cardinality of tmp graph: " + str(tmp_path_set_size))
        if int(tmp_path_set_size) > int(option["max_cardinality"]):
            mode = 1
            priority_check = True
            path_score = pb_entire_priority_mode(protein_id, query_path, search_graph, search_features, weights,
                                                 query_features, seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                                                 protein_lengths, option)

    if not priority_check:
        mode = 0
        path_score = pb_entire_graphtraversal(search_graph, query_path, search_features, weights, query_features,
                                              a_s_f, a_q_f, clan_dict, tmp_timelimit, option)
    return path_score, mode


def pb_entire_priority_mode(protein, query_path, search_graph, search_features, weights, query_features, seed_proteome,
                            a_s_f, a_q_f, clan_dict, query_clans, protein_lengths, option):
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
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param option: dictionary that contains the main option variables of FAS
    :return: best_path
    """
    logging.info("pb_entire_priority_mode")

    # best_path (Path, (MS_score(float), PS_score(float), CS_score(float), final_score(float), path_weight(float),
    # common_feature(bool)), QPath)
    best_path = ("NULL", (0.0, 0.0, 0.0, 0.0, 0.0, False), "NULL")
    priority_list = []
    for domain_type in seed_proteome[protein]:
        if seed_proteome[protein][domain_type][0]:
            priority_list.append(domain_type)
    priority_list.append("NONE")
    for domain_type in priority_list:
        path_score_w = pb_entire_graphtraversal_priority(search_graph, domain_type, query_path, 0, search_features,
                                                         weights, query_features, seed_proteome, a_s_f, a_q_f,
                                                         clan_dict, query_clans, protein_lengths, option)
        if path_score_w[1][5]:
            if path_score_w[1][3] >= best_path[1][3]:
                best_path = (path_score_w[0], path_score_w[1], path_score_w[2])
        else:
            # if so far no path with common features found AND the weight is higher
            if (not best_path[1][5]) and (path_score_w[1][4] >= best_path[1][4]):
                best_path = (path_score_w[0], path_score_w[1], path_score_w[2])
    return best_path


def pb_entire_graphtraversal(search_graph, query_path, search_features, weights, query_features, a_s_f, a_q_f,
                             clan_dict, timelimit, option):
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
    :param timelimit: timelimit for exhaustive mode
    :return: best_path
    """
    logging.info("pb_entire_graphtraversal")
    logging.debug("search graph: " + str(search_graph))

    calcstart = time.time()
    v_stack = ["START"]
    p_stack = []
    best_path = ([], (0.0, 0.0, 0.0, 0.0, 0.0, False, 0.0), [])
    first = 1
    timecheck = time.time()
    if timelimit == "OFF":
        timelimit_off = True
        timelimit = 0.0
    else:
        timelimit_off = False
    if len(query_path) > 0:
        while v_stack and (timecheck - calcstart < timelimit or first or timelimit_off):
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
                    logging.info("search path " + str(path_ad) + " in pb_entire_graphtraversal")
                    logging.debug("Score info: " + str(score_w) + " for query_path " + str(
                        query_path_ad) + " in pb_entire_graphtraversal")
                    if (score_w[4] >= best_path[1][4] and score_w[3] == best_path[1][3]) or \
                       score_w[3] >= best_path[1][3]:
                        best_path = (path_ad, score_w, query_path_ad)
                        logging.debug("new best_path by weight: " + str(best_path))
                    first = 0
                else:
                    v_stack.append(next_vertex)
                    p_stack.append(path + [next_vertex])
            timecheck = time.time()
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
                    logging.info("search path " + str(path_ad) + " in pb_entire_graphtraversal")
                    logging.debug("Score info: " + str(score_w) + " for query_path " + str(
                        query_path_ad) + " in pb_entire_graphtraversal")
                    if (score_w[4] >= best_path[1][4] and score_w[3] == best_path[1][3]) or \
                       score_w[3] >= best_path[1][3]:
                        best_path = (path_ad, score_w, query_path_ad)
                        logging.debug("new best_path by weight: " + str(best_path))
                    first = 0
                else:
                    v_stack.append(next_vertex)
                    p_stack.append(path + [next_vertex])
    logging.info("return: " + str(best_path))
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
    logging.debug(graph)
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
                    logging.info("cardinality: " + str(len(paths)))
                    if len(paths) > option["max_cardinality"]:
                        return paths
                else:
                    return path, v_stack, p_stack
            else:
                v_stack.append(next_vertex)
                p_stack.append(path + [next_vertex])
    logging.debug("returning paths: " + str(paths))
    logging.info("cardinality: " + str(len(paths)))
    return paths


def pb_entire_graphtraversal_priority(search_graph, priority, query_path, mode, search_features, weights,
                                      query_features, seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                                      protein_lengths, option):
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
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param option: dictionary that contains the main option variables of FAS
    :return: paths or best_path
    """
    logging.info("pb_entire_graphtraversal_priority")
    logging.debug(
        "search_features: " + str(search_features) + "\nsearch_graph: " + str(search_graph) + "\tpriority: " + str(
            priority) + "\ta_s_f: " + str(a_s_f) + "\tquery_path: " + str(query_path) + "\ta_q_f: " + str(a_q_f))
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
                    logging.debug("path_ad: " + str(path_ad))
                    logging.debug("query_path_ad: " + str(query_path_ad))

                    score_w = sf_entire_calc_score(path_ad, query_path_ad, weights, search_features, a_s_f,
                                                   query_features, a_q_f, clan_dict, option)
                    if (score_w[4] >= best_path[1][4] and score_w[3] == best_path[1][3]) or \
                       score_w[3] >= best_path[1][3]:
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
                        # baustelle: fixed: path elongation without additives
                        score = sf_calc_score(path + [next_vertex], protein, weights, search_features,
                                              query_features, seed_proteome, clan_dict, query_clans,
                                              protein_lengths, option)
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
                                          seed_proteome, clan_dict, query_clans, protein_lengths, option)
                    if score[3] >= best_priority_bridger[1][3]:
                        best_priority_bridger = (next_vertex, score)
                v_stack.append(best_priority_bridger[0])
                p_stack.append(path + [best_priority_bridger[0]])

            elif mode == 0:
                # some kind of greedy strategy: if feature type (p priority) not found
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
        logging.debug("returning: " + str(paths))
        logging.debug("cardinality: " + str(len(paths)))
        return paths
    elif mode == 0:
        logging.debug(best_path)
        return best_path


# start #
def main():
    version = "1.0.0"
#    greedyfas_path = inspect.getfile(inspect.currentframe())
    expath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

    parser = argparse.ArgumentParser(description="You are running greedyFAS.py version " + str(version) + ".")
    parser.add_argument("-q", "--query", default=None, type=str, required=True,
                        help="Path to the folder containing the xml-files of the query protein set")
    parser.add_argument("-s", "--seed", default=None, type=str, required=True,
                        help="Path to the folder containing the xml-files of the template protein set")
    parser.add_argument("--query_id", default=None, nargs='*', type=str,
                        help="Choose a single protein from the query input for calculation, by default this is off "
                             "(all proteins are used for calculation)")
    parser.add_argument("--seed_id", default=None, nargs='*', type=str,
                        help="Choose a single protein from the seed input for calculation, by default this is off (all "
                             "proteins are used for calculation)")
    parser.add_argument("-w", "--score_weights", nargs=3, default=[0.7, 0.0, 0.3], type=float,
                        help="Defines how the three scores MS, CS and PS are weighted, takes three float arguments, sum"
                             " should be 1.0, the default is 0.7, 0.0, 0.3")
    parser.add_argument("-j", "--jobname", default="out", type=str,
                        help="Defines the name for the outputfile, can also be used to define the output path, if no "
                             "path is given the output will be named out.xml and out_architecture.xml")
    parser.add_argument("-r", "--ref_proteome", default=None, type=str,
                        help="Path to a reference proteome which can be used for the weighting of features, "
                             "by default there is no reference proteome used, the weighting will be uniform")
    parser.add_argument("-a", "--raw_output", default=0, type=int,
                        help="If set to 1, the FAS score will be printed to STDOUT. If 0, scores will be printed into "
                             "output file (XML format). If 2, both output variants are conducted.")
    parser.add_argument("--no_priority_check", action="store_false",
                        help="deactivate priority checks prior to the calculation, skipping all priority checks "
                             "and going through all proteins exhaustively until time-limit is reached where it will "
                             "switch to priority mode, NOT RECOMMENDED as this may result in a very long runtime")
    parser.add_argument("-t", "--priority_threshold", type=int, default=50,
                        help="Change to define the feature number threshold for activating priority mode in the path "
                             "evaluation.")
    parser.add_argument("-m", "--max_cardinality", default=5000, type=int,
                        help="Change to define the threshold for the maximal cardinality of feature paths in a graph. "
                             "If max. cardinality is exceeded the priority mode will be used to for the path "
                             "evaluation.")
    parser.add_argument("-l", "--loglevel", default="ERROR", type=str,
                        help="Change the verbosity of the standard output to the screen. "
                             "Levels: DEBUG, INFO, WARNING, ERROR, CRITICAL")
    parser.add_argument("-f", "--efilter", default="0.001", type=float,
                        help="E-value filter for hmm based search methods (feature based/complete sequence).")
    parser.add_argument("-i", "--inst_efilter", default="0.01", type=float,
                        help="E-value filter for hmm based search methods (instance based/complete sequence).")
    parser.add_argument("-g", "--weightcorrection", default="loge", type=str,
                        help="Function applied to the frequency of feature types during weighting, options are "
                             "linear(no function), loge(natural logarithm[Default]), log10(base-10 logarithm), "
                             "root4(4th root) and root8(8th root).")
    parser.add_argument("-x", "--weight_constraints", default=None, type=str,
                        help="Apply weight constraints via constraints file, by default there are no constraints.")
    parser.add_argument("-d", "--featuretypes", default=None, type=str,
                        help="inputfile that contains the tools/databases used to predict features")
    parser.add_argument("-e", "--no_arch", action="store_false",
                        help="deactivate creation of _architecture file")
    parser.add_argument("--feature_info", action="store_true",
                        help="create a file with information on the abundance of all seed and query features in the "
                             "reference")
    parser.add_argument("--bidirectional", action="store_true",
                        help="calculate both scoring directions (separate files), creates csv file with combined "
                             "scores")
    parser.add_argument("-c", "--max_overlap", dest="max_overlap", default=0, type=int,
                        help="maximum size overlap allowed, default is 0 amino acids")
    parser.add_argument("--max_overlap_percentage", dest="max_overlap_percentage", default=0.4, type=float,
                        help="defines how much percent of a feature the overlap is allowed to cover, default "
                             "is 0.4 (40%%)")
    parser.add_argument("--classicMS", action="store_true",
                        help="activate classic MS score")
    parser.add_argument("--timelimit", dest="timelimit", default=7200, type=int,
                        help="Sets a maximum time-limit in seconds for the calculation between a pair of proteins,"
                             "default is 2 hours after which it will stop, set to 0 to deactivate; As FAS divides this "
                             "time among multiple processes, this limit does not necessarily represent the actual "
                             "runtime,especially if multiple cores are used")
    parser.add_argument("--cores", dest="cores", default=1, type=int,
                        help="Number of cores available for calculation, only useful when not using priority_mode")
    parser.add_argument("--ref_2", dest="ref_2", default=None, type=str,
                        help="Give a second reference for bidirectional mode, does not do anything if bidirectional "
                             "mode is not active or if no main reference was given")
    parser.add_argument("--phyloprofile", dest="phyloprofile", default=None, type=str,
                        help="activate phyloprofile output, needs mapping file for all query proteins, single seed "
                             "only, will run with more but output won't work without editing")
    parser.add_argument("--domain", dest="domain", action="store_true",
                        help="activate domain tabular output")
    parser.add_argument("--pairwise", dest="pairwise", default=None, type=str,
                        help="deactivate all against all comparison, needs a pairing file with the ids that should be "
                             "compared (one pair per line tab seperated)")
    options = parser.parse_args()

    # READ OPTIONS #
    option_dict = {"p_path": options.seed, "s_path": options.query, "ref_proteome": options.ref_proteome,
                   "weight_const": False, "version": version, "seed_id": options.seed_id, "query_id": options.query_id,
                   "priority_mode": options.no_priority_check, "priority_threshold": options.priority_threshold,
                   "max_cardinality": options.max_cardinality, "efilter": options.efilter, "cores": options.cores,
                   "inst_efilter": options.inst_efilter, "e_output": options.no_arch,
                   "feature_info": options.feature_info, "bidirectional": options.bidirectional,
                   "max_overlap": options.max_overlap, "classicMS": options.classicMS, "timelimit": options.timelimit,
                   "ref_2": options.ref_2, "phyloprofile": options.phyloprofile, "domain": options.domain}
    loglevel = options.loglevel.upper()
    if options.phyloprofile and not os.path.exists(options.phyloprofile):
        raise Exception(str(options.phyloprofile) + " does not exist")
    try:
        option_dict["score_weights"] = []
        test = 0.0
        for weight in options.score_weights:
            option_dict["score_weights"].append(weight)
            test += float(weight)
        if test != 1.0:
            raise Exception("The sum of all values in -w [--score_weights] should be 1.0")
        option_dict["score_weights"] = tuple(option_dict["score_weights"])
    except ValueError:
        raise Exception(str(options.score_weights) + " is not a valid input for -w [--score_weights]")
    if option_dict["cores"] < 1:
        raise Exception(str(options.cores) + " is not a valid input for [--cores], must be 1 or more")
    if not (option_dict["timelimit"] >= 0):
        raise Exception(str(options.timelimit) + " is not a valid input for [--timelimit], must be integer higher than"
                                                 " 0, or 0 to deactivate")
    if not option_dict["ref_proteome"]:
        option_dict["MS_uni"] = 1
        option_dict["weight_correction"] = None
    else:
        option_dict["MS_uni"] = 0
        if options.weightcorrection == "log10":
            option_dict["weight_correction"] = "log10"
        elif options.weightcorrection == "loge":
            option_dict["weight_correction"] = "loge"
        elif options.weightcorrection == "root4":
            option_dict["weight_correction"] = "root4"
        elif options.weightcorrection == "root8":
            option_dict["weight_correction"] = "root8"
        elif options.weightcorrection == "linear":
            option_dict["weight_correction"] = None
        else:
            raise Exception(str(options.weightcorrection) + "is not a valid input for -g [--weightcorrection]")
    if 0 <= options.raw_output <= 2:
        option_dict["output"] = options.raw_output
    else:
        raise Exception(str(options.raw_output) + "is not a valid input for -a [--raw_output]")
    if options.ref_proteome and options.weight_constraints:
        option_dict["weight_const"] = True
    elif options.weight_constraints:
        raise Exception("[--weight_constraints] only works with a reference proteome")
    if 0.0 <= float(options.max_overlap_percentage) <= 1.0:
        option_dict["max_overlap_percentage"] = float(options.max_overlap_percentage)
    else:
        raise Exception("[--max_overlap_percentage] should be between 0.0 and 1.0")

    # SETUP LOGGING OPTIONS
    # possible logging configurations
    # logging into stdout with time stamp:
    logging.basicConfig(level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')
    # logging into file with line number:
    # logging.basicConfig(filename='testlog.log', filemode='w', level=loglevel,
    #                    format='%(lineno)s - %(levelname)s - %(message)s')
    # logging into stdout with line number:
    # logging.basicConfig(level=loglevel, format='%(lineno)s - %(levelname)s - %(message)s')

    logging.info(
        'greedyFAS.py started with options: priority_threshold=' + str(option_dict["priority_threshold"]) +
        ', log_level=' + str(loglevel))
    logging.info(
        'score_weights are set to: ' + str(option_dict["score_weights"][0]) + " " + str(option_dict["score_weights"][1])
        + " " + str(option_dict["score_weights"][2]))
    logging.info('ref_proteome is set to: ' + str(option_dict["ref_proteome"]))

    if options.jobname.find("/") != -1:
        option_dict["outpath"] = options.jobname

    else:
        option_dict["outpath"] = expath + "/out/" + su_set_path(options.jobname, expath)

    if option_dict["weight_const"]:
        option_dict["constraints"] = constraints_in(options.weight_constraints)
        option_dict["constname"] = options.weight_constraints.rstrip("/").split("/")[-1]

    if options.featuretypes is not None:
        option_dict = featuretypes(options.featuretypes, option_dict)
    else:
        option_dict["input_linearized"] = ["pfam", "smart"]
        option_dict["input_normal"] = ["flps", "coils", "seg", "signalp", "tmhmm"]
    if options.pairwise:
        option_dict["pairwise"] = read_pairwise(options.pairwise)
    else:
        option_dict["pairwise"] = None

    fc_start(option_dict)


if __name__ == '__main__':
    main()
