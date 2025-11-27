#!/bin/env python

#######################################################################
# Copyright (C) 2025 Vinh Tran
#
# This file is part of FAS.
#
#  It is used to convert InterPro TSV output to annotation JSON file
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
import os
import re
import sys
from Bio import SeqIO
from importlib.metadata import version, PackageNotFoundError
import greedyFAS.annoFAS.annoModules as annoModules


def get_seq_len(seq_file):
    id_to_len = {rec.id: len(rec.seq) for rec in SeqIO.parse(seq_file, "fasta")}
    return(id_to_len)

def readDatFile(toolPath):
    """
    Reads the Pfam-A.hmm.dat file to get clan, ID and domain len
    """
    datFile = f'{toolPath}/Pfam/Pfam-hmms/Pfam-A.hmm.dat'
    if not os.path.exists(datFile):
        print(f'{datFile} does not exist')
        return {}

    with open(datFile, 'r') as file:
        datFileR = file.read()
    blocks = datFileR.split('//')
    pfamDict = {}
    for bl in blocks:
        # Extract accession (MUST exist, otherwise skip)
        m_ac = re.search(r'^\s*#=GF AC\s+(\S+)', bl, re.MULTILINE)
        if not m_ac:
            continue
        acc_full = m_ac.group(1)        # e.g. PF00569.20
        acc = acc_full.split('.')[0]    # e.g. PF00569  ← remove version
        # Extract other fields safely
        m_id  = re.search(r'^\s*#=GF ID\s+(\S+)', bl, re.MULTILINE)
        m_cl  = re.search(r'^\s*#=GF CL\s+(\S+)', bl, re.MULTILINE)
        m_len = re.search(r'^\s*#=GF ML\s+(\S+)', bl, re.MULTILINE)
        pfamDict[acc] = {
            'id':    m_id.group(1)  if m_id else acc,
            'clan':  m_cl.group(1)  if m_cl else 'NA',
            'length': m_len.group(1) if m_len else 'NA'
        }
    return pfamDict


def read_tsv_to_features(ignore_lines, columns, tool_names, tool_column, path, seq_file, toolPath):
    """
    Reads a TSV file and generates the full feature dictionary with extended metrics.
    """
    features = {}
    interpro = {}
    clan_dict = {}
    domain_len_dict = {}
    pfamDict = readDatFile(toolPath)

    if 'pfam' in tool_names:
        pfam_len = annoModules.read_file_to_dict(f'{toolPath}/Pfam/Pfam-hmms/Pfam-A.hmm.length')
    if 'smart' in tool_names:
        smart_len = annoModules.read_file_to_dict(f'{toolPath}/SMART/SMART-hmms/SMART.hmm.length')
    with open(path, 'r') as f:
        for _ in range(ignore_lines):
            f.readline()
        for line in f:
            cells = line.rstrip('\n').split('\t')
            p_id = cells[columns['protein_id']]
            length = int(cells[columns['length']])
            f_name = cells[columns['feature_id']]
            start = int(cells[columns['start']])
            stop = int(cells[columns['stop']])
            evalue = cells[columns.get('evalue', -1)] if 'evalue' in columns else "NA"
            if evalue == "-":
                evalue = "NA"
            score = float(cells[columns.get('score', -1)]) if 'score' in columns and cells[columns['score']] != '-' else None
            tool_name = cells[tool_column].lower()

            if tool_name == 'smart':
                if f_name in smart_len:
                    domain_len_dict[f"smart_{f_name}"] = pfam_len[f_name]

            if tool_name == 'pfam':
                if f_name in pfamDict:
                    clan_id = 'NA'
                    dom_len = 'NA'
                    if pfamDict[f_name]['clan'] != 'NA':
                        clan_id = pfamDict[f_name]['clan']
                    f_name = pfamDict[f_name]['id']
                    clan_dict[f"pfam_{f_name}"] = clan_id
                    if f_name in pfam_len:
                        domain_len_dict[f"pfam_{f_name}"] = pfam_len[f_name]

            f_id = f"{tool_name}_{f_name}"

            if tool_name not in tool_names:
                raise Exception(f"Found tool '{tool_name}' not in {tool_names}")

            if p_id not in features:
                features[p_id] = {
                    'length': length,
                    'pfam': {},
                    'smart': {},
                    'tmhmm': {},
                    'signalp': {},
                    'coils2': {}
                }

            if f_id not in features[p_id][tool_name]:
                features[p_id][tool_name][f_id] = {'evalue': 'NA', 'instance': []}

            # Construct instance: [start, stop, evalue, score, start_index, stop_index]
            instance = [start, stop, evalue]
            if score is not None:
                instance.append(score)
                instance.extend([start, stop])  # replicates start_index, stop_index in your JSON example
            features[p_id][tool_name][f_id]['instance'].append(instance)

            # InterPro ID
            if tool_name in ['pfam', 'smart']:
                interpro_id = cells[columns.get('interpro', 5)] if 'interpro' in columns or len(cells) > 5 else '-'
                if interpro_id != '-':
                    interpro[f_id] = interpro_id

    if seq_file and os.path.exists(seq_file):
        all_seqs = get_seq_len(seq_file)
        for seq_id in all_seqs:
            if not seq_id in features:
                features[seq_id] = {
                    'length': all_seqs[seq_id],
                    'pfam': {},
                    'smart': {},
                    'tmhmm': {},
                    'signalp': {},
                    'coils2': {}
                }

    return features, interpro, clan_dict, domain_len_dict

def compute_count(features):
    count = {}
    for p_id, tools in features.items():
        for tool_name, feats in tools.items():
            if tool_name != 'length':
                for f_id, f_data in feats.items():
                    count[f_id] = count.get(f_id, 0) + len(f_data['instance'])
    return count

def write_json(data, path):
    with open(path, 'w') as f:
        json.dump(data, f, ensure_ascii=False, separators=(',', ':'))

def main():
    fas_version = version("greedyFAS")
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(fas_version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-i", "--input", type=str, default=None, required=True,
                          help="Path to input file, needs to be in tsv format")
    optional.add_argument("-o", "--output", type=str, default=None,
                          help="Path to output file, should be a .json file")
    required.add_argument("-f", "--feature_columns", type=int, nargs=6, required=True,
                          help="6 integer values that point to the columns that contain: "
                               "(1) the protein id, (2) the protein length, (3) the feature/interpro id, "
                               "(4) the start position of the feature, (5) the end position of the feature, "
                               "(6) the e-value of the feature. NOTE: the column "
                               "indices start at 0 so the first column has index 0, the second index 1, etc. "
                               "Example: -f 0 2 4 6 7 8")
    required.add_argument("-c", "--tool_column", type=int, default=None,
                          help="The index (integer) of column that contains the annotation tool names")
    optional.add_argument("--seed_fasta", help="Path to fasta file of the seed sequences (RECOMMENDED!)", default='')
    optional.add_argument("--ignore_lines", type=int, help="skip the first n lines (for example headers)", default=0)
    optional.add_argument('--toolPath', help='Path to annotation tools', action='store', default='')
    optional.add_argument('--annoToolFile', help='Path to file containing list of annotation tools', action='store', default='')
    args = parser.parse_args()


    toolPath = args.toolPath
    if toolPath == '':
        pathconfigFile = os.path.realpath(__file__).replace('annoFAS/annoParserInterPro.py', 'pathconfig.txt')
        if not os.path.exists(pathconfigFile):
            sys.exit('No pathconfig.txt found. Please run fas.setup (https://github.com/BIONF/FAS/wiki/setup).')
        with open(pathconfigFile) as f:
            toolPath = f.readline().strip()
    else:
        toolPath = os.path.abspath(args.toolPath)

    # Map column names for extended fields
    columns = {
        'protein_id': args.feature_columns[0],
        'length': args.feature_columns[1],
        'feature_id': args.feature_columns[2],
        'start': args.feature_columns[3],
        'stop': args.feature_columns[4],
        'evalue': args.feature_columns[5],
        'interpro': args.feature_columns[2]
    }

    annoToolFile = args.annoToolFile
    toolList = annoModules.getAnnoTools(annoToolFile, toolPath)

    input = os.path.abspath(args.input) if os.path.exists(args.input) else sys.exit(f'{args.input} not found')
    features, interpro, clan, length = read_tsv_to_features(args.ignore_lines, columns, [t.lower() for t in toolList],
                                              args.tool_column, input, args.seed_fasta, toolPath)
    count = compute_count(features)

    data_out = {
        'feature': features,
        'interproID': interpro,
        'clan': clan,
        'count': count,
        'length': length,
        'version': {
            'pfam': {"version": "interpro", "evalue": ['NA', 'NA']},
            'smart': {"version": "interpro", "evalue": ['NA', 'NA']},
            'coils2': {"version": "interpro"},
            'signalp': {"version": "4.1", "org": "euk"},
            'tmhmm': {"version": "TMHMM2.0c"}
        }
    }

    output = args.output if args.output else f'{input}.json'
    write_json(data_out, output)
    print(f'Finished! Output {output}')

if __name__ == "__main__":
    main()
