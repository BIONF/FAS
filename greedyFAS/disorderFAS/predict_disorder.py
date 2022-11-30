#!/bin/env python

#######################################################################
# Copyright (C) 2022 Julian Dosch
#
# This file is part of greedyFAS.
#
#  fas_disorder_addon is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  greedyFAS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with greedyFAS.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################


from Bio import SeqIO
import subprocess
import os
import sys
from pathlib import Path
import multiprocessing as mp
from tqdm import tqdm
import argparse
import json
from pkg_resources import get_distribution


def save2json(dict2save, outpath):
    jsonOut = json.dumps(dict2save, ensure_ascii=False)
    f = open(outpath, 'w')
    f.write(jsonOut)
    f.close()


def parse_aucpred(infile):
    struc = []
    with open(infile, 'r') as toparse:
        for line in toparse.readlines():
            if not line[0] == '#':
                if line.split()[2] == '*':
                    struc.append(1)
                elif line.split()[2] == '.':
                    struc.append(0)
                else:
                    raise Exception('Disordered annotation seems to be faulty!')
    dis = []
    switch = False
    start = None
    for i in range(len(struc)):
        if struc[i] and not switch:
            start = i
            switch = True
        elif not struc[i] and switch:
            dis.append([start, i, 'NA'])
            switch = False
    return {'aucpred': {'aucpred_disordered_region': {'evalue': 'NA', 'instance': dis}}}


def prepare_annojobs(infile, tmppath, aucpred):
    annojobs = []
    ids = []
    i = 0
    for seq in SeqIO.parse(infile, 'fasta'):
        if seq.id in ids:
            raise Exception('Headers in the fasta are not unique.')
        else:
            ids.append(seq.id)
            annojobs.append([seq.id, str(seq.seq), tmppath, aucpred, 'seq_' + str(i)])
        i += 1
    return annojobs


def remove_tmp(filepath):
    if os.path.exists(filepath):
        os.remove(filepath)


def make_tmp_fasta(header, seq, path):
    with open(path + '/' + header + '.fasta', 'w') as out:
        out.write('>sequence\n' + seq)


def run_anno(inpath, outpath, tmppath, cpus, aucpred):
    name = ''.join(inpath.split('/')[-1].split('.')[0:-1])
    Path(outpath).mkdir(parents=True, exist_ok=True)
    Path(tmppath + name + '/').mkdir(parents=True, exist_ok=True)
    joblist = prepare_annojobs(inpath, tmppath + name + '/', aucpred)
    out = []
    pool = mp.Pool(cpus)
    try:
        for _ in tqdm(pool.imap_unordered(run_anno_single, joblist), total=len(joblist), mininterval=5.0):
            out.append(_)
    except subprocess.CalledProcessError:
        raise 'Something went wrong during annotation.'
    pool.close()
    write_output(out, outpath + '/' + name)
    Path(tmppath + name + '/').rmdir()


def write_output(out, outpath):
    failed = []
    count = 0
    fdict = {}
    for seq in out:
        if not seq[1]:
            failed.append(seq[0])
        else:
            fdict[seq[0]] = seq[1]
            count += len(seq[1]['aucpred']['aucpred_disordered_region']['instance'])
    if not len(failed) == len(out):
        outdict = {'feature': fdict, 'clan': {}, 'count': {'aucpred_disordered_region': count}}
        save2json(outdict, outpath + '_disorder.json')
    if failed:
        print("Errors occured for the annotation of the following sequences:")
        for i in failed:
            print(i)


def run_anno_single(args):
    header, seq, tmppath, aucpred, tmpname = args
    make_tmp_fasta(tmpname, seq, tmppath)
    try:
        disorder = run_aucpred(tmpname, tmppath, aucpred)
        if os.path.exists(tmppath + tmpname + '.fasta'):
            os.remove(tmppath + tmpname + '.fasta')
        return header, disorder
    except:
        if os.path.exists(tmppath + tmpname + '.fasta'):
            os.remove(tmppath + tmpname + '.fasta')
        return header, None


def run_aucpred(header, tmppath, aucpred):
    cmd = aucpred + ' -i "' + tmppath + header + '.fasta" -o ' + tmppath
    subprocess.run([cmd], shell=True, capture_output=True, check=True)
    disorder = parse_aucpred(tmppath + header + '.diso_noprof')
    for i in ['.diso_noprof', '.diso_prev']:
        remove_tmp(tmppath + header + i)
    return disorder


def main():
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-i", "--input", default='.', type=str, required=True,
                          help="path to input fasta")
    required.add_argument("-o", "--outPath", default='.', type=str, required=True,
                          help="path to output directory. File will be named after fasta file")
    optional.add_argument("-t", "--tmp", default='.', type=str, required=False,
                          help="Path to a temporary directory.")
    optional.add_argument("--cpus", default=1, type=int, required=False,
                          help="number of cpus used")
    optional.add_argument("-a", "--aucpred", default=None, type=str, required=False,
                          help="Alternative path to AUCpreD.sh.")
    args = parser.parse_args()
    if args.aucpred:
        toolpath = args.aucpred
    else:
        pathconfigFile = os.path.realpath(__file__).replace('disorderFAS/predict_disorder.py', 'pathconfig.txt')
        if not os.path.exists(pathconfigFile):
            sys.exit('No pathconfig.txt found. Please run fas.setup (https://github.com/BIONF/FAS/wiki/setup).')
        with open(pathconfigFile, 'r') as input:
            toolpath = input.readline().strip() + '/Predict_Property/AUCpreD.sh'
    if not os.path.exists(toolpath):
        sys.exit('AUCpreD.sh not found either give a path to the AUCpreD.sh with [-a|--aucpred] or install via '
                 '"fas.disorder.setup".')
    run_anno(args.input, args.outPath, args.tmp, args.cpus, toolpath)


if __name__ == '__main__':
    main()
