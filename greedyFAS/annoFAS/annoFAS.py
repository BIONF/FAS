#!/bin/env python

#######################################################################
# Copyright (C) 2020 Vinh Tran
#
#  This file is part of FAS.
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

import os
import sys
import argparse
import time
from os.path import expanduser
from pathlib import Path
import multiprocessing as mp
import greedyFAS.annoFAS.annoModules as annoModules
from tqdm import tqdm
from shutil import copyfile

home = expanduser('~')

def runAnnoFas(args):
    (seqFile, outPath, toolPath, force, outName, eFlps, signalpOrg, eFeature, eInstance, hmmCores, redo, extract,
     annoFile, cpus) = args
    # do annotation
    outFile = outPath+'/'+outName+'.json'
    if annoModules.checkFileEmpty(outFile) == True or force:
        if extract == '':
            print('Doing annotation for %s...' % seqFile)
            annoJobs = annoModules.createAnnoJobs([outName, outPath, seqFile, toolPath, annoModules.getAnnoTools(toolPath),
                                                   eFlps, signalpOrg, eFeature, eInstance, hmmCores])
            # do annotation and save to json output
            pool = mp.Pool(cpus)
            annoOut = []
            for _ in tqdm(pool.imap_unordered(annoModules.doAnno, annoJobs), total=len(annoJobs)):
                annoOut.append(_)
            annoDict = annoModules.mergeNestedDic(annoOut)
            annoDict['clan'] = annoModules.getClans(toolPath, annoDict['feature'])
            annoDict['count'] = annoModules.countFeatures(annoDict['feature'])
            annoModules.save2json(annoDict, outName, outPath)
            pool.close()
        else:
            if annoFile == '':
                print('No reference annotaion given! Please specify with --annoFile <path to exising annotation file>')
            else:
                print('Extracting annotations...')
                annoDict = annoModules.extractAnno(seqFile, annoFile)
                annoDict['clan'] = annoModules.getClans(toolPath, annoDict['feature'])
                annoDict['count'] = annoModules.countFeatures(annoDict['feature'])
                annoModules.save2json(annoDict, outName, outPath)
    else:
        if not redo == '':
            print('Redoing annotation for %s...' % redo)
            redoJobs = annoModules.createAnnoJobs([outName, outPath, seqFile, toolPath, [redo], eFlps, signalpOrg, eFeature,
                                                   eInstance, hmmCores])
            # redo annotation
            pool = mp.Pool(mp.cpu_count()-1)
            annoOut = []
            for _ in tqdm(pool.imap_unordered(annoModules.doAnno, redoJobs), total=len(redoJobs)):
                annoOut.append(_)
            redoAnnoDict = annoModules.mergeNestedDic(annoOut)
            # replace old annotations and save to json output
            annoDict = annoModules.replaceAnno(outFile, redoAnnoDict, redo)
            annoDict['clan'] = annoModules.getClans(toolPath, annoDict['feature'])
            annoDict['count'] = annoModules.countFeatures(annoDict['feature'])
            annoModules.save2json(annoDict, outName, outPath)
            pool.close()
        else:
            print(outFile + ' already exists!')


def main():
    version = '1.3.3'
    parser = argparse.ArgumentParser(description='You are running annoFAS version ' + str(version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--fasta', help='Input sequence(s) in fasta format', action='store', default='',
                          required=True)
    required.add_argument('-o', '--outPath', help='Output directory', action='store', default='', required=True)
    optional.add_argument('--force', help='Force override annotations', action='store_true')
    optional.add_argument('-n', '--name', help='Name of annotation file', action='store', default='')
    optional.add_argument('-t', '--toolPath', help='Path to annotation tools', action='store', default='')
    optional.add_argument('--eFeature', help='eValue cutoff for PFAM/SMART domain. Default = 0.001', action='store',
                          default=0.001, type=float)
    optional.add_argument('--eInstance', help='eValue cutoff for PFAM/SMART instance. Default = 0.01', action='store',
                          default=0.01, type=float)
    optional.add_argument('--cpus', help='Number of CPUs used for annotation. Default = available cores - 1',
                          action='store', default=0, type=int)
    optional.add_argument('--hmmCores', help='Number of CPUs used for hmm search', action='store', default=1, type=int)
    optional.add_argument('--eFlps', help='eValue cutoff for fLPS. Default = 0.0000001', action='store',
                          default=0.0000001, type=float)
    optional.add_argument('--org', help='Organism of input sequence(s) for SignalP search. Default = "euk"',
                          choices=['euk', 'gram+', 'gram-'], action='store', default='euk', type=str)
    optional.add_argument('--redo', help='Re-annotation the sequence with flps|coils2|seg|pfam|signalp|smart|tmhmm. '
                                         'Only one selection allowed!',
                          choices=['flps', 'tmhmm', 'signalp', 'coils2', 'seg', 'smart', 'pfam'],
                          action='store', default='', type=str)
    optional.add_argument('-e', '--extract', help='Path to save the extracted annotation for input sequence(s). --annoFile required for specifying existing annotation file!',
                          action='store_true', default='')
    optional.add_argument('-a', '--annoFile', help='Path to existing annotation JSON file',
                          action='store', default='')

    args = parser.parse_args()

    # options for doing annotation
    toolPath = args.toolPath
    if toolPath == '':
        pathconfigFile = os.path.realpath(__file__).replace('annoFAS/annoFAS.py', 'pathconfig.txt')
        if not os.path.exists(pathconfigFile):
            sys.exit('No pathconfig.txt found. Please run setupFAS (https://github.com/BIONF/FAS/wiki/setupFAS).')
        with open(pathconfigFile) as f:
            toolPath = f.readline().strip()
    else:
        toolPath = os.path.abspath(args.toolPath)

    seqFile = args.fasta
    if seqFile == 'test_annofas.fa':
        seqFile = os.path.realpath(__file__).replace('annoFAS.py', 'test_annofas.fa')
    else:
        seqFile = os.path.abspath(seqFile)
    annoModules.checkFileExist(seqFile)

    eFlps = args.eFlps
    signalpOrg = args.org
    eFeature = args.eFeature
    eInstance = args.eInstance
    hmmCores = args.hmmCores
    cpus = args.cpus
    if cpus == 0:
        cpus = mp.cpu_count()-1

    # option for saving json file
    outPath = os.path.abspath(args.outPath)
    Path(outPath).mkdir(parents=True, exist_ok=True)
    outName = args.name
    if len(outName) == 0:
        outName = seqFile.split('/')[-1].split('.')[0]

    # other options
    force = args.force
    redo = args.redo
    extract = args.extract
    annoFile = args.annoFile
    if extract:
        annoModules.checkFileExist(annoFile)
        if annoFile == '':
            sys.exit('No existing annotation file given for extraction! Please specify it using --annoFile')

    # run annoFAS
    start = time.time()
    print('PID ' + str(os.getpid()))
    runAnnoFas([seqFile, outPath, toolPath, force, outName, eFlps, signalpOrg, eFeature, eInstance, hmmCores, redo,
                extract, annoFile, cpus])
    ende = time.time()
    print('Finished in ' + '{:5.3f}s'.format(ende-start))


if __name__ == '__main__':
    main()
