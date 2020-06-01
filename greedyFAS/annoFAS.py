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
import greedyFAS.annoModules as annoModules

home = expanduser('~')

def runAnnoFas(args):
    (seqFile, outPath, toolPath, force, outName, eFlps, signalpOrg, eFeature, eInstance, hmmCores, redo, extract, oldName, cpus) = args
    # do annotation
    outFile = outPath+'/'+outName+'.json'
    if annoModules.checkFileEmpty(outFile) == True or force:
        if extract == '':
            print('Doing annotation...')
            annoJobs = annoModules.createAnnoJobs([outName, seqFile, toolPath, annoModules.getAnnoTools(toolPath), eFlps, signalpOrg, eFeature, eInstance, hmmCores])
            # do annotation and save to json output
            pool = mp.Pool(cpus)
            annoOut = pool.map(annoModules.doAnno, annoJobs)
            annoDict = annoModules.mergeNestedDic(annoOut)
            annoDict['clan'] = annoModules.getClans(toolPath, annoDict['feature'])
            annoDict['count'] = annoModules.countFeatures(annoDict['feature'])
            annoModules.save2json(annoDict, outName, outPath)
            # remove tmp sequence files
            annoModules.removeTmpFasta(outName)
            pool.close()
        else:
            if oldName == '':
                print('No reference annotaion given! Use --name <name of annotation file>')
            else:
                print('Extracting annotations...')
                annoDict = annoModules.extractAnno(seqFile, outPath+'/'+oldName+'.json')
                annoDict['clan'] = annoModules.getClans(toolPath, annoDict['feature'])
                annoDict['count'] = annoModules.countFeatures(annoDict['feature'])
                annoModules.save2json(annoDict, outName, outPath)
    else:
        if not redo == '':
            print('Redoing annotation for %s...' % redo)
            redoJobs = annoModules.createAnnoJobs([outName, seqFile, toolPath, [redo], eFlps, signalpOrg, eFeature, eInstance, hmmCores])
            # redo annotation
            pool = mp.Pool(mp.cpu_count()-1)
            annoOut = pool.map(annoModules.doAnno, redoJobs)
            redoAnnoDict = annoModules.mergeNestedDic(annoOut)
            # replace old annotations and save to json output
            annoDict = annoModules.replaceAnno(outFile, redoAnnoDict, redo)
            annoDict['clan'] = annoModules.getClans(toolPath, annoDict['feature'])
            annoDict['count'] = annoModules.countFeatures(annoDict['feature'])
            annoModules.save2json(annoDict, outName, outPath)
            # remove tmp sequence files
            annoModules.removeTmpFasta(outName)
            pool.close()
        else:
            print(outFile + ' already exists!')


def main():
    version = '1.1.0'
    parser = argparse.ArgumentParser(description='You are running annoFAS version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--fasta', help='Input sequence(s) in fasta format', action='store', default='',
                          required=True)
    required.add_argument('-o', '--outPath', help='Output directory', action='store', default='', required=True)
    required.add_argument('-t', '--toolPath', help='Path to annotation tools', action='store', default='',
                          required=True)
    optional.add_argument('--force', help='Force override annotations', action='store_true')
    optional.add_argument('-n', '--name', help='Name of annotation file', action='store', default='')
    optional.add_argument('--eFeature', help='eValue cutoff for PFAM/SMART domain', action='store', default=0.001, type=float)
    optional.add_argument('--eInstance', help='eValue cutoff for PFAM/SMART instance', action='store', default=0.01, type=float)
    optional.add_argument('--cpus', help='Number of CPUs used for annotation', action='store', default=0, type=int)
    optional.add_argument('--hmmCores', help='Number of CPUs used for hmm search', action='store', default=1, type=int)
    optional.add_argument('--eFlps', help='eValue cutoff for fLPS', action='store', default=0.0000001, type=float)
    optional.add_argument('--org', help='Organism of input sequence(s) for SignalP search', choices=['euk', 'gram+', 'gram-'], action='store', default='euk', type=str)
    optional.add_argument('--redo', help='Re-annotation the sequence with flps|coils2|seg|pfam|signalp|smart|tmhmm. '
                                         'Only one selection allowed!', choices=['flps', 'tmhmm', 'signalp', 'coils2', 'seg', 'smart', 'pfam'],
                                         action='store', default='', type=str)
    optional.add_argument('-e', '--extract', help='Path to save the extracted annotation for input sequence',
                          action='store_true', default='')

    args = parser.parse_args()

    # options for doing annotation
    seqFile = os.path.abspath(args.fasta)
    annoModules.checkFileExist(seqFile)
    toolPath = os.path.abspath(args.toolPath)
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
    try:
        my_abs_path = Path(outPath).resolve(strict=True)
    except FileNotFoundError:
        Path(outPath).mkdir(parents = True, exist_ok = True)
    if args.extract == '':
        oldName = ''
        outName = args.name
        if len(outName) == 0:
            outName = seqFile.split('/')[-1].split('.')[0]
    else:
        oldName = args.name
        outName = seqFile.split('/')[-1].split('.')[0]

    # other options
    force = args.force
    redo = args.redo
    extract = args.extract

    # run annoFAS
    start = time.time()
    print('PID ' + str(os.getpid()))
    runAnnoFas([seqFile, outPath, toolPath, force, outName, eFlps, signalpOrg, eFeature, eInstance, hmmCores, redo, extract, oldName, cpus])
    ende = time.time()
    print('Finished in ' + '{:5.3f}s'.format(ende-start))


if __name__ == '__main__':
    main()
