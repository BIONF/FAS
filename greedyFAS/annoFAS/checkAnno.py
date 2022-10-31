# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2021 Vinh Tran
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

import sys
import os
import argparse
from pathlib import Path
import subprocess
from Bio import SeqIO
import shutil
import multiprocessing as mp
import json
import greedyFAS.annoFAS.annoModules as annoModules

def getFasToolPath():
    cmd = 'fas.setup -t ~/ -c'
    try:
        flpsOut = subprocess.run([cmd], shell=True, capture_output=True, check=True)
    except:
        sys.exit('Error running\n%s' % cmd)
    lines = flpsOut.stdout.decode().split('\n')
    if 'FAS is ready to run' in lines[0]:
        return(lines[0].replace('Annotation tools can be found at ','').replace('. FAS is ready to run!',''))
    else:
        sys.exit('FAS has not been setup!')

def checkCompleteAnno(seqFile, jsonFile):
    allSeq = []
    allSeqDict = SeqIO.to_dict((SeqIO.parse(open(seqFile), 'fasta')))
    allSeq = list(allSeqDict.keys())
    with open(jsonFile) as jf:
        dt = json.load(jf)
        annotated = dt['feature']
        missingSeqs = list(set(allSeq).difference(annotated))
        if len(missingSeqs) > 0:
            out = []
            for seq in missingSeqs:
                out.append('>%s\n%s\n' % (seq, allSeqDict[seq].seq))
            return(out)
        else:
            return()

def doAnnoForMissing(taxon, missingAnno, jsonFile, outPath, cpus, silent, annoToolFile):
    toolPath = getFasToolPath()
    # do annotation for missing proteins
    faFile = '%s/%s_tmp.fa' % (outPath, taxon)
    with open(faFile, 'w') as f:
        f.write(''.join(missingAnno))

    annoCmd = 'fas.doAnno -i %s -o %s --cpus %s --annoToolFile %s' % (faFile, outPath, cpus, annoToolFile)
    if silent:
        annoCmd = annoCmd + ' > /dev/null 2>&1'
    try:
        subprocess.call([annoCmd], shell = True)
    except:
        print('\033[91mProblem with running fas.doAnno. You can check it with this command:\n%s\033[0m' % annoCmd)
    # merge with old annotation json file
    with open(jsonFile) as oldJson:
        oldAnno = json.load(oldJson)
    with open('%s/%s_tmp.json' % (outPath, taxon)) as newJson:
        newAnno = json.load(newJson)
    annoOut = [oldAnno, newAnno]
    annoDict = annoModules.mergeNestedDic(annoOut)
    annoDict['clan'] = annoModules.getClans(toolPath, annoDict['feature'])
    annoDict['count'] = annoModules.countFeatures(annoDict['feature'])
    annoModules.save2json(annoDict, taxon, outPath)
    # remove tmp files
    os.remove(faFile)
    os.remove('%s/%s_tmp.json' % (outPath, taxon))
    shutil.rmtree('%s/tmp' % outPath)

def main():
    version = '1.14.2'
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-s', '--seqFile', help='Input sequence file in fasta format', action='store', default='',
                          required=True)
    required.add_argument('-a', '--annoFile', help='Input annotation file in json format', action='store', default='',
                          required=True)
    required.add_argument('-o', '--outPath', help='Output directory', action='store', default='', required=True)
    optional.add_argument('-n', '--noAnno', help='Do not annotate missing proteins', action='store_true')
    optional.add_argument('--silent', help='Turn off terminal output', action='store_true')
    optional.add_argument('--cpus', help='Number of CPUs used for annotation. Default = available cores - 1',
                          action='store', default=0, type=int)
    optional.add_argument('--annoToolFile', help='Path to files contains annotation tool names',
                          action='store', default='')

    args = parser.parse_args()

    seqFile = args.seqFile
    annoModules.checkFileExist(seqFile)
    seqFile = os.path.abspath(seqFile)
    annoFile = args.annoFile
    annoModules.checkFileExist(annoFile)
    annoFile = os.path.abspath(annoFile)
    outPath = os.path.abspath(args.outPath)
    Path(outPath).mkdir(parents=True, exist_ok=True)
    cpus = args.cpus
    if cpus == 0:
        cpus = mp.cpu_count()-1
    noAnno = args.noAnno
    silent = args.silent
    annoToolFile = args.annoToolFile
    annoModules.checkFileExist(annoToolFile)

    taxon = annoFile.split('/')[-1].replace('.json', '')
    missingAnno = checkCompleteAnno(seqFile, annoFile)
    if len(missingAnno) > 0:
        if noAnno == False:
            doAnnoForMissing(taxon, missingAnno, annoFile, outPath, cpus, silent, annoToolFile)
            if not silent:
                print('%s missing proteins of %s has been annotated!' % (len(missingAnno), taxon))
        else:
            if not silent:
                print('Annotation for %s proteins of %s are missing!\n%s' % (len(missingAnno), taxon, missingAnno))
    else:
        if not silent:
            print('Annotations found for all %s sequences!' % (taxon))

if __name__ == '__main__':
    main()
