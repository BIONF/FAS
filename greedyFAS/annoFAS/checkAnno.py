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
from pkg_resources import get_distribution


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


def checkOutdatedAnno(jsonFile):
    with open(jsonFile) as jf:
        annoDict = json.load(jf)
        if not 'interproID' in annoDict:
            return('interproID')
        elif not 'version' in annoDict:
            return('version')
        elif 'version' in annoDict:
            tools = list(annoDict['version'].keys())
            toolPath = annoModules.getFasToolPath()
            cutoffs = (0.001, 0.01, 0.0000001, 'euk')
            currentVer = annoModules.getVersions(tools, toolPath, cutoffs)
            outdatedTools = []
            for tool in tools:
                if 'version' in currentVer[tool]:
                    if 'version' in annoDict['version'][tool]:
                        if not annoDict['version'][tool]['version'] == currentVer[tool]['version']:
                            outdatedTools.append('%s (%s vs %s)' % (tool, annoDict['version'][tool]['version'], currentVer[tool]['version']))
                    else:
                        return('version')
            if len(outdatedTools) > 0:
                return(', '.join(outdatedTools))
            else:
                return('updated')
        else:
            return('updated')


def doAnnoForMissing(taxon, missingAnno, jsonFile, outPath, cpus, silent, annoToolFile):
    toolPath = annoModules.getFasToolPath()
    # do annotation for missing proteins
    faFile = '%s/%s_tmp.fa' % (outPath, taxon)
    with open(faFile, 'w') as f:
        f.write(''.join(missingAnno))

    annoCmd = 'fas.doAnno -i %s -o %s --cpus %s' % (faFile, outPath, cpus)
    if annoToolFile:
        annoCmd += ' --annoToolFile %s' % annoToolFile
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
    version = get_distribution('greedyFAS').version
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
    optional.add_argument('-u', '--update', help='Update data structure of annotation file (not the annotations themselves)', action='store_true')
    optional.add_argument('--keep', help='Keep old annotaion file', action='store_true')
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
    update = args.update
    keep = args.keep
    silent = args.silent
    annoToolFile = args.annoToolFile
    annoModules.checkFileExist(annoToolFile)

    taxon = annoFile.split('/')[-1].replace('.json', '')

    ### check for missing annotation
    missingAnno = checkCompleteAnno(seqFile, annoFile)
    if len(missingAnno) > 0:
        if noAnno == False:
            doAnnoForMissing(taxon, missingAnno, annoFile, outPath, cpus, silent, annoToolFile)
            annoModules.printMsg(silent, '%s missing proteins of %s has been annotated!' % (len(missingAnno), taxon))
        else:
            print('WARNING: Annotation for %s proteins of %s are missing!\n%s' % (len(missingAnno), taxon, missingAnno))
    else:
        annoModules.printMsg(silent, 'Annotations found for all %s sequences!' % (taxon))


    ### check for outdated annotation file
    outdatedCheck = checkOutdatedAnno(annoFile)
    if outdatedCheck == 'updated':
        annoModules.printMsg(silent, 'Annotation file updated!')
    elif outdatedCheck == 'interproID' or outdatedCheck == 'version':
        if update == True:
            oldAnnoFile = annoModules.updateAnnoFile(annoFile)
            annoModules.printMsg(silent, '%s updated!' % annoFile)
            if keep == False:
                os.remove(oldAnnoFile)
        else:
            print('WARNING: %s is not updated! Please consider update it with --update option' % annoFile)
    else:
        print('WARNING: Annotations from the tools below are outdated. Please consider to re-annotate them!\n%s' % outdatedCheck)



if __name__ == '__main__':
    main()
