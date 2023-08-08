#!/bin/env python

#######################################################################
# Copyright (C) 2020 Vinh Tran
#
#  This file is part of FAS.
#  This script is used for preparing the annotation tools for annoFAS
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
from Bio import SeqIO
from pathlib import Path
import subprocess
import re
import collections
import json
import shutil

# general functions
def mergeNestedDic(dictList):
    out = collections.defaultdict(list)
    out.update(dictList.pop(0))
    for dd in dictList:
        for key, value in dd.items():
            out[key].update(value)
    return out


def save2json(dict2save, outName, outDir):
    Path(outDir).mkdir(parents=True, exist_ok=True)
    jsonOut = json.dumps(dict2save, ensure_ascii=False)
    f = open(outDir+'/'+outName+'.json', 'w')
    f.write(jsonOut)
    f.close()


def checkFileExist(file):
    if not os.path.exists(os.path.abspath(file)):
        sys.exit('%s not found' % file)


def checkFileEmpty(file):
    flag = False
    try:
        if os.path.getsize(file) == 0:
            flag = True
    except OSError as e:
        flag = True
    return flag


def read_file_to_dict(file):
    """ Read a file with 2 columns {id}\t{val} and return a dictionary"""
    if os.path.exists(file):
        dict = {}
        with open(file, 'r') as f:
            for l in [line.rstrip() for line in f]:
                dict[l.split()[0]] = l.split()[1].strip()
        return(dict)
    else:
        sys.exit('%s not found' % file)


def get_property(prop, project):
    result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop), open(project + '/__init__.py').read())
    return result.group(1)


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


def printMsg(silent, msg):
    if silent == False:
        print(msg)


# functions for doing annotation with single tool
def doFlps(args):
    (seqFile, toolPath, threshold) = args
    # load fasta seq
    inSeq = SeqIO.to_dict((SeqIO.parse(open(seqFile), 'fasta')))
    # do fLPS
    cmd = '%s/fLPS/fLPS -s -t %s "%s"' % (toolPath, threshold, seqFile)
    try:
        flpsOut = subprocess.run([cmd], shell=True, capture_output=True, check=True)
    except:
        sys.exit('Error running\n%s' % cmd)
    lines = flpsOut.stdout.decode().split('\n')
    # save to dict
    annoOut = {}
    annotatedSeq = {}
    if len(lines) > 0:
        for line in lines:
            tmp = line.strip().split('\t')
            if tmp[0]:
                if not tmp[0] in annoOut:
                    annoOut[tmp[0]] = {}
                    annoOut[tmp[0]]['length'] = len(inSeq[tmp[0]])
                    annoOut[tmp[0]]['flps'] = {}
                flps_id = 'flps_' + tmp[1] + '_' + tmp[7]
                if len(tmp) > 7:
                    flps_id = 'flps_' + tmp[2] + '_' + tmp[8]
                if not flps_id in annoOut[tmp[0]]['flps']:
                    annoOut[tmp[0]]['flps'][flps_id] = {}
                    annoOut[tmp[0]]['flps'][flps_id]['evalue'] = float(threshold)
                    annoOut[tmp[0]]['flps'][flps_id]['instance'] = []
                if len(tmp) > 7:
                    annoOut[tmp[0]]['flps'][flps_id]['instance'].append((int(tmp[4]), int(tmp[5]), float(tmp[7])))
                else:
                    annoOut[tmp[0]]['flps'][flps_id]['instance'].append((int(tmp[3]), int(tmp[4]), float(tmp[6])))
                annotatedSeq[tmp[0]] = 1

    for id in inSeq:
        if not id in annotatedSeq:
            annoOut[id] = {}
            annoOut[id]['length'] = len(inSeq[id])
            annoOut[id]['flps'] = {}
    return annoOut


def doTmhmm(args):
    (seqFile, toolPath) = args
    # load fasta seq
    inSeq = SeqIO.to_dict((SeqIO.parse(open(seqFile),'fasta')))
    # do TMHMM
    cmd = 'cat "%s" | %s/TMHMM/decodeanhmm -f %s/TMHMM//lib/TMHMM2.0.options -modelfile %s/TMHMM/lib/TMHMM2.0.model' % (
        seqFile, toolPath, toolPath, toolPath)
    try:
        tmhmmOut = subprocess.run([cmd], shell=True, capture_output=True) #check=True not works
    except:
        sys.exit('Error running\n%s' % cmd)
    lines = tmhmmOut.stdout.decode().split('\n')
    # save to dict
    annoOut = {}
    annotatedSeq = {}
    if len(lines) > 0:
        for line in lines:
            if line.startswith('>'):
                id = line.replace('>', '').strip()
                annoOut[id] = {}
                annoOut[id]['length'] = len(inSeq[id])
                annoOut[id]['tmhmm'] = {}
                annoOut[id]['tmhmm']['tmhmm_transmembrane'] = {}
                annoOut[id]['tmhmm']['tmhmm_transmembrane']['evalue'] = 'NA'
                annoOut[id]['tmhmm']['tmhmm_transmembrane']['instance'] = []
            if '%pred' in line:
                tmp = line.strip().split(',')
                for item in tmp:
                    if item.startswith(' M'):
                        pos = item.split()
                        annoOut[id]['tmhmm']['tmhmm_transmembrane']['instance'].append((int(pos[1]), int(pos[2]), 'NA'))
                annotatedSeq[id] = 1
    for id in inSeq:
        id = id
        if not id in annotatedSeq:
            annoOut[id] = {}
            annoOut[id]['length'] = len(inSeq[id])
            annoOut[id]['tmhmm'] = {}
        else:
            if len(annoOut[id]['tmhmm']['tmhmm_transmembrane']['instance']) == 0:
                annoOut[id]['tmhmm'].pop('tmhmm_transmembrane', None)
    return annoOut


def doSignalp(args):
    (seqFile, toolPath, org, outPath) = args
    # create tmp folder for signalp output
    seqFileName = seqFile.split('/')[-1].split('.')[0]
    Path(outPath+'/tmp/signalp/'+seqFileName).mkdir(parents = True, exist_ok = True)
    # load fasta seq
    inSeq = SeqIO.to_dict((SeqIO.parse(open(seqFile), 'fasta')))
    # do signalp
    cmd = '%s/SignalP/signalp -T %s/tmp/signalp/%s -t %s "%s"' % (toolPath, outPath, seqFileName, org, seqFile)
    try:
        signalpOut = subprocess.run([cmd], shell=True, capture_output=True, check=True)
    except:
        sys.exit('Error running\n%s' % cmd)
    lines = signalpOut.stdout.decode().split('\n')
    # print(signalpOut.stdout.decode())
    # save to dict
    annoOut = {}
    if len(lines) > 0:
        for line in lines:
            if not line.startswith('#') and not 'Temporary files in' in line:
                    if len(line) > 0:
                        tmp = line.strip().split()
                        if tmp[9] == 'Y':
                            annoOut[tmp[0]] = {}
                            annoOut[tmp[0]]['length'] = len(inSeq[tmp[0]])
                            annoOut[tmp[0]]['signalp'] = {}
                            annoOut[tmp[0]]['signalp']['signalp_SIGNALP'] = {}
                            annoOut[tmp[0]]['signalp']['signalp_SIGNALP']['evalue'] = 'NA'
                            annoOut[tmp[0]]['signalp']['signalp_SIGNALP']['instance'] = [[1, int(tmp[4])-1, 'NA']]
                        else:
                            annoOut[tmp[0]] = {}
                            annoOut[tmp[0]]['length'] = len(inSeq[tmp[0]])
                            annoOut[tmp[0]]['signalp'] = {}
    return annoOut


def doCoils(args):
    (seqFile, toolPath) = args
    currentDir = os.path.abspath(os.getcwd())
    # load fasta seq
    inSeq = SeqIO.to_dict((SeqIO.parse(open(seqFile),'fasta')))
    # do COILS2
    cmd = '%s/COILS2/COILS2 -f < "%s"' % (toolPath, os.path.abspath(seqFile))
    os.chdir(toolPath + '/COILS2')
    try:
        coilsOut = subprocess.run([cmd], shell=True, capture_output=True, check=True)
    except:
        sys.exit('Error running\n%s' % cmd)
    results = coilsOut.stdout.decode().split('>')
    # save to dict
    annoOut = {}
    if len(results) > 0:
        for outSeq in results:
            if len(outSeq) > 0:
                tmp = outSeq.split('\n')
                id = tmp.pop(0).strip()
                concatSeq = ''.join(tmp).strip()
                annoOut[id] = {}
                annoOut[id]['length'] = len(inSeq[id])
                annoOut[id]['coils2'] = {}
                annoOut[id]['coils2']['coils_coiled_coil'] = {}
                annoOut[id]['coils2']['coils_coiled_coil']['evalue'] = 'NA'
                annoOut[id]['coils2']['coils_coiled_coil']['instance'] = []
                finished = False
                while not finished:
                    match = re.search(r'x+', concatSeq)
                    if match:
                        start = concatSeq.find(match.group()) + 1
                        end = start + len(match.group()) - 1
                        annoOut[id]['coils2']['coils_coiled_coil']['instance'].append((int(start), int(end), 'NA'))
                        concatSeq = concatSeq.replace(match.group(), 'N' * len(match.group()), 1)
                    else:
                        finished = True
    os.chdir(currentDir)
    for id in inSeq:
        if len(annoOut[id]['coils2']['coils_coiled_coil']['instance']) == 0:
            annoOut[id]['coils2'].pop('coils_coiled_coil', None)
    return annoOut


def doSeg(args):
    (seqFile, toolPath) = args
    # load fasta seq
    inSeq = SeqIO.to_dict((SeqIO.parse(open(seqFile),'fasta')))
    # do signalp
    cmd = '%s/SEG/seg "%s" -l -n -p ' % (toolPath, seqFile)
    try:
        signalpOut = subprocess.run([cmd], shell=True, capture_output=True, check=True)
    except:
        sys.exit('Error running\n%s' % cmd)
    results = signalpOut.stdout.decode().split('>')
    # save to dict
    annoOut = {}
    if len(results) > 0:
        for result in results:
            if len(result) > 0:
                lines = result.split('\n')
                id = lines.pop(0).strip().replace('>','')
                annoOut[id] = {}
                annoOut[id]['length'] = len(inSeq[id])
                annoOut[id]['seg'] = {}
                annoOut[id]['seg']['seg_low_complexity_regions'] = {}
                annoOut[id]['seg']['seg_low_complexity_regions']['evalue'] = 'NA'
                annoOut[id]['seg']['seg_low_complexity_regions']['instance'] = []
                for line in lines:
                    tmp = line.strip().split()
                    if (len(tmp) == 2) and ('-' in tmp[1]):
                        pos = tmp[1].split('-')
                        annoOut[id]['seg']['seg_low_complexity_regions']['instance'].append((int(pos[0]), int(pos[1]),
                                                                                             'NA'))
    for id in inSeq:
        if len(annoOut[id]['seg']['seg_low_complexity_regions']['instance']) == 0:
            annoOut[id]['seg'].pop('seg_low_complexity_regions', None)
    return(annoOut)


def doSmart(args):
    (seqFile, toolPath, cpus, eFeature, eInstance) = args
    currentDir = os.path.abspath(os.getcwd())
    checkDB(toolPath, 'SMART')
    smartOut = hmmScan(seqFile, toolPath, 'SMART', cpus)
    annoOut = parseHmmscan(smartOut, 'SMART', eFeature, eInstance)
    return(annoOut)


def doPfam(args):
    (seqFile, toolPath, cpus, eFeature, eInstance) = args
    currentDir = os.path.abspath(os.getcwd())
    checkDB(toolPath, 'Pfam')
    pfamOut = hmmScan(seqFile, toolPath, 'Pfam', cpus)
    annoOut = parseHmmscan(pfamOut, 'Pfam', eFeature, eInstance)
    return(annoOut)


def checkDB(toolPath, toolName):
    dbDir = '%s/%s/%s-hmms/' % (toolPath, toolName, toolName)
    ext = '.hmm'
    if toolName == 'Pfam':
        ext = '-A.hmm'
    checkFileExist(dbDir + toolName + ext)
    if not (os.path.exists(dbDir + toolName + ext + '.h3f') or os.path.exists(dbDir + toolName + ext + '.h3i') or
            os.path.exists(dbDir + toolName + ext + '.h3m') or os.path.exists(dbDir + toolName + ext + '.h3p')):
        cmd = 'hmmpress -f %s/%s%s' % (dbDir, toolName, ext)
        try:
            subprocess.run([cmd], shell=True, check=True)
        except:
            print('Problem occurred while creating binary files for %s/%s%s' % (dbDir, toolName, ext))


def hmmScan(seqFile, toolPath, toolName, cpus):
    ext = '.hmm'
    if toolName == 'Pfam':
        ext = '-A.hmm'
    scanCmd = 'hmmscan -E 0.01 --domE 0.1 --noali --cpu %s %s/%s/%s-hmms/%s%s "%s"' % (cpus, toolPath, toolName,
                                                                                     toolName, toolName, ext, seqFile)
    flag = False
    try:
        FNULL = open(os.devnull, 'w')
        hmmOut = subprocess.run([scanCmd], shell=True, capture_output=True, check=True)
        return(hmmOut.stdout.decode())
    except:
        print('Error running hmmscan! %s' % (scanCmd))
        return(False)


def parseHmmscan(hmmOut, toolName, eFeature, eInstance):
    outDict = {}
    toolName = toolName.lower()
    hmmResults = hmmOut.split('//')
    for result in hmmResults:
        if len(result.split('\n')) > 3:
            query = re.search(r'Query:(.)+', result).group().split()[1]
            seqLen = re.search(r'Query:(.)+', result).group().split()[2].replace('[L=', '').replace(']', '')
            outDict[query] ={}
            outDict[query]['length'] = int(seqLen)
            outDict[query][toolName] = {}
            if not 'No hits detected that satisfy reporting thresholds' in result:
                tmp = result.split('Domain annotation for each model:')
                flag_instance = {} # check if any instance exceeds the cutoff
                for line in tmp[0].split('\n'):
                    if re.search('^\d+', line.lstrip()):
                        items = line.lstrip().split()
                        if float(items[0]) <= eFeature:
                            outDict[query][toolName][toolName+'_'+items[8]] = {}
                            outDict[query][toolName][toolName+'_'+items[8]]['evalue'] = float(items[0])
                            outDict[query][toolName][toolName+'_'+items[8]]['instance'] = []
                            flag_instance[toolName+'_'+items[8]] = 0
                for line in tmp[1].split('\n'):
                    if line.startswith('>>'):
                        dom = line.strip().split()[1]
                    if re.search('^\d+', line.lstrip()):
                        items2 = line.lstrip().split()
                        if toolName + '_' + dom in outDict[query][toolName]:
                            if float(items2[4]) <= eInstance:
                                flag_instance[toolName+'_'+dom] = 1
                                outDict[query][toolName][toolName+'_'+dom]['instance'].append((int(items2[12]), # envfrom
                                                                                               int(items2[13]), # envto
                                                                                               float(items2[4]), # c-Evalue
                                                                                               float(items2[2]), # bit-score
                                                                                               int(items2[6]), # hmmfrom
                                                                                               int(items2[7]))) # hmmto
                for i in flag_instance:
                    if flag_instance[i] == 0:
                        del outDict[query][toolName][i]
    return outDict

# get clan and acc from PFAM dat file
def readDatFile(toolPath):
    datFile = '%s/Pfam/Pfam-hmms/Pfam-A.hmm.dat' % toolPath
    try:
        with open(datFile, 'r') as file:
            datFileR = file.read()
            file.close()
            blocks = datFileR.split('//')
            pfamDict = {}
            for bl in blocks:
                if '#=GF ID' in bl:
                    dom = re.search(r'#=GF ID(.)+', bl).group().split()[-1]
                    pfamDict[dom] = {}
                    if '#=GF CL' in bl:
                        clan = re.search(r'#=GF CL(.)+', bl).group().split()[-1]
                    else:
                        clan = 'NA'
                    acc = re.search(r'#=GF AC(.)+', bl).group().split()[-1]
                    len = re.search(r'#=GF ML(.)+', bl).group().split()[-1]
                    pfamDict[dom]['clan'] = clan
                    pfamDict[dom]['acc'] = acc
                    pfamDict[dom]['length'] = len
            return pfamDict
    except:
        print('%s does not exist or no accession ID can be found' % datFile)
        return

# functions for doing annotation for multiple tools
def getAnnoTools(annoToolFile, toolPath):
    if annoToolFile == "" or annoToolFile == "None":
        checkFileExist(toolPath+'/annoTools.txt')
        toolList = []
        with open(toolPath+'/annoTools.txt') as f:
            file =  f.readlines()
            if '#checked' not in ''.join(file):
                sys.exit('Annotation tools not ready. Please run fas.setup first!')
            else:
                for tool in file:
                    if ('#' not in tool) and (len(tool) > 1):
                        toolList.append(tool.strip().lower())
    else:
        checkFileExist(annoToolFile)
        toolList = []
        with open(annoToolFile) as f:
            file =  f.readlines()
            if '#checked' not in ''.join(file):
                sys.exit('Annotation tools not ready. Please run fas.setup first!')
            else:
                for tool in file:
                    if ('#' not in tool) and (len(tool) > 1):
                        toolList.append(tool.strip().lower())
    return toolList


def createAnnoJobs(args):
    (outName, outPath, seqFile, toolPath, toolList, eFlps, signalpOrg, eFeature, eInstance, hmmCores, pid) = args
    annoJobs = []
    for s in SeqIO.parse(seqFile, 'fasta'):
        annoJobs.append([s.id, s.seq, outName, outPath, toolPath, toolList, eFlps, signalpOrg, eFeature,
                         eInstance, hmmCores, pid])
    return annoJobs

def doAnno(args):
    (seqId, seq, outName, outPath, toolPath, toolList, eFlps, signalpOrg, eFeature, eInstance, hmmCores, pid) = args
    # create temp fasta file
    Path(f'{outPath}/tmp/{pid}').mkdir(parents = True, exist_ok = True)
    outNameTmp = outName.replace("|","_")
    seqIdTmp = seqId.replace("|","_").replace(".","_")
    seqFile = f'{outPath}/tmp/{pid}/{outNameTmp}_{seqIdTmp}.fa'
    tmpFile = open(seqFile, 'w')
    tmpFile.write(str('>' + seqId + '\n' + seq))
    tmpFile.close()
    # run annotation
    annoList = []
    final = {}
    if 'flps' in toolList:
        anno = doFlps([seqFile, toolPath, eFlps])
        annoList.append(anno)
    if 'tmhmm' in toolList:
        anno = doTmhmm([seqFile, toolPath])
        annoList.append(anno)
    if 'signalp' in toolList:
        anno = doSignalp([seqFile, toolPath, signalpOrg, outPath])
        annoList.append(anno)
    if 'coils2' in toolList:
        anno = doCoils([seqFile, toolPath])
        annoList.append(anno)
    if 'seg' in toolList:
        anno = doSeg([seqFile, toolPath])
        annoList.append(anno)
    if 'smart' in toolList:
        anno = doSmart([seqFile, toolPath, hmmCores, eFeature, eInstance])
        annoList.append(anno)
    if 'pfam' in toolList:
        anno = doPfam([seqFile, toolPath, hmmCores, eFeature, eInstance])
        annoList.append(anno)
    final['feature'] = mergeNestedDic(annoList)
    # remove tmp fasta file
    seqFileTmp = seqFile.replace("|","\|")
    rmCmd = 'rm %s' % (seqFileTmp)
    try:
        subprocess.run([rmCmd], shell=True, check=True)
    except:
        sys.exit('Error running\n%s' % rmCmd)
    if 'signalp' in toolList:
        try:
            seqFileName = seqFile.split('/')[-1].split('.')[0]
            shutil.rmtree('%s/tmp/signalp/%s' % (outPath, seqFileName))
        except:
            sys.exit('Error deleting %s/tmp/signalp/%s_%s' % (outPath, outNameTmp, seqIdTmp))
    return final

# functions for posprocessing annotation dictionary
def countFeatures(annoDict):
    out = []
    for prot in list(annoDict.keys()):
        for toolName in list(annoDict[prot].keys()):
            if not toolName == 'length':
                for feat, value in annoDict[prot][toolName].items():
                    out.extend([feat] * len(annoDict[prot][toolName][feat]['instance']))
    out.sort()
    count = collections.Counter(out)
    return dict(count)


def getClans(toolPath, annoDict):
    pfamDict = readDatFile(toolPath)
    outDict = {}
    for prot in list(annoDict.keys()):
        for dom, value in annoDict[prot]['pfam'].items():
            if dom.replace('pfam_','') in pfamDict:
                if not pfamDict[dom.replace('pfam_','')]['clan'] == 'NA':
                    outDict[dom] = pfamDict[dom.replace('pfam_','')]['clan']
    return outDict


def getPfamAcc(toolPath, annoDict):
    pfamDict = readDatFile(toolPath)
    outDict = {}
    for prot in list(annoDict.keys()):
        for dom, value in annoDict[prot]['pfam'].items():
            if dom.replace('pfam_','') in pfamDict:
                outDict[dom] = pfamDict[dom.replace('pfam_','')]['acc']
    return outDict


def getPhmmLength(toolPath, annoDict, toolList):
    outDict = {}
    if 'pfam' in toolList:
        pfam_len = read_file_to_dict(f'{toolPath}/Pfam/Pfam-hmms/Pfam-A.hmm.length')
    if 'smart' in toolList:
        smart_len = read_file_to_dict(f'{toolPath}/SMART/SMART-hmms/SMART.hmm.length')
    for prot in list(annoDict.keys()):
        if 'pfam' in annoDict[prot]:
            for dom, value in annoDict[prot]['pfam'].items():
                if dom.replace('pfam_','') in pfam_len:
                    outDict[dom] = pfam_len[dom.replace('pfam_','')]
        if 'smart' in annoDict[prot]:
            for dom, value in annoDict[prot]['smart'].items():
                if dom.replace('smart_','') in smart_len:
                    outDict[dom] = smart_len[dom.replace('smart_','')]
    return outDict


# get version and used cutoffs of annotation tools
def getVersions(tools, toolPath, cutoffs):
    eFeature = cutoffs[0]
    eInstance = cutoffs[1]
    eFlps = cutoffs[2]
    signalpOrg = cutoffs[3]
    versionDict = {}
    for tool in tools:
        versionDict[tool] = {}
    if 'pfam' in tools:
        pfamCmd = 'grep "  RELEASE " %s/Pfam/Pfam-hmms/relnotes.txt | rev | cut -d " " -f1 | rev' % toolPath
        try:
            pfamVersion = subprocess.run([pfamCmd], shell=True, capture_output=True, check=True)
            versionDict['pfam']['version'] = pfamVersion.stdout.decode().strip()
            versionDict['pfam']['evalue'] = (eFeature, eInstance)
        except subprocess.CalledProcessError as e:
            print(e.output.decode(sys.stdout.encoding))
            print('Error getting Pfam version! %s' % (pfamVersion))
    if 'smart' in tools:
        versionDict['smart']['version'] = 'NA'
        versionDict['smart']['evalue'] = (eFeature, eInstance)
        smartCmd = 'head %s/SMART/version.txt' % toolPath
        try:
            smartVersion = subprocess.run([smartCmd], shell=True, capture_output=True, check=True)
            if (len(smartVersion.stdout.decode().strip()) > 1):
                versionDict['smart']['version'] = smartVersion.stdout.decode().strip()
        except subprocess.CalledProcessError as e:
            print(e.output.decode(sys.stdout.encoding))
            print('SMART version cannot be identified! %s' % (smartVersion))
    if 'flps' in tools:
        versionDict['flps']['version'] = '1.0'
        versionDict['flps']['evalue'] = eFlps
        flpsCmd = 'grep README %s/fLPS/README.first | grep version | cut -d " " -f 4 | sed "s/)//"' % toolPath
        try:
            fLPSVersion = subprocess.run([flpsCmd], shell=True, capture_output=True, check=True)
            if (len(fLPSVersion.stdout.decode().strip()) > 1):
                versionDict['flps']['version'] = fLPSVersion.stdout.decode().strip()
        except subprocess.CalledProcessError as e:
            print(e.output.decode(sys.stdout.encoding))
            print('Error getting fLPS version! %s' % (fLPSVersion))
    if 'signalp' in tools:
        signalpCmd = 'head %s/SignalP/signalp-*.readme | grep "INSTALLATION INSTRUCTIONS" | cut -f1 | cut -d " " -f2' % toolPath
        try:
            signalpVersion = subprocess.run([signalpCmd], shell=True, capture_output=True, check=True)
            versionDict['signalp']['version'] = signalpVersion.stdout.decode().strip()
            versionDict['signalp']['org'] = (signalpOrg)
        except subprocess.CalledProcessError as e:
            print(e.output.decode(sys.stdout.encoding))
            print('Error getting SignalP version! %s' % (signalpVersion))
    if 'tmhmm' in tools:
        tmhmmCmd = 'head -n1 %s/TMHMM/README' % toolPath
        try:
            tmhmmVersion = subprocess.run([tmhmmCmd], shell=True, capture_output=True, check=True)
            versionDict['tmhmm']['version'] = tmhmmVersion.stdout.decode().strip()
        except subprocess.CalledProcessError as e:
            print(e.output.decode(sys.stdout.encoding))
            print('Error getting TMHMM version! %s' % (tmhmmVersion))
    return(versionDict)

# functions for replacing and extrating annotation from json file
def replaceAnno(oldAnnoFile, newAnnoDict, redo):
    with open(oldAnnoFile) as f:
        annoDict = json.load(f)
        f.close()
        for prot in list(annoDict['feature'].keys()):
            annoDict['feature'][prot][redo] = dict(newAnnoDict['feature'])[prot][redo]
        if 'clan' in newAnnoDict:
            annoDict['clan'] = newAnnoDict['clan']
        return annoDict


def extractAnno(seqFile, existingAnno):
    inSeq = SeqIO.to_dict((SeqIO.parse(open(seqFile),'fasta')))
    with open(existingAnno) as f:
        existingDict = json.load(f)
        f.close()
        annoDict = {}
        annoDict['feature'] = dict((prot, existingDict['feature'][prot]) for prot in list(inSeq.keys()))
        return annoDict

# update old json file to add interproIDs and tool versions
def updateAnnoFile(jsonFile):
    with open(jsonFile) as f:
        annoDict = json.load(f)
        toolPath = getFasToolPath()
        taxon = jsonFile.split('/')[-1].replace('.json', '')
        annoPath = jsonFile.replace('%s.json' % taxon, '')
        toolList = getAnnoTools("", toolPath)
        shutil.move('%s' % jsonFile, '%s.old' % jsonFile)
        if not "interproID" in annoDict:
            annoDict['interproID'] = getPfamAcc(toolPath, annoDict['feature'])
        if not 'version' in annoDict:
            cutoffs = (0.001, 0.01, 0.0000001, 'euk')
            annoDict['version'] = getVersions(getAnnoTools('', toolPath), toolPath, cutoffs)
        if not 'length' in annoDict:
            annoDict['length'] = getPhmmLength(toolPath, annoDict['feature'], toolList)
        save2json(annoDict, taxon, annoPath)
        return('%s.old' % jsonFile)
