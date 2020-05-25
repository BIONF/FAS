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
from collections import defaultdict
import json

def mergeNestedDic(dictList):
    out = defaultdict(list)
    out.update(dictList.pop(0))
    for dd in dictList:
        for key, value in dd.items():
            out[key].update(value)
    return out

def save2json(outDict, toolName, outDir):
    Path(outDir).mkdir(parents = True, exist_ok = True)
    jsonOut = json.dumps(outDict, ensure_ascii = False)
    f = open(outDir+"/"+toolName+'.json', 'w')
    f.write(jsonOut)
    f.close()


def doFlps(args):
    (seqFile, toolPath, threshold) = args
    # load fasta seq
    inSeq = SeqIO.to_dict((SeqIO.parse(open(seqFile),'fasta')))
    # do fLPS
    cmd = '%s/fLPS/fLPS -s -t %s %s' % (toolPath, threshold, seqFile)
    flpsOut = subprocess.run([cmd], shell=True, capture_output=True)
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
                    annoOut[tmp[0]]["length"] = len(inSeq[tmp[0]])
                    annoOut[tmp[0]]["flps"] = {}
                if not "flps_" + tmp[1] + "_" + tmp[7] in annoOut[tmp[0]]["flps"]:
                    annoOut[tmp[0]]["flps"]["flps_" + tmp[1] + "_" + tmp[7]] = []
                annoOut[tmp[0]]["flps"]["flps_" + tmp[1] + "_" + tmp[7]].append((tmp[3],tmp[4]))
                annotatedSeq[tmp[0]] = 1

    for id in inSeq:
        if not id in annotatedSeq:
            annoOut[id] = {}
            annoOut[id]["length"] = len(inSeq[id])
            annoOut[id]["flps"] = {}
    return(annoOut)


def doTmhmm(args):
    (seqFile, toolPath) = args
    # load fasta seq
    inSeq = SeqIO.to_dict((SeqIO.parse(open(seqFile),'fasta')))
    # do TMHMM
    cmd = 'cat %s | %s/TMHMM/decodeanhmm -f %s/TMHMM//lib/TMHMM2.0.options -modelfile %s/TMHMM/lib/TMHMM2.0.model' % (seqFile, toolPath, toolPath, toolPath)
    tmhmmOut = subprocess.run([cmd], shell=True, capture_output=True)
    lines = tmhmmOut.stdout.decode().split('\n')
    # save to dict
    annoOut = {}
    annotatedSeq = {}
    if len(lines) > 0:
        for line in lines:
            if line.startswith(">"):
                id = line.replace(">","").strip()
                annoOut[id] = {}
                annoOut[id]["length"] = len(inSeq[id])
                annoOut[id]["tmhmm"] = {}
                annoOut[id]["tmhmm"]["tmhmm_transmembrane"] = []
            if "%pred" in line:
                tmp = line.strip().split(',')
                for item in tmp:
                    if item.startswith(" M"):
                        pos = item.split()
                        annoOut[id]["tmhmm"]["tmhmm_transmembrane"].append((pos[1], pos[2]))
                annotatedSeq[id] = 1
    for id in inSeq:
        if not id in annotatedSeq:
            annoOut[id] = {}
            annoOut[id]["length"] = len(inSeq[id])
            annoOut[id]["tmhmm"] = {}
    return(annoOut)


def doSignalp(args):
    (seqFile, toolPath, org) = args
    # load fasta seq
    inSeq = SeqIO.to_dict((SeqIO.parse(open(seqFile),'fasta')))
    # do signalp
    cmd = '%s/SignalP/signalp -t %s %s' % (toolPath, org, seqFile)
    signalpOut = subprocess.run([cmd], shell=True, capture_output=True)
    lines = signalpOut.stdout.decode().split('\n')
    # save to dict
    annoOut = {}
    if len(lines) > 0:
        for line in lines:
            if not line.startswith("#"):
                if len(line) > 0:
                    tmp = line.strip().split()
                    # annoOut[tmp[0] + ";" + str(len(inSeq[tmp[0]]))] = {}
                    if tmp[9] == "Y":
                        annoOut[tmp[0]] = {}
                        annoOut[tmp[0]]["length"] = len(inSeq[tmp[0]])
                        annoOut[tmp[0]]["signalp"] = {}
                        annoOut[tmp[0]]["signalp"]["signalp_SIGNALP"] = (1,int(tmp[4])-1)
                    else:
                        annoOut[tmp[0]] = {}
                        annoOut[tmp[0]]["length"] = len(inSeq[tmp[0]])
                        annoOut[tmp[0]]["signalp"] = {}
    return(annoOut)


def doCoils(args):
    (seqFile, toolPath) = args
    currentDir = os.path.abspath(os.getcwd())
    # load fasta seq
    inSeq = SeqIO.to_dict((SeqIO.parse(open(seqFile),'fasta')))
    # do COILS2
    cmd = '%s/COILS2/COILS2 -f < %s' % (toolPath, os.path.abspath(seqFile))
    os.chdir(toolPath + "/COILS2")
    coilsOut = subprocess.run([cmd], shell=True, capture_output=True)
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
                annoOut[id]["length"] = len(inSeq[id])
                annoOut[id]["coils"] = {}
                annoOut[id]["coils"]["coils_coiled_coil"] = []
                finished = False
                while not finished:
                    match = re.search(r'x+', concatSeq)
                    if match:
                        start = concatSeq.find(match.group()) + 1
                        end = start + len(match.group()) - 1
                        annoOut[id]["coils"]["coils_coiled_coil"].append((start, end))
                        concatSeq = concatSeq.replace(match.group(), "N" * len(match.group()), 1)
                    else:
                        finished = True
    os.chdir(currentDir)
    return(annoOut)


def doSeg(args):
    (seqFile, toolPath) = args
    # load fasta seq
    inSeq = SeqIO.to_dict((SeqIO.parse(open(seqFile),'fasta')))
    # do signalp
    cmd = '%s/SEG/seg %s -l -n -p ' % (toolPath, seqFile)
    signalpOut = subprocess.run([cmd], shell=True, capture_output=True)
    results = signalpOut.stdout.decode().split('>')
    # save to dict
    annoOut = {}
    if len(results) > 0:
        for result in results:
            if len(result) > 0:
                lines = result.split('\n')
                id = lines.pop(0).strip().replace(">","")
                annoOut[id] = {}
                annoOut[id]["length"] = len(inSeq[id])
                annoOut[id]["seg"] = {}
                annoOut[id]["seg"]["seg_low_complexity_regions"] = []
                for line in lines:
                    tmp = line.strip().split()
                    if (len(tmp) == 2) and ("-" in tmp[1]):
                        pos = tmp[1].split('-')
                        annoOut[id]["seg"]["seg_low_complexity_regions"].append(pos)
    return(annoOut)


def doSmart(args):
    (seqFile, toolPath, cpus) = args
    currentDir = os.path.abspath(os.getcwd())
    # load fasta seq
    inSeq = SeqIO.to_dict((SeqIO.parse(open(seqFile),'fasta')))
    # do smart
    cmd = 'perl %s/SMART/smart_scan_v5.pl %s "" %s' % (toolPath, os.path.abspath(seqFile), cpus)
    os.chdir(toolPath + "/SMART")
    smartOut = subprocess.run([cmd], shell=True, capture_output=True)
    lines = smartOut.stdout.decode().split('\n')
    os.chdir(currentDir)
    # save to dict
    if ".out" in lines[-1]:
        smartOutFile = toolPath + "/SMART/output_files/" + lines[-1]
        return(parseHmmscan(inSeq, smartOutFile,"smart"))
    else:
        return("No output for SMART search")


def doPfam(args):
    (seqFile, toolPath, cpus) = args
    currentDir = os.path.abspath(os.getcwd())
    # load fasta seq
    inSeq = SeqIO.to_dict((SeqIO.parse(open(seqFile),'fasta')))
    # do smart
    cmd = 'perl %s/Pfam/pfam_scan_v2.pl %s "" %s' % (toolPath, os.path.abspath(seqFile), cpus)
    os.chdir(toolPath + "/Pfam")
    pfamOut = subprocess.run([cmd], shell=True, capture_output=True)
    lines = pfamOut.stdout.decode().split('\n')
    os.chdir(currentDir)
    # save to dict
    if ".out" in lines[-1]:
        pfamOutFile = toolPath + "/Pfam/output_files/" + lines[-1]
        return(parseHmmscan(inSeq, pfamOutFile, "pfam"))
    else:
        return("No output for PFAM search")


def parseHmmscan(inSeq, hmmOut, toolName):
    annoOut = {}
    annotatedSeq = {}
    with open(hmmOut, 'r') as file:
        outFile = file.read()
        outSeqs = outFile.split('>')
        for outSeq in outSeqs:
            if (len(outSeq) > 3) and (not (outSeq.startswith(" queryID"))) and (not ("No hits detected" in outSeq)):
                lines = outSeq.strip().split('\n')
                id = lines.pop(0).strip()
                annoOut[id] = {}
                annoOut[id]["length"] = len(inSeq[id])
                annoOut[id][toolName] = {}
                for line in lines:
                    items = line.strip().split('|')
                    if "family:" in line:
                        type = items[0].replace("# family: ","").strip()
                        clan = items[2].replace(" clan: ","").strip()
                        type_evalue = items[-1].replace(" E-value: ","").strip()
                        annoOut[id][toolName][toolName+"_"+type] = {}
                        annoOut[id][toolName][toolName+"_"+type]["clan"] = clan
                        annoOut[id][toolName][toolName+"_"+type]["evalue"] = type_evalue
                        annoOut[id][toolName][toolName+"_"+type]["instance"] = []
                        annotatedSeq[id] = 1
                    else:
                        start = items[3].replace("# family: ","").strip()
                        end = items[4].replace("# family: ","").strip()
                        instance_evalue = items[10].replace("# family: ","").strip()
                        annoOut[id][toolName][toolName+"_"+type]["instance"].append((start,end,instance_evalue))
    for id in inSeq:
        if not id in annotatedSeq:
            annoOut[id] = {}
            annoOut[id]["length"] = len(inSeq[id])
            annoOut[id][toolName] = {}
    subprocess.run(['rm', hmmOut])
    return(annoOut)
