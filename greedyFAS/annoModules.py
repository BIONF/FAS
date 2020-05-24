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

import time # remove later

def appendToDict(key,value,aDict):
    if not key in aDict.keys():
        aDict[key] = []
        aDict[key].append(value)
    else:
        if not value in aDict[key]:
            aDict[key].append(value)

# def doPfam(args):
#     (seq, toolPath, cpus) = args
#     pfamPath = toolPath + "/Pfam"
#     outPfam = pfamPath + "/output_files"

def doFlps(args):
    (seqFile, toolPath, threshold) = args
    # load fasta seq
    inSeq = SeqIO.to_dict((SeqIO.parse(open(seqFile),'fasta')))
    # do fLPS
    cmd = '%s/fLPS/fLPS -s -t %s %s' % (toolPath, threshold, seqFile)
    flpsOut = subprocess.run([cmd], shell=True, capture_output=True)
    lines = flpsOut.stdout.decode().split('\n')
    # save to dict
    annoOut = {} # annoOut[seqID,seqLen][(type,instance)] = [(start,end),(start,end),(start,end)]
    annotatedSeq = {}
    forCounting = {}
    if len(lines) > 0:
        for line in lines:
            tmp = line.strip().split('\t')
            if tmp[0]:
                if not (tmp[0], len(inSeq[tmp[0]])) in forCounting:
                    forCounting[(tmp[0], len(inSeq[tmp[0]]))] = {}
                appendToDict(tmp[1] + "_" + tmp[7], (tmp[3], tmp[4]), forCounting[(tmp[0], len(inSeq[tmp[0]]))])
                annotatedSeq[tmp[0]] = 1
    for key in forCounting:
        annoOut[key] = {}
        for feature in forCounting[key]:
            appendToDict((feature, len(forCounting[key][feature])), forCounting[key][feature], annoOut[key])

    for id in inSeq:
        if not id in annotatedSeq:
            annoOut[(id, len(inSeq[id]))] = {}
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
    annoOut = {} # annoOut[seqID,seqLen][(type,instance)] = [(start,end),(start,end),(start,end)]
    annotatedSeq = {}
    if len(lines) > 0:
        for line in lines:
            if line.startswith(">"):
                id = line.replace(">","").strip()
                annoOut[(id, len(inSeq[id]))] = {}
            posList = []
            if "%pred" in line:
                tmp = line.strip().split(',')
                for item in tmp:
                    if item.startswith(" M"):
                        pos = item.split()
                        posList.append((pos[1], pos[2]))
            if len(posList) > 0:
                appendToDict(("transmembrane", len(posList)), posList, annoOut[(id, len(inSeq[id]))])
                annotatedSeq[id] = 1
    for id in inSeq:
        if not id in annotatedSeq:
            annoOut[(id, len(inSeq[id]))][("transmembrane", 0)] = {}
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
    annoOut = {} # annoOut[seqID,seqLen][(type,instance)] = [(start,end),(start,end),(start,end)]
    annotatedSeq = {}
    if len(lines) > 0:
        for line in lines:
            if not line.startswith("#"):
                if len(line) > 0:
                    tmp = line.strip().split()
                    annoOut[(tmp[0], len(inSeq[tmp[0]]))] = {}
                    if tmp[9] == "Y":
                        appendToDict(("SIGNALP", 1), (1, int(tmp[4]) - 1), annoOut[(tmp[0], len(inSeq[tmp[0]]))])
                    else:
                        annoOut[(tmp[0], len(inSeq[tmp[0]]))][("SIGNALP", 0)] = {}
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
    annoOut = {} # annoOut[seqID,seqLen][(type,instance)] = [(start,end),(start,end),(start,end)]
    if len(results) > 0:
        for outSeq in results:
            if len(outSeq) > 0:
                tmp = outSeq.split('\n')
                id = tmp.pop(0).strip()
                concatSeq = ''.join(tmp).strip()
                annoOut[(id, len(inSeq[id]))] = {}
                posList = []
                finished = False
                while not finished:
                    match = re.search(r'x+', concatSeq)
                    if match:
                        start = concatSeq.find(match.group()) + 1
                        end = start + len(match.group()) - 1
                        posList.append((start, end))
                        concatSeq = concatSeq.replace(match.group(), "N" * len(match.group()), 1)
                    else:
                        finished = True
                if len(posList) > 0:
                    appendToDict(("coiled_coil", len(posList)), posList, annoOut[(id, len(inSeq[id]))])
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
    annoOut = {} # annoOut[seqID,seqLen][(type,instance)] = [(start,end),(start,end),(start,end)]
    if len(results) > 0:
        for result in results:
            if len(result) > 0:
                lines = result.split('\n')
                id = lines.pop(0).strip().replace(">","")
                annoOut[(id, len(inSeq[id]))] = {}
                posList = []
                for line in lines:
                    tmp = line.strip().split()
                    if (len(tmp) == 2) and ("-" in tmp[1]):
                        pos = tmp[1].split('-')
                        posList.append(pos)
                if len(posList) > 0:
                    appendToDict(("low complexity regions", len(posList)), posList, annoOut[(id, len(inSeq[id]))])
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
    smartOutFile = toolPath + "/SMART/output_files/" + lines[-1]
    # save to dict
    annoOut = {} # annoOut[seqID,seqLen][(type,instance,clan,evalue)] = [(evalue,start,end),(evalue,start,end),(evalue,start,end)]
    annotatedSeq = {}
    if ".out" in lines[-1]:
        with open(smartOutFile, 'r') as file:
            outFile = file.read()
            outSeqs = outFile.split('>')
            for outSeq in outSeqs:
                if (len(outSeq) > 3) and (not (outSeq.startswith(" queryID"))) and (not ("No hits detected" in outSeq)):
                    lines = outSeq.strip().split('\n')
                    id = lines.pop(0).strip()
                    annoOut[id, len(inSeq[id])] = {}
                    tmpDict = {}
                    for line in lines:
                        items = line.strip().split('|')
                        if "family:" in line:
                            type = items[0].replace("# family: ","").strip()
                            clan = items[2].replace(" clan: ","").strip()
                            type_evalue = items[-1].replace(" E-value: ","").strip()
                            tmpDict[(type,clan,type_evalue)] = []
                        else:
                            start = items[3].replace("# family: ","").strip()
                            end = items[4].replace("# family: ","").strip()
                            instance_evalue = items[10].replace("# family: ","").strip()
                            tmpDict[(type,clan,type_evalue)].append((instance_evalue,start,end))
                    for key in tmpDict:
                        appendToDict((key[0], len(tmpDict[key]), key[1], key[2]), tmpDict[key], annoOut[id, len(inSeq[id])])
                        annotatedSeq[id] = 1

    for id in inSeq:
        if not id in annotatedSeq:
            annoOut[(id, len(inSeq[id]))] = {}
    subprocess.run(['rm', smartOutFile])
    os.chdir(currentDir)
    return(annoOut)

# def write_xml(outpath, proteome, prot_lengths):
#     with open(outpath, 'w') as out:
#         out.write('<?xml version="1.0"?>\n<tool name="fLPS">\n')
#         for protein in proteome:
#             out.write('\t<protein id="' + protein + '" length="' + str(prot_lengths[protein]) + '">\n')
#             for feature in proteome[protein]:
#                 out.write('\t\t<feature type="' + feature + '" instance="' + str(len(proteome[protein][feature])) +
#                           '">\n')
#                 for instance in proteome[protein][feature]:
#                     out.write('\t\t\t<start start="' + instance[0] + '">\n\t\t\t<end end="' + instance[1] + '">\n')
#                 out.write('\t\t</feature>\n')
#             out.write('\t</protein>\n')
#         out.write('</tool>')
