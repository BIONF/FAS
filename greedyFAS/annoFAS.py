#!/bin/env python

#######################################################################
# Copyright (C) 2019 Vinh Tran
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

import os
import sys
import re
import subprocess
import argparse
import readline
import glob
from os.path import expanduser


def getPath():
    fasInfo = subprocess.check_output(['pip', 'show', 'greedyFAS']).decode(sys.stdout.encoding).strip()
    for line in fasInfo.split('\n'):
        if 'Location' in line:
            return(line.replace('Location: ', '').rstrip() + '/greedyFAS')

def downloadData(file, checksum):
    url = 'https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo' + '/' + file
    parameter = '--no-check-certificate'
    subprocess.call(['wget', url, parameter])
    if os.path.isfile(file):
        checksumFile = subprocess.check_output(['cksum', file]).decode(sys.stdout.encoding).strip()
        if checksumFile == checksum:
            print('Extracting %s ...' % file)
            subprocess.call(['tar', 'xf', file])
        else:
            sys.exit('Downloaded file corrupted!')
    else:
        sys.exit('Cannot download annotation tools!')

def query_yes_no(question, default = "yes"):
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

def complete(text, state):
    return (glob.glob(os.path.expanduser(text)+'*')+[None])[state]

def installPath():
    annoPath = expanduser("~") + "/annotation_fas"
    print('Annotation dir:', end = ' ')
    if query_yes_no(annoPath) == False:
        readline.set_completer_delims(' \t\n;')
        readline.parse_and_bind("tab: complete")
        readline.set_completer(complete)
        annoPath = input('Enter annotation dir: ')
    return(annoPath)

def main():
    currentDir = os.getcwd()

    # get arguments
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('--fasta',help = 'Input sequence in fasta format', action = 'store', default = '', required = True)
    required.add_argument('--path', help = 'Output directory', action = 'store', default = '', required = True)
    required.add_argument('--name', help = 'Name of output annotation folder', action = 'store', default = '', required = True)
    optional.add_argument('--extract', help = 'Path to save the extracted annotation for input sequence', action = 'store', default = '')
    optional.add_argument('--redo', help = 'Re-annotation the sequence with cast|coils|seg|pfam|signalp|smart|tmhmm. Only one selection allowed!', action = 'store', default = 0)
    optional.add_argument('--force', help = 'Force override annotations [y/n, default = n]', action = 'store', default = 'n')
    optional.add_argument('--prepare', help = 'Download annotation tools and do configuration [y/n, default = n]', action = 'store', default = 'n')
    optional.add_argument('--annoPath', help = 'Path to annotation dir', action = 'store', default = '')
    args = parser.parse_args()

    # get config status and annoPath (if availible)
    perlScript = getPath() + '/annoFAS.pl'
    status = subprocess.check_output("grep 'my $config' %s" % perlScript, shell = True).decode(sys.stdout.encoding).strip()
    flag = 0
    if status == 'my $config = 0;':
        print('Annotation tools need to be downloaded!')
        flag = 1
    else:
        annoPathTmp = subprocess.check_output("grep 'my $annotationPath' %s" % perlScript, shell = True).decode(sys.stdout.encoding).strip()
        annoPathTmp = annoPathTmp.replace('my $annotationPath =', '')
        annoPathTmp = re.sub(r'[;"\s]', '', annoPathTmp)
        if not os.path.isfile(annoPathTmp + '/Pfam/Pfam-hmms/Pfam-A.hmm'):
            flag = 1
        else:
            print('Annotation tools found in %s!' % annoPathTmp)

    if flag == 1:
        # annotation directory
        if not args.annoPath == '':
            if os.path.isdir(os.path.abspath(args.annoPath)):
                annoPath = os.path.abspath(args.annoPath)
            else:
                annoPath = installPath()
        else:
            annoPath = installPath()
        if not os.path.isdir(annoPath):
            os.mkdir(annoPath)
        annoPath = os.path.abspath(annoPath)

        # create folders for annotation tools
        folders = ['CAST', 'COILS2', 'Pfam', 'Pfam/Pfam-hmms', 'Pfam/output_files', 'SEG', 'SignalP', 'SMART', 'TMHMM']
        for folder in folders:
            if not os.path.isdir(annoPath + '/' + folder):
                os.mkdir(annoPath + '/' + folder)

        # download annotation tools
        os.chdir(annoPath)
        print('Annotation tools will be saved in', end = ' ')
        print(os.getcwd())
        if not os.path.isfile('Pfam/Pfam-hmms/Pfam-A.hmm'):
            file = 'annotation_FAS2018b.tar.gz'
            checksum = '1548260242 1055860344 ' + file
            if os.path.isfile(file):
                checksumFile = subprocess.check_output(['cksum', file]).decode(sys.stdout.encoding).strip()
                if checksumFile == checksum:
                    print('Extracting %s ...' % file)
                    subprocess.call(['tar', 'xf', file])
                else:
                    subprocess.call(['rm', file])
                    downloadData(file, checksum)
            else:
                downloadData(file, checksum)

            # copy tools to their folders
            tools = ['Pfam', 'SMART', 'CAST', 'COILS2', 'SEG', 'SignalP', 'TMHMM']
            for tool in tools:
                print('Moving %s ...' % tool)
                sourceDir = 'annotation_FAS/' + tool + '/'
                targetDir = tool + '/'
                subprocess.call(['rsync', '-rva', '--include=*', sourceDir, targetDir])
                print('---------------------')

            # remove temp files
            subprocess.call(['rm', '-rf', annoPath + '/annotation_FAS'])
            subprocess.call(['rm', annoPath + '/' + file])

            # add path to annotation dir to annFASpl script
            modAnnoPath = annoPath.replace('/', '\/')
            sedCMD1 = 'sed -i -e \'s/my $annotationPath = .*/my $annotationPath = \"%s\";/\' %s' % (modAnnoPath, perlScript)
            sedCMD2 = 'sed -i -e \'s/$config = 0/$config = 1/\' %s' % (perlScript)
            subprocess.call([sedCMD1], shell = True)
            subprocess.call([sedCMD2], shell = True)
            print('Annotation tools downloaded!')
        else:
            print('Annotation tools found!')
            # add path to annotation dir to annFASpl script
            modAnnoPath = annoPath.replace('/', '\/')
            sedCMD1 = 'sed -i -e \'s/my $annotationPath = .*/my $annotationPath = \"%s\";/\' %s' % (modAnnoPath, perlScript)
            sedCMD2 = 'sed -i -e \'s/$config = 0/$config = 1/\' %s' % (perlScript)
            subprocess.call([sedCMD1], shell = True)
            subprocess.call([sedCMD2], shell = True)

    # run annotation.pl script
    if args.prepare == 'y':
        sys.exit('Config done!')

    os.chdir(currentDir)
    requiredArgs = '--fasta %s --path %s --name %s' % (os.path.abspath(args.fasta), os.path.abspath(args.path), args.name)
    optionalArgs = '--force %s --extract %s --redo %s' % (args.force, args.extract, args.redo)
    cmd = 'perl ' + perlScript + ' ' + requiredArgs
    if args.force == 'y':
        cmd = cmd + ' --force'
    if not args.extract == '':
        if os.path.isdir(os.path.abspath(args.extract)):
            cmd = cmd + ' --extract ' + os.path.abspath(args.extract)
    if args.redo in ['cast', 'coils', 'seg', 'pfam', 'signalp','smart', 'tmhmm']:
        cmd = cmd + ' --redo ' + args.redo
    # print(cmd)
    subprocess.call([cmd], shell = True)

if __name__ == '__main__':
    main()
