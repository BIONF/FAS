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

import getopt
import os
import sys
import subprocess
import argparse

def subprocess_cmd(commands):
	for cmd in commands:
		subprocess.call(cmd, shell = True)

def getPath():
    fasInfo = subprocess.check_output(['pip', 'show', 'greedyFAS']).decode(sys.stdout.encoding).strip()
    for line in fasInfo.split('\n'):
        if 'Location' in line:
            return(line.replace('Location: ', '').rstrip() + '/greedyFAS')

def downloadData(file, checksum, output):
    url = 'https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo' + '/' + file
    parameter = '--no-check-certificate'
    subprocess.call(['wget', '-O', output, url, parameter])
    if os.path.isfile(file):
        checksumFile = subprocess.check_output(['cksum', file]).decode(sys.stdout.encoding).strip()
        if checksumFile == checksum:
            subprocess.call(['tar', 'xfv', file])
        else:
            sys.exit('Downloaded file corrupted!')
    else:
        sys.exit('Cannot download annotation tools!')

def parse_arguments():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('--fasta',help = 'Input sequence in fasta format', action = 'store', default = '', required = True)
    required.add_argument('--path', help = 'Output directory', action = 'store', default = '', required = True)
    required.add_argument('--name', help = 'Query or Proteom [q or p]', action = 'store', default = '', required = True)
    optional.add_argument('--extract', help = 'Path to save the extracted annotation for input sequence', action = 'store', default = '')
    optional.add_argument('--redo', help = 'Re-annotation the sequence with cast|coils|seg|pfam|signalp|smart|tmhmm. Only one selection allowed!', action = 'store', default = 0)
    optional.add_argument('--force', help = 'Force override annotations [1, default = 0]', action = 'store', default = 0)
    args = parser.parse_args()
    return args

def main(args):
    currentDir = os.getcwd()
    # create folders for annotation tools
    fasPath = getPath()
    folders = ['CAST', 'COILS2', 'Pfam', 'Pfam/Pfam-hmms', 'Pfam/output_files', 'SEG', 'SignalP', 'SMART', 'TMHMM']
    for folder in folders:
        if not os.path.isdir(fasPath + '/' + folder):
            os.mkdir(fasPath + '/' + folder)

    # download annotation tools
    os.chdir(fasPath)
    print('Installed FAS folder:')
    print(os.getcwd())
    if not os.path.isfile('Pfam/Pfam-hmms/Pfam-A.hmm'):
        file = 'data_HaMStR.tar'
        checksum = '4100986910 5840435200 ' + file
        if os.path.isfile(file):
            checksumFile = subprocess.check_output(['cksum', file]).decode(sys.stdout.encoding).strip()
            if checksumFile == checksum:
                subprocess.call(['tar', 'xfv', file])
            else:
                subprocess.call(['rm', file])
                downloadData(file, checksum)
        else:
            downloadData(file, checksum)

        # copy tools to their folders
        tools = ['Pfam', 'SMART', 'CAST', 'COILS2', 'SEG', 'SignalP', 'TMHMM']
        for tool in tools:
            print('Moving %s ...' % tool)
            sourceDir = 'data_HaMStR/' + tool + '/'
            targetDir = tool + '/'
            subprocess.call(['rsync', '-rva', '--include=*', sourceDir, targetDir])
            print('---------------------')

        # remove temp files
        subprocess.call('rm', '-rf', 'data_HaMStR')
        subprocess.call('rm', file)
        print('Annotation tools downloaded!')
    else:
        print('Annotation tools found!')

    # run annotation.pl script
    os.chdir(currentDir)
    print('Current working dir:')
    print(os.getcwd())
    perlScript = fasPath + '/annoFAS.pl'
    requiredArgs = '--fasta %s --path %s --name %s' % (args.fasta, args.path, args.name)
    optionalArgs = '--force %s --extract %s --redo %s' % (args.force, args.extract, args.redo)
    subprocess.call(['perl', perlScript, requiredArgs, optionalArgs])

if __name__ == '__main__':
    arguments = parse_arguments()
    main(arguments)
