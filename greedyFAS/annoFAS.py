#!/bin/env python

#######################################################################
# Copyright (C) 2019 Vinh Tran
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
from sys import platform
from pathlib import Path
import re
import subprocess
import argparse
import readline
import glob
import inspect
import greedyFAS
from os.path import expanduser

home = expanduser("~")


def get_path():
    return(os.path.dirname(inspect.getfile(greedyFAS)).strip())

def search_string_in_file(file_name, string_to_search):
    with open(file_name, 'r') as read_obj:
        for line in read_obj:
            if string_to_search in line:
                return(line.rstrip())

def check_status(perl_script):
    status = search_string_in_file(perl_script, "my $config")
    flag = 0
    if status == 'my $config = 0;':
        flag = 1
    else:
        anno_path_tmp = search_string_in_file(perl_script, "my $annotationPath")
        anno_path_tmp = anno_path_tmp.replace('my $annotationPath =', '')
        anno_path_tmp = re.sub(r'[;"\s]', '', anno_path_tmp)
        if not os.path.isfile(anno_path_tmp + '/Pfam/Pfam-hmms/Pfam-A.hmm'):
            flag = 1
    return flag

def get_toolPath(perl_script):
    status = search_string_in_file(perl_script, "my $config")
    flag = 0
    if status == 'my $config = 0;':
        sys.exit('annoFAS has not been run yet!')
    else:
        anno_path_tmp = search_string_in_file(perl_script, "my $annotationPath")
        anno_path_tmp = anno_path_tmp.replace('my $annotationPath =', '')
        anno_path_tmp = re.sub(r'[;"\s]', '', anno_path_tmp)
        if not os.path.isfile(anno_path_tmp + '/Pfam/Pfam-hmms/Pfam-A.hmm'):
            sys.exit('Annotation tools found at %s, but it seems to be not complete (e.g. PfamDB is missing)!' % anno_path_tmp)
        else:
            sys.exit('Annotation tools found at %s!' % anno_path_tmp)

def call_annoFAS_perl(options):
    # check status
    perl_script = get_path() + '/annoFAS.pl'
    if check_status(perl_script) == 1:
        sys.exit("Please run prepareFAS first!")

    # run annoFAS.pl
    current_dir = os.getcwd()
    os.chdir(current_dir)
    inputFullPath = os.path.abspath(options['fasta'])
    if "~" in options['fasta']:
        inputFullPath = options['fasta'].replace("~", home)
    outputFullPath = os.path.abspath(options['path'])
    if "~" in options['path']:
        outputFullPath = options['path'].replace("~", home)
    required_args = '--fasta %s --path %s --name %s' % (inputFullPath, outputFullPath,
                                                        options['name'])
    optional_args = '--force %s --extract %s --redo %s' % (options['force'], options['extract'], options['redo'])
    cmd = 'perl ' + perl_script + ' ' + required_args
    if options['cores']:
        cmd = cmd + ' --cores ' + options['cores']
    if options['force']:
        cmd = cmd + ' --force'
    if not options['extract'] == '':
        if not os.path.isdir(os.path.abspath(options['extract'])):
            os.mkdir(os.path.abspath(options['extract']))
        cmd = cmd + ' --extract ' + os.path.abspath(options['extract'])
    if options['redo'] in ['flps', 'coils', 'seg', 'pfam', 'signalp', 'smart', 'tmhmm']:
        cmd = cmd + ' --redo ' + options['redo']
    # print(cmd)
    subprocess.call([cmd], shell=True)

def main():
    version = "1.0.3"
    parser = argparse.ArgumentParser(description="You are running annoFAS version " + str(version) + ".")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--fasta', help='Input sequence in fasta format', action='store', default='',
                          required=True)
    required.add_argument('-o', '--path', help='Output directory', action='store', default='', required=True)
    required.add_argument('-n', '--name', help='Name of output annotation folder', action='store', default='',
                          required=True)
    optional.add_argument('-e', '--extract', help='Path to save the extracted annotation for input sequence',
                          action='store', default='')
    optional.add_argument('-r', '--redo', help='Re-annotation the sequence with flps|coils|seg|pfam|signalp|smart|tmhmm. '
                                         'Only one selection allowed!', action='store', default=0)
    optional.add_argument('-f', '--force', help='Force override annotations', action='store_true')
    optional.add_argument('-c', '--cores', help='number of cores', action='store', default='')
    optional.add_argument('--getToolPath', help='Get path to annotation tools', action='store_true')

    args = parser.parse_args()

    options = {
        'fasta': args.fasta,
        'path': args.path,
        'name': args.name,
        'extract': args.extract,
        'redo': args.redo,
        'force': args.force,
        'cores': args.cores,
    }
    if args.getToolPath:
        perl_script = get_path() + '/annoFAS.pl'
        get_toolPath(perl_script)
    else:
        call_annoFAS_perl(options)


if __name__ == '__main__':
    main()
