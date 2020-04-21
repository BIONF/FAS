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
from sys import platform
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


def complete(text, state):
    return(glob.glob(os.path.expanduser(text)+'*')+[None])[state]


def subprocess_cmd(commands):
    for cmd in commands:
        subprocess.call(cmd, shell = True)


def download_data(file, checksum):
    url = 'https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo' + '/' + file
    parameter = '--no-check-certificate'
    subprocess.call(['wget', url, parameter])
    if os.path.isfile(file):
        checksum_file = subprocess.check_output(['cksum', file]).decode(sys.stdout.encoding).strip()
        if checksum_file == checksum:
            print('Extracting %s ...' % file)
            subprocess.call(['tar', 'xf', file])
        else:
            sys.exit('Downloaded file corrupted!')
    else:
        sys.exit('Cannot download annotation tools!')


def query_yes_no(question, default="yes"):
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
        # sys.stdout.write(question + prompt)
        choice = sys.stdin.readline().rstrip().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def install_path():
    anno_path = expanduser("~") + "/annotation_fas"
    print("Annotation dir: %s (y/n)" % anno_path)
    if not query_yes_no("anno_path?"):
        readline.set_completer_delims(' \t\n;')
        readline.parse_and_bind("tab: complete")
        readline.set_completer(complete)
        anno_path = input('Enter annotation dir: ')
    return(anno_path)


def get_dtu_path():
    print("Do you want to use TMHMM and SignalP? (y/n)")
    if query_yes_no("Source path?"):
        readline.set_completer_delims(' \t\n;')
        readline.parse_and_bind("tab: complete")
        readline.set_completer(complete)
        dtu_path = input('Enter path to TMHMM and SignalP tar files: ')
        cmd = 'ls %s | grep \"signalp\|tmhmm\"' % dtu_path
        dtu_tool = subprocess.check_output([cmd], shell=True).decode(sys.stdout.encoding).strip()
        if not "tmhmm" in dtu_tool:
            sys.exit('TMHMM not found in %s' % dtu_path)
        if not "signalp" in dtu_tool:
            sys.exit('SignalP not found in %s' % dtu_path)
        return dtu_path, dtu_tool
    else:
        return 0, 0


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
        # else:
        #     print('Annotation tools found in %s!' % anno_path_tmp)
    return flag


def prepare_annoTool(annoPathIn):
    current_dir = os.getcwd()
    perl_script = get_path() + '/annoFAS.pl'
    if check_status(perl_script) == 1:
        print('Annotation tools need to be downloaded!')
        # get DTU tools path
        (dtu_path, dtu_tool) = get_dtu_path()

        # get annotation directory
        if not annoPathIn == '':
            if os.path.isdir(os.path.abspath(annoPathIn)):
                anno_path = os.path.abspath(annoPathIn)
            else:
                anno_path = install_path()
        else:
            anno_path = install_path()
        if not os.path.isdir(anno_path):
            os.mkdir(anno_path)
        anno_path = os.path.abspath(anno_path)
        os.chdir(anno_path)
        print('Annotation tools will be saved in ')
        print(os.getcwd())
        print('----------------------------------')
        tools = ['fLPS', 'Pfam', 'SMART', 'COILS2', 'SEG'] #, 'SignalP', 'TMHMM']
        with open('annoTools.txt', mode = 'wt') as tool_file:
            tool_file.write("#linearized\nPfam\nSMART\n#normal\nfLPS\nCOILS2\nSEG\n")

        # install TMHMM and SignalP
        if not dtu_path == 0:
            print('Installing SignalP and TMHMM...')
            for tool in dtu_tool.split('\n'):
                print(dtu_path + "/" + tool)
                tar_cmd = 'tar xf %s/%s --directory %s' % (dtu_path, tool, anno_path)
                subprocess.call([tar_cmd], shell=True)

            mv_cmd1 = 'mv signalp* SignalP'
            subprocess.call([mv_cmd1], shell=True)
            os.chdir("SignalP")
            makelink_signalp = 'ln -s -f bin/signalp signalp'
            subprocess.call([makelink_signalp], shell=True)
            os.chdir(anno_path)

            mv_cmd2 = 'mv tmhmm* TMHMM'
            subprocess.call([mv_cmd2], shell=True)
            machine = os.uname()[4]
            os.chdir("TMHMM")
            makelink_tmhmm = 'ln -s -f bin/decodeanhmm.Linux_%s decodeanhmm' % machine
            subprocess.call([makelink_tmhmm], shell=True)
            os.chdir(anno_path)

            dtu_tool_list = ["TMHMM", "SignalP"]
            with open('annoTools.txt', mode = 'a') as tool_file:
                tool_file.write("TMHMM\nSignalP\n")
        tool_file.close()

        print('Installing other tools...')
        # create folders for other annotation tools
        folders = ['fLPS', 'COILS2', 'Pfam', 'Pfam/Pfam-hmms', 'Pfam/output_files', 'SEG', 'SMART'] #, 'SignalP', 'TMHMM']
        for folder in folders:
            if not os.path.isdir(anno_path + '/' + folder):
                os.mkdir(anno_path + '/' + folder)

        # download annotation tools
        if not os.path.isfile('Pfam/Pfam-hmms/Pfam-A.hmm'):
            file = 'annotation_FAS2020b.tar.gz'
            checksum = '4256933429 1119115794 ' + file
            if os.path.isfile(file):
                checksum_file = subprocess.check_output(['cksum', file]).decode(sys.stdout.encoding).strip()
                if checksum_file == checksum:
                    print('Extracting %s ...' % file)
                    subprocess.call(['tar', 'xf', file])
                else:
                    subprocess.call(['rm', file])
                    download_data(file, checksum)
            else:
                download_data(file, checksum)

            # copy tools to their folders
            for tool in tools:
                if tool == "fLPS":
                    print('Downloading fLPS ...')
                    fLPS_file = 'fLPS.tar.gz'
                    fLPS_url = 'http://biology.mcgill.ca/faculty/harrison/' + fLPS_file
                    subprocess.call(['wget', fLPS_url, '--no-check-certificate'])
                    subprocess.call(['tar', 'xf', fLPS_file])
                    subprocess.call(['rm', fLPS_file])
                else:
                    print('Moving %s ...' % tool)
                    source_dir = 'annotation_FAS/' + tool + '/'
                    target_dir = tool + '/'
                    subprocess.call(['rsync', '-ra', '--include=*', source_dir, target_dir])
                print('---------------------')

            # make symlink for fLPS (depend on OS system)
            source = os.getcwd() + "/fLPS/bin"
            target = os.getcwd() + "/fLPS/"
            if platform == "darwin":
                source = source + "/mac64/fLPS"
                subprocess.call(['ln', '-fs', source, target])
            else:
                source = source + "/linux/fLPS"
                subprocess.call(['ln', '-fs', source, target])
            # re-compile COILS2 for mac OS
            if platform == "darwin":
                coils_path = anno_path + "/COILS2"
                os.chdir(coils_path)
                subprocess.call(['tar', 'xf', 'ncoils.tar.gz'])
                coils_bin = coils_path + "/coils"
                os.chdir(coils_bin)
                compile_cmd = 'cc -O2 -I. -o ncoils-osf ncoils.c read_matrix.c -lm'
                COILSDIR = 'echo \"export COILSDIR=%s\" >> ~/.bash_profile' % coils_bin
                subprocess_cmd((compile_cmd, COILSDIR))
                os.chdir(coils_path)
                makelink_cmd = 'ln -s -f coils/ncoils-osf ./COILS2'
                source_cmd = 'source ~/.bash_profile'
                subprocess_cmd((makelink_cmd, source_cmd))
                os.chdir(anno_path)

            # remove temp files
            subprocess.call(['rm', '-rf', anno_path + '/annotation_FAS'])
            subprocess.call(['rm', anno_path + '/' + file])

            # add path to annotation dir to annFASpl script
            mod_anno_path = anno_path.replace('/', '\/')
            sed_cmd1 = 'sed -i -e \'s/my $annotationPath = .*/my $annotationPath = \"%s\";/\' %s' % (mod_anno_path,
                                                                                                     perl_script)
            sed_cmd2 = 'sed -i -e \'s/$config = 0/$config = 1/\' %s' % perl_script
            subprocess.call([sed_cmd1], shell=True)
            subprocess.call([sed_cmd2], shell=True)
            print('Annotation tools downloaded!')
        else:
            print('Annotation tools found!')
            # add path to annotation dir to annFAS.pl script
            mod_anno_path = anno_path.replace('/', '\/')
            sed_cmd1 = 'sed -i -e \'s/my $annotationPath = .*/my $annotationPath = \"%s\";/\' %s' % (mod_anno_path,
                                                                                                     perl_script)
            sed_cmd2 = 'sed -i -e \'s/$config = 0/$config = 1/\' %s' % perl_script
            subprocess.call([sed_cmd1], shell=True)
            subprocess.call([sed_cmd2], shell=True)
    os.chdir(current_dir)


def call_annoFAS_perl(options):
    # check status
    perl_script = get_path() + '/annoFAS.pl'
    if check_status(perl_script) == 1:
        prepare_annoTool(options['annoPath'])

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
    optional.add_argument('--prepare', help='Download annotation tools and do configuration', action='store_true')
    optional.add_argument('-p', '--annoPath', help='Path to annotation dir', action='store', default='')
    optional.add_argument('-c', '--cores', help='number of cores', action='store', default='')

    args = parser.parse_args()

    options = {
        'fasta': args.fasta,
        'path': args.path,
        'name': args.name,
        'extract': args.extract,
        'redo': args.redo,
        'force': args.force,
        'annoPath': args.annoPath,
        'cores': args.cores
    }
    if args.prepare:
        prepare_annoTool(args.annoPath)
    else:
        call_annoFAS_perl(options)


if __name__ == '__main__':
    main()
