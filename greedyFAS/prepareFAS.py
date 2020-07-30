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
from os.path import expanduser

home = expanduser('~')

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

def query_yes_no(question, default='yes'):
    valid = {'yes': True, 'y': True, 'ye': True,
             'no': False, 'n': False}
    if default is None:
        prompt = ' [y/n] '
    elif default == 'yes':
        prompt = ' [Y/n] '
    elif default == 'no':
        prompt = ' [y/N] '
    else:
        raise ValueError('invalid default answer: "%s"' % default)
    while True:
        # sys.stdout.write(question + prompt)
        choice = sys.stdin.readline().rstrip().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write('Please respond with "yes" or "no" '
                             '(or "y" or "n").\n')


def get_dtu_path(dtuPathIn):
    if not dtuPathIn == '':
        if not os.path.isdir(dtuPathIn):
            sys.exit(dtuPathIn + ' not found!')
        else:
            dtu_path = os.path.abspath(dtuPathIn)
    else:
        print('Do you want to use TMHMM (v2.0c) and SignalP (v4.1g)? (y/n)')
        if query_yes_no('Source path?'):
            readline.set_completer_delims(' \t\n;')
            readline.parse_and_bind('tab: complete')
            readline.set_completer(complete)
            dtu_path = input('Please download TMHMM 2.0c and SignalP 4.1g at https://services.healthtech.dtu.dk/ and enter path to the downloaded tar files:')
        else:
            dtu_path = ''

    if not dtu_path == '':
        cmd = 'ls %s | grep \'signalp-4.1g\|tmhmm-2.0c\'' % dtu_path
        dtu_tool = subprocess.check_output([cmd], shell=True).decode(sys.stdout.encoding).strip()
        if not 'tmhmm' in dtu_tool:
            sys.exit('TMHMM not found in %s' % dtu_path)
        if not 'signalp' in dtu_tool:
            sys.exit('signalp-4.1g not found in %s' % dtu_path)
        return dtu_path, dtu_tool
    else:
        return 0, 0


def check_status(toolPath, force, tarfile):
    flag = 1
    if os.path.isfile(toolPath+'/annoTools.txt'):
        with open(toolPath+'/annoTools.txt') as f:
            if '#checked' in f.read():
                if force:
                    print('Annotation tools found in %s will be deleted and reinstalled! Enter to continue.' % toolPath)
                    if query_yes_no(''):
                        backupCmd = 'mv %s/%s ../' % (toolPath, tarfile)
                        rm_annotools = 'rm -rf %s/*' % (toolPath)
                        mvCmd = 'mv ../%s %s/' % (tarfile, toolPath)
                        subprocess_cmd([backupCmd, rm_annotools, mvCmd])
                        flag = 1
                    else:
                        flag = toolPath
                else:
                    flag = toolPath
        f.close()
    return flag

def install_signalp():
    mv_cmd1 = 'mv signalp-4.1* SignalP'
    subprocess.call([mv_cmd1], shell=True)
    os.chdir('SignalP')
    signalp_path = os.getcwd().replace('/','\/')
    addPath_signalp = 'sed -i -e \'s/$ENV{SIGNALP} = .*/$ENV{SIGNALP} = \"%s\";/\' %s' % (signalp_path, signalp_path+'/signalp')
    subprocess.call([addPath_signalp], shell=True)
    # makelink_signalp = 'ln -s -f bin/signalp signalp' # for signalp 5.0
    # subprocess.call([makelink_signalp], shell=True)

def install_tmhmm():
    mv_cmd2 = 'mv tmhmm* TMHMM'
    subprocess.call([mv_cmd2], shell=True)
    machine = os.uname()[4]
    os.chdir('TMHMM')
    makelink_tmhmm = 'ln -s -f bin/decodeanhmm.Linux_%s decodeanhmm' % machine
    subprocess.call([makelink_tmhmm], shell=True)

def prepare_annoTool(options):
    anno_path = options['toolPath']
    dtuPathIn = options['dtuPath']
    force = options['force']

    file = 'annotation_FAS2020d.tar.gz'
    checksum = '1818703744 1108970315 ' + file

    if os.path.exists(os.path.abspath(anno_path+'/annoTools.txt')):
        with open(anno_path+'/annoTools.txt') as anno_path:
            if '#checked' in anno_path.read():
                print('Annotation tools already found at %s' % anno_path)
                if not os.path.exists(os.path.abspath(options['greedyFasPath']+'/pathconfig.txt')):
                    print('If you want to add %s to config file of FAS, please rerun this function with --savePath!' % anno_path)
                print('If you want to re-install them, rerun this function with --force!')
                sys.exit()

    current_dir = os.getcwd()
    if check_status(anno_path, force, file) == 1:
        Path(anno_path).mkdir(parents = True, exist_ok = True)
        anno_path = os.path.abspath(anno_path)
        os.chdir(anno_path)

        print('Annotation tools will be installed in ')
        print(os.getcwd())
        print('----------------------------------')
        tools = ['fLPS', 'Pfam', 'SMART', 'COILS2', 'SEG'] #, 'SignalP', 'TMHMM']
        with open('annoTools.txt', mode = 'wt') as tool_file:
            if platform == 'darwin':
                tool_file.write('#linearized\nPfam\nSMART\n#normal\nfLPS\nCOILS2\n')
            else:
                tool_file.write('#linearized\nPfam\nSMART\n#normal\nfLPS\nCOILS2\nSEG\n')

        # get DTU tools path
        (dtu_path, dtu_tool) = get_dtu_path(dtuPathIn)

        # install TMHMM and SignalP
        if not dtu_path == 0:
            print('Installing SignalP and TMHMM...')
            for tool in dtu_tool.split('\n'):
                print(dtu_path + '/' + tool)
                tar_cmd = 'tar xf %s/%s --directory %s 2> /dev/null' % (dtu_path, tool, anno_path)
                subprocess.call([tar_cmd], shell=True)

            if not os.path.isdir(anno_path + '/SignalP'):
                install_signalp()
            else:
                rm_signalp = 'rm -rf signalp*'
                subprocess.call([rm_signalp], shell=True)
            os.chdir(anno_path)

            if not os.path.isdir(anno_path + '/TMHMM'):
                install_tmhmm()
            else:
                rm_tmhmm = 'rm -rf tmhmm*'
                subprocess.call([rm_tmhmm], shell=True)
            os.chdir(anno_path)

            with open('annoTools.txt', mode = 'a') as tool_file:
                if platform == 'darwin':
                    tool_file.write('SignalP\n')
                else:
                    tool_file.write('TMHMM\nSignalP\n')
        tool_file.close()

        # create folders for other annotation tools
        folders = ['fLPS', 'COILS2', 'Pfam', 'Pfam/Pfam-hmms', 'SEG', 'SMART'] #, 'SignalP', 'TMHMM']
        emptyFlag = 0
        for folder in folders:
            Path(anno_path + '/' + folder).mkdir(parents = True, exist_ok = True)
            if len(os.listdir(anno_path + '/' + folder)) == 0:
                print(anno_path + '/' + folder)
                emptyFlag = 1

        # download annotation tools
        if emptyFlag == 1:
            print('Installing other tools...')
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
                if tool == 'fLPS':
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
            source = os.getcwd() + '/fLPS/bin'
            target = os.getcwd() + '/fLPS/'
            if platform == 'darwin':
                source = source + '/mac64/fLPS'
                subprocess.call(['ln', '-fs', source, target])
            else:
                source = source + '/linux/fLPS'
                subprocess.call(['ln', '-fs', source, target])
            # re-compile COILS2
            coils_path = anno_path + '/COILS2'
            os.chdir(coils_path)
            subprocess.call(['tar', 'xf', 'ncoils.tar.gz'])
            coils_bin = coils_path + '/coils'
            os.chdir(coils_bin)
            compile_cmd = 'cc -O2 -I. -o ncoils-osf ncoils.c read_matrix.c -lm'
            if platform == 'darwin':
                COILSDIR = 'echo \'export COILSDIR=%s\' >> ~/.bash_profile' % coils_bin
            else:
                COILSDIR = 'echo \'export COILSDIR=%s\' >> ~/.bashrc' % coils_bin
            subprocess_cmd((compile_cmd, COILSDIR))
            os.chdir(coils_path)
            makelink_cmd = 'ln -s -f coils/ncoils-osf ./COILS2'
            if platform == 'darwin':
                source_cmd = 'source ~/.bash_profile'
            else:
                source_cmd = 'source ~/.bashrc'
            subprocess_cmd((makelink_cmd, source_cmd))
            os.chdir(anno_path)

            # remove temp files
            subprocess.call(['rm', '-rf', anno_path + '/annotation_FAS'])
            if not options['keep']:
                subprocess.call(['rm', anno_path + '/' + file])
        os.chdir(current_dir)
        return(anno_path)
    else:
        os.chdir(current_dir)
        return(check_status(anno_path, force, file))


def checkExecutable(anno_path):
    with open(anno_path+'/annoTools.txt') as file:
        availTool = [line.strip() for line in file]
    sedCmd ='sed -i \'s/#checked//\' %s/annoTools.txt' % anno_path
    subprocess.call([sedCmd], shell=True)

    print('Checking if annotation tools are excutable...')
    # test pfam and smart
    if not os.path.isfile(anno_path + '/Pfam/Pfam-hmms/Pfam-A.hmm'):
        sys.exit('Pfam hmm file not found. Please run prepareFAS with --force!')
    if not os.path.isfile(anno_path + '/SMART/SMART-hmms/SMART.hmm'):
        sys.exit('SMART hmm file not found. Please run prepareFAS with --force!')
    # test seg
    if 'SEG' in availTool:
        segCmd = '%s/SEG/seg' % anno_path
        try:
            p1 = subprocess.Popen([segCmd], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output1, err1 = p1.communicate()
            if not 'Usage' in err1.decode('UTF-8').strip():
                sys.exit('Error with SEG. You can reinstall it by running prepareFAS with --force!')
        except:
            sys.exit('Error with SEG. You can reinstall it by running prepareFAS with --force!')
    # test fLPS
    if 'fLPS' in availTool:
        flpsCmd = '%s/fLPS/fLPS' % anno_path
        try:
            p2 = subprocess.Popen([flpsCmd], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output2, err2 = p2.communicate()
            if not err2.decode('UTF-8').strip() == 'There is no sequence file. Please supply one.':
                sys.exit('Error with fLPS. You can reinstall it by running prepareFAS with --force!')
        except:
            sys.exit('Error with fLPS. You can reinstall it by running prepareFAS with --force!')
    # test COILS2
    if 'COILS2' in availTool:
        coilsCmd = '%s/COILS2/COILS2' % anno_path
        try:
            flag = 1
            p3 = subprocess.Popen([coilsCmd], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output3, err3 = p3.communicate()
            if err3.decode('UTF-8').strip() == 'Error reading '+os.getcwd()+'/new.mat':
                flag = 0
            if err3.decode('UTF-8').strip() == 'error: environment variable COILSDIR must be set':
                print('NOTE: THE TERMINAL MUST BE RESTARTED BEFORE USING FAS!!!')
                flag = 0
            if '0 sequences' in err3.decode('UTF-8').strip():
                flag = 0
            if flag == 1:
                sys.exit('Error with COILS2. You can reinstall it by running prepareFAS with --force!')
        except:
            sys.exit('Error with COILS2. You can reinstall it by running prepareFAS with --force!')
    # test tmhmm
    if 'TMHMM' in availTool:
        tmhmmCmd = '%s/TMHMM/decodeanhmm' % anno_path
        try:
            p4 = subprocess.Popen([tmhmmCmd], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output4, err4 = p4.communicate()
            if not 'Error: No modelfile given' in err4.decode('UTF-8').strip():
                sys.exit('Error with TMHMM. You can reinstall it by running prepareFAS with --force!')
        except:
            sys.exit('Error with TMHMM. You can reinstall it by running prepareFAS with --force!')
    # test signalp
    if 'SignalP' in availTool:
        signalpCmd = '%s/SignalP/signalp' % anno_path
        try:
            p5 = subprocess.Popen([signalpCmd, '-V'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output5, err5 = p5.communicate()
            if not err5.decode('UTF-8').strip() == '':
                sys.exit('Error with SignalP. You can reinstall it by running prepareFAS with --force!')
        except:
            sys.exit('Error with SignalP. You can reinstall it by running prepareFAS with --force!')
    return(True)

def checkAnnoToolsFile(toolPath):
    if not os.path.exists(os.path.abspath(toolPath+'/annoTools.txt')):
        sys.exit('ERROR: %s not found' % (toolPath+'/annoTools.txt'))
    else:
        with open(toolPath+'/annoTools.txt') as f:
            if not '#checked' in f.read():
                sys.exit('ERROR: Some errors occur with annotation tools. Please install them again!')

def saveConfigFile(checkResult, anno_path, greedyFasPath):
    if checkResult:
        with open(anno_path+'/annoTools.txt') as f:
            if not '#checked' in f.read():
                with open(anno_path+'/annoTools.txt','a') as file:
                    file.write('#checked')
                    file.close()
        f.close()
        with open(greedyFasPath+'/pathconfig.txt','w') as config:
            config.write(os.path.abspath(anno_path))
            config.close()
        sys.exit('Done! Annotation tools can be found in %s' % anno_path)
    else:
        sys.exit('Some errors occur with annotation tools. Please check if they can be excuted at %s' % anno_path)

def main():
    version = '1.2.3'
    parser = argparse.ArgumentParser(description='You are running prepareFAS version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-t', '--toolPath', help='Set path to save annotation tools', action='store', default='', required=True)
    optional.add_argument('-d', '--dtuPath', help='Set path to DTU tools (SignalP and TMHMM)', action='store', default='')
    optional.add_argument('-f', '--force', help='Overwrite old annotation tools if exist', action='store_true')
    optional.add_argument('-k', '--keep', help='Keep downloaded source file', action='store_true')
    optional.add_argument('-s', '--savePath', help='Save annotation tool path to config file for FAS', action='store_true')
    optional.add_argument('-c', '--check', help='Check if FAS ready to run. NO real tool path need to be given!', action='store_true')
    optional.add_argument('--checkExecutable', help='Check if annotation tools are executable!', action='store_true')

    args = parser.parse_args()
    greedyFasPath = os.path.realpath(__file__).replace('/prepareFAS.py','')
    options = {
        'toolPath': args.toolPath,
        'dtuPath': args.dtuPath,
        'force': args.force,
        'keep': args.keep,
        'greedyFasPath': greedyFasPath
    }

    if args.checkExecutable:
        if not os.path.exists(os.path.abspath(args.toolPath+'/annoTools.txt')):
            sys.exit('ERROR: %s not found' % (args.toolPath+'/annoTools.txt'))
        else:
            allRun = checkExecutable(args.toolPath)
            saveConfigFile(allRun, args.toolPath, greedyFasPath)

    if args.check:
        if not os.path.exists(os.path.abspath(greedyFasPath+'/pathconfig.txt')):
            sys.exit('ERROR: %s not found' % (greedyFasPath+'/pathconfig.txt'))
        else:
            with open(greedyFasPath+'/pathconfig.txt', 'r') as file:
                savedPath = file.read().strip()
                checkAnnoToolsFile(savedPath)
                sys.exit('Annotation tools can be found at %s. FAS is ready to run!' % savedPath)

    if args.savePath:
        checkAnnoToolsFile(args.toolPath)
        with open(greedyFasPath+'/pathconfig.txt','w') as config:
            config.write(os.path.abspath(args.toolPath))
            config.close()
        sys.exit('Annotation tools can be found at %s. FAS is ready to run!' % args.toolPath)

    anno_path = prepare_annoTool(options)
    allRun = checkExecutable(anno_path)
    saveConfigFile(allRun, anno_path, greedyFasPath)

if __name__ == '__main__':
    main()
