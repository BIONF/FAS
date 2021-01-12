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
import shutil
import subprocess
import argparse
import readline
import glob
from os.path import expanduser
import ssl
import urllib.request
import time

home = expanduser('~')


def complete(text, state):
    return(glob.glob(os.path.expanduser(text)+'*')+[None])[state]


def subprocess_cmd(commands):
    for cmd in commands:
        subprocess.call(cmd, shell=True)


def download_progress(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    if percent > 100:
        percent = 100
    sys.stdout.write("\r...%d%%, %d MB, %d KB/s, %d seconds passed" %
                     (percent, progress_size / (1024 * 1024), speed, duration))
    sys.stdout.flush()


def download_file(url, file):
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    download_file = urllib.request.URLopener(context=ctx)
    print('Downloading %s' % (url + '/' + file))
    download_file.retrieve(url + '/' + file, file, download_progress)
    print(' ... done!')


def download_data(file, checksum, toolPath):
    url = 'https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo'
    download_file(url, file)
    if os.path.isfile(file):
        checksum_file = subprocess.check_output(['cksum', file]).decode(sys.stdout.encoding).strip()
        if checksum_file == checksum:
            print('Extracting %s ...' % file)
            shutil.unpack_archive(file, toolPath, 'gztar')
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
            dtu_path = input('Please download TMHMM 2.0c and SignalP 4.1g at https://services.healthtech.dtu.dk/ and '
                             'enter path to the downloaded tar files:')
        else:
            dtu_path = ''

    if not dtu_path == '':
        files = [f for f in os.listdir(dtu_path) if os.path.isfile(os.path.join(dtu_path, f))]
        signalp_source = [i for i in files if 'signalp' in i]
        tmhmm_source = [i for i in files if 'tmhmm' in i]
        if not tmhmm_source:
            sys.exit('TMHMM not found in %s' % dtu_path)
        if not signalp_source:
            sys.exit('signalp-4.1g not found in %s' % dtu_path)
        return dtu_path, signalp_source[0], tmhmm_source[0]
    else:
        return 0, 0, 0


def check_status(toolPath, force, tarfile):
    flag = 1
    cwd = os.getcwd()
    if os.path.isfile(toolPath+'/annoTools.txt'):
        with open(toolPath+'/annoTools.txt') as f:
            if '#checked' in f.read():
                if force:
                    print('Annotation tools found in %s will be deleted and reinstalled! Enter to continue.' % toolPath)
                    if query_yes_no(''):
                        if os.path.exists(os.path.abspath(toolPath + '/' + tarfile)):
                            shutil.move(toolPath + '/' + tarfile, cwd + '/' + tarfile)
                        try:
                            shutil.rmtree(toolPath)
                        except:
                            sys.exit('Failed to delete %s. Please manually remove it and run setupFAS again!'
                                     % toolPath)
                        if os.path.exists(os.path.abspath(cwd + '/' + tarfile)):
                            shutil.move(cwd + '/' + tarfile, toolPath + '/' + tarfile)
                        flag = 1
                    else:
                        flag = toolPath
                else:
                    flag = toolPath
        f.close()
    return flag


def install_signalp():
    shutil.move('signalp-4.1', 'SignalP')
    subprocess.call(['chmod', '-R', '+w', './SignalP'])
    os.chdir('SignalP')
    signalp_path = os.getcwd()
    with open(signalp_path+'/signalp', "r") as infile, open(signalp_path+'/signalp.mod', "w") as outfile:
        for line in infile:
            if 'ENV{SIGNALP} = ' in line:
                line = '    $ENV{SIGNALP} = "%s";\n' % signalp_path
            outfile.write(line)
    shutil.move(signalp_path+'/signalp.mod', signalp_path+'/signalp')
    subprocess.call(['chmod', '0755', signalp_path+'/signalp'])
    # makelink_signalp = 'ln -s -f bin/signalp signalp' # for signalp 5.0
    # subprocess.call([makelink_signalp], shell=True)


def install_tmhmm():
    shutil.move('tmhmm-2.0c', 'TMHMM')
    machine = os.uname()[4]
    os.chdir('TMHMM')
    makelink_tmhmm = 'ln -s -f bin/decodeanhmm.Linux_%s decodeanhmm' % machine
    subprocess.call([makelink_tmhmm], shell=True)


def write_coilsdir(coilsdir):
    home = str(Path.home())
    filename = home + '/.bashrc'
    if platform == 'darwin':
        filename = home + '/.bash_profile'
    with open(filename, "r+") as file:
        for line in file:
            if coilsdir in line:
                break
        else:
            file.write(coilsdir)


def prepare_annoTool(options):
    anno_path = options['toolPath']
    dtuPathIn = options['dtuPath']
    force = options['force']

    file = 'annotation_FAS2020d.tar.gz'
    checksum = '1818703744 1108970315 ' + file

    if os.path.exists(os.path.abspath(anno_path+'/annoTools.txt')):
        with open(anno_path+'/annoTools.txt') as checkfile:
            if '#checked' in checkfile.read():
                if not force:
                    print('Annotation tools already found at %s' % anno_path)
                    if not os.path.exists(os.path.abspath(options['greedyFasPath']+'/pathconfig.txt')):
                        print('If you want to add %s to config file of FAS, please rerun this function with --savePath!'
                              % anno_path)
                    print('If you want to re-install them, rerun this function with --force!')
                    sys.exit()

    current_dir = os.getcwd()
    if check_status(anno_path, force, file) == 1:
        Path(anno_path).mkdir(parents=True, exist_ok=True)
        anno_path = os.path.abspath(anno_path)
        os.chdir(anno_path)

        print('Annotation tools will be installed in ')
        print(os.getcwd())
        print('----------------------------------')
        tools = ['fLPS', 'Pfam', 'SMART', 'COILS2', 'SEG'] #, 'SignalP', 'TMHMM']
        with open('annoTools.txt', mode='wt') as tool_file:
            if platform == 'darwin':
                tool_file.write('#linearized\nPfam\nSMART\n#normal\nfLPS\nCOILS2\n')
            else:
                tool_file.write('#linearized\nPfam\nSMART\n#normal\nfLPS\nCOILS2\nSEG\n')

        # get DTU tools path
        (dtu_path, signalp_source, tmhmm_source) = get_dtu_path(dtuPathIn)

        # install TMHMM and SignalP
        if not dtu_path == 0:
            print('Installing SignalP and TMHMM...')
            if not os.path.isdir(anno_path + '/SignalP'):
                shutil.unpack_archive(dtu_path + '/' + signalp_source, anno_path, 'gztar')
                install_signalp()
            os.chdir(anno_path)
            if not os.path.isdir(anno_path + '/TMHMM'):
                shutil.unpack_archive(dtu_path + '/' + tmhmm_source, anno_path, 'gztar')
                install_tmhmm()
            os.chdir(anno_path)
            with open('annoTools.txt', mode='a') as tool_file:
                if platform == 'darwin':
                    tool_file.write('SignalP\n')
                else:
                    tool_file.write('TMHMM\nSignalP\n')
        tool_file.close()

        # download other annotation tools
        print('Installing other tools...')
        if os.path.isfile(file):
            checksum_file = subprocess.check_output(['cksum', file]).decode(sys.stdout.encoding).strip()
            if checksum_file == checksum:
                print('Extracting %s ...' % file)
                shutil.unpack_archive(anno_path + '/' + file, anno_path, 'gztar')
            else:
                os.remove(file)
                download_data(file, checksum, anno_path)
        else:
            download_data(file, checksum, anno_path)

        # copy tools to their folders
        for tool in tools:
            if tool == 'fLPS':
                print('Downloading fLPS ...')
                fLPS_file = 'fLPS.tar.gz'
                fLPS_url = 'http://biology.mcgill.ca/faculty/harrison/'
                download_file(fLPS_url, fLPS_file)
                shutil.unpack_archive(fLPS_file, anno_path, 'gztar')
                os.remove(fLPS_file)
            else:
                print('Moving %s ...' % tool)
                source_dir = 'annotation_FAS/' + tool + '/'
                if os.path.exists(os.path.abspath(anno_path + '/' + tool)):
                    shutil.rmtree(anno_path + '/' + tool)
                shutil.move(source_dir, anno_path, copy_function=shutil.copytree)

        # make symlink for fLPS (depend on OS system)
        source = os.getcwd() + '/fLPS/bin'
        target = os.getcwd() + '/fLPS/fLPS'
        if platform == 'darwin':
            source = source + '/mac64/fLPS'
            try:
                os.symlink(source, target)
            except FileExistsError:
                os.remove(target)
                os.symlink(source, target)
        else:
            source = source + '/linux/fLPS'
            try:
                os.symlink(source, target)
            except FileExistsError:
                os.remove(target)
                os.symlink(source, target)

        # re-compile COILS2
        coils_path = anno_path + '/COILS2'
        os.chdir(coils_path)
        shutil.unpack_archive('ncoils.tar.gz', coils_path, 'gztar')
        coils_bin = coils_path + '/coils'
        os.chdir(coils_bin)
        compile_cmd = 'cc -O2 -I. -o ncoils-osf ncoils.c read_matrix.c -lm' + ' > /dev/null 2>&1'
        try:
            subprocess.call(compile_cmd, shell=True)
        except:
            print('ERROR: Failed to compile COILS2.\nPlease read instruction at %s and do it manually!' % coils_path)
        COILSDIR = 'export COILSDIR=%s\n' % coils_bin
        write_coilsdir(COILSDIR)
        os.chdir(coils_path)
        if os.path.exists('COILS2'):
            os.remove('COILS2')
        os.symlink('coils/ncoils-osf', './COILS2')
        home = str(Path.home())
        if platform == 'darwin':
            os.system('. %s/.bash_profile' % home)
        else:
            os.system('. %s/.bashrc' % home)
        os.chdir(anno_path)
        print('----------------------------------')
        # remove temp files
        shutil.rmtree(anno_path + '/annotation_FAS')
        if not options['keep']:
            os.remove(anno_path + '/' + file)
        os.chdir(current_dir)
        return anno_path
    else:
        os.chdir(current_dir)
        return check_status(anno_path, force, file)


def checkExecutable(anno_path):
    with open(anno_path+'/annoTools.txt') as file:
        availTool = [line.strip() for line in file]
    sedCmd ='perl -pi -e \'s/#checked\\n//\' %s/annoTools.txt' % anno_path
    subprocess.call([sedCmd], shell=True)

    print('Checking if annotation tools are excutable...')
    # test pfam and smart
    if not os.path.isfile(anno_path + '/Pfam/Pfam-hmms/Pfam-A.hmm'):
        sys.exit('Pfam hmm file not found. Please run setupFAS with --force!')
    if not os.path.isfile(anno_path + '/SMART/SMART-hmms/SMART.hmm'):
        sys.exit('SMART hmm file not found. Please run setupFAS with --force!')
    # test seg
    if 'SEG' in availTool:
        segCmd = '%s/SEG/seg' % anno_path
        try:
            p1 = subprocess.Popen([segCmd], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output1, err1 = p1.communicate()
            if not 'Usage' in err1.decode('UTF-8').strip():
                sys.exit('Error with SEG. You can reinstall it by running setupFAS with --force!')
        except:
            sys.exit('Error with SEG. You can reinstall it by running setupFAS with --force!')
    # test fLPS
    if 'fLPS' in availTool:
        flpsCmd = '%s/fLPS/fLPS' % anno_path
        try:
            p2 = subprocess.Popen([flpsCmd], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output2, err2 = p2.communicate()
            if not err2.decode('UTF-8').strip() == 'There is no sequence file. Please supply one.':
                sys.exit('Error with fLPS. You can reinstall it by running setupFAS with --force!')
        except:
            sys.exit('Error with fLPS. You can reinstall it by running setupFAS with --force!')
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
                home = str(Path.home())
                bash_file = '.bashrc'
                if platform == 'darwin':
                    bash_file = '.bash_profile'
                with open(home + '/' + bash_file) as file:
                    if not 'export COILSDIR' in file.read():
                        print("PLEASE PUT THE FOLLOWING LINE TO %s/%s:\nexport COILSDIR=%s" %
                              (home, bash_file, anno_path + '/COILS2/coils'))
                print('NOTE: THE TERMINAL MUST BE RESTARTED BEFORE USING FAS!!!')
                flag = 0
            if '0 sequences' in err3.decode('UTF-8').strip():
                flag = 0
            if flag == 1:
                sys.exit('Error with COILS2. Please restart the terminal and run setupFAS again with --checkExecutable!')
        except:
            sys.exit('Error with COILS2. Please check https://github.com/BIONF/FAS/wiki/FAQ#Error-with-COILS2!')
    # test tmhmm
    if 'TMHMM' in availTool:
        tmhmmCmd = '%s/TMHMM/decodeanhmm' % anno_path
        try:
            p4 = subprocess.Popen([tmhmmCmd], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output4, err4 = p4.communicate()
            if not 'Error: No modelfile given' in err4.decode('UTF-8').strip():
                sys.exit('Error with TMHMM. You can reinstall it by running setupFAS with --force!')
        except:
            sys.exit('Error with TMHMM. You can reinstall it by running setupFAS with --force!')
    # test signalp
    if 'SignalP' in availTool:
        signalpCmd = '%s/SignalP/signalp' % anno_path
        try:
            p5 = subprocess.Popen([signalpCmd, '-V'], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
            output5, err5 = p5.communicate()
            if not err5.decode('UTF-8').strip() == '':
                sys.exit('Error with SignalP. Please check https://github.com/BIONF/FAS/wiki/FAQ#Error-with-SignalP!')
        except:
            sys.exit('Error with SignalP. Please check https://github.com/BIONF/FAS/wiki/FAQ#Error-with-SignalP!')
    return True


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
            if '#checked' not in f.read():
                with open(anno_path+'/annoTools.txt', 'a') as file:
                    file.write('#checked')
                    file.close()
        f.close()
        with open(greedyFasPath+'/pathconfig.txt', 'w') as config:
            config.write(os.path.abspath(anno_path))
            config.close()
        print('Done! Annotation tools can be found in %s' % anno_path)
        print('You should test annoFAS with this command:')
        print('annoFAS -i test_annofas.fa -o testFas_output')
        sys.exit()
    else:
        sys.exit('Some errors occur with annotation tools. Please check if they can be excuted at %s' % anno_path)


def main():
    version = '1.5.3'
    parser = argparse.ArgumentParser(description='You are running setupFAS version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-t', '--toolPath', help='Set path to save annotation tools', action='store', default='',
                          required=True)
    optional.add_argument('-d', '--dtuPath', help='Set path to DTU tools (SignalP and TMHMM)', action='store',
                          default='')
    optional.add_argument('-f', '--force', help='Overwrite old annotation tools if exist', action='store_true')
    optional.add_argument('-k', '--keep', help='Keep downloaded source file', action='store_true')
    optional.add_argument('-s', '--savePath', help='Save annotation tool path to config file for FAS',
                          action='store_true')
    optional.add_argument('-c', '--check', help='Check if FAS ready to run. NO real tool path need to be given!',
                          action='store_true')
    optional.add_argument('--checkExecutable', help='Check if annotation tools are executable!', action='store_true')

    args = parser.parse_args()
    greedyFasPath = os.path.realpath(__file__).replace('/setupFAS.py','')
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
                print('Annotation tools can be found at %s. FAS is ready to run!' % savedPath)
                print('You should test annoFAS with this command:')
                print('annoFAS -i test_annofas.fa -o testFas_output')
                sys.exit()

    if args.savePath:
        checkAnnoToolsFile(args.toolPath)
        with open(greedyFasPath+'/pathconfig.txt', 'w') as config:
            config.write(os.path.abspath(args.toolPath))
            config.close()
        print('Annotation tools can be found at %s. FAS is ready to run!' % args.toolPath)
        print('You should test annoFAS with this command:')
        print('annoFAS -i test_annofas.fa -o testFas_output')
        sys.exit()

    anno_path = prepare_annoTool(options)
    allRun = checkExecutable(anno_path)
    saveConfigFile(allRun, anno_path, greedyFasPath)


if __name__ == '__main__':
    main()
