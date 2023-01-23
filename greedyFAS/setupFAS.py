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
import gnureadline
import glob
from os.path import expanduser
import ssl
import urllib.request
import time
from greedyFAS.disorderFAS import install_aucpred
from pkg_resources import get_distribution

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
            sys.exit('ERROR: Downloaded file corrupted!')
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


def check_conda_env():
    """ Return if a conda env is currently using """
    if 'CONDA_DEFAULT_ENV' in os.environ:
        if not os.environ['CONDA_DEFAULT_ENV'] == 'base':
            return(True)
    return(False)


def get_dtu_path(dtuPathIn):
    if not dtuPathIn == '':
        if not os.path.isdir(dtuPathIn):
            sys.exit('ERROR: ' + dtuPathIn + ' not found!')
        else:
            dtu_path = os.path.abspath(dtuPathIn)
    else:
        print('Do you want to use TMHMM (v2.0c) and SignalP (v4.1g)? (y/n)')
        if query_yes_no('Source path?'):
            gnureadline.set_completer_delims(' \t\n;')
            gnureadline.parse_and_bind('tab: complete')
            gnureadline.set_completer(complete)
            dtu_path = input('Please download TMHMM 2.0c and SignalP 4.1g at https://services.healthtech.dtu.dk/ and '
                             'enter path to the downloaded tar files:')
        else:
            dtu_path = ''

    if not dtu_path == '':
        files = [f for f in os.listdir(dtu_path) if os.path.isfile(os.path.join(dtu_path, f))]
        signalp_source = [i for i in files if 'signalp' in i]
        tmhmm_source = [i for i in files if 'tmhmm' in i]
        if not tmhmm_source:
            sys.exit('ERROR: TMHMM not found in %s' % dtu_path)
        if not signalp_source:
            sys.exit('ERROR: signalp-4.1g not found in %s' % dtu_path)
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
                            sys.exit('ERROR: Failed to delete %s. Please manually remove it and run fas.setup again!'
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


def install_tmhmm():
    shutil.move('tmhmm-2.0c', 'TMHMM')
    machine = os.uname()[4]
    os.chdir('TMHMM')
    makelink_tmhmm = 'ln -s -f bin/decodeanhmm.Linux_%s decodeanhmm' % machine
    subprocess.call([makelink_tmhmm], shell=True)


def install_smart(smart_path, anno_path):
    check = glob.glob('%s/hmm/*.HMM' % smart_path)
    if len(check) < 1:
        return(0)
    else:
        Path(anno_path + '/SMART/SMART-hmms').mkdir(parents=True, exist_ok=True)
        catCmd = 'cat %s/hmm/*.HMM > %s/SMART/SMART-hmms/SMART.hmm' % (smart_path, anno_path)
        subprocess.call([catCmd], shell=True)
        hmmpressCmd = 'hmmpress -f %s/SMART/SMART-hmms/SMART.hmm' % (anno_path)
        try:
            subprocess.run([hmmpressCmd], shell=True, check=True)
        except:
            print('ERROR: Problem occurred while creating binary files for SMART at %s/SMART/SMART-hmms' % (anno_path))
        getVersionCmd = 'tail -n 1 %s/README | cut -d " " -f2 | sed "s/\//_/g" > %s/SMART/version.txt' % (smart_path, anno_path)
        subprocess.call([getVersionCmd], shell=True)
        return(1)

def install_pfam(pfam_version, anno_path):
    url = 'https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam%s' % pfam_version
    file = 'Pfam-A.hmm.gz'
    dat_file = 'Pfam-A.hmm.dat.gz'
    relnotes_file = 'relnotes.txt'
    download_file(url, file)
    download_file(url, dat_file)
    download_file(url, relnotes_file)
    if os.path.isfile(file):
        Path(anno_path + '/Pfam/Pfam-hmms').mkdir(parents=True, exist_ok=True)
        # move donwloaded file to anno_path/Pfam/Pfam-hmms
        try:
            shutil.move(file, '%s/Pfam/Pfam-hmms/%s' % (anno_path, file))
            shutil.move(dat_file, '%s/Pfam/Pfam-hmms/%s' % (anno_path, dat_file))
            shutil.move(relnotes_file, '%s/Pfam/Pfam-hmms/%s' % (anno_path, relnotes_file))
        except:
            sys.exit('ERROR: Cannot move %s, %s & %s to %s/Pfam/Pfam-hmms' % (file, dat_file, relnotes_file, anno_path))
        # unzip downloaded files
        unzipCmd = 'gzip -d %s/Pfam/Pfam-hmms/*.gz' % anno_path
        try:
            subprocess.run([unzipCmd], shell=True, check=True)
        except:
            sys.exit('ERROR: Cannot upzip gz files in %s/Pfam/Pfam-hmms' % anno_path)
        # create bin files for Pfam-A.hmm
        hmmpressCmd = 'hmmpress -f %s/Pfam/Pfam-hmms/Pfam-A.hmm' % (anno_path)
        try:
            subprocess.run([hmmpressCmd], shell=True, check=True)
        except:
            sys.exit('ERROR: Problem occurred while creating binary files for PFAM at %s/Pfam/Pfam-hmms' % (anno_path))
    else:
        sys.exit('ERROR: No Pfam-A.hmm.gz found at %s' % url)


def check_hmmer():
    try:
        subprocess.check_output(['which hmmsearch'], shell = True, stderr = subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        if check_conda_env() == True:
            conda_install_cmd = 'conda install -c bioconda hmmer -y'
            try:
                subprocess.call([conda_install_cmd], shell = True)
            except:
                sys.exit('\033[91mERROR: Cannot install hmmer using this command\n%s\033[0m' % conda_install_cmd)
        else:
            sys.exit('\033[91mERROR: Please install HMMER before using FAS (http://hmmer.org/documentation.html)\033[0m')


def install_annoTool(options):
    anno_path = options['toolPath']
    dtuPathIn = options['dtuPath']
    smart_path = options['smartPath']
    force = options['force']
    ignoreList = options['ignore']
    pfam_version = options['pfamVersion']
    check_hmmer()

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
        defaultTools = ['fLPS', 'Pfam', 'COILS2', 'SEG']
        if platform == 'darwin':
            ignoreList.extend(('SEG', 'TMHMM'))
        tools = [tool for tool in defaultTools if tool not in ignoreList]

        # get DTU tools path
        if not ('TMHMM' in ignoreList and 'SignalP' in ignoreList):
            (dtu_path, signalp_source, tmhmm_source) = get_dtu_path(dtuPathIn)
        else:
            dtu_path = 0

        # install TMHMM and SignalP
        if not dtu_path == 0:
            if not 'SignalP' in ignoreList:
                print('Installing SignalP...')
                if not os.path.isdir(anno_path + '/SignalP'):
                    shutil.unpack_archive(dtu_path + '/' + signalp_source, anno_path, 'gztar')
                    install_signalp()
                os.chdir(anno_path)
                tools.append('SignalP')
            if not 'TMHMM' in ignoreList:
                if not platform == 'darwin':
                    print('Installing TMHMM...')
                    if not os.path.isdir(anno_path + '/TMHMM'):
                        shutil.unpack_archive(dtu_path + '/' + tmhmm_source, anno_path, 'gztar')
                        install_tmhmm()
                    os.chdir(anno_path)
                    tools.append('TMHMM')

        # install SMART
        smart = 0
        if not smart_path == '':
            print('Installing SMART...')
            smart = install_smart(smart_path, anno_path)
            if smart == 0:
                print('ERROR: Failed to install SMART from %s' % smart_path)

        # close annoTools file
        with open('annoTools.txt', mode='wt') as tool_file:
            tool_file.write('#linearized\nPfam\n#normal\n')
            if (smart == 1):
                tool_file.write('#linearized\nPfam\nSMART\n#normal\n')
            for t in tools:
                if not t in ['Pfam','SMART']:
                    tool_file.write('%s\n' % t)
        tool_file.close()

        # download other annotation tools
        print('Installing other tools...')
        if options['disorder']:
            install_aucpred.install_auc(options['toolPath'])

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
                fLPS_file = 'fLPS2programs.tar' #'fLPS.tar.gz'
                fLPS_url = 'http://biology.mcgill.ca/faculty/harrison/'
                download_file(fLPS_url, fLPS_file)
                shutil.unpack_archive(fLPS_file, anno_path, 'gztar')
                os.remove(fLPS_file)
                shutil.move('fLPS2programs', 'fLPS')
            else:
                if not tool in ('SignalP', 'TMHMM'):
                    print('Moving %s ...' % tool)
                    source_dir = 'annotation_FAS/' + tool + '/'
                    if os.path.exists(os.path.abspath(anno_path + '/' + tool)):
                        shutil.rmtree(anno_path + '/' + tool)
                    shutil.move(source_dir, anno_path, copy_function=shutil.copytree)

        # replace PFAM by user-defined version
        if not pfam_version == '':
            print('Downloading Pfam version %s...' % pfam_version)
            shutil.rmtree(anno_path + '/Pfam')
            install_pfam(pfam_version, anno_path)

        # compile and make symlink for fLPS
        flps_path = anno_path + '/fLPS'
        os.chdir(flps_path + '/src')
        makeCmd = 'make'
        try:
            subprocess.call([makeCmd], shell=True)
            source = flps_path + '/src/fLPS2'
            target = flps_path + '/fLPS'
            try:
                os.symlink(source, target)
            except FileExistsError:
                os.remove(target)
                os.symlink(source, target)
        except:
            print('ERROR: Failed to compile fLPS.\nPlease try to do it manually at' % flps_path)

        if not 'COILS2' in ignoreList:
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
            os.chdir(coils_path)
            if os.path.exists('COILS2'):
                os.remove('COILS2')
            os.symlink('coils/ncoils-osf', './COILS2')
            os.chdir(anno_path)
            COILSDIR = 'export COILSDIR=%s\n' % coils_bin
            with open(anno_path + '/fas.profile', 'w') as outfile:
                outfile.write(COILSDIR)
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
        sys.exit('Pfam hmm file not found. Please run fas.setup with --force!')
    if 'SMART' in availTool:
        if not os.path.isfile(anno_path + '/SMART/SMART-hmms/SMART.hmm'):
            sys.exit('SMART hmm file not found. Please run fas.setup with --force!')
    # test seg
    if 'SEG' in availTool:
        segCmd = '%s/SEG/seg' % anno_path
        try:
            p1 = subprocess.Popen([segCmd], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output1, err1 = p1.communicate()
            if not 'Usage' in err1.decode('UTF-8').strip():
                sys.exit('Error with SEG. You can reinstall it by running fas.setup with --force!')
        except:
            sys.exit('Error with SEG. You can reinstall it by running fas.setup with --force!')
    # test fLPS
    if 'fLPS' in availTool:
        flpsCmd = '%s/fLPS/fLPS' % anno_path
        try:
            p2 = subprocess.Popen([flpsCmd], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output2, err2 = p2.communicate()
            if not err2.decode('UTF-8').strip() == 'There is no sequence file. Please supply one in FASTA format.': # 'There is no sequence file. Please supply one.':
                sys.exit('Error with fLPS. You can reinstall it by running fas.setup with --force!')
        except:
            sys.exit('Error with fLPS. You can reinstall it by running fas.setup with --force!')
    # test COILS2
    if 'COILS2' in availTool:
        profileFile = '%s/fas.profile' % anno_path
        if not os.path.isfile(profileFile):
            sys.exit('No config file for COILSDIR found. Please check https://github.com/BIONF/FAS/wiki/FAQ#Error-with-COILS2!')
    # test tmhmm
    if 'TMHMM' in availTool:
        tmhmmCmd = '%s/TMHMM/decodeanhmm' % anno_path
        try:
            p4 = subprocess.Popen([tmhmmCmd], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output4, err4 = p4.communicate()
            if not 'Error: No modelfile given' in err4.decode('UTF-8').strip():
                sys.exit('Error with TMHMM. You can reinstall it by running fas.setup with --force!')
        except:
            sys.exit('Error with TMHMM. You can reinstall it by running fas.setup with --force!')
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
        print('You should test fas.doAnno with this command:')
        print('fas.doAnno -i test_annofas.fa -o testFas_output')
        print('NOTE: YOU NEED TO source %s/fas.profile BEFORE USING FAS!' % anno_path)
        print('Check https://github.com/BIONF/FAS/wiki/setup#add-COILSDIR-to-bashrc-or-bash_profile for more details')
        sys.exit()
    else:
        sys.exit('Some errors occur with annotation tools. Please check if they can be excuted at %s' % anno_path)


def main():
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-t', '--toolPath', help='Set path to save annotation tools', action='store', default='',
                          required=True)
    optional.add_argument('--smartPath', help='Set path to your downloaded SMART folder', action='store', default='')
    optional.add_argument('-d', '--dtuPath', help='Set path to DTU tools (SignalP and TMHMM)', action='store',
                          default='')
    optional.add_argument('--pfamVersion', help='Specify your own PFAM version. E.g. 28.0', action='store', default='')
    optional.add_argument('-i', '--ignore', help='List of tools should be ignored', nargs="*",
                        choices=['SignalP', 'TMHMM', 'COILS2', 'fLPS', 'SEG'],
                        action='store', default=[])
    optional.add_argument('-f', '--force', help='Overwrite old annotation tools if exist', action='store_true')
    optional.add_argument('-k', '--keep', help='Keep downloaded source file', action='store_true')
    optional.add_argument('-s', '--savePath', help='Save annotation tool path to config file for FAS',
                          action='store_true')
    optional.add_argument('-c', '--check', help='Check if FAS ready to run. NO real tool path need to be given!',
                          action='store_true')
    optional.add_argument('--checkExecutable', help='Check if annotation tools are executable!', action='store_true')
    optional.add_argument('--disorder', help='Install Disorder annotation, can be installed separately later as well',
                          action='store_true')

    args = parser.parse_args()
    greedyFasPath = os.path.realpath(__file__).replace('/setupFAS.py','')
    options = {
        'toolPath': args.toolPath,
        'smartPath': args.smartPath,
        'dtuPath': args.dtuPath,
        'pfamVersion': args.pfamVersion,
        'ignore': args.ignore,
        'force': args.force,
        'keep': args.keep,
        'greedyFasPath': greedyFasPath,
        'disorder': args.disorder
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
                print('You should test fas.doAnno with this command:')
                print('fas.doAnno -i test_annofas.fa -o testFas_output')
                sys.exit()

    if args.savePath:
        checkAnnoToolsFile(args.toolPath)
        with open(greedyFasPath+'/pathconfig.txt', 'w') as config:
            config.write(os.path.abspath(args.toolPath))
            config.close()
        print('Annotation tools can be found at %s. FAS is ready to run!' % args.toolPath)
        print('You should test fas.doAnno with this command:')
        print('fas.doAnno -i test_annofas.fa -o testFas_output')
        print('NOTE: For using FAS you need to source %s/fas.profile first!' % args.toolPath)
        sys.exit()

    anno_path = install_annoTool(options)
    allRun = checkExecutable(anno_path)
    saveConfigFile(allRun, anno_path, greedyFasPath)


if __name__ == '__main__':
    main()
