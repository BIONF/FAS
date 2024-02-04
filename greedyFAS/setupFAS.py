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
import multiprocessing as mp
from tqdm import tqdm
from greedyFAS.disorderFAS import install_aucpred
import greedyFAS.annoFAS.annoModules as annoModules
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


def check_status(toolPath, force, reinstall):
    flag = 1
    cwd = os.getcwd()
    if os.path.isfile(toolPath+'/annoTools.txt'):
        with open(toolPath+'/annoTools.txt') as f:
            if '#checked' in f.read():
                if force:
                    print('Annotation tools found in %s will be deleted and reinstalled! Enter to continue.' % toolPath)
                    if query_yes_no(''):
                        try:
                            shutil.rmtree(toolPath)
                        except:
                            sys.exit('ERROR: Failed to delete %s. Please manually remove it and run fas.setup again!'
                                     % toolPath)
                        flag = 1
                    else:
                        flag = toolPath
                else:
                    if len(reinstall) > 0:
                        flag = 1
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
            if 'uname -m' in line:
                if 'Debian' in os.uname()[3]:
                    line = 'my $architecture = \"i686\";\n'
            outfile.write(line)
    shutil.move(signalp_path+'/signalp.mod', signalp_path+'/signalp')
    subprocess.call(['chmod', '0755', signalp_path+'/signalp'])


def install_tmhmm():
    shutil.move('tmhmm-2.0c', 'TMHMM')
    machine = os.uname()[4]
    if 'Debian' in os.uname()[3]:
        machine = 'i686'
    os.chdir('TMHMM')
    makelink_tmhmm = 'ln -s -f bin/decodeanhmm.Linux_%s decodeanhmm' % machine
    subprocess.call([makelink_tmhmm], shell=True)


def correct_aln(args):
    (aln_file, anno_path) = args
    tmp = {}
    file_name = aln_file.split('/')[-1].split('.')[0]
    Path(anno_path + '/SMART/aln').mkdir(parents = True, exist_ok = True)
    with open(aln_file, 'r') as f:
        for l in [line.rstrip() for line in f]:
            if not l.startswith('CLUSTAL') and not l.startswith(' ') and len(l) > 0:
                if not l.startswith('>'):
                    id = l.split()[0].replace('.', '_')
                    if len(l.split()) > 1:
                        if not id in tmp:
                            tmp[id] = l.split()[1].strip()
                        else:
                            tmp[id] += l.split()[1].strip()
                else:
                    break
    if not len(tmp) > 0:
        if not os.path.exists(f'{anno_path}/SMART/aln/{file_name}.aln'):
            os.symlink(aln_file, f'{anno_path}/SMART/aln/{file_name}.aln')
    else:
        with open(f'{anno_path}/SMART/aln/{file_name}.aln', 'w') as o:
            for i in tmp:
                o.write(f'>{i}\n{tmp[i]}\n')


def create_phmm(args):
    (aln_file, anno_path) = args
    file_name = aln_file.split('/')[-1].split('.')[0]
    Path(anno_path + '/SMART/hmm').mkdir(parents = True, exist_ok = True)
    if os.path.exists(f'{anno_path}/SMART/aln/{file_name}.aln'):
        try:
            hmmbuildCmd = f'hmmbuild {anno_path}/SMART/hmm/{file_name}.hmm {anno_path}/SMART/aln/{file_name}.aln > /dev/null 2>&1'
            subprocess.run([hmmbuildCmd], shell=True, check=True)
        except:
            print(f'ERROR: Problem occurred while creating HMM profile for {aln_file}!')


def install_smart(smart_path, anno_path):
    aln_files = glob.glob(f'{smart_path}/aln/*.aln')
    if len(aln_files) < 1:
        sys.exit(f'ERROR: Cannot find aln files in {smart_path}')
    if not os.path.exists('%s/SMART/SMART-hmms/SMART.hmm.length' % anno_path):
        # create profile hmms from aln files
        aln_jobs = []
        hmm_jobs = []
        for aln in aln_files:
            aln_jobs.append([aln, anno_path])
            hmm_jobs.append([aln, anno_path])
        cpus = mp.cpu_count()-1
        pool = mp.Pool(cpus)
        for _ in tqdm(pool.imap_unordered(correct_aln, aln_jobs), total = len(aln_jobs)):
            pass
        for _ in tqdm(pool.imap_unordered(create_phmm, hmm_jobs), total = len(hmm_jobs)):
            pass
        pool.close()
        # create hmm database for hmmscan
        Path(anno_path + '/SMART/SMART-hmms').mkdir(parents=True, exist_ok=True)
        catCmd = f'cat {anno_path}/SMART/hmm/*.hmm > {anno_path}/SMART/SMART-hmms/SMART.hmm'
        subprocess.call([catCmd], shell=True)
        hmmpressCmd = 'hmmpress -f %s/SMART/SMART-hmms/SMART.hmm' % (anno_path)
        try:
            subprocess.run([hmmpressCmd], shell=True, check=True)
        except:
            sys.exit('ERROR: Problem occurred while creating binary files for SMART at %s/SMART/SMART-hmms' % (anno_path))
        getVersionCmd = 'tail -n 1 %s/README | cut -d " " -f2 | sed "s/\//_/g" > %s/SMART/version.txt' % (smart_path, anno_path)
        subprocess.call([getVersionCmd], shell=True)
        # write length file
        write_phmm_length(anno_path, 'SMART')
        # remove temp folders (aln and hmm)
        shutil.rmtree(f'{anno_path}/SMART/aln')
        shutil.rmtree(f'{anno_path}/SMART/hmm')
    else:
        print('WARNING: SMART is already installed!')


def install_pfam(pfam_version, pfam_path, anno_path):
    if not os.path.exists('%s/Pfam/Pfam-hmms/Pfam-A.hmm.length' % anno_path):
        file = 'Pfam-A.hmm.gz'
        dat_file = 'Pfam-A.hmm.dat.gz'
        relnotes_file = 'relnotes.txt'
        if not pfam_path:
            url = 'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/'
            if pfam_version:
                url = 'https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam%s' % pfam_version
            download_file(url, file)
            download_file(url, dat_file)
            download_file(url, relnotes_file)
        else:
            annoModules.checkFileExist(f'{pfam_path}/{file}')
            os.symlink(f'{pfam_path}/{file}', f'{anno_path}/{file}')
            annoModules.checkFileExist(f'{pfam_path}/{dat_file}')
            os.symlink(f'{pfam_path}/{dat_file}', f'{anno_path}/{dat_file}')
            annoModules.checkFileExist(f'{pfam_path}/{relnotes_file}')
            os.symlink(f'{pfam_path}/{relnotes_file}', f'{anno_path}/{relnotes_file}')
        if os.path.isfile(file):
            Path(anno_path + '/Pfam/Pfam-hmms').mkdir(parents=True, exist_ok=True)
            # move donwloaded file to anno_path/Pfam/Pfam-hmms
            try:
                if not pfam_path:
                    shutil.move(file, '%s/Pfam/Pfam-hmms/%s' % (anno_path, file))
                    shutil.move(dat_file, '%s/Pfam/Pfam-hmms/%s' % (anno_path, dat_file))
                    shutil.move(relnotes_file, '%s/Pfam/Pfam-hmms/%s' % (anno_path, relnotes_file))
                else:
                    shutil.copy(f'{pfam_path}/{file}', f'{anno_path}/Pfam/Pfam-hmms/{file}')
                    shutil.copy(f'{pfam_path}/{dat_file}', f'{anno_path}/Pfam/Pfam-hmms/{dat_file}')
                    shutil.copy(f'{pfam_path}/{relnotes_file}', f'{anno_path}/Pfam/Pfam-hmms/{relnotes_file}')
            except:
                sys.exit('ERROR: Cannot move/copy %s, %s & %s to %s/Pfam/Pfam-hmms' % (file, dat_file, relnotes_file, anno_path))
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
            # write length file for PFAM
            write_phmm_length(anno_path, 'Pfam')
            # remove symlinks
            if os.path.isfile(f'{anno_path}/{file}'):
                os.remove(f'{anno_path}/{file}')
                os.remove(f'{anno_path}/{dat_file}')
                os.remove(f'{anno_path}/{relnotes_file}')
        else:
            sys.exit('ERROR: No Pfam-A.hmm.gz found at %s' % url)
    else:
        print('WARNING: PFAM is already installed!')


def write_phmm_length(anno_path, tool_name):
    hmm_file_name = f'{tool_name}-A.hmm'
    if tool_name == 'SMART':
        hmm_file_name = f'{tool_name}.hmm'
    hmm_file = f'{anno_path}/{tool_name}/{tool_name}-hmms/{hmm_file_name}'
    with open(f'{anno_path}/{tool_name}/{tool_name}-hmms/{hmm_file_name}.length', 'w') as len_f:
        with open(hmm_file, 'r') as file:
            for l in [line.rstrip() for line in file]:
                if l.startswith('NAME'):
                    id = l.split()[1].strip()
                if l.startswith('LENG'):
                    len_f.write(f'{id}\t{l.split()[1].strip()}\n')


def install_coils(anno_path):
    if not os.path.islink('%s/COILS2/COILS2' % anno_path):
        # url = 'http://ftp.ebi.ac.uk/pub/software/unix/coils-2.2/'
        # file = 'ncoils.tar.gz'
        url = 'http://www.russelllab.org/coils/'
        file = 'coils.tar.gz'
        download_file(url, file)
        if os.path.isfile(file):
            coils_path = anno_path + '/COILS2'
            Path(coils_path).mkdir(parents=True, exist_ok=True)
            # move donwloaded file to anno_path/COILS2
            try:
                shutil.move(file, '%s/%s' % (coils_path, file))
            except:
                sys.exit('ERROR: Cannot move %s to %s' % (file, coils_path))
            # re-compile COILS2
            os.chdir(coils_path)
            shutil.unpack_archive(file, coils_path, 'gztar')
            coils_bin = coils_path + '/coils'
            os.chdir(coils_bin)
            compile_cmd = 'cc -O2 -I. -o ncoils-osf ncoils.c read_matrix.c -lm' + ' > /dev/null 2>&1'
            try:
                subprocess.call(compile_cmd, shell=True)
            except:
                os.chdir(anno_path)
                sys.exit('ERROR: Failed to compile COILS2.\nPlease read instruction at %s and do it manually!' % coils_path)
            os.chdir(coils_path)
            if os.path.exists('COILS2'):
                os.remove('COILS2')
            os.symlink('coils/ncoils-osf', './COILS2')
            os.symlink('coils/new.mat', './new.mat')
            os.symlink('coils/old.mat', './old.mat')
            os.chdir(anno_path)
            COILSDIR = 'export COILSDIR=%s\n' % coils_bin
            with open(anno_path + '/fas.profile', 'w') as outfile:
                outfile.write(COILSDIR)
    else:
        print('WARNING: COILS2 is already installed!')


def install_flps(anno_path):
    if not os.path.islink('%s/fLPS/fLPS' % anno_path):
        # download flps
        fLPS_file = 'fLPS2programs.tar' #'fLPS.tar.gz'
        fLPS_url = 'http://biology.mcgill.ca/faculty/harrison/'
        download_file(fLPS_url, fLPS_file)
        shutil.unpack_archive(fLPS_file, anno_path, 'gztar')
        os.remove(fLPS_file)
        os.chdir(anno_path)
        shutil.move('fLPS2programs', 'fLPS')
        # compile and make symlink for fLPS
        flps_path = anno_path + '/fLPS'
        shutil.rmtree(anno_path + '/fLPS/bin')
        os.chdir(flps_path + '/src')
        makeCmd = 'make'
        os.mkdir(anno_path + '/fLPS/bin')
        try:
            subprocess.call([makeCmd], shell=True)
            for i in ['CompositionMaker', 'DomainFilter', 'fLPS2']:
                shutil.move(flps_path + '/src/' + i, flps_path + '/bin/' + i)
            source = flps_path + '/bin/fLPS2'
            target = flps_path + '/fLPS'
            try:
                os.symlink(source, target)
            except FileExistsError:
                os.remove(target)
                os.symlink(source, target)
        except:
            os.chdir(anno_path)
            sys.exit('ERROR: Failed to compile fLPS.\nPlease try to do it manually at' % flps_path)
        os.chdir(anno_path)
    else:
        print('WARNING: fLPS is already installed!')


def install_seg(anno_path):
    if not os.path.islink('%s/SEG/seg' % anno_path):
        url = 'https://raw.githubusercontent.com/BIONF/data4travis/main'
        file = 'SEG.tar.gz'
        download_file(url, file)
        shutil.unpack_archive(file, anno_path, 'gztar')
        os.remove(file)
    else:
        print('WARNING: SEG is already installed!')


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


def check_installed_tools(anno_path, ignoreList, reinstall, force, greedyFasPath):
        # list of annotation tools
        defaultTools = ['Pfam', 'fLPS', 'COILS2', 'SEG']
        if platform == 'darwin':
            ignoreList.extend(('SEG', 'TMHMM'))
        tools = [tool for tool in defaultTools if tool not in ignoreList]

        # check if anno tools exist
        if os.path.exists(os.path.abspath(anno_path+'/annoTools.txt')):
            with open(anno_path+'/annoTools.txt') as checkfile:
                if '#checked' in checkfile.read():
                    if not force:
                        if len(reinstall) > 0:
                            # get existing tools
                            tools = []
                            with open(anno_path+'/annoTools.txt') as f:
                                file =  f.readlines()
                                for line in file:
                                    if ('#' not in line) and (len(line) > 1):
                                        tools.append(line.strip())
                            # delete installed tools to reinstall
                            for r_tool in reinstall:
                                if os.path.exists('%s/%s' % (anno_path, r_tool)):
                                    print('Deleting %s/%s' % (anno_path, r_tool))
                                    shutil.rmtree('%s/%s' % (anno_path, r_tool))
                                    if not r_tool in tools:
                                        tools.append(r_tool)
                        else:
                            print('Annotation tools already found at %s' % anno_path)
                            if not os.path.exists(os.path.abspath(greedyFasPath + '/pathconfig.txt')):
                                print('If you want to add %s to config file of FAS, please rerun this function with --savePath!'
                                      % anno_path)
                            print('If you want to re-install them, rerun this function with --force!')
                            sys.exit()
        return(tools)


def install_annoTool(options):
    anno_path = options['toolPath']
    dtuPathIn = options['dtuPath']
    smart_path = options['smartPath']
    force = options['force']
    ignoreList = options['ignore']
    reinstall = options['reinstall']
    pfam_version = options['pfamVersion']
    pfam_path = options['pfamPath']
    check_hmmer()


    tools = check_installed_tools(anno_path, ignoreList, reinstall, force, options['greedyFasPath'])

    # do install anno tools
    current_dir = os.getcwd()
    if check_status(anno_path, force, reinstall) == 1:
        Path(anno_path).mkdir(parents=True, exist_ok=True)
        anno_path = os.path.abspath(anno_path)
        os.chdir(anno_path)

        print('Annotation tools will be installed in ')
        print(os.getcwd())
        print('----------------------------------')
        # install TMHMM and SignalP
        if len(reinstall) == 0 or 'SignalP' in reinstall or 'TMHMM' in reinstall:
            if not ('TMHMM' in ignoreList and 'SignalP' in ignoreList):
                (dtu_path, signalp_source, tmhmm_source) = get_dtu_path(dtuPathIn)
            else:
                dtu_path = 0
            if not dtu_path == 0:
                if not 'SignalP' in ignoreList:
                    print('==> Installing SignalP...')
                    if not os.path.isdir(anno_path + '/SignalP'):
                        shutil.unpack_archive(dtu_path + '/' + signalp_source, anno_path, 'gztar')
                        install_signalp()
                    os.chdir(anno_path)
                    if not 'SignalP' in tools:
                        tools.append('SignalP')
                if not 'TMHMM' in ignoreList:
                    if not platform == 'darwin':
                        print('==> Installing TMHMM...')
                        if not os.path.isdir(anno_path + '/TMHMM'):
                            shutil.unpack_archive(dtu_path + '/' + tmhmm_source, anno_path, 'gztar')
                            install_tmhmm()
                        os.chdir(anno_path)
                        if not 'TMHMM' in tools:
                            tools.append('TMHMM')

        # install SMART
        if len(reinstall) == 0 or 'SMART' in reinstall:
            if not smart_path == '':
                print('==> Installing SMART...')
                install_smart(smart_path, anno_path)
                if not 'SMART' in tools:
                    tools.append('SMART')
            else:
                if 'SMART' in reinstall:
                    sys.exit('ERROR: You must provide your downloaded SMART data. Check https://github.com/BIONF/FAS/wiki/setup for more details!')

        # install PFAM
        if len(reinstall) == 0 or 'Pfam' in reinstall:
            if pfam_path:
                print('==> Installing PFAM from %s...' % pfam_path)
            elif pfam_version:
                print('==> Installing PFAM version %s...' % pfam_version)
            else:
                print('==> Installing PFAM latest version...')
            install_pfam(pfam_version, pfam_path, anno_path)

        # install fLPS, COILS2 and SEG
        if len(reinstall) == 0 or 'fLPS' in reinstall:
            if not 'fLPS' in ignoreList:
                print('==> Installing fLPS ...')
                install_flps(anno_path)

        if len(reinstall) == 0 or 'COILS2' in reinstall:
            if 'COILS2' not in ignoreList:
                print('==> Installing COILS2...')
                install_coils(anno_path)

        if len(reinstall) == 0 or 'SEG' in reinstall:
            if 'SEG' not in ignoreList:
                print('==> Installing SEG...')
                install_seg(anno_path)

        # install disorder
        if options['disorder']:
            install_aucpred.install_auc(options['toolPath'])

        # write annoTools file
        with open('annoTools.txt', mode='wt') as tool_file:
            tool_file.write('#linearized\nPfam\n')
            if 'SMART' in tools:
                tool_file.write('SMART\n')
            tool_file.write('#normal\n')
            for t in tools:
                if not t in ['Pfam','SMART']:
                    tool_file.write('%s\n' % t)
        tool_file.close()

        print('----------------------------------')
        os.chdir(current_dir)
        return anno_path
    else:
        os.chdir(current_dir)
        return check_status(anno_path, force, reinstall)


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
        print('==> fas.doAnno -i test_annofas.fa -o testFas_output <==')
        print('*** NOTE: YOU NEED TO source %s/fas.profile BEFORE USING FAS!' % anno_path)
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
    optional.add_argument('--pfamVersion', help='Specify PFAM version. E.g.: 35.0', action='store', default='')
    optional.add_argument('--pfamPath', help='Set path to your downloaded PFAM folder', action='store', default='')
    optional.add_argument('-i', '--ignore', help='List of tools should be ignored, e.g. SEG COILS2', nargs="*",
                        choices=['SignalP', 'TMHMM', 'COILS2', 'fLPS', 'SEG'],
                        action='store', default=[])
    optional.add_argument('--noAnno', help='Run setup without installing annotation tools',
                        action='store_true')
    optional.add_argument('-f', '--force', help='Overwrite old annotation tools if exist', action='store_true')
    optional.add_argument('-r', '--reinstall', help='List of tools should be reinstalled, e.g. SMART SignalP', nargs="*",
                        choices=['SMART','Pfam','SignalP', 'TMHMM', 'COILS2', 'fLPS', 'SEG'],
                        action='store', default=[])
    optional.add_argument('-s', '--savePath', help='Save annotation tool path to config file for FAS',
                          action='store_true')
    optional.add_argument('-c', '--check', help='Check if FAS ready to run. NO real tool path need to be given!',
                          action='store_true')
    optional.add_argument('--checkExecutable', help='Check if annotation tools are executable!', action='store_true')
    optional.add_argument('--addLength', help='Create length files for PFAM and SMART', action='store_true')
    optional.add_argument('--disorder', help='Install Disorder annotation, can be installed separately later as well',
                          action='store_true')

    args = parser.parse_args()
    greedyFasPath = os.path.realpath(__file__).replace('/setupFAS.py','')
    options = {
        'toolPath': args.toolPath,
        'smartPath': args.smartPath,
        'dtuPath': args.dtuPath,
        'pfamVersion': args.pfamVersion,
        'pfamPath': args.pfamPath,
        'ignore': args.ignore,
        'force': args.force,
        'reinstall': args.reinstall,
        'greedyFasPath': greedyFasPath,
        'disorder': args.disorder,
        'noAnno': args.noAnno
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
                print('==> fas.doAnno -i test_annofas.fa -o testFas_output <==')
                sys.exit()

    if args.savePath:
        checkAnnoToolsFile(args.toolPath)
        with open(greedyFasPath+'/pathconfig.txt', 'w') as config:
            config.write(os.path.abspath(args.toolPath))
        print('Annotation tools can be found at %s. FAS is ready to run!' % args.toolPath)
        print('You should test fas.doAnno with this command:')
        print('==> fas.doAnno -i test_annofas.fa -o testFas_output <==')
        print('*** NOTE: For using FAS you need to source %s/fas.profile first!' % args.toolPath)
        sys.exit()

    if args.addLength:
        if not os.path.exists(os.path.abspath(args.toolPath+'/annoTools.txt')):
            sys.exit('ERROR: %s not found' % (args.toolPath+'/annoTools.txt'))
        else:
            if not os.path.exists(os.path.abspath(f'{args.toolPath}/Pfam/Pfam-hmms/Pfam-A.hmm.length')):
                print(f'Creating length file for PFAM...')
                write_phmm_length(args.toolPath, 'Pfam')
            if not os.path.exists(os.path.abspath(f'{args.toolPath}/SMART/SMART-hmms/SMART.hmm.length')):
                print(f'Creating length file for SMART...')
                write_phmm_length(args.toolPath, 'SMART')
            sys.exit()

    if args.noAnno:
        with open(greedyFasPath+'/pathconfig.txt', 'w') as config:
            config.write(os.path.abspath(greedyFasPath))
        with open(greedyFasPath+'/annoTools.txt', 'w') as annoFile:
            annoFile.write("#linearized\nPfam\nSMART\n#normal\nfLPS\nCOILS2\nSEG\nTMHMM\nSignalP\n#checked")
        print(f'*** NOTE: This FAS version works by default only with complete annotations from all tools.')
        print(f'Modify {greedyFasPath}/annoTools.txt if needed!')
        sys.exit()


    anno_path = install_annoTool(options)
    allRun = checkExecutable(anno_path)
    saveConfigFile(allRun, anno_path, greedyFasPath)


if __name__ == '__main__':
    main()
