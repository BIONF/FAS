#!/bin/env python

#######################################################################
# Copyright (C) 2022 Julian Dosch
#
# This file is part of greedyFAS.
#
#  greedyFAS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  greedyFAS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with greedyFAS.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################


import os
import sys
import subprocess
import argparse
from os.path import expanduser
from git import Repo
from git import RemoteProgress
from tqdm import tqdm
from pkg_resources import get_distribution


home = expanduser('~')


class CloneProgress(RemoteProgress):
    def __init__(self):
        super().__init__()
        self.pbar = tqdm()

    def update(self, op_code, cur_count, max_count=None, message=''):
        self.pbar.total = max_count
        self.pbar.n = cur_count
        self.pbar.refresh()


def compile_AUCpreD(path):
    auc_src = path + '/source_code'
    os.chdir(auc_src)
    compile_cmd = 'make' + ' > /dev/null 2>&1'
    try:
        subprocess.call(compile_cmd, shell=True)
    except:
        raise Exception('Failed to compile PredictProperty.'
                        '\nPlease read instruction at %s and do it manually!' % auc_src)


def install_auc(path):
    print('Cloning https://github.com/realbigws/Predict_Property.git...')
    try:
        Repo.clone_from('https://github.com/realbigws/Predict_Property.git', path.rstrip('/') +
                        '/Predict_Property/', branch='master', progress=CloneProgress())
    except:
        raise Exception('Could not clone github repository.'
                        '\nMake sure that you have a working internet connection and check if the install directory '
                        'does not already have a Predict_Property installation.')
    print('done!\ncompiling PredictProperty...')
    compile_AUCpreD(path.rstrip('/') + '/Predict_Property/')
    print('Finished installing PredictProperty.')


def main():
    version = get_distribution('greedyFAS').version
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    optional = parser.add_argument_group('required arguments')
    optional.add_argument("-p", "--path", default=None, type=str, required=False,
                          help="install path for PredictProperty, per default it will be installed in the same "
                               "directory as the other annotation tools")
    args = parser.parse_args()
    if args.path:
        path = args.path
    else:
        pathconfigFile = os.path.realpath(__file__).replace('disorderFAS/install_aucpred.py', 'pathconfig.txt')
        print(pathconfigFile)
        if not os.path.exists(pathconfigFile):
            sys.exit('No pathconfig.txt found. Please run fas.setup (https://github.com/BIONF/FAS/wiki/setup).')
        with open(pathconfigFile, 'r') as input:
            path = input.readline().strip()
    install_auc(path)


if __name__ == '__main__':
    main()
