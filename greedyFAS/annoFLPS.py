#!/bin/env python

#######################################################################
# Copyright (C) 2020 Julian Dosch
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


import subprocess
from sys import argv


def main(fastapath, outpath, flpspath, threshold):
    flps_out = subprocess.run([flpspath + ' ' + fastapath + ' -s -t ' + threshold], shell=True, capture_output=True)
    print(flps_out.stdout.decode().split('\n'))
    ### to be added: parse stdout, somehow get all prot ids from sequence files stop proteins from disappearing


if __name__ == '__main__':
    main(argv[1], argv[2], argv[3], argv[4])
