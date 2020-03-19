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


def main(fastapath, outpath, flpspath, threshold, prot_lengths):
    flps_out = subprocess.run([flpspath + ' ' + fastapath + ' -s -t ' + threshold], shell=True, capture_output=True)
    lines = flps_out.stdout.decode().split('\n')
    proteome = {}
    for i in prot_lengths:
        proteome[i] = {}
    for line in lines:
        cells = line.split('\t')
        feature = cells[1] + '_' + cells[7]
        if feature in proteome[cells[0]]:
            proteome[cells[0]][feature].append((cells[3], cells[4]))
        else:
            proteome[cells[0]][feature] = [(cells[3], cells[4])]
    write_xml(outpath, proteome, prot_lengths)
    ### to be added: somehow get all prot ids from sequence files to stop proteins from disappearing (should be in the main anno script)


def write_xml(outpath, proteome, prot_lengths):
    with open(outpath, 'w') as out:
        out.write('<?xml version="1.0"?>\n<tool name="fLPS">\n')
        for protein in proteome:
            out.write('\t<protein id="' + protein + '" length="' + str(prot_lengths[protein]) + '">\n')
            for feature in proteome[protein]:
                out.write('\t\t<feature type="' + feature + '" instance="' + str(len(proteome[protein][feature])) +
                          '">\n')
                for instance in proteome[protein][feature]:
                    out.write('\t\t\t<start start="' + instance[0] + '">\n\t\t\t<end end="' + instance[1] + '">\n')
                out.write('\t\t</feature>\n')
            out.write('\t</protein>\n')
        out.write('</tool>')


if __name__ == '__main__':
    main(argv[1], argv[2], argv[3], argv[4], argv[5])
