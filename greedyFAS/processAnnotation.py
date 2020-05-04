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


from sys import version_info
import json
import argparse

if version_info.major == 3:
    from greedyFAS.fasInput import xmlreader
    from greedyFAS.fasWeighting import w_count_ref
    from greedyFAS.fasWeighting import w_weight_correction
elif version_info.major == 2:
    from fasInput import xmlreader
    from fasWeighting import w_count_ref
    from fasWeighting import w_weight_correction


def main():
    parser = argparse.ArgumentParser(description="InterPro table output parser to FAS XML input")
    parser.add_argument("-i", "--input", type=str, default=None, required=True)
    parser.add_argument("-o", "--output", type=str, default=None, required=True)
    args = parser.parse_args()
    databases = [["pfam", "smart"], ["flps", "coils", "seg", "signalp", "tmhmm"]]
    outdict = manage_xml(databases, args.input, "loge")
    dump_data(outdict, args.output)


def manage_xml(databases, path, weight_correction):
    option = {'seed_id': None, 'query_id': None, 'efilter': 0.001, 'inst_efilter': 0.01}
    for ftype in databases[0]:
        proteome, protein_lengths, clan_dict = xmlreader(path + "/" + ftype + ".xml", 2, ftype, True,
                                                         proteome, protein_lengths, clan_dict, option)
    for ftype in databases[1]:
        proteome, protein_lengths, clan_dict = xmlreader(path + "/" + ftype + ".xml", 2, ftype, False,
                                                         proteome, protein_lengths, clan_dict, option)
    prot_count, domain_count = w_count_ref(proteome)
    if weight_correction:
        domain_count = w_weight_correction(weight_correction, domain_count)
    out_dict = {'proteome': proteome, 'protein_lengths': protein_lengths, 'clan_dict': clan_dict,
                'domain_count': domain_count, 'databases': databases}
    return out_dict


def dump_data(out_dict, path):
    with open(path, 'w') as out:
        out.write(json.dumps(out_dict, indent=4, sort_keys=True))


def get_data(path):
    in_dict = json.loads(path)
    return in_dict


if __name__ == '__main__':
    main()
