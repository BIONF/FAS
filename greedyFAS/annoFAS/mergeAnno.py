#!/bin/env python

#######################################################################
# Copyright (C) 2025 Vinh Tran
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


import argparse
import json
import os
from collections import defaultdict
from importlib.metadata import version, PackageNotFoundError


def merge_json_files(input_folder, output_path):

    # Collect all JSON files in the folder
    filenames = sorted(
        f for f in os.listdir(input_folder)
        if f.endswith(".json")
    )

    if not filenames:
        raise RuntimeError(f"No JSON files found in: {input_folder}")

    # Initialize merged structure
    merged = {
        "feature": {},
        "interproID": {},
        "clan": {},
        "count": defaultdict(int),
        "length": {},
        "version": {}
    }

    print(f"Found {len(filenames)} JSON files. Merging...")

    # Process each file
    for fname in filenames:
        fullpath = os.path.join(input_folder, fname)
        print(f" → Reading {fname}")

        with open(fullpath) as f:
            data = json.load(f)

        # Merge feature dictionaries
        for k, v in data.get("feature", {}).items():
            if k in merged["feature"]:
                raise ValueError(f"Duplicate feature ID encountered: {k}")
            merged["feature"][k] = v

        # Merge interproID
        merged["interproID"].update(data.get("interproID", {}))

        # Merge clan
        merged["clan"].update(data.get("clan", {}))

        # Merge length
        merged["length"].update(data.get("length", {}))

        # Merge version – warn on conflicts
        for vkey, vdict in data.get("version", {}).items():
            if vkey not in merged["version"]:
                merged["version"][vkey] = vdict
            else:
                if merged["version"][vkey] != vdict:
                    print(f"WARNING: Version mismatch for '{vkey}' "
                          f"between files. Keeping the first value.")

        # Sum counts
        for k, val in data.get("count", {}).items():
            merged["count"][k] += val

    # Convert defaultdict → dict
    merged["count"] = dict(merged["count"])

    # Write merged output as a single line JSON
    with open(output_path, "w") as out:
        json.dump(merged, out, separators=(",", ":"))

    print(f"\nMerge complete! Output saved to: {output_path}")


def main():
    fas_version = version("greedyFAS")
    parser = argparse.ArgumentParser(description='You are running FAS version ' + str(fas_version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/FAS/wiki")
    parser.add_argument("-i", "--input_folder", required=True, help="Folder containing JSON files to merge.")
    parser.add_argument("-o", "--output", required=True, help="Name/path of the output merged JSON file.")

    args = parser.parse_args()
    merge_json_files(args.input_folder, args.output)


if __name__ == "__main__":
    main()
