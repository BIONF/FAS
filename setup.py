#!/bin/env python

#######################################################################
# Copyright (C) 2019 Julian Dosch
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

from setuptools import setup, find_packages

with open("README.md", "r") as input:
    long_description = input.read()

setup(
    name="greedyFAS",
    version="1.5.3",
    python_requires='>=3.7.0',
    description="A tool to compare protein feature architectures",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Julian Dosch",
    author_email="Dosch@bio.uni-frankfurt.de",
    url="https://github.com/BIONF/FAS",
    packages=find_packages(),
    package_data={'': ['*']},
    install_requires=[
        'biopython',
        'tqdm',
        'graphviz'
    ],
    entry_points={
        'console_scripts': ["annoParserFAS = greedyFAS.annoFAS.annoParserFAS:main",
                            "annoFAS = greedyFAS.annoFAS.annoFAS:main",
                            "setupFAS = greedyFAS.setupFAS:main",
                            "calcFAS = greedyFAS.calcFAS:main",
                            "fdogFAS = greedyFAS.fdogFAS:main",
                            "complexityFAS = greedyFAS.complexityFAS:main",
                            "domainFAS = greedyFAS.domainFAS:main"],
    },
    license="GPL-3.0",
    classifiers=[
        "Environment :: Console",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: End Users/Desktop",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],
)
