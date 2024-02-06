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
import re

with open('README.md', 'r') as input:
    long_description = input.read()

setup(
    name='greedyFAS',
    version='1.18.10',
    python_requires='>=3.7.0',
    description='A tool to compare protein feature architectures',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Julian Dosch',
    author_email='Dosch@bio.uni-frankfurt.de',
    url='https://github.com/BIONF/FAS',
    packages=find_packages(),
    package_data={'': ['*']},
    install_requires=[
        'biopython',
        'tqdm',
        'graphviz',
        'gnureadline',
        'GitPython',
        'pathlib'
    ],
    entry_points={
        'console_scripts': ['fas.parseAnno = greedyFAS.annoFAS.annoParserFAS:main',
                            'fas.doAnno = greedyFAS.annoFAS.annoFAS:main',
                            'fas.checkAnno = greedyFAS.annoFAS.checkAnno:main',
                            'fas.updateAnnoInfo = greedyFAS.annoFAS.updateAnnoFile:main',
                            'fas.getProtByAnno = greedyFAS.annoFAS.getProtByAnno:main',
                            'fas.setup = greedyFAS.setupFAS:main',
                            'fas.run = greedyFAS.calcFAS:main',
                            'fas.runMultiTaxa = greedyFAS.calcFASmulti:main',
                            'fas.runFdogFas = greedyFAS.fdogFAS:main',
                            'fas.calcComplexity = greedyFAS.complexityFAS:main',
                            'fas.getDomains = greedyFAS.domainFAS:main',
                            'fas.splitAnno = greedyFAS.extractAnnoFAS:main',
                            'fas.disorder.predict = greedyFAS.disorderFAS.predict_disorder:main',
                            'fas.disorder.setup = greedyFAS.disorderFAS.install_aucpred:main',
                            'fas.mergeAnno = greedyFAS.mergeAnno:main',
                            'fas.updateAnno = greedyFAS.updateAnno:main',
                            'fas.mergeJson = greedyFAS.mergeJson:main',
                            'fas.updateJson = greedyFAS.updateJson:main',
                            'fas.tsv2json = greedyFAS.tsv2json:main',
                            'fas.getAnnoVersion = greedyFAS.annoFAS.getAnnoVersion:main',
                            'fas.overlapStatistics = greedyFAS.getOverlapStats:main',
                            'fad.run = greedyFAS.calcFAD:main'
                            ],
    },
    license='GPL-3.0',
    classifiers=[
        'Environment :: Console',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: End Users/Desktop',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
)
