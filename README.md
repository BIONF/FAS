# FAS - Feature Architecture Similarity
FAS is a new release of the original [FACT](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-417) algorithm. It calculates the so called FAS-score which is a measure of how similar the feature architectures of two proteins are. This is done by combining the Multiplicity Score (MS) and the Positional Score (PS) from FACT. Unlike the original FACT, FAS can resolve feature architectures that have overlapping features by searching for the best overlap-free path. This can be done either extensively or by using the priority mode, a greedy approach. FAS also allows for more options in the weighting of features.

The main FAS script is written in Python and should run on both Python 2 and Python 3. The additional annotation script that generates the standart input for FAS is written in Perl

# Table of Contents
* [Installation](#installation)
* [Usage](#usage)
* [Additional Information](#additional-information)
* [Contributors](#contributors)
* [Contact](#contact)



# Installation
FAS can be cloned from github here:

```
git clone --depth=1 https://github.com/BIONF/FAS
```

Then go to FAS directory and run pip install:
```
cd FAS
pip install .
```


Alternatively, you can use Anaconda:
```
conda install -c BIONF fas
```

# Usage
FAS comes with three main scripts: annoFAS and parseInterPro, which generate the standart input for FAS, and the actual FAS script greedyFAS.
To get started using FAS you need the protein sequence of the two (or more) proteins you want to compare. You should have two file in (Multi-)Fasta format, one for the seed protein(s) and one for the ortholog(s). Begin by using the annoFAS command:

```
annoFAS --fasta seed.fasta --path PATH --name seed
annoFAS --fasta orthologs.fasta --path PATH --name ortholog
```

This should give you an output folder of the chosen name containing seven xml files, one for each feature type used in the default set from FACT. Once you have annotated the features of both, the seed and ortholog proteins, you are ready to use the actual FAS algorithm with the two output folders of annotation script. The -j variable allows you to set an outputname and output path. If no path is given the output will be named out:

```
greedyFAS -q PATH/ortholog -s PATH/seed -j PATH/JOBNAME
```

This should give two xml files as output, out.xml and out_architecture.xml.
The third command parseInterPro can be used to parse the interPro tsv format and create an input for FAS:

```
parseInterPro -i INPUT.tsv -s PATH/seed -j output
```

# Additional Information

A thorough guide to all FAS commands and options can be found at [our WIKI page](https://github.com/BIONF/FAS/wiki).

# Contributors
- [Ingo Ebersberger](https://github.com/ebersber)
- [Julian Dosch](https://github.com/JuRuDo)
- [Vinh Tran](https://github.com/trvinh)
- [Holger Bergmann](https://github.com/holgerbgm)

# Contact
Julian Dosch dosch@bio.uni-frankfurt.de

Ingo Ebersberger ebersberger@bio.uni-frankfurt.de
