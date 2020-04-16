# FAS - Feature Architecture Similarity

FAS is a new release of the original [FACT](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-417) algorithm. It calculates the so called FAS-score which is a measure of how similar the feature architectures of two proteins are. This is done by combining the Multiplicity Score (MS) and the Positional Score (PS) from FACT. Unlike the original FACT, FAS can resolve feature architectures that have overlapping features by searching for the best overlap-free path. This can be done either extensively or by using the priority mode, a greedy approach. FAS also allows for more options in the weighting of features.

# Table of Contents
* [Installation](#installation)
* [Usage](#usage)
  * [Annotate protein features](#annotate-protein-features)
  * [Compare protein feature architectures](#compare-protein-feature-architectures)
* [Additional Information](#additional-information)
* [Contributors](#contributors)
* [Contact](#contact)

# Installation

FAS is provided as a python package and compatible with both Python2 and Python3.

First, get the source code of FAS from our github:

```
git clone --depth=1 https://github.com/BIONF/FAS
```

Then go to FAS directory and run pip install:
```
cd FAS
pip install .
```

In case you do not have admin rights you need to use the --user option:
```
pip install --user .
```

and then add the following line to the end of your `.bashrc` or `.bash_profile` file, restart the current terminal to apply the change:
```
export PATH=$HOME/.local/bin:$PATH
```

Alternatively, you can install FAS directly from a Conda environment without the need of downloading the source code:
```
conda install -c BIONF fas
```
You may have to add bioconda to your channels beforehand as FAS requires the hmmer package:
```
conda config --add channels bioconda
```


# Usage

FAS comes with three main functions: **annoFAS** and **parseInterPro**, which generate the standard input for FAS, and the main FAS function **greedyFAS**.

## Annotate protein features
To compare the feature architecture of two proteins, first we need to have the feature annotation of those sequences. We provide the `annoFAS` function to do this task by assigning the features to your proteins based on 7 databases/annotation tools: [PFAM](https://pfam.xfam.org/) and [SMART](http://smart.embl-heidelberg.de/), [fLPS](http://biology.mcgill.ca/faculty/harrison/flps.html), [SEG](http://www.biology.wustl.edu/gcg/seg.html), [COILS](https://embnet.vital-it.ch/software/COILS_form.html), [THMHH](http://www.cbs.dtu.dk/services/TMHMM/) and [SignalP](http://www.cbs.dtu.dk/services/SignalP/). *NOTE: we provide compiled code for Pfam, SMART, COILS and SEG. fLPS will be automatically downloaded and installed. For TMHMM and SignalP, you can decide if you want to include those two tools to the annotation step (recommended) or ignore them. For using TMHMM and SignalP, you need to request a license from the authors at https://services.healthtech.dtu.dk, and save the downloaded files in the same directory. FAS will do the rest for you ;-)*

First, you need to download the annotation tools to use FAS.
```
annoFAS --fasta seed.fa --path PATH --name anno --prepare
```
The annotation tools will be download and saved in your selected directory. Inside this directory you will find a file called *annoTools.txt* that contains all installed annotation tools. If you wish to discard any of them from the annotation process, you can just edit that file.

Your two input proteins (seed and ortholog) must be in FASTA format. Using the following commands to do the annotations for the two sequences:

```
annoFAS --fasta seed.fa --path PATH --name seed
annoFAS --fasta ortholog.fa --path PATH --name ortholog
```

This will output two folders `seed` and `ortholog` (as being defined using the `--name` parameter), each contains 7 XML files corresponding for 7 reference databases/annotation tools. These folders will be the input for FAS.


Alternatively, you can do the annotation using [InterProScan](https://www.ebi.ac.uk/interpro/about/interproscan/) and use the function `parseInterPro` to convert the InterProScan's *tsv* output into *XML format* for using with FAS

```
parseInterPro -i INPUT.tsv -s PATH/seed -j output
```

## Compare protein feature architectures

Once you have annotated the features of the seed and ortholog proteins, you are ready to use FAS algorithm with the two output annotation folders

```
greedyFAS -q PATH/ortholog -s PATH/seed -j PATH/OUTPUTNAME
```

The `-j` variable allows you to set an *output name* and *output path*. If no path is given, the default outputs will be `out.xml` and `out_architecture.xml`.

# Additional Information

A thorough guide to all FAS commands and options can be found at [our WIKI page](https://github.com/BIONF/FAS/wiki).

# Contributors
- [Ingo Ebersberger](https://github.com/ebersber)
- [Julian Dosch](https://github.com/JuRuDo)
- [Holger Bergmann](https://github.com/holgerbgm)
- [Vinh Tran](https://github.com/trvinh)

# Contact
Julian Dosch dosch@bio.uni-frankfurt.de

Ingo Ebersberger ebersberger@bio.uni-frankfurt.de
