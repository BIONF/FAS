# FAS - Feature Architecture Similarity
[![PyPI version](https://badge.fury.io/py/greedyFAS.svg)](https://badge.fury.io/py/greedyFAS)
[![GPLv3-license](https://anaconda.org/bionf/hamstr/badges/license.svg)](https://www.gnu.org/licenses/gpl-3.0.de.html)
![Github Build](https://github.com/BIONF/FAS/workflows/build/badge.svg)

FAS is a new release of the original [FACT](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-417) algorithm. It calculates the so called FAS-score which is a measure of how similar the feature architectures of two proteins are. This is done by combining the Multiplicity Score (MS) and the Positional Score (PS) from FACT. Unlike the original FACT, FAS can resolve feature architectures that have overlapping features by searching for the best overlap-free path. This can be done either extensively or by using the priority mode, a greedy approach. FAS also allows for more options in the weighting of features.

# Table of Contents
* [Installation](#installation)
* [Usage](#usage)
  * [Annotate protein features](#annotate-protein-features)
  * [Compare protein feature architectures](#compare-protein-feature-architectures)
* [Additional Information](#additional-information)
* [How-To Cite](#how-to-cite)
* [Contributors](#contributors)
* [Contact](#contact)

# Installation

FAS is provided as a python package and compatible with **Python3**.

You can install FAS with pip:
```
python3 -m pip install greedyFAS
```

(\*) In case you **do not have admin rights**, and don't use package systems like Anaconda to manage environments you need to use the **--user** option (not recommended):
```
python3 -m pip install --user greedyFAS
```

and then add the following line to the end of your `.bashrc` or `.bash_profile` file, restart the current terminal to apply the change:
```
export PATH=$HOME/.local/bin:$PATH
```

# Usage

## Download and install annotation tools
Before using FAS, some annotation tools and databases need to be installed. FAS' standard databases/annotation tools are: [PFAM](https://www.ebi.ac.uk/interpro/download/Pfam/), [SMART](https://software.embl-em.de/software/18), [fLPS](http://biology.mcgill.ca/faculty/harrison/flps.html), [SEG](https://mendel.imp.ac.at/METHODS/seg.server.html), [COILS](https://mybiosoftware.com/coils-2-2-prediction-coiled-coil-regions-proteins.html), [THMHH 2.0c](http://www.cbs.dtu.dk/services/TMHMM/) and [SignalP 4.1g](http://www.cbs.dtu.dk/services/SignalP/). To get these tools and make a configuration file for FAS, please use the `setupFAS` function:
```
fas.setup -t /directory/where/you/want/to/save/annotation/tools
```
Inside the output directory you will find a file called *annoTools.txt* that contains all installed annotation tools. If you wish to discard any of them from the annotation process, you can just remove the unneeded tools from that file.

*Please read our [wiki page of setupFAS](https://github.com/BIONF/FAS/wiki/setupFAS) for other use-cases, such as how to use your old annotation tools with the new FAS, etc.*

__*NOTE: we provide compiled code only for PFAM, COILS and SEG. fLPS will be automatically downloaded and installed. For SMART, you need to download it from [EMBLEM](https://software.embl-em.de/software/18) and give the path to `fas.setup`. For TMHMM and SignalP, you can decide if you want to include those two tools to the annotation step (recommended) or ignore them. For using TMHMM version 2.0c and SignalP version 4.1g, you need to request a license from the authors at https://services.healthtech.dtu.dk, and save the downloaded files in the same directory. FAS will do the rest for you ;-)*__

__*NOTE2: SignalP 5.0b is not supported yet!!!*__

We suggest you test the annotation tools by running this command:
```
fas.doAnno -i test_annofas.fa -o test_output
```
*`test_annofas.fa` is a demo multiple fasta file, which is saved in the installed greedyFAS directory.*

## Perform feature annotation

If you only want to annotate your protein sequences without calculating the FAS scores, you can use the `doAnno` function.

```
fas.doAnno --fasta your_proteins.fa --outPath /annotation/path/
```

The annotation output (`your_proteins.json` by default) will be saved in `/annotation/path/`.

Alternatively, you can do the annotation using [InterProScan](https://www.ebi.ac.uk/interpro/about/interproscan/) and use the function `parseAnno` to convert the InterProScan's *tsv* output into *json format* for using with FAS.

```
fas.parseAnno -i INPUT.tsv -o /annotation/path/INPUT.json -t <tool_name> -f <feature columns> ...
```

Please check the usage of `parseAnno` for more info (using `fas.parseAnno -h`).

## Compare protein feature architectures

The main purpose of FAS is to calculate the similarity score between 2 given proteins (or two list of proteins). This can be done using the `run` function.

```
fas.run -s seed.fa -q query.fa -a /annotation/path/ -o /output/path/
```
If the annotations of *seed* and *query* protein(s) already exist in `/annotation/path/` (*seed.json* and *query.json*, respectively), `run` will use these annotations for calculating the FAS scores. Otherwise, it will first annotate the proteins and then compare the feature architectures of those two protein sets.

# Additional Information

A thorough guide to all FAS commands and options can be found at [our WIKI page](https://github.com/BIONF/FAS/wiki).

# How-To Cite

Julian Dosch, Holger Bergmann, Vinh Tran, Ingo Ebersberger, FAS: assessing the similarity between proteins using multi-layered feature architectures, Bioinformatics, Volume 39, Issue 5, May 2023, btad226, https://doi.org/10.1093/bioinformatics/btad226

# Contributors
- [Ingo Ebersberger](https://github.com/ebersber)
- [Julian Dosch](https://github.com/JuRuDo)
- [Holger Bergmann](https://github.com/holgerbgm)
- [Vinh Tran](https://github.com/trvinh)

# Contact
Julian Dosch dosch@bio.uni-frankfurt.de

Ingo Ebersberger ebersberger@bio.uni-frankfurt.de
