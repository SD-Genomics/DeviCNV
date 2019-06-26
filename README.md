# DeviCNV
Detection and Visualization of Exon-Level Copy Number Variants in Targeted Next Generation Sequencing Data

# Environment 
DeviCNV runs on Python 2.7 and R 3.2.0.

# Python dependencies
- sys
- intervaltree
- intervaltree_bio
- numpy
- operator
- random
- pysam
- pyvcf
- scipy

# R dependencies
- ggplot2
- PSCBS

# Installation
To install DeviCNV, simply download 9 scripts in “Code” directory.

# Documentation
 PDF documentation is included in the package. 
 - DeviCNV1.5 Manual20171101.pdf

# Version description
We uploaded DeviCNV_v1.5.1 in 26/06/2019
We fixed some bugs in code.
1. Delete codes for running with Slurm Workload Manager in "DeviCNV_Example.runningScript.sh"
2. Fix codes for selecting large segments in "python.scoreCNV.py"
