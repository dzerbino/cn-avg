CN-AVG code
===========

This package was used to create the simulations and figures for the CN-AVG theory publication (In preparation).

Pre-requisites:
---------------

- DNAcopy package from BioConductor/R. Inside R:
'''
source("http://bioconductor.org/biocLite.R")
biocLite("DNAcopy")
'''

- Python libraries
'''
pip install pyvcf 
'''

- Scipy: 
I suffered many times to install Scipy till
I discovered Anaconda: http://docs.continuum.io/anaconda/index.html

- Cactus library
'''
git clone git://github.com/benedictpaten/sonLib.git
git clone git://github.com/benedictpaten/cactus.git
'''

- For the automated pipeline:
'''
git clone git://github.com/benedictpaten/jobTree.git
'''
bedGraphToBigWig

Installation:
-------------

'''
make
make install
'''

Point PYTHONPATH to this directory.
