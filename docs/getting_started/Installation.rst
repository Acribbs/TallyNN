.. _getting_started-Installation:

============
Installation
============

The following sections describe how to install the TallyNN analysis pipelines.


We recommend installing [miniconda](https://docs.conda.io/en/latest/miniconda.html), then creating
a new environment and install mamba::


  conda create -n tallynn
  conda install mamba -c conda-forge


Next install the required software::


  mamba install cgatcore samtools minimap2 subread


Then, at the moment, you will need to manually install the fork of umi tools::


  git clone https://github.com/Acribbs/UMI-tools
  python setup.py install


To install tallynn code::


  python setup.py install
