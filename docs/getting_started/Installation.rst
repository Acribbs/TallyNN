.. _getting_started-Installation:

============
Installation
============

The following sections describe how to install the TallyNN analysis pipelines.

First clone the Github repo and activate the UMI-tools submodule::

   git clone https://github.com/Acribbs/TallyNN.git
   git submodule init
   git submodule update


We recommend installing [miniconda](https://docs.conda.io/en/latest/miniconda.html), then creating
a new environment and install mamba::


  conda create -n tallynn
  conda install mamba -c conda-forge


Next install the required software::


  mamba install cgatcore samtools minimap2 subread


Then, at the moment, you will need to manually install the fork of umi tools.
The fork is provided as a submodule within TallyNN. You can install it as follows::


  cd tallynn/UMI-tools
  python setup.py install


To install tallynn code, navigate to the tallynn directory and run::


  python setup.py develop
