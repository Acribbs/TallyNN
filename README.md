
<img src="img/Nanopore-workflow.png" height=300>

Installation
============

We reccomend installing [miniconda](https://docs.conda.io/en/latest/miniconda.html), then creating
a new environment and install mamba

  ```
  conda create -n tallynn
  conda install mamba -c conda-forge
  ```
  
Next install the required software

  ```
  mamba install cgatcore samtools minimap2 subread
  ```

Then, at the moment, you will need to manually install the fork of umi tools

  ```
  git clone https://github.com/Acribbs/UMI-tools
  python setup.py install
  ```
  
To install tallynn code

  ```
  python setup.py install
  ```

Documentation
=============

Further information how you can run the pipelines can be found at [read the docs](https://tallynn.readthedocs.io/en/latest/)

Usage
=====

Run the ``tallynn --help`` command to see what workflows are available and ``tallynn nanopore -help`` to see how to use them.


For example, to generate a configuration file run

   ```
   tallynn nanopore config
   ```

To set up the configuration file please refer to [read the docs]().

Then configure your directory as directed by read the docs [set up]().

To run the pipeline with all tasks then run
   
   ```
   tallynn nanopore make full -v5 
   ```
