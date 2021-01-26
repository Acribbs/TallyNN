.. _manual-main:

=====================
TallyNN documentation
=====================

TallyNN is a collection of single-cell workflows that allow users to perform barcode and UMI correction
for oligonucleotide sequences that are synthesised using double phosphoramidites for droplet based
single-cell sequencing.

The workflow management systems that underpins all of our pipelines is `CGAT-core <>https://github.com/cgat-developers/cgat-core>`_.
We demonstrate the functionality of CGAT-core using a simple RNA-seq pipeline in `cgat-showcase <https://github.com/cgat-developers/cgat-showcase>`_.
Therefore, if you are not familiar with how we build our workflows I suggest that you look at these two pipelines first.


.. _manual-quick_example:

--------
Citation
--------

Watch this space....

.. _manual-support:

-------
Support
-------

- Please refer to our :ref:`FAQ` section
- For bugs and issues, please raise an issue on `github <https://github.com/Acribbs/TallyNN/issues>`_
- For contributions, please refer to our contributor section and `<project_info/Contributingt>`_ source code.


.. toctree::
   :caption: Getting started
   :name: getting_started
   :maxdepth: 1
   :hidden:

   getting_started/Installation.rst
   getting_started/Cluster_config.rst
   getting_started/Tutorial.rst

.. toctree::
   :caption: Pipelines
   :name: build
   :maxdepth: 1
   :hidden:

   pipelines/nanopore.rst
   pipelines/illumina.rst

.. toctree::
   :caption: Project Info
   :name: project-info
   :maxdepth: 1
   :hidden:

   project_info/Contributing.rst
   project_info/how_to_contribute.rst
   project_info/FAQ.rst
   project_info/Licence.rst
