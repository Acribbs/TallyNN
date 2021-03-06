.. _getting_started-tutorial:

==============
Tutorial guide
==============

In the following page, we describe the steps involved in running TallyNN to
error correct barcodes and UMIs. The pipeline can be executed using a test dataset
that contains 100,000 Nanopore sequencing reads (https://www.cgat.org/downloads/public/adam/TallyNN/test_nanopore.fastq.gz)
or full data (https://www.cgat.org/downloads/public/adam/TallyNN/E1.fastq.gz).

The fastq files are derived from an experiment in which an equal mix of JJN3, NCI-H929 and DF15 myeloma cells.

Download the files
------------------

The demo data is available at the following link - https://www.cgat.org/downloads/public/adam/TallyNN/::

    wget https://www.cgat.org/downloads/public/adam/TallyNN/test_nanopore.fastq.gz
    wget https://www.cgat.org/downloads/public/adam/TallyNN/Homo_sapiens-GRCh38.fa
    wget https://www.cgat.org/downloads/public/adam/TallyNN/Homo_sapiens.GRCh38.101.bed
    wget https://www.cgat.org/downloads/public/adam/TallyNN/hg38.fasta
    wget https://www.cgat.org/downloads/public/adam/TallyNN/Homo_sapiens.GRCh38.101.gtf
    mkdir data
    mv test_nanopore.fastq.gz data/

To run the demo you will need to download the fasta file of human transcripts, the fastq file of
nanopore sequenced reads, a genome fasta file and a junction bed file generated according to minimap2 [documentation](https://lh3.github.io/minimap2/minimap2.html).

Add the downloaded data to a directory called `data.dir/`

Generate the config file
------------------------

In order to run the nanopore pipeline you will need to generate a configuration
file that can be used to modify the execution of the pipeline.

This can be performed by::

    tallynn nanopore config

This will generate a pipeline.yml file in the current directory. However for the
purpose of this tutorial please download the pipeline.yml as follows::

    wget https://www.cgat.org/downloads/public/adam/TallyNN/pipeline.yml

Modify the config file
----------------------

Open up the pipeline.yml file. You can see a set of key value pairs that can be modified
to change the running of the pipeline. For the current analysis, no changes need to be made to this file.

Run the pipeline
----------------

The defult behaviour is for the pipeline to execute across a cluster. This is why it is important to set up
your `.cgat.yml` file correctly. In order to run the pipeline to completion run the following in the commandline::

    tallynn nanopore make full -v5

Otherwise, the pipeline can also be ran locally without a cluster::

    tallynn nanopore make full -v5 --no-cluster
