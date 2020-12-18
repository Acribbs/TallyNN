========
Nanopore
========

Overview
========

This pipeline performs alignment free based quantification of scBUC-seq for
single-cell sequencing long-read nanopore data.

The pipeline performs the following analyses:
* Barcode and UMI correction
* Alignment using kallisto
* QC of reads and generation of a barnyard plot for species mixing experiments

Input files
-----------

The pipeline is ran using fastq files added to a directory called data/

 * a fastq file
 * a single cDNA of all transcripts. Downloaded from Gencode. e.g. Homo_sapiens.GRCh38.cdna.all.fa.gz. If you
 have a species mixed experiment then also download Mus_musculus.GRCm38.cdna.all.fa.gz and cat both cdna files
 into one and modify the pipeline.yml file.
 * Genome fasta file, downloaded from Gencode

Configuring the pipeline
------------------------

To generate a configuration file run::

  tallynn nanopore config -v5

Running the pipeline
--------------------

To run the pipeline you will need to set up the cluster configuration according
to the cluster documentation.

However the pipeline can also be run locally without the cluster using the
commandline flag `--no-cluster`.

The following command will run the pipeline::

   tallynn nanopore make full -v5

output
------

The output of the pipeline is count matrix in market matrix format. This can be found within the
directory mtx.dir/
