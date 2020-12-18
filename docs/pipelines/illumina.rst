
========
Illumina
========

Overview
========

This pipeline performs alignment free based quantification of scBUC-seq for
single-cell sequencing short-read illumina data.

The pipeline performs the following analyses:
* Barcode and UMI correction
* Alignment using kallisto
* QC of reads and generation of a barnyard plot for species mixing experiments

Input files
-----------

The pipeline is ran using fastq files that follow the naming convention Read1: Name.fastq.1.gz
and read2: Name.fastq.2.gz. The fastq files are added to a directory called data.dir/

 * a fastq file (paired end following the naming convention below)
 * a cDNA of all transcripts. Downloaded from Gencode. e.g. Homo_sapiens.GRCh38.cdna.all.fa.gz. If you
 have a species mixed experiment then also download Mus_musculus.GRCm38.cdna.all.fa.gz.

The default file format assumes the following convention:
fastq.1.gz and fastq.2.gz for paired data, where fastq.1.gz contains
UMI/cellular barcode data and fastq.2.gz contains sequencing reads.

Configuring the pipeline
------------------------

To generate a configuration file run::

  tallynn illumina config -v5

Running the pipeline
--------------------

To run the pipeline you will need to set up the cluster configuration according
to the cluster documentation.

However the pipeline can also be run locally without the cluster using the
commandline flag `--no-cluster`.

The following command will run the pipeline::

   tallynn illumina make full -v5

output
------

The output of the pipeline is count matrix in market matrix format. This can be found within the
directory kallistao.dir/<name_of_fastq_file>/bus/genecount/

 
