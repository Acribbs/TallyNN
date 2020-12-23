"""
====================
Pipeline fusiontrans
====================


Overview
==================

The pipeline takes as an input a bam file that was ran using featurecounts from pipeline_nanopore.py and generates a bed file of
fusion transcripts.

Usage
=====

Configuration
-------------

The pipeline uses CGAT-core and CGAT-apps throughout the pipeline. Please see installation
and setup and installation instructions at `cgat-core documentation <>`_


Input files
-----------



Pipeline output
==================


Code
==================

"""
from ruffus import *

import sys
import os
import re
import sqlite3
import glob
import pysam

import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.database as database

# Load options from the config file

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


# Determine the location of the input fastq files

try:
    PARAMS['data']
except NameError:
    DATADIR = "."
else:
    if PARAMS['data'] == 0:
        DATADIR = "."
    elif PARAMS['data'] == 1:
        DATADIR = "data"
    else:
        DATADIR = PARAMS['data']


SEQUENCESUFFIXES = ("*.bam")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

@transform(SEQUENCEFILES,
           regex("data/(\S+).bam"),
           r"\1_SA.bam")
def make_sabam(infile, outfile):
    ''' Generate a bam file containing all reads listed as
        chimeric'''

    bamfile = pysam.AlignmentFile(infile,"rb")
    split = pysam.AlignmentFile(outfile, "wb", template=bamfile)

    for line in bamfile:
        if line.has_tag("SA"):
            split.write(line)
        else:
            pass
    
    bamfile.close()
    split.close()


@transform(PARAMS['bed'],
           suffix(".bed"),
           ".bed.gz.tbi")
def tabix_bed(infile, outfile):
    ''' convert bed file to tabix'''

    statement = '''cat %(infile)s | bgzip -c > %(infile)s.gz &&
                   tabix -p bed %(infile)s.gz'''

    P.run(statement)


@transform(make_sabam,
           regex("(\S+).bam"),
           add_inputs(tabix_bed),
           r"\1_FusionAnnotate.bam")
def fusion_annotate(infiles, outfile):
    '''Annotate fusion and original gene'''

    
    infile, bed = infiles
    bed = bed.replace(".tbi","")

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python  %(PYTHON_ROOT)s/fusion_annotate.py --infile=%(infile)s --outfile=%(outfile)s --bedfile=%(bed)s'''

    P.run(statement)


@transform(fusion_annotate,
           regex("(\S+)_FusionAnnotate.bam"),
           add_inputs(tabix_bed),
           r"\1_FinalAnnotate.bam")
def gene_annotate(infiles, outfile):
    '''Annotate gene and original gene'''

    
    infile, bed = infiles
    bed = bed.replace(".tbi","")
    outfile_tmp = outfile +  ".tmp"

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python  %(PYTHON_ROOT)s/gene_annotate.py --infile=%(infile)s --outfile=%(outfile_tmp)s --bedfile=%(bed)s &&
                   samtools sort %(outfile_tmp)s -o %(outfile)s &&
                   rm -rf %(outfile_tmp)s'''

    P.run(statement)


@transform(gene_annotate,
          regex("(\S+)_FinalAnnotate.bam"),
          [r"\1_fusion1.bed",r"\1_fusion2.bed"])
def generate_bedout(infile, outfiles):
    '''Generate output bed file '''

    outfile1, outfile2 = outfiles
    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/bed_fusion.py --infile=%(infile)s --bed1=%(outfile1)s --bed2=%(outfile2)s'''

    P.run(statement)


@transform(generate_bedout,
           suffix("_fusion1.bed"),
           "_counts.txt")
def generate_counts(infiles, outfile):
    '''Generate counts for each fusion gene'''
    
    bed1, bed2 = infiles

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/generate_counts.py --bed1=%(bed1)s --bed2=%(bed2)s --outfile=%(outfile)s'''

    P.run(statement)


@follows(generate_counts)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))

