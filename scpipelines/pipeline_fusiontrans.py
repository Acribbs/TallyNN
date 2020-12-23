"""
====================
Pipeline fusiontrans
====================


Overview
==================

The pipeline takes as an input a bam file and generates a bed file of
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

    statement = '''cat %(infile)s | bgzip -c > coding_gene_region.bed.gz &&
                   tabix -p bed coding_gene_region.bed.gz'''

    P.run(statement)


@transform(make_sabam,
           regex("(\S+).bam"),
           add_inputs(tabix_bed),
           r"\1_FusionAnnotate.bam")
def fusion_annotate(infiles, outfile):
    '''Annotate fusion and original gene'''

    infile, bed = infiles
    bed = bed.replace(".tbi","")

    bamfile = pysam.AlignmentFile(infile)
    tabixfile = pysam.TabixFile(bed)

    out_bam = pysam.AlignmentFile("test.bam", "wb", template=bamfile)

    for line in bamfile:

        a = line.get_tag("SA").split(";")[0]

        chrom, start, strand, cigar, mapq, nm = a.split(",")
    
        match = re.findall(r'[0-9]+[A-Za-z]', cigar)
        for m in match:
            if "M" in m:
                m = m.replace("M", "")
                if "chrUn" in chrom or "alt" in chrom or "KI" in chrom or "GL" in chrom:
                    pass
                else:
                    for gtf in tabixfile.fetch(chrom, int(start), int(start)+int(m)):
                        gtf = gtf.split("\t")
                        if gtf[0] == "chrM":
                            pass
                        else:
                            line.tags += [("Ta", gtf[0])]
                            line.tags += [("Tb", gtf[1])]
                            line.tags += [("Tc", gtf[2])]
                            line.tags += [("Td", gtf[3])]

                            out_bam.write(line)
        out_bam.close()
        bamfile.close()

@follows(make_sabam)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))

