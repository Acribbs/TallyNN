##############################################################################
#
#   Botnar Resaerch Centre
#
#   $Id$
#
#   Copyright (C) 2020 Adam Cribbs
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################

"""
====================
Pipeline single cell
====================


Overview
==================


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

import cgatcore.pipeline as P
import cgatcore.experiment as E

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
        DATADIR = "data.dir"
    else:
        DATADIR = PARAMS['data']


SEQUENCESUFFIXES = ("*.fastq.1.gz")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])




@follows(mkdir("split_tmp.dir"))
@split('data.dir/*.fastq.1.gz',"split_tmp.dir/read1*")
def split_fastq_perfect(infile, outfiles):
    '''
    Split the fastq file before identifying perfect barcodes
    '''
    infile = "".join(infile)
    second_read = infile.replace(".fastq.1.gz", ".fastq.2.gz")

    statement = '''zcat %(infile)s | split -l %(split)s - read1. &&
                   zcat %(second_read)s | split -l %(split)s - read2. &&
                   mv read*.* split_tmp.dir/'''

    P.run(statement)


@follows(mkdir("perfect_reads_split.dir"))
@transform(split_fastq_perfect,
           regex("split_tmp.dir/read1.(\S+)"),
           [r"perfect_reads_split.dir/\1_perfect.fastq.1.gz",r"perfect_reads_split.dir/\1_error.fastq.1.gz"])
def perfect_reads(infile, outfiles):
    '''calulate perfect and error reads in fastq'''
    
    read1 = infile
    read2 = infile.replace("read1.", "read2.")

    name = infile.replace("read1.","")
    name = name.replace("split_tmp.dir/","")

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)sidentify_perfect.py --read1=%(read1)s --read2=%(read2)s --outname=%(name)s'''

    P.run(statement)


@merge(perfect_reads, "perfect.fastq.1.gz")
def merge_fastq_perfect(infiles, outfile):
    '''merge read1 of fastq in perparation for running umi-tools whitelist'''

    infile = []
    infile2 = []

    for i in infiles:
        infile2.append(i[0].replace(".fastq.1.gz",".fastq.2.gz"))
        infile.append(str(i[0]))

    outfile2 = outfile.replace(".fastq.1.gz",".fastq.2.gz")

    infiles = " ".join(infile)
    infiles2 = " ".join(infile2)

    statement = '''cat %(infiles)s > %(outfile)s &&
                   cat %(infiles2)s > %(outfile2)s'''

    P.run(statement)

@follows(mkdir("whitelist.dir"))
@transform(merge_fastq_perfect,
           regex("(\S+).fastq.1.gz"),
           r"whitelist.dir/\1.txt")
def whitelist(infile, outfile):
    ''' '''

    cell_num = PARAMS['whitelist']
    statement = '''umi_tools whitelist --stdin=%(infile)s --bc-pattern=CCCCCCCCCCCCCCCCCCCCCCCCNNNNNNNNNNNNNNNN
                   --set-cell-number=%(cell_num)s -L extract.log > %(outfile)s
 '''

    P.run(statement)


@follows(mkdir("corrected_reads.dir"))
@transform(perfect_reads,
           regex("perfect_reads_split.dir/(\S+)_perfect.fastq.1.gz"),
           add_inputs(whitelist),
           r"corrected_reads.dir/\1_corrected.fastq.1.gz")
def correct_reads(infiles, outfile):
    '''Correct barcodes that are from the error files '''

    infile1, infile2 = infiles

    read1 = infile1[1].replace("_perfect.fastq.1.gz","_error.fastq.1.gz")
    read2 = infile1[0].replace("_perfect.fastq.1.gz","_error.fastq.2.gz")

    name = infile1[0].replace("_perfect.fastq.1.gz","")
    name = name.replace("perfect_reads_split.dir/","")

    distance = PARAMS['distance']

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")
    statement = '''python %(PYTHON_ROOT)s/correct_barcode.py --whitelist=%(infile2)s --read1=%(read1)s --read2=%(read2)s --outname=%(name)s --distance=%(distance)s'''

    P.run(statement)


@merge(correct_reads, "corrected.fastq.1.gz")
def merge_corrected(infiles, outfile):
    '''merge the corrected reads '''

    infile = []
    infile2 = []

    for i in infiles:
        infile2.append(i.replace(".fastq.1.gz",".fastq.2.gz"))
        infile.append(str(i))
 
    infiles = " ".join(infile)
    infiles2 = " ".join(infile2)

    outfile2 = outfile.replace(".fastq.1.gz",".fastq.2.gz")

    statement = '''cat %(infiles)s > %(outfile)s &&
                   cat %(infiles2)s > %(outfile2)s'''

    P.run(statement)


@transform(merge_corrected,
           suffix(".fastq.1.gz"),
           add_inputs(merge_fastq_perfect),
           r"final")
def combine_perfect_corrected(infiles, outfile):
    '''Combine the the perfect and the corrected reads together'''

    infile_1_R1, infile_2_R1 = infiles

    infile_1_R2 = infile_1_R1.replace(".fastq.1.gz",".fastq.2.gz") 
    infile_2_R2 = infile_2_R1.replace(".fastq.1.gz",".fastq.2.gz")

    statement = '''cat %(infile_1_R1)s %(infile_2_R1)s > final.fastq.1.gz &&
                  cat %(infile_1_R2)s %(infile_2_R2)s > final.fastq.2.gz '''

    P.run(statement)

@follows(combine_perfect_corrected)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
