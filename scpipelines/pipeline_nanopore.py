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


SEQUENCESUFFIXES = ("*.fastq")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])




@follows(mkdir("split_tmp.dir"))
@split('data.dir/*.fastq.gz',"split_tmp.dir/out*")
def split_fastq(infile, outfiles):
    '''
    Split the fastq file before identifying perfect barcodes
    '''

    infile = "".join(infile)

    statement = '''zcat %(infile)s | split -l %(split)s - out. &&
                   mv out*.* split_tmp.dir/'''

    P.run(statement)


@follows(mkdir("polyA_correct.dir"))
@transform(split_fastq,
           regex("split_tmp.dir/out.(\S+)"),
           r"polyA_correct.dir/\1_correct_polya.fastq")
def correct_polyA(infile, outfile):
    '''
    Look across all fastq files and correct the polyA so its on same strand
    '''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/complement_polyA.py --infile=%(infile)s --outname=%(outfile)s'''

    P.run(statement)


@follows(mkdir("perfect_reads.dir"))
@transform(correct_polyA,
           regex("polyA_correct.dir/(\S+)_correct_polya.fastq"),
           r"perfect_reads.dir/\1_ambiguous_barcode_R1.fastq")
def identify_perfect(infile, outfile):
    '''
    Identify ambigous and unabiguous reads
    '''

    name  = outfile.replace("_ambiguous_barcode_R1.fastq","")

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/identify_perfect_nano.py --outname=%(name)s --infile=%(infile)s --whitelist=%(name)s.whitelist.txt'''

    P.run(statement)


@merge(identify_perfect, "unambiguous.fastq.1")
def merge_unambiguous(infiles, outfile):
    '''
    in order to identiify the true whitelist barcodes the next step is to merge the fastq files together
    '''

    infile = []
    infile2 = []

    for i in infiles:
        infile2.append(i.replace("_ambiguous_barcode_R1.fastq","_unambiguous_barcode_R2.fastq"))
        infile.append(str(i.replace("_ambiguous_barcode_R1.fastq","_unambiguous_barcode_R1.fastq")))

    outfile2 = outfile.replace(".fastq.1",".fastq.2")
    infiles = " ".join(infile)
    infiles2 = " ".join(infile2)

    statement = '''cat %(infiles)s > %(outfile)s &&
                   cat %(infiles2)s > %(outfile2)s'''

    P.run(statement)


@merge(identify_perfect, "whitelist.txt")
def merge_whitelist(infiles, outfile):
    '''
    merge whitelists
    '''

    whitelists = []

    for i in infiles:
        print(i)
        whitelists.append(i.replace("_ambiguous_barcode_R1.fastq",".whitelist.txt"))


    whitelist_files = " ".join(whitelists)

    statement = '''cat %(whitelist_files)s | sort | uniq > %(outfile)s'''

    P.run(statement)


@merge(identify_perfect, "ambiguous.fastq.1")
def merge_ambiguous(infiles, outfile):
    '''
    in order to identiify the true whitelist barcodes the next step is to merge the fastq files together
    '''

    infile = []
    infile2 = []

    for i in infiles:
        infile2.append(i.replace("_ambiguous_barcode_R1.fastq","_ambiguous_barcode_R2.fastq"))
        infile.append(str(i))

    outfile2 = outfile.replace(".fastq.1",".fastq.2")
    infiles = " ".join(infile)
    infiles2 = " ".join(infile2)

    statement = '''cat %(infiles)s > %(outfile)s &&
                   cat %(infiles2)s > %(outfile2)s'''

    P.run(statement)


@follows(mkdir("whitelist.dir"))
@transform(merge_unambiguous,
           regex("unambiguous.fastq.1"),
           r"whitelist.dir/whitelist.txt")
def whitelist_umitools(infile, outfile):
    ''' '''

    
    cell_num = PARAMS['whitelist']
    statement = '''umi_tools whitelist --stdin=%(infile)s --bc-pattern=CCCCCCCCCCCCCCCCCCCCCCCCNNNNNNNNNNNNNNNN
                   --set-cell-number=%(cell_num)s -L extract.log > whitelist.dir/whitelist.txt
 '''

    P.run(statement)


@follows(mkdir("correct_reads.dir"))
@transform(identify_perfect,
         regex("perfect_reads.dir/(\S+)_ambiguous_barcode_R1.fastq"),
         add_inputs(merge_whitelist),
         r"correct_reads.dir/\1_unambiguous_fixed_barcode_R1.fastq")
def correct_reads(infiles, outfile):
    '''Use levenshtein distance and correct the barcodes '''

    infile, whitelist_file = infiles

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    infile2 = infile.replace("_ambiguous_barcode_R1.fastq","_ambiguous_barcode_R2.fastq")
    name = outfile.replace("_unambiguous_fixed_barcode_R1.fastq","")

    ld = PARAMS['distance']

    statement = '''python %(PYTHON_ROOT)s/correct_barcode_nano.py --whitelist=%(whitelist_file)s --read1=%(infile)s --read2=%(infile2)s
                   --outname=%(name)s --distance=%(ld)s'''

    P.run(statement)


@merge(correct_reads, "merge_corrected.fastq.1")
def merge_correct_reads(infiles, outfile):
    '''Merge the corrected reads '''

    infile = []
    infile2 = []

    for i in infiles:
        infile2.append(i.replace("_R1.fastq","_R2.fastq"))
        infile.append(str(i))
 
    infiles = " ".join(infile)
    infiles2 = " ".join(infile2)

    outfile2 = outfile.replace(".fastq.1",".fastq.2")

    statement = '''cat %(infiles)s > %(outfile)s &&
                   cat %(infiles2)s > %(outfile2)s'''

    P.run(statement)


@transform(merge_correct_reads,
           regex("merge_corrected.fastq.1"),
           add_inputs(merge_unambiguous),
           r"final.fastq.1.gz" )
def merge_full(infiles, outfile):
    ''' '''

    infile1_R1, infile2_R1 = infiles

    infile1_R2 = infile1_R1.replace(".fastq.1",".fastq.2") 
    infile2_R2 = infile2_R1.replace(".fastq.1",".fastq.2")

    statement = '''cat %(infile1_R1)s %(infile2_R1)s| gzip  > final.fastq.1.gz &&
                  cat %(infile1_R2)s %(infile2_R2)s | gzip  > final.fastq.2.gz'''

    P.run(statement)


@transform(merge_full,
           regex("final.fastq.1.gz"),
           r"final_extract.fastq.1.gz")
def extract_barcodeumi(infile, outfile):
    ''' '''

    infile2 = infile.replace(".fastq.1.gz",".fastq.2.gz")
    name = outfile.replace(".fastq.1.gz","")
    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/extract_umibc_readname.py --read1=%(infile)s --read2=%(infile2)s --outname=%(name)s'''

    P.run(statement)


@follows(mkdir("whitelist.dir/"))
@transform(merge_full,
           regex("final.fastq.1.gz"),
           r"whitelist.dir/whitelist2.txt")
def whitelist_umitools2(infile, outfile):
    ''' '''

    
    cell_num = PARAMS['whitelist']
    statement = '''umi_tools whitelist --stdin=%(infile)s --bc-pattern=CCCCCCCCCCCCCCCCCCCCCCCCNNNNNNNNNNNNNNNN
                   --set-cell-number=%(cell_num)s -L extract.log > whitelist.dir/whitelist2.txt
 '''

    P.run(statement)


@transform(merge_full,
           regex("final.fastq.1.gz"),
           add_inputs(whitelist_umitools2),
           r"final_extract_umitools.fastq.1.gz")
def extract_barcodeumitools(infiles, outfile):
    ''' '''
    infile, whitelist = infiles

    whitelist = "".join(whitelist)

    infile2 = infile.replace(".fastq.1.gz",".fastq.2.gz")
    outfile2 = outfile.replace(".fastq.1.gz",".fastq.2.gz")

    statement = '''umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCCCCCCCCCNNNNNNNNNNNNNNNN --stdin %(infile)s 
                   --stdout=%(outfile)s --read2-in %(infile2)s --read2-out=%(outfile2)s --whitelist=%(whitelist)s'''

    P.run(statement)

@transform(extract_barcodeumitools,
           regex("final_extract_umitools.fastq.1.gz"),
           r"final.sam")
def mapping(infile, outfile):
    '''Run minimap2 to map the fastq files'''

    infile = infile.replace(".fastq.1.gz",".fastq.2.gz")
    
    cdna = PARAMS['minimap2_fasta_cdna']
    options = PARAMS['minimap2_options']

    statement = '''minimap2  %(options)s %(cdna)s  %(infile)s > %(outfile)s 2> %(outfile)s.log'''

    P.run(statement)


@transform(mapping,
           regex("final.sam"),
           r"final_sorted.bam")
def run_samtools(infile, outfile):
    '''convert sam to bam and sort'''

    statement = '''samtools view -bS %(infile)s > final.bam &&
                   samtools sort final.bam -o final_sorted.bam &&
                   samtools index final_sorted.bam'''

    P.run(statement)


@transform(run_samtools,
           regex("final_sorted.bam"),
           r"final_XT.bam")
def add_xt_tag(infile, outfile):
    '''Add trancript name to XT tag in bam file so umi-tools counts can be  perfromed'''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")


    statement = '''python %(PYTHON_ROOT)s/xt_tag_nano.py --infile=%(infile)s --outfile=%(outfile)s &&
                   samtools index %(outfile)s'''

    P.run(statement)


@transform(add_xt_tag,
         regex("final_XT.bam"),
         r"counts.tsv.gz")
def count(infile, outfile):
    '''use umi_tools to count the reads - need to adapt umi tools to double oligo'''

    statement = '''umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I %(infile)s -S counts.tsv.gz'''

    P.run(statement)

# Need to merge correct reads with the unambiguous reads
# Then run umi-tools whitelist
# Extract umi and barcode and add to read name
# Then map using minimap2
# Run UMI toos to demultiplex
# Add tag XT to bam file and run umi-tools counts


@follows(count)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
