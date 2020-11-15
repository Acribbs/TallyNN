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


SEQUENCESUFFIXES = ("*.fastq.gz")
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


@follows(correct_polyA)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
