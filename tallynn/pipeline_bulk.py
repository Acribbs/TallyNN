"""===========================
Pipeline bulk
===========================

Overview
========

This code is a CGAT pipeline for generating a count matrix for
downstream differential expression analysis using nanopore reads.
The pipeline takes an input fastq file, processes it, and outputs 
a count matrix with samples as columns and rows as either transcripts
or genes. The pipeline makes use of multiple Python libraries and tools
like Minimap2, Samtools, UMI-tools, and mclumi.

Pipeline tasks
==============

The pipeline consists of the following steps:

* Filtering reads with incorrect polyA orientation.
* Identifying polyA and TSO UMIs for each read.
* Mapping reads to transcripts using Minimap2.
* Processing and sorting BAM files using Samtools.
* Adding XT tags to SAM files.
* Generating counts tables for transcripts using UMI-tools and mclumi.
* Merging counts tables for transcripts.
* Mapping reads to genes using Minimap2.
* Sorting BAM files using Samtools.
* Running featureCounts on BAM files.
* Generating counts tables for genes using UMI-tools and mclumi.
* Merging counts tables for genes.

Usage
=====

The pipeline requires a configured pipeline.yml file,
which contains various settings and parameters required for
the pipeline to run. 

To generate the config file to change the running of the pipeline you need to
run:

tallynn bulk config

This will generate a pipeline.yml file that the user can modify to change the
output of the pipeline. Once the user has modified the pipeline.yml file the
pipeline can then be ran using the following commandline command:

tallynn bulk make full -v5

You can run the pipeline locally (without a cluster) using --local

tallynn bulk make full -v5 --local

Input files
===========

Input files should be fastq.gz files of nanopore reads
that have been sequenced with trimers at the 5' and 3' ends,
and should be added to the data.dir folder.

Pipeline output
===============

The pipeline outputs a counts matrix with sample as columns and
rows as either transcripts or genes.


Code
====

"""
import sys
import os
import pysam
import glob
import pandas as pd
from ruffus import *
import cgatcore.iotools as iotools
import cgatcore.pipeline as P
import cgatcore.experiment as E
from cgatcore.pipeline import cluster_runnable

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


SEQUENCESUFFIXES = ("*.fastq.gz")

FASTQTARGET = tuple([os.path.join("data.dir/", suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])


def merge_trans_noumi_data(infiles):
    """
    This function merges all the dataframes across samples
    for transcripts. It takes a list of input files and
    returns the merged dataframe.
    """

    final_df = pd.DataFrame()
    for infile in infiles:
        name = infile.replace(".counts_noumis.tsv.gz", "")
        tmp_df = pd.read_table(infile, sep="\t", header=1, index_col=0, skiprows = 1)
        tmp_df = tmp_df.iloc[:,-1:]
        tmp_df.columns = [name]
        final_df = final_df.merge(tmp_df, how="outer", left_index=True, right_index=True, suffixes=("","_drop"))

    return final_df

def merge_feature_data(infiles):
    '''
    This function merges all the input files provided in
    the list of input files. It returns a merged dataframe.
    '''

    final_df = pd.DataFrame()
    for infile in infiles:
        name = infile.replace(".counts.tsv.gz", "")
        tmp_df = pd.read_table(infile, sep="\t", header=0, names=["transcript_name", name], index_col=0)
        final_df = final_df.merge(tmp_df, how="outer", left_index=True, right_index=True)

    return final_df

def merge_featurecounts_data(infiles):
    '''
    This function merges all of the input files from featurecounts
    count output. It returns a merged dataframe with the
    column names replaced.
    '''

    final_df = pd.DataFrame()
    for infile in infiles:
    
        tmp_df = pd.read_table(infile, sep="\t", header=0, index_col=0, skiprows = 1)
        tmp_df = tmp_df.iloc[:,-1:]
        tmp_df.columns = ["count"]
        final_df = final_df.merge(tmp_df, how="outer", left_index=True, right_index=True, suffixes=("","_drop"))

    names = [x.replace("_gene_assigned.txt", "") for x in infiles]
    final_df.columns = names
    return final_df

# fastqsplitter input.fastq.gz -n 3 --prefix split

@follows(mkdir("processed_fastq.dir"))
@transform(FASTQTARGET,
         regex("data.dir/(\S+).fastq.gz"),
         r"processed_fastq.dir/\1_polyA.fastq.gz")
def polya_correct(infile, outfile):
    '''
    This function filters reads less than 300 bp and ensures the
    polyA is in the correct orientation. It takes an input file
    and an output file as arguments.
    '''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/complement_polyA_bulk.py --infile=%(infile)s --outname=%(outfile)s'''

    P.run(statement)

@transform(polya_correct,
         regex("processed_fastq.dir/(\S+)_polyA.fastq.gz"),
         r"processed_fastq.dir/\1_tso_UMI.fastq.gz")
def polya_umi(infile, outfile):
    '''
    This function identifies the polyA UMI for each read
    in the input file and writes the result to the output file.
    '''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = """python %(PYTHON_ROOT)s/polya_umi_nocorrect_bulk.py --infile=%(infile)s --outname=%(outfile)s """

    P.run(statement)


@transform(polya_umi,
         regex("processed_fastq.dir/(\S+)_tso_UMI.fastq.gz"),
         r"processed_fastq.dir/\1_polya_tso_UMI.fastq.gz")
def tso_umi(infile, outfile):
    '''
    This function identifies the TSO UMI for each read in
    the input file and writes the result to the output file.
    '''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")


    statement = """cp %(infile)s %(outfile)s"""


    P.run(statement)


@follows(mkdir("mapped_files.dir"))
@transform(tso_umi,
           regex("processed_fastq.dir/(\S+)_polya_tso_UMI.fastq.gz"),
           r"mapped_files.dir/\1_tso_polya_UMI.sam")
def mapping_trans(infile, outfile):
    '''
    This function maps the transcripts using minimap2 and
    writes the output to the specified outfile.
    '''

    statement = '''minimap2 -ax map-ont -p 0.9 --end-bonus 10 -N 3 %(cdna_fasta)s %(infile)s  > %(outfile)s 2> %(outfile)s.log'''

    P.run(statement)


@transform(mapping_trans,
           regex("mapped_files.dir/(\S+)_tso_polya_UMI.sam"),
           r"mapped_files.dir/\1_final_sorted.bam")
def samtools(infile, outfile):
    '''
    This function runs samtools on the input file
    and writes the output to the specified outfile.
    '''

    name = outfile.replace("_final_sorted.bam", "")

    statement = '''samtools view -bS %(infile)s > %(name)s_final.bam && 
                   samtools sort %(name)s_final.bam -o %(name)s_final_sorted.bam && 
                   samtools index %(name)s_final_sorted.bam '''

    P.run(statement)


@active_if(PARAMS['correct'])
@transform(samtools,
           regex("mapped_files.dir/(\S+)_final_sorted.bam"),
           r"mapped_files.dir/\1_XT.bam")
def xt_tag(infile, outfile):
    '''
    This function adds the XT tag to the input SAM file and
    writes the output to the specified outfile.
    '''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/add_XT.py --infile=%(infile)s --outname=%(outfile)s && samtools index %(outfile)s'''

    P.run(statement)


@active_if(PARAMS['correct'])
@follows(mkdir('counts_trans.dir'))
@transform(xt_tag,
           regex("mapped_files.dir/(\S+)_XT.bam"),
           r"counts_trans.dir/\1.counts.tsv.gz")
def count_trans(infile, outfile):
    '''Use umi-tools to collapse UMIs and generate counts table'''

    statement = '''umi_tools count --per-gene --gene-tag=XT --dual-nucleotide -I %(infile)s -S %(outfile)s'''

    P.run(statement)


@active_if(PARAMS['correct'])
@follows(mkdir('counts_trans_mclumi.dir'))
@transform(xt_tag,
           regex("mapped_files.dir/(\S+)_XT.bam"),
           r"counts_trans_mclumi.dir/\1.counts_mclumi.txt")
def count_trans_mclumi(infile, outfile):
    '''
    This function uses umi-tools to collapse UMIs and generate a counts
    table for transcripts. It takes an input BAM file and writes the
    output to the specified outfile.
    '''

    output_bam = infile.replace("_XT.bam", "_mclumi_XT.bam")

    statement = '''mclumi dedup_gene -m mcl_ed -gt XT -gist XS -ed %(mclumi_editdistance)s -ibam %(infile)s -obam %(output_bam)s -odsum %(outfile)s'''

    P.run(statement, job_memory=PARAMS['mclumi_memory'])


@active_if(PARAMS['correct'])
@transform(xt_tag,
           regex("mapped_files.dir/(\S+)_XT.bam"),
           r"counts_trans.dir/\1.counts_unique.tsv.gz")
def count_trans_unique(infile, outfile):
    '''
    This function uses mclumi to collapse UMIs and generate a
    counts table for transcripts. It takes an input BAM file
    and writes the output to the specified outfile.
    '''

    statement = '''umi_tools count --per-gene --method=unique --gene-tag=XT -I %(infile)s -S %(outfile)s'''

    P.run(statement)


@active_if(PARAMS['correct'])
@transform(xt_tag,
           regex("mapped_files.dir/(\S+)_XT.bam"),
           r"counts_trans.dir/\1.counts_noumis.tsv.gz")
def count_trans_noumis(infile, outfile):
    '''
    This function uses umi-tools to count unique UMIs for each
    transcript in the input BAM file and writes the output to
    the specified outfile.
    '''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/trans_count.py --infile=%(infile)s --outfile=%(outfile)s'''

    P.run(statement)


@active_if(PARAMS['correct'])
@merge(count_trans, "counts_trans.dir/counts.tsv.gz")
def merge_count(infiles, outfile):
    '''
    This function merges counts from each sample into one file.
    It takes a list of input files and writes the merged output
    to the specified outfile.
    '''

    df = merge_feature_data(infiles)
    df = df.fillna(0)
    df.to_csv(outfile, sep="\t", compression="gzip")


@active_if(PARAMS['correct'])
@merge(count_trans_unique, "counts_trans.dir/counts_unique.tsv.gz")
def merge_count_unique(infiles, outfile):
    '''
    This function merges unique counts from each sample into one
    file. It takes a list of input files and writes the merged
    output to the specified outfile.
    '''

    df = merge_feature_data(infiles)
    df = df.fillna(0)
    df.to_csv(outfile, sep="\t", compression="gzip")


@active_if(PARAMS['correct'])
@merge(count_trans_noumis, "counts_trans.dir/counts_noumis.tsv.gz")
def merge_trans_noumi(infiles, outfile):
    '''
    This function merges counts from each sample into one file without
    taking into account the UMI. It takes a list of input files
    and writes the merged output to the specified outfile.
    '''

    df = merge_trans_noumi_data(infiles)
    df = df.fillna(0)
    df.to_csv(outfile, sep="\t", compression="gzip")



#############################
## Gene level analysis ######
#############################

@active_if(PARAMS['correct'])
@transform(tso_umi,
           regex("processed_fastq.dir/(\S+)_polya_tso_UMI.fastq.gz"),
           r"mapped_files.dir/\1_gene.sam")
def mapping_gene(infile, outfile):
    '''
    This function maps the genes using minimap2
    and writes the output to the specified outfile.
    '''

    name = infile.replace('processed_fastq.dir/','')
    name = name.replace('_polya_tso_UMI.fastq.gz','_tmp_sam')

    if PARAMS['minimap2_splitprefix']:
        statement = """minimap2 -ax splice  -k 14 --split-prefix %(name)s --sam-hit-only --secondary=no --junc-bed %(junc_bed)s %(genome_fasta)s %(infile)s > %(outfile)s  2> %(outfile)s.log"""

    else:
        statement = '''minimap2 -ax splice  -k 14 --sam-hit-only --secondary=no --junc-bed %(junc_bed)s %(genome_fasta)s %(infile)s > %(outfile)s  2> %(outfile)s.log'''

    P.run(statement, job_memory="60G")


@active_if(PARAMS['correct'])
@transform(mapping_gene,
           regex("mapped_files.dir/(\S+)_gene.sam"),
           r"mapped_files.dir/\1_gene_sorted.bam")
def samtools_sort(infile, outfile):
    '''
    This function sorts the input BAM file using samtools and
    writes the output to the specified outfile.
    '''

    name = infile.replace("_gene.sam", "")

    statement = '''cgat bam2bam --method=strip-sequence -L strip.log -I %(infile)s -S %(name)s_gene_strip.sam &&
                   samtools view -bh %(name)s_gene_strip.sam > %(name)s_gene.bam &&
                   samtools sort %(name)s_gene.bam -o %(name)s_gene_sorted.bam &&
                   samtools index %(name)s_gene_sorted.bam'''

    P.run(statement)


@active_if(PARAMS['correct'])
@transform(samtools_sort,
           regex("mapped_files.dir/(\S+)_gene_sorted.bam"),
           r"mapped_files.dir/\1_featurecounts_gene_sorted.bam")
def featurecounts(infile, outfile):
    '''
    This function runs featurecounts on the input BAM file
    and writes the output to the specified outfile.
    '''

    name = infile.replace("_gene_sorted.bam", "")
    statement = '''featureCounts -a %(gtf)s -o %(name)s_gene_assigned.txt -R BAM %(infile)s &&
                   samtools sort %(infile)s.featureCounts.bam -o %(name)s_featurecounts_gene_sorted.bam &&
                   samtools index %(name)s_featurecounts_gene_sorted.bam'''

    P.run(statement)


@active_if(PARAMS['correct'])
@follows(mkdir('counts_genes.dir'))
@transform(featurecounts,
           regex("mapped_files.dir/(\S+)_featurecounts_gene_sorted.bam"),
           r"counts_genes.dir/\1.counts_gene.tsv.gz")
def count_gene(infile, outfile):
    '''
    This function uses umi-tools to collapse UMIs and generate a counts
    table for genes. It takes an input BAM file and writes the output
    to the specified outfile.
    '''

    statement = '''umi_tools count --per-gene --gene-tag=XT -I %(infile)s -S %(outfile)s'''

    P.run(statement)


@active_if(PARAMS['correct'])
@follows(mkdir('counts_genes_mclumi.dir'))
@transform(featurecounts,
           regex("mapped_files.dir/(\S+)_featurecounts_gene_sorted.bam"),
           r"counts_genes_mclumi.dir/\1.counts_gene.txt")
def count_gene_mclumi(infile, outfile):
    '''Use mclumi to count genes'''

    output_bam = infile.replace("featurecounts_gene_sorted.bam", "_mclumi_featurecounts_gene_sorted.bam")

    statement = '''mclumi dedup_gene -m mcl_ed -gt XT -gist XS -ed %(mclumi_editdistance)s -ibam %(infile)s -obam %(output_bam)s -odsum %(outfile)s'''

    P.run(statement, job_memory=PARAMS['mclumi_memory'])


@active_if(PARAMS['correct'])
@merge(count_gene, "counts_genes.dir/counts_gene.tsv.gz")
def merge_count_gene(infiles, outfile):
    '''merge counts from ech sample into one'''

    df = merge_feature_data(infiles)
    df = df.fillna(0)
    df.to_csv(outfile, sep="\t", compression="gzip")


@active_if(PARAMS['correct'])
@transform(featurecounts,
           regex("mapped_files.dir/(\S+)_featurecounts_gene_sorted.bam"),
           r"counts_genes.dir/\1.count_gene_unique.tsv.gz")
def count_gene_unique(infile, outfile):
    '''Use umi-tools to collapse UMIs and generate counts table'''

    statement = '''umi_tools count --per-gene --gene-tag=XT --method=unique -I %(infile)s -S %(outfile)s'''

    P.run(statement)


@active_if(PARAMS['correct'])
@merge(count_gene_unique, "counts_genes.dir/gene_counts_unique.tsv.gz")
def merge_count_gene_unique(infiles, outfile):
    '''merge counts from ech sample into one'''

    df = merge_feature_data(infiles)
    df = df.fillna(0)
    df.to_csv(outfile, sep="\t", compression="gzip")


########
### analyse without UMI sequences
#########

@active_if(PARAMS['correct'])
@follows(featurecounts)
@originate("counts_gene_noumis.tsv.gz")
def merge_featurecounts(outfile):
    ''' '''

    infiles = glob.glob("mapped_files.dir/*_gene_assigned.txt")
    final_df = merge_featurecounts_data(infiles)
    names = [x.replace("_gene_assigned.txt", "") for x in infiles]
    final_df.columns = names
    df = final_df.fillna(0)
    df.to_csv(outfile, sep="\t", compression="gzip")



########
### Correct using greedy algorithm
#########


@active_if(PARAMS['correct'])
@follows(mkdir('counts_trans_greedy.dir'))
@transform(xt_tag,
           regex("mapped_files.dir/(\S+)_XT.bam"),
           r"counts_trans_greedy.dir/\1.counts_greedy.txt")
def count_trans_greedy(infile, outfile):
    '''Use greedy to collapse UMIs'''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")
    
    statement = '''python %(PYTHON_ROOT)s/greedy_bulk.py count -i %(infile)s
                   -t XT -o %(outfile)s'''

    P.run(statement)


########
### Analyse with no UMI added to the read (ONTs workflow)
########

@active_if(PARAMS['no_umi'])
@follows(mkdir('counts_trans.dir'))
@transform(samtools,
           regex("mapped_files.dir/(\S+)_final_sorted.bam"),
           r"counts_trans.dir/\1_counts_noumi.tsv")
def count_trans_noumi(infile, outfile):
    '''Count the reads for every gene'''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python  %(PYTHON_ROOT)s/trans_count.py --infile=%(infile)s --outfile=%(outfile)s'''

    P.run(statement)


@active_if(PARAMS['no_umi'])
@follows(count_trans_noumi)
@originate("counts_trans_noumis.tsv.gz")
def merge_noumis(outfile):
    ''' '''

    infiles = glob.glob("counts_trans.dir/*_counts_noumi.tsv")
    final_df = merge_featurecounts_data(infiles)
    names = [x.replace("_counts_noumi.tsv", "") for x in infiles]
    names = [x.replace("counts_trans.dir/", "") for x in names]
    final_df.columns = names
    df = final_df.fillna(0)
    df.to_csv(outfile, sep="\t", compression="gzip")


@follows(merge_count, merge_count_unique, merge_count_gene, merge_count_gene_unique, merge_featurecounts, count_trans_mclumi, count_gene_mclumi, count_trans_greedy, merge_noumis)
def full():
    '''
    A placeholder function that serves as a checkpoint
    to run all previous ruffus tasks and ensure that all
    previous tasks are completed.
    '''
    pass


def main(argv=None):
    '''
    The main function that runs the pipeline using the cgatcore.pipeline module.
    Takes an optional argument list (default is sys.argv).

    Please note that some of these functions use external Python scripts or
    tools. For a complete understanding of their functionality, it is
    necessary to examine the code of those scripts as well.
    '''
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))    
