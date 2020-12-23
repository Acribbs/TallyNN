import sys
import regex
import cgatcore.iotools as iotools
import pysam
import logging
import argparse
import re



# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("complement_ployA.py")

# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--infile", default=None, type=str,
                    help='input bam  file')
parser.add_argument("--bedfile", default=None, type=str,
                    help='bedfile  file')
parser.add_argument("--outfile", default=None, type=str,
                    help='name for output bam file')

args = parser.parse_args()

L.info("args:")
print(args)



# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #


bamfile = pysam.AlignmentFile(args.infile)
tabixfile = pysam.TabixFile(args.bedfile)
out_bam = pysam.AlignmentFile(args.outfile, "wb", template=bamfile)

for line in bamfile:
    
    chrom = line.reference_name
    start = line.reference_start
    end = line.reference_end
    if "chrUn" in chrom or "alt" in chrom or "KI" in chrom or "GL" in chrom:
                pass
    else:
        for gtf in tabixfile.fetch(chrom, int(start), int(end)):
            gtf = gtf.split("\t")
            line.tags += [("XT", gtf[3])]
            out_bam.write(line)

out_bam.close()
bamfile.close() 
