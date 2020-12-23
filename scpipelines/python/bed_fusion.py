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
parser.add_argument("--bed1", default=None, type=str,
                    help='bedfile1 out  file')
parser.add_argument("--bed2", default=None, type=str,
                    help='bedfile2 out  file')

args = parser.parse_args()

L.info("args:")
print(args)



# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #


bamfile = pysam.AlignmentFile(args.infile)
bed1 = open(args.bed1, "w")
bed2 = open(args.bed2, "w")

for line in bamfile:

    barcode = line.query_name.split("_")[1]
    umi = line.query_name.split("_")[1]
    barcode_umi = barcode+umi
    
    chrom1 = line.get_tag("Ta")
    start1 = line.get_tag("Tb")
    end1 = line.get_tag("Tc")
    name1 = line.get_tag("Td")
    
    chrom2 = line.reference_name
    start2 = line.reference_start
    end2 = line.reference_end
    name2 = line.get_tag("XT")
    

    bed1.write("%s\t%s\t%s\t%s\t%s\n" % (chrom1, start1, end1, name1, barcode_umi))
    bed2.write("%s\t%s\t%s\t%s\t%s\n" % (chrom2, start2, end2, name2, barcode_umi))
