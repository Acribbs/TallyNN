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
