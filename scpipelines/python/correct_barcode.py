import sys
import regex
import cgatcore.iotools as iotools
import pysam
import Levenshtein
import logging
import argparse

# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("correct_barcode.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--whitelist", default=None, type=str,
                    help='a file barcodes extracted using umi whitelist')
parser.add_argument("--read1", default=None, type=str,
                    help='read1 fastq  file')
parser.add_argument("--read2", default=None, type=str,
                    help='read2 fastq file')
parser.add_argument("--outname", default=None, type=str,
                    help='name for output fastq files')
parser.add_argument("--distance", default=None, type=str,
                    help='levenshtein distance')

args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ######################## Import whitelist    ############################## #
# ########################################################################### #

whitelist = iotools.open_file(args.whitelist, "r")
barcodes = []

for line in whitelist:
    barcodes.append(line.split("\t")[0])


# ########################################################################### #
# ######################## Apply Levenstein distance ######################## #
# ########################################################################### #

outf = iotools.open_file("corrected_reads.dir/" + args.outname + "_corrected.fastq.1.gz","w")
outf2 = iotools.open_file("corrected_reads.dir/" + args.outname + "_corrected.fastq.2.gz","w")


with pysam.FastxFile(args.read1) as fh, pysam.FastxFile(args.read2) as fh2:
    
    n = 0
    y = 0
    for record_fh, record_fh2  in zip(fh, fh2):
        barcode = record_fh.sequence[0:24]

        y += 1

        for b in barcodes:

            

            if Levenshtein.distance(barcode, b) <= int(args.distance):
                n +=1
                b = b + record_fh.sequence[24:]

                outf.write("@%s\n%s\n+\n%s\n" % (record_fh.name, b, record_fh.quality))
                outf2.write("@%s\n%s\n+\n%s\n" % (record_fh2.name, record_fh2.sequence, record_fh2.quality))
                break
                
            else:
                pass

outf.close()
outf2.close()
whitelist.close()
