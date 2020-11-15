import sys
import regex
import cgatcore.iotools as iotools
import pysam
import logging
import argparse
import Levenshtein


# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("correct_barcode_nano.py")


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
# ######################## Code                ############################## #
# ########################################################################### #


outf =  iotools.open_file(args.outname + "unambiguous_fixed_barcode_R1.fastq","w")
outf2 = iotools.open_file(args.outname + "unambiguous_fixed_barcode_R2.fastq","w")
log =  iotools.open_file(args.outname + ".log","w")


with pysam.FastxFile("ambiguous_barcode_R1.fastq") as fh, pysam.FastxFile("ambiguous_barcode_R2.fastq") as fh2:
    
    n = 0
    y = 0
    for record_fh, record_fh2  in zip(fh, fh2):
        barcode = record_fh.sequence[0:24]
        
        odd_bc = barcode[1::2]
        even_bc = barcode[::2]
        y += 1
        for b in barcodes:

            

            if Levenshtein.distance(odd_bc, b) <= args.distance:
                n +=1
                b = b + record_fh.sequence[24:]
                
                outf.write("@%s\n%s\n+\n%s\n" % (record_fh.name, b, record_fh.quality[12:]))
                outf2.write("@%s\n%s\n+\n%s\n" % (record_fh2.name, record_fh2.sequence, record_fh2.quality))
                break
                
            elif Levenshtein.distance(even_bc, b) <= args.distance:
                n +=1
                b = b + record_fh.sequence[24:]
                outf.write("@%s\n%s\n+\n%s\n" % (record_fh.name, b, record_fh.quality[12:]))
                outf2.write("@%s\n%s\n+\n%s\n" % (record_fh2.name, record_fh2.sequence, record_fh2.quality))
                break
                
            else:
                pass


log.write("The number of total reads with levenshtein distance less than %s: %s\n" %(args.distance, n))
log.write("The number of total reads is: %s\n" %(y))
log.write("The number of total recovered percent is: %s\n" %((y/n)*100))

log.close()
outf.close()
outf2.close()
