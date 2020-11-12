import sys
import regex
import cgatcore.iotools as iotools
import pysam
import logging
import argparse

# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("collapse_bcumi.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--read1", default=None, type=str,
                    help='read1 fastq  file')
parser.add_argument("--read2", default=None, type=str,
                    help='read2 fastq file')
parser.add_argument("--outname", default=None, type=str,
                    help='name for output fastq files')
args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ######################## Identify perfect barcodde    ##################### #
# ########################################################################### #


outf = iotools.open_file(args.outname + ".fastq.1.gz","w")
outf2 = iotools.open_file(args.outname + ".fastq.2.gz","w")


with pysam.FastxFile(args.read1) as fh, pysam.FastxFile(args.read2) as fh2:

    for record_fh, record_fh2  in zip(fh, fh2):
        
        
        sequence = record_fh.sequence[::2]
        quality = record_fh.quality[::2]
        
        outf.write("@%s\n%s\n+\n%s\n" % (record_fh.name, sequence, quality))
        outf2.write("@%s\n%s\n+\n%s\n" % (record_fh2.name, record_fh2.sequence, record_fh2.quality))
        
        
outf.close()
outf2.close()
