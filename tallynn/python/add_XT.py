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
L = logging.getLogger("tso_umi.py")

# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--infile", default=None, type=str,
                    help='infile fastq  file')
parser.add_argument("--outname", default=None, type=str,
                    help='name for output fastq files')

args = parser.parse_args()

L.info("args:")
print(args)



# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #

samfile = pysam.AlignmentFile(args.infile, "rb")
outfile = pysam.AlignmentFile(args.outname, "wb", template=samfile)

for read in samfile:

    if read.reference_name is not None:
        read.tags += [('XT',read.reference_name)]
    
    else:
        pass

    outfile.write(read)

samfile.close()
outfile.close()

