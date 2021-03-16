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
L = logging.getLogger("naive_collapse.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--infile", default=None, type=str,
                    help='nanopore infile fastq  file')
parser.add_argument("--outname", default=None, type=str,
                    help='name for output fastq files')

args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #



log =  iotools.open_file(args.outname + "_naive_collapse" + ".log","w")





read1 = iotools.open_file(args.outname + "_R1.fastq","w")
read2 = iotools.open_file(args.outname + "_R2.fastq","w")



a = 0
ua = 0
total = 0
# generate set of barcodes for whitelist
barcodes = []
with pysam.FastxFile(args.infile) as fh:
    for record in fh:
        

        
        
        seq_nano = record.sequence
        total += 1

        m=regex.finditer("(AAGCAGTGGTATCAACGCAGAGT){e<=4}", str(seq_nano))


        for i in m:
            a += 1
            
            new_seq = seq_nano[i.end():]
            new_seq_quality = record.quality[i.end():]
            
            m=regex.finditer("(TTTTTTTTTTTTTTTTTTTT){e<=3}", str(seq_nano))
            for i in m:
                read2_seq = seq_nano[i.start():]
                read2_seq_quality = record.quality[i.start():]

            barcode = new_seq[2:26] 
            barcode_quality = new_seq_quality[2:26] 
            umi = new_seq[26:42]
            umi_quality = new_seq_quality[26:42]
            
            barcode_umi = barcode + umi
            barcode_umi_quality = barcode_quality + umi_quality
            barcode_umi = barcode_umi[::2]
            barcode_umi_quality = barcode_umi_quality[::2]
            
            umi = umi[::2]

            barcode = barcode[::2]

            record_name = record.name + "_" + barcode + "_" + umi
            
            if len(barcode) == 12 and len(umi) == 8:
                read1.write("@%s\n%s\n+\n%s\n" % (record_name, barcode_umi, barcode_umi_quality))
                read2.write("@%s\n%s\n+\n%s\n" % (record_name, read2_seq, read2_seq_quality))
            else:
                pass
            
read1.close()
read2.close()


log.write("Total reads: %s\n" %(total))
log.write("Total reads with barcode in: %s\n" %(a))

log.close()
