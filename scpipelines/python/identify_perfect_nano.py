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
L = logging.getLogger("identify_perfect_nano.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--whitelist", default=None, type=str,
                    help='a file naming the outfile for the whitelist of barcodes.')
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



log =  iotools.open_file(args.outname + "_perfect_nano" + ".log","w")



def count_pairs(s):
    pairs_cnt = 0
    unique_chars = set(s)
    for char in unique_chars:
        pairs_cnt += s.count(char + char)
    return pairs_cnt

read1 = iotools.open_file(args.outname + "_unambiguous_barcode_R1.fastq","w")
read2 = iotools.open_file(args.outname + "_unambiguous_barcode_R2.fastq","w")

read1_no = iotools.open_file(args.outname + "_ambiguous_barcode_R1.fastq","w")
read2_no = iotools.open_file(args.outname + "_ambiguous_barcode_R2.fastq","w")


a = 0
ua = 0
total = 0
# generate set of barcodes for whitelist
barcodes = []
with pysam.FastxFile(args.infile) as fh:
    for record in fh:
        

        
        
        seq_nano = record.sequence
        total += 1

        m=regex.finditer("(AAGCAGTGGTATCAACGCAGAGT){e<=3}", str(seq_nano))


        for i in m:
            
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

            if count_pairs(barcode) == 12:
                barcode_umi = barcode + umi
                barcodes.append(barcode)
                
                barcode_umi_quality = barcode_quality + umi_quality

                if len(barcode_umi) == 40:
                    ua += 1
                    read1.write("@%s\n%s\n+\n%s\n" % (record.name, barcode_umi, barcode_umi_quality))
                    read2.write("@%s\n%s\n+\n%s\n" % (record.name, read2_seq, read2_seq_quality))
                else:
                    pass
            else:
                if len(barcode_umi) == 40: 
                    a += 1
                    read1_no.write("@%s\n%s\n+\n%s\n" % (record.name, barcode_umi, barcode_umi_quality))
                    read2_no.write("@%s\n%s\n+\n%s\n" % (record.name, read2_seq, read2_seq_quality))
                else:
                    pass

                
    # Write out a list of whitelist barcodes
    out_barcodes = open(args.whitelist,"w")
    y = 0
    for i in set(barcodes):
        y += 1
        out_barcodes.write("%s\n" % (i))
    out_barcodes.close()
        
read1.close()
read2.close()
read1_no.close()
read2_no.close()


log.write("The number of unambiguous barcodes in whitelist identified is: %s\n" %(y))
log.write("The number of unambiguous barcode reads identified is: %s\n" %(ua))
log.write("The number of ambiguous barcode reads identified is: %s\n" %(a))
log.write("The total number of barcode reads identified is: %s\n" %(total))
log.write("The percent of unambiguous barcode reads : %s\n" %((ua/total)*100))
log.write("The percent of ambiguous barcode reads : %s\n" %((a/total)*100))

log.close()
