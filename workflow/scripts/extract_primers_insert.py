import pysam
import sys

read_dict = dict()
bam_file = pysam.AlignmentFile(snakemake.input[0], "rb")
with open(snakemake.output[0], "w") as out:
    for read in bam_file.fetch():
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            read_dict[read.query_name] = read.query_sequence
        else:
            if read.is_read1:
                primer1 = read.query_sequence
                primer2 = read_dict.pop(qname, "")
            else:
                primer1 = read_dict.pop(qname, "")
                primer2 = read.query_sequence
            insert_size = abs(read.template_length) + len(primer2)
            print("{}\t{}\t{}".format(primer1, primer2, insert_size), file=out)

for read in read_dict.keys():
    print("{} has no mate".format(read.query_name), file=sys.stderr)
