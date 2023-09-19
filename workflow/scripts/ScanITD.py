#!/usr/bin/env python3
#-*- coding: utf-8 -*-
#===============================================================================
# the modifications of this script (original ScanITD https://github.com/ylab-hi/ScanITD)
# can be found in this repository: https://github.com/dawidkrzeciesa/ScanITD
#===============================================================================
__version__ = 0.3

import sys
import os
import time
import subprocess
import re
import argparse
from pyfaidx import Fasta

try:
    import pysam
except:
    sys.exit('pysam module not found.\nPlease install it.')
try:
    import numpy as np
except:
    sys.exit('numpy module not found.\nPlease install it.')
try:
    import skbio
except:
    sys.exit('scikit-bio module not found.\nPlease install it.')


def detect_softclip_mode(cigarstring):
    # if cigarstring is 40M25N5M then cigartuple is [('40', 'M'), ('25', 'N'), ('5', 'M')] with the statement below
    # mode 1 = MS; 2=SM; 3=other
    if 'I' in cigarstring or 'D' in cigarstring or 'H' in cigarstring:
        return 3, -1, -1
    cigartuple = re.findall(r'(\d+)(\w)', cigarstring)
    if cigartuple[0][1] == 'S':
        if cigartuple[-1][1] == 'S':
            if int(cigartuple[0][0]) > int(cigartuple[-1][0]):
                return 2, int(cigartuple[1][0]), int(cigartuple[-1][0])
            else:
                return 1, int(cigartuple[1][0]), int(cigartuple[0][0])
        else:
            return 2, int(cigartuple[1][0]), 0
    else:
        return 1, int(cigartuple[0][0]), 0


def detect_itd_from_cigar(chr, read, mapq_cutoff):
    strand = '-' if read.is_reverse else '+'
    try:
        chimeric_aln = read.get_tag('SA').split(';')[:-1]
    except KeyError:
        return 'NA', []
    if len(chimeric_aln) == 1:
        chr_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = read.get_tag(
            'SA').split(',')
        if int(mapq_sa) < mapq_cutoff:
            return 'NA', []
    else:
        return 'NA', []
    read_mode, read_match, read_clipped = detect_softclip_mode(
        read.cigarstring)
    sa_mode, sa_match, sa_clipped = detect_softclip_mode(cigar_sa)
    if read_mode == 3 or sa_mode == 3:
        return 'NA', []
    target_start = 0
    target_end = 0
    newcigar = []
    newpos = -1
    if chr == chr_sa:
        if strand_sa == strand:  # deletion, insertion, duplication
            if read_mode == 2 and sa_mode == 1:
                target_start = int(pos_sa) - 1
                target_end = read.aend
                if sa_clipped != 0:
                    newcigar.append([4, sa_clipped])
                newcigar.append([0, sa_match])
                target_offset = target_end - target_start
                query_offset = read.rlen - sa_clipped - read_clipped
                indel_size = query_offset - target_offset
                newpos = target_start
                if indel_size <= 0:  # deletion
                    return 'NA', []
                elif indel_size >= query_offset:  # large ITD
                    return 'TDUP', [read.pos, indel_size, 2, 1]
                else:  # small ITDs
                    newcigar.append([1, indel_size])
                    last_match_size = read.aend - (int(pos_sa) - 1 + sa_match)
                    if last_match_size <= 0:
                        for each in newcigar:
                            if each[0] == 0:
                                each[1] = each[1] + last_match_size
                                break
                    else:
                        newcigar.append([0, last_match_size])
                if read_clipped != 0:
                    newcigar.append([4, read_clipped])
            elif read_mode == 1 and sa_mode == 2:
                target_start = read.pos
                target_end = int(pos_sa) - 1 + sa_match
                if read_clipped != 0:
                    newcigar.append([4, read_clipped])
                newcigar.append([0, read_match])
                target_offset = target_end - target_start
                query_offset = read.rlen - sa_clipped - read_clipped
                indel_size = query_offset - target_offset
                newpos = target_start
                if indel_size <= 0:
                    return 'NA', []
                elif indel_size >= query_offset:
                    return 'TDUP', [int(pos_sa) - 1, indel_size, 1, 2]
                else:
                    newcigar.append([1, indel_size])
                    last_match_size = int(pos_sa) - 1 + sa_match - read.aend
                    if last_match_size <= 0:
                        for each in newcigar:
                            if each[0] == 0:
                                each[1] = each[1] + last_match_size
                                break
                    else:
                        newcigar.append([0, last_match_size])
                if sa_clipped != 0:
                    newcigar.append([4, sa_clipped])
            else:
                return 'NA', []
        else:
            return 'NA', []
    else:
        return 'NA', []
    return list(map(tuple, newcigar)), newpos


def softclipping_realignment(mapq_cutoff, input, output):
    bwa_bam = pysam.AlignmentFile(input, 'rb')
    output_bam = pysam.AlignmentFile(output + '.temp.bam', 'wb', template=bwa_bam)
    try:
        for read in bwa_bam.fetch(until_eof=True):
            if read.mapq >= mapq_cutoff and not read.is_secondary and not read.has_tag(
                    'XA'):
                chr = bwa_bam.getrname(read.rname)
                newcigar, newpos = detect_itd_from_cigar(
                    chr, read, mapq_cutoff)
                if newcigar == 'TDUP':
                    read.setTag('SV', newcigar +
                                ',' +
                                str(newpos[0] +
                                    1) +
                                ',' +
                                str(newpos[1]) +
                                ',' +
                                str(newpos[2]) +
                                ',' +
                                str(newpos[3]))
                elif newcigar != 'NA' and newcigar != read.cigar:
                    read.setTag('OA', str(read.pos + 1) +
                                ',' + read.cigarstring)
                    read.cigar, read.pos = newcigar, newpos
            output_bam.write(read)
    except ValueError as e:
        print('Bam index file is not found!', e, file=sys.stderr)
        sys.exit(1)
    bwa_bam.close()
    output_bam.close()
    pysam.sort('-o', output, '{0}.temp.bam'.format(output))
    pysam.index(output)

    '''
    try:
        subprocess.check_call(
            "samtools sort {0}.temp.bam -o {0}".format(output),
            shell=True)
    except subprocess.CalledProcessError as e:
        print('Execution failed for samtools:', e, file=sys.stderr)
        sys.exit(1)

    subprocess.check_call("samtools index {}".format(output), shell=True)
    '''


def version_checking(samtools_path):
    '''HTSlib-based Samtools is needed!
    '''
    try:
        response = subprocess.check_output('{} --version-only'.format(samtools_path), shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print >> sys.stderr, "HTSlib-based Samtools is not available, please install latest Samtools at http://www.htslib.org/"
        print >> sys.stderr, "Exiting."
        sys.exit(0)


def external_tool_checking():
    """checking dependencies are installed"""
    software = ['samtools']
    cmd = "which"
    for each in software:
        try:
            path = subprocess.check_output([cmd, each], stderr=subprocess.STDOUT)
            path = str(path, 'utf-8')
        except subprocess.CalledProcessError:
            print("Checking for {0}: ERROR - could not find {0}".format(each), file=sys.stderr)
            print("Exiting.", file=sys.stderr)
            sys.exit(0)
        print("Checking for '" + each + "': found " + path)
    return path.rstrip()


def remove(infile):
    if os.path.isfile(infile):
        os.remove(infile)


def vcf_header(output_prefix):
    header = ['##fileformat=VCFv4.1']
    header.append('##source=ScanITD {}'.format(__version__))
    header.append('##ALT=<ID=TDUP,Description="Tandem duplication">')
    header.append('##ALT=<ID=INS,Description="Insertion">')
    header.append('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">')
    header.append('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">')
    header.append('##INFO=<ID=AO,Number=1,Type=Integer,Description="Alternate allele observations, with partial observations recorded fractionally">')
    header.append('##INFO=<ID=AB,Number=1,Type=Float,Description="Estimated allele frequency in the range (0,1], representing the ratio of reads showing the alternative allele to all reads">')
    header.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="The type of event. The type was restricted to INS only to use Varlociraptor variant caller">')
    header.append('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">')
    header.append('##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">')
    header.append('##INFO=<ID=END,Number=1,Type=Integer,Description="nd position of the structural variant">')
    header.append('##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">')
    header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+output_prefix)
    return '\n'.join(header)


def mismatch_count(sr_pos, sr_seq, itd_id, itd_seq,
                   mismatch_cutoff, ref_fa, bamfile):
    items = itd_id.split(',')
    itd_type = items[0]
    itd_pos = int(items[1])
    itd_len = int(items[2])
    sr_mode = sr_pos[0]
    ref_file = pysam.Fastafile(ref_fa)
    if sr_mode == 2:
        seq_target_from_read = skbio.DNA(itd_seq[1])
        seq_target_start = sr_pos[-1]
        seq_target_end = sr_pos[-1] + itd_len
        cov_pos = seq_target_start
    else:
        seq_target_from_read = skbio.DNA(itd_seq[0])
        seq_target_start = sr_pos[-1] - itd_len + 1
        if seq_target_start < 0:
            seq_target_start = 0
        seq_target_end = sr_pos[-1] + 1
        cov_pos = seq_target_end - 1

    genome_seq = Fasta(ref_fa, sequence_always_upper=True, as_raw=True)
    seq_target_from_genome = skbio.DNA( str(genome_seq[sr_pos[1]][seq_target_start:seq_target_end]) )

    ao = 0

    dp_line = pysam.depth('{0}'.format(bamfile), '-q', '15', '-r', '{0}:{1}-{1}'.format(sr_pos[1], cov_pos+1))
    dp = int(dp_line.rstrip().split('\t')[-1])

    for each_seq in sr_seq:
        seq_query = skbio.DNA(each_seq)
        try:
            alignment1, score1, start_end_pos1 = skbio.alignment.local_pairwise_align_ssw(
                seq_query, seq_target_from_genome)
            alignment2, score2, start_end_pos2 = skbio.alignment.local_pairwise_align_ssw(
                seq_query, seq_target_from_read)
        except IndexError:
            continue
        if sr_mode == 1:  # MS mode
            if (start_end_pos1[0][0] == 0 and start_end_pos1[1][0] == 0 and sum(alignment1[0].mismatches(alignment1[1])) <= mismatch_cutoff) or (
                    start_end_pos2[0][0] == 0 and start_end_pos2[1][0] == 0 and sum(alignment2[0].mismatches(alignment2[1])) <= mismatch_cutoff):
                ao += 1
        else:  # SM mode
            if (start_end_pos1[0][1] == len(seq_query) - 1 and start_end_pos1[1][1] == len(seq_target_from_genome) - 1 and sum(alignment1[0].mismatches(alignment1[1])) <= mismatch_cutoff) or (
                    start_end_pos2[0][1] == len(seq_query) - 1 and start_end_pos2[1][1] == len(seq_target_from_read) - 1 and sum(alignment2[0].mismatches(alignment2[1])) <= mismatch_cutoff):
                ao += 1
    return ao, dp, seq_target_start


def is_itd(read, ref_site, mapq_cutoff):
    if read.alignment.has_tag('SV'):
        chr_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = read.alignment.get_tag(
            'SA').split(',')
        if int(mapq_sa) < mapq_cutoff:
            return False
        itd_type, itd_pos, itd_length, read_mode, sa_mode = read.alignment.get_tag(
            'SV').split(',')
        if itd_type != 'TDUP':
            return False
        if int(itd_pos) - 2 + (int(read_mode) - 1) == ref_site:
            return True
        elif itd_type != 'TRA' and int(itd_pos) - 2 + int(itd_length) + (int(read_mode) - 1) == ref_site:
            return True
    return False


def get_softclip_info(read):
    if read.cigar[0][0] == 4:
        if read.cigar[-1][0] == 4:
            if read.cigar[0][1] > read.cigar[-1][1]:
                return 2, read.seq[:read.cigar[0][1]], read.pos
            else:
                return 1, read.seq[read.rlen -
                                   read.cigar[-1][1]:], read.aend - 1
        else:
            return 2, read.seq[:read.cigar[0][1]], read.pos
    elif read.cigar[-1][0] == 4:
        return 1, read.seq[read.rlen - read.cigar[-1][1]:], read.aend - 1
    else:
        return 0, '', -1

def mismatch_seq(s1, s2):
    if len(s1) != len(s2):
        sys.exit('Error in fetching insertion sequence!')
    else:
        count = 0
        for i,j in zip(s1, s2):
            if i != j:
                count +=1
        return count

def self_loop_checker(insertion_seq, left_seq, right_seq):
    ins_seq = insertion_seq
    ins_len = len(insertion_seq)
    if ins_len % 2 == 0:
        steps = int(ins_len/2)
    else:
        steps = int(ins_len/2) + 1
    # left rolling
    count = 1
    for i in range(steps):
        ins_seq = ins_seq[-1:] + ins_seq[0:len(ins_seq)-1]
        combo_seq = left_seq[-count:] + right_seq[:ins_len-count]
        count += 1
        if mismatch_seq(ins_seq, combo_seq) <= 5:
            return True
    # right rolling
    ins_seq = insertion_seq
    count = 1
    for i in range(steps):
        ins_seq = ins_seq[1:] + ins_seq[:1]
        combo_seq = left_seq[-(ins_len-count):] + right_seq[:count]
        count += 1
        if mismatch_seq(ins_seq, combo_seq) <= 5:
            return True
    # is a insertion of novel sequence
    return False

def itd_scan(input_bam, output_prefix, ao_cutoff, dp_cutoff, vaf_cutoff,
             itd_len_cutoff, target, mapq_cutoff, mismatch_cutoff, fasta_file):
    itd_set = set()
    samfile = pysam.AlignmentFile(input_bam, 'rb')

    genome_seq = Fasta(fasta_file, sequence_always_upper=True, as_raw=True)

    vcffile = open(output_prefix + '.ITD.vcf', 'w')
    print(vcf_header(output_prefix), file=vcffile)
    vcf_field_gt = 'GT\t0/1'
    regions = []
    try:
        ifp = open(target)
        for line in ifp:
            items = line.rstrip().split()
            regions.append(items[0] + ':' + items[1] + '-' + items[2])
    except IOError:
        if target == '':
            regions = [None]
        else:
            regions = [target]
    for each_region in regions:
        try:
            for col in samfile.pileup(region=each_region, stepper='all', truncate=True):
                dp = col.n
                itd_dict = {}
                sr_dict = {}
                itd_seq_dict = {}
                for read in col.pileups:
                    if read.alignment.mapq >= mapq_cutoff:
                        if 'S' in read.alignment.cigarstring and not read.alignment.has_tag(
                                'SV'):
                            soft_mode, soft_seq, soft_pos = get_softclip_info(
                                read.alignment)
                            if 'N' not in soft_seq:
                                try:
                                    sr_dict[(soft_mode, col.reference_name, soft_pos)].append(
                                        soft_seq)
                                except BaseException:
                                    sr_dict[(soft_mode, col.reference_name, soft_pos)] = [
                                        soft_seq]
                        if read.indel >= itd_len_cutoff:
                            saving_flag = False
                            if re.search(r'\d+M\d+I\d+M', str(read.alignment.cigarstring)):
                                left_seq_from_genome = genome_seq[col.reference_name][(col.pos-read.indel+2):col.pos+1]
                                right_seq_from_genome = genome_seq[col.reference_name][col.pos+1:(col.pos+read.indel)]
                                insertion_seq_in_read = read.alignment.query_sequence[read.query_position+1:(read.query_position + read.indel+1)]
                                if self_loop_checker(insertion_seq_in_read, left_seq_from_genome, right_seq_from_genome):
                                    saving_flag = True
                            else:
                                saving_flag = True
                            if saving_flag:
                                itd_id = 'INS,' + \
                                    str(col.pos + 1) + ',' + str(read.indel) + ',1'
                                itd_seq_dict[itd_id] = [read.alignment.query_sequence[read.query_position + 1:],
                                          read.alignment.query_sequence[:read.query_position + read.indel]]
                                try:
                                    itd_dict[itd_id][-1] = itd_dict[itd_id][-1] + 1
                                except KeyError:
                                    ref_seq = read.alignment.query_sequence[read.query_position]
                                    alt_seq = read.alignment.query_sequence[
                                        read.query_position:read.query_position + read.indel + 1]
                                    itd_dict[itd_id] = [ref_seq, alt_seq, 1]
                        elif is_itd(read, col.pos, mapq_cutoff):
                            itd_id = ','.join(
                                read.alignment.get_tag('SV').split(',')[:4])
                            itd_seq_dict[itd_id] = [
                                read.alignment.query_sequence[read.query_position + 1:], read.alignment.query_sequence[:read.query_position]]
                            try:
                                itd_dict[itd_id][-1] = itd_dict[itd_id][-1] + 1
                            except KeyError:
                                itd_dict[itd_id] = [
                                    '.', '<' + itd_id.split(',')[0] + '>', 1]
                for itd_id in itd_dict:  # only consider one event per genomic position
                    ao = itd_dict[itd_id][-1]
                    candidate_itd = {}
                    for sr_id in sr_dict:
                        new_ao, new_dp, new_pos = mismatch_count(
                            sr_id, sr_dict[sr_id], itd_id, itd_seq_dict[itd_id], mismatch_cutoff, fasta_file, input_bam)
                        if new_pos not in candidate_itd or new_ao > candidate_itd[new_pos][0]:
                            candidate_itd[new_pos] = [new_ao, new_dp]
                    try:
                        seed_itd_pos1 = int(itd_id.split(',')[1]) - 1
                        candidate_itd[seed_itd_pos1][0] += ao
                    except BaseException:
                        candidate_itd[seed_itd_pos1] = [ao, dp]
                    try:
                        seed_itd_pos2 = int(itd_id.split(
                            ',')[1]) - int(itd_id.split(',')[2])
                        candidate_itd[seed_itd_pos2][0] += ao
                    except BaseException:
                        candidate_itd[seed_itd_pos2] = [ao, dp]
                    sorted_itd = sorted(
                        candidate_itd.items(),
                        key=lambda kv: kv[1][0],
                        reverse=True)
                    ao = sorted_itd[0][1][0]
                    dp = sorted_itd[0][1][1]
                    vaf = float(0) if dp <= 0 else float(ao)/float(dp)
                    itd_pos = sorted_itd[0][0]
                    if dp >= dp_cutoff and ao >= ao_cutoff and vaf >= vaf_cutoff:
                        itd_type, itd_pos, itd_length, itd_mode = itd_id.split(
                            ',')
                        end = sorted_itd[0][0] + int(itd_length)
                        chr2 = col.reference_name
                        itd_length = itd_length
                        
                        # get ref and alt allele sequence
                        itd_length = int(itd_id.split(',')[2])
                        ref_allele_seq = str(genome_seq[col.reference_name][sorted_itd[0][0]])
                        alt_allele_seq = str(genome_seq[col.reference_name][sorted_itd[0][0] : sorted_itd[0][0] + itd_length])
                        
                        # varlociraptor call the tandem duplication events if SVTYPE=INS, therefore SVTYPE was changed from TDUP to INS
                        vcf_field_info = 'NS=1;AO=' + str(ao) + ';DP=' + str(dp) + ';AB=' + \
                            '{:.2g}'.format(vaf) + ';SVLEN=' + str(itd_length) + ';SVTYPE=INS;SVMETHOD=ScanITD_ALN;CHR2=' + chr2 + ';END=' + str(end)
                        
                        # fetch ref and alt allele sequence
                        out_itd = col.reference_name + '\t' + \
                            str(sorted_itd[0][0] + 1) + \
                            '\t.\t' + ref_allele_seq + '\t' + ref_allele_seq + alt_allele_seq + '\t.\t.\t' + vcf_field_info
                        if out_itd in itd_set:
                            continue
                        print(col.reference_name +
                              '\t' +
                              str(sorted_itd[0][0] + 1) +
                              '\t.\t' + ref_allele_seq + '\t' + ref_allele_seq + alt_allele_seq + '\t.\t.\t' +
                              vcf_field_info +
                              '\t' +
                              vcf_field_gt, file=vcffile)
                        itd_set.add(out_itd)
        except ValueError:
            continue
    samfile.close()
    return


def itd_len_type(x):
    x = int(x)
    if x == 0:
        raise argparse.ArgumentTypeError("Minimum ITD length has to be 1 base")
    return x


def parse_args():
    parser = argparse.ArgumentParser(description = "%(prog)s -i input_bam_file -r indexed_refenence_genome_fasta -o output_vcf_filename_prefix [opts]", epilog="ScanITD: detecting internal tandem duplication with robust variant allele frequency estimation")
    parser.add_argument('-i', '--input', action='store', dest='input', help="BWA-MEM BAM file", required=True)
    parser.add_argument('-r', '--ref', action='store', dest='ref', help="reference genome in FASTA format (with fai index)", required=True)
    parser.add_argument('-o', '--output', action='store', dest='output', help="output prefix", required=True)
    parser.add_argument('-m','--mapq', action='store', dest='mapq', type=int, help="minimal MAPQ in BAM for calling ITD (default: %(default)s)", default=15)
    parser.add_argument('-c','--ao', action='store', dest='ao', type=int, help="minimal observation count for ITD (default: %(default)s)", default=4)
    parser.add_argument('-d','--depth', action='store', dest='dp', type=int, help="minimal depth to call ITD (default: %(default)s)", default=10)
    parser.add_argument('-f', '--vaf', action='store', dest='vaf', type=float, help="minimal variant allele frequency (default: %(default)s)", default=0.1)
    parser.add_argument('-l', '--len', action='store', dest='itd_len', type=itd_len_type, help="minimal ITD length to report (default: %(default)s)", default=10)
    parser.add_argument('-n', action='store', dest='mismatch', type=int, help="maximum allowed mismatch bases of pairwise local alignment (default: %(default)s)", default=3)
    parser.add_argument('-t', '--target', action='store', dest='target', help="Limit analysis to targets listed in the BED-format file or a samtools region string", default='')
    parser.add_argument('-k','--keep', action='store_true', dest='keep', help="Kepp the ITD build BAM file")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__))

    return parser


def main():

    if sys.version_info < (3,4):
        sys.exit('Sorry, this code need Python 3.4 or higher. Please update. Aborting...')
    parser = parse_args()

    if len(sys.argv[1:]) < 1:
        parser.print_help()
        sys.exit(1)
    else:
        options = parser.parse_args()

    print('ScanITD build starts running: ' + time.strftime("%Y-%m-%d %H:%M:%S"))
    start = time.time()

    #check external tools used
    path = external_tool_checking()
    version_checking(path)

    # CIGAR string refinement or add SV tag 
    #itd_build = '{}'.format(os.path.join(os.path.dirname(options.input), '{}.itd_build.bam'.format(os.path.splitext(os.path.basename(options.input))[0])))
    itd_build = '{}'.format(os.path.join(os.path.dirname(options.output), '{}.itd_build.bam'.format(os.path.splitext(os.path.basename(options.input))[0])))

    softclipping_realignment(mapq_cutoff=options.mapq, input=options.input, output=itd_build)

    remove('{}.temp.bam'.format(itd_build))

    print("ScanITD build running done: " + time.strftime("%Y-%m-%d %H:%M:%S"))
    end = time.time()
    print('ScanITD build takes ' + str(end - start) + ' seconds.')

    ##########################################################################
    print('ScanITD call starts running: ' + time.strftime("%Y-%m-%d %H:%M:%S"))
    start = time.time()
    # predicint SV events: indel, duplication, inversion and translocation
    itd_scan(input_bam=itd_build, output_prefix=options.output, ao_cutoff=options.ao, dp_cutoff=options.dp, vaf_cutoff=options.vaf, itd_len_cutoff=options.itd_len, target=options.target, mapq_cutoff=options.mapq, mismatch_cutoff=options.mismatch, fasta_file=options.ref)
    print("ScanITD call running done: " + time.strftime("%Y-%m-%d %H:%M:%S"))
    end = time.time()

    # remove ITD build BAM by default
    if not options.keep:
        remove(itd_build)
        remove('{}.bai'.format(itd_build))

    print('ScanITD call takes ' + str(end - start) + ' seconds.')



if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(1)
