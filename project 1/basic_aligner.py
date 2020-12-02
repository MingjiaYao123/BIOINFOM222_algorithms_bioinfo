import sys
import argparse
import numpy as np
import time
import zipfile


def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                if count % 1000 == 0:
                    print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None


def parse_ref_file(ref_fn):
    """
    :param ref_fn: the file containing the reference genome
    :return: a string containing the reference genome
    """
    try:
        with open(ref_fn, 'r') as gFile:
            print("Parsing Ref")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_aligner.py takes in data for homework assignment 1 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPs based on this alignment')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the '
                             'online submission system recognizes which leaderboard this file should be submitted to.'
                             'This HAS to be practice_W_1_chr_1 for the practice data and hw1_W_2_chr_1 for the '
                             'for-credit assignment!')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
    reference = parse_ref_file(reference_fn)
    if reference is None:
        sys.exit(1)
        
    hash_table={}
    for i in range(len(reference)-25):
        seq=reference[i:i+25]
        if seq not in hash_table:
            hash_table[seq]=[i]
        else:
            hash_table[seq]=hash_table[seq]+[i]  #build the hash table for checking
    
    seq_table=[]
    for line in input_reads:
        st1=line[0]
        seq_table=seq_table+[[st1[0:25],st1[25:50]]]
        st2=line[1]
        seq_table=seq_table+[[st2[0:25],st2[25:50]]]  #build the sequences for checking
    
    cand_snp=[]
    for seqs in seq_table:
        seq1=seqs[0]
        seq2=seqs[1]
        seq=seq1+seq2
        if seq1 in hash_table:
            inds=hash_table[seq1]
            for ind in inds:
                SNP_ind_list=[]
                if ind<=9950:
                    for i in range(25,50):
                        compare_seq1=seq[i]
                        compare_seq2=reference[ind+i]
                        if compare_seq1!=compare_seq2:
                            SNP_ind_list=SNP_ind_list+[[compare_seq2,compare_seq1,ind+i]]
                    if len(SNP_ind_list)<=3 and len(SNP_ind_list)>0:
                        cand_snp=cand_snp+SNP_ind_list
        if seq2 in hash_table:
            inds=hash_table[seq2]
            for ind in inds:
                SNP_ind_list=[]
                if ind>=25:
                    for i in range(25):
                        compare_seq1=seq1[-i]
                        compare_seq2=reference[ind-i]
                        if compare_seq1!=compare_seq2:
                            SNP_ind_list=SNP_ind_list+[[compare_seq2,compare_seq1,ind-i]] #find the candidate SNPs
                    if len(SNP_ind_list)<=3 and len(SNP_ind_list)>0:
                        cand_snp=cand_snp+SNP_ind_list
    pos_count={}
    for cand in cand_snp:
        pos=cand[2]
        if pos not in pos_count:   #find the appearance for each position
            pos_count[pos]=1
        else:
            pos_count[pos]+=1
    ret=[]
    for cand in cand_snp:
        pos=cand[2]
        if pos_count[pos]>=5:
            if cand not in ret:  #find the returned list of SNPs
                ret=ret+[cand]
            
    snps = ret
    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        header = '>' + args.output_header + '\n>SNP\n'
        output_file.write(header)
        for x in snps:
            line = ','.join([str(u) for u in x]) + '\n'
            output_file.write(line)

        tails = ('>' + x for x in ('STR', 'CNV', 'ALU', 'INV', 'INS', 'DEL'))
        output_file.write('\n'.join(tails))

    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
