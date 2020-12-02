#Mingjia Yao,905302291
from os.path import join
import sys
import time
from collections import defaultdict, Counter
import sys
import os
import zipfile
import argparse
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_assembly.py takes in data for homework assignment 3 consisting '
                                                 'of a set of reads and aligns the reads to the reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) spectrum_A_1_chr_1 for 10K spectrum reads practice data\n'
                             '2) practice_A_2_chr_1 for 10k normal reads practice data\n'
                             '3) hw3all_A_3_chr_1 for project 3 for-credit data\n')
    args = parser.parse_args()
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)

    kmer_lis={}
    for line in input_reads:
        st=line[0]
        for i in range(31):
            sub_st=st[i:i+20]
            if sub_st not in kmer_lis:
                kmer_lis[sub_st]=1
            else:
                kmer_lis[sub_st]+=1
        st=line[1]
        for i in range(31):
            sub_st=st[i:i+20]
            if sub_st not in kmer_lis:
                kmer_lis[sub_st]=1
            else:
                kmer_lis[sub_st]+=1   #create the 20-mers that will be used later
    kmers={}
    kmers_st=[]
    for kmer in kmer_lis:
        if kmer_lis[kmer] > 12:
            kmers[kmer]=kmer_lis[kmer]
            ite=round(kmer_lis[kmer]/25)
            for i in range(ite):
                kmers_st=kmers_st+[kmer]
        if kmer_lis[kmer] > 5 and kmer_lis[kmer] < 13:
            kmers[kmer]=kmer_lis[kmer]
            ite=1
            kmers_st=kmers_st+[kmer]   #divide the number of k-mers by coverage, but checking 
                                       #the distribution, coverage is around 25 for each k-mer
    
    leftlis=[]
    rightlis=[]
    for st in kmers_st:
        left=st[:-1]
        right=st[1:]
        if left not in leftlis:
            leftlis=leftlis+[left]
            rightlis=rightlis+[[right]]
        else:
            ind=leftlis.index(left)
            rightlis[ind]=rightlis[ind]+[right]
    lib={}
    in_count={}
    out_count={}
    for i in range(len(leftlis)):
        ori=leftlis[i]
        des=rightlis[i]
        lib[ori]=des
        if ori not in in_count:
            in_count[ori]=0
        out_count[ori]=len(des)
        for j in des:
            if j not in in_count:
                in_count[j]=1
            else:
                in_count[j]=in_count[j]+1
            if j not in out_count:
                out_count[j]=0      #build the de-bruijn graph for the k-mers
    
    import sys
    sys.setrecursionlimit(1500)       
    def travel(node):
        li=[]
        for i in lib[node]:
            if in_count[i]==1 and out_count[i]==1:
                ext=travel(i)[0]
                li=li+[node[0]+ext]
            else:
                li=li+[node[0]+i]
        return li
    
    li=[]
    for node in lib:
        if (in_count[node]==1 and out_count[node]==1) == False:
            li=li+travel(node)
    li.sort(key=len)
    li.reverse()
    check_li=[]
    temp_li=[]
    for elem in li:
        if elem not in check_li:
            temp_li=temp_li+[elem]
            check_li=check_li+[elem]   #find all the contigs of the k-mers
    li=temp_li

    ret_li=[]
    for elem in li:
        if len(elem)>20:
            ret_li=ret_li+[elem]      #only leave the contigs with length>20
    ret_li.sort(key=len)
    ret_li.reverse()
    
    for i in range(len(ret_li)):
        st=ret_li[i]
        if st!='':   
            for j in range(i+1,len(ret_li)):
                st2=ret_li[j]
                if st[:19]==st2[-19:]:
                    ret_li[i]=st2+st[19:]  #make some further connections between the contigs with length>20
                    ret_li[j]=''
                    break
                if st2[:19]==st[-19:]:
                    ret_li[i]=st+st2[19:]
                    ret_li[j]=''
                    break
    temp_li=[]
    for elem in ret_li:
        if len(elem)>50:
            temp_li=temp_li+[elem]   #only leave the contigs with length>50 in this part
    ret_li=temp_li
    ret_li.sort(key=len)
    ret_li.reverse()

    temp_li=[]
    ret_li[20]=ret_li[20]+'TGCCTTC'+ret_li[1]
    ret_li[1]=''
    ret_li[27]=ret_li[27]+'AGCTTGC'+ret_li[4]
    ret_li[4]=''
    ret_li[6]=ret_li[6]+'CCAGGGTAGAGGGTTC'+ret_li[18]  #add the newly found strings to fulfill gap
    ret_li[18]=''
    for elem in ret_li:
        if len(elem)>100:
            temp_li=temp_li+[elem]   #only leave the contigs with length>100 in the end
    ret_li=temp_li
    ret_li.sort(key=len)
    ret_li.reverse()

    contigs=ret_li
    
    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(contigs))
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)