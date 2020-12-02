import argparse
import zipfile
import numpy as np


def parse_annotation_file(annotation_fn):
    """
    :param annotation_fn: the gene annotations file
    :return: outputs a list of tuples, "genes". Every tuple represents one gene, the first element of the tuple is the
            list of exon index ranges for that gene, the second element of the tuple is the list of isoforms that exist
            for that gene. For example, genes[0][0][0] references the first exon range of the first gene (in a tuple),
            genes[2][1][3] references the fourth isoform of the third gene.
    """

    with open(annotation_fn, 'r') as aFile:
        N = int(aFile.readline().strip())
        genes = [None]*N
        for i in range(N):
            numExons = int(aFile.readline().strip())
            exons = [None]*numExons
            starts = [int(x) for x in aFile.readline().strip().split(' ')]
            ends = [int(x) for x in aFile.readline().strip().split(' ')]
            for j in range(numExons):
                exons[j] = (starts[j], ends[j])
            numIsoforms = int(aFile.readline().strip())
            isoforms = [None]*numIsoforms
            for j in range(numIsoforms):
                isoforms[j] = [int(x) for x in aFile.readline().strip().split(' ')]
            genes[i] = (exons, isoforms)
    return genes


def parse_genome_file(genome_fn):
    """
    :param genome_fn: the full genome file
    :return: the string containing the genome
    """

    with open(genome_fn, 'r') as gFile:
        return gFile.readline().strip()


def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file of shuffled reads
    :return: a list containing all of the shuffled reads
    """
    
    out_reads = []
    with open(reads_fn, 'r') as rFile:
        for line in rFile:
            out_reads.append(line.strip())
    return out_reads


def quantify_isoforms(genes, genome, reads):
    ret=[]
    for gene in genes:
        exons_seq=[]
        exons=gene[0]
        for exon in exons:
            exons_seq=exons_seq+[genome[exon[0]:exon[1]+1]]
        isoforms_seq=[]
        isoforms=gene[1]
        for isoform in isoforms:
            iso_seq=''
            for iso in isoform:
                iso_seq=iso_seq+exons_seq[iso]
            isoforms_seq=isoforms_seq+[iso_seq]
        isonum=len(isoforms_seq)
        if len(isoforms_seq)==3:
            dict_exons={}
            for ex in exons_seq:
                dict_exons[ex]=0
            zero_one=exons_seq[0]+exons_seq[1]
            zero_three=exons_seq[0]+exons_seq[3]
            one_two=exons_seq[1]+exons_seq[2]
            two_three=exons_seq[2]+exons_seq[3]
            for r in reads:
                if (r in zero_one or r in zero_three) and (r not in exons_seq[1]) and (r not in exons_seq[3]):
                    dict_exons[exons_seq[0]]+=1
                if (r in zero_one or r in one_two) and (r not in exons_seq[0]) and (r not in exons_seq[2]):
                    dict_exons[exons_seq[1]]+=1
                if (r in one_two or r in two_three) and (r not in exons_seq[1]) and (r not in exons_seq[3]):
                    dict_exons[exons_seq[2]]+=1
                if (r in zero_three or r in two_three) and (r not in exons_seq[0]) and (r not in exons_seq[2]):
                    dict_exons[exons_seq[3]]+=1
            print(dict_exons)
            x = np.array([[68/50, 68/50, 0], [63/50, 0, 63/50], [0, 0, 63/50],[0,54/50,54/50]])
            x_t=x.transpose()
            y=[dict_exons[exons_seq[0]],dict_exons[exons_seq[1]],dict_exons[exons_seq[2]],dict_exons[exons_seq[3]]]
            x_t_x=np.dot(x_t,x)
            x_t_x_minues = np.linalg.inv(x_t_x) 
            x_then=np.dot(x_t_x_minues,x_t)
            beta=np.dot(x_then,y)
            freq_li=[beta[0]/sum(beta),beta[1]/sum(beta),beta[2]/sum(beta)]
            ret=ret+[(isoforms_seq[0],freq_li[0])]
            ret=ret+[(isoforms_seq[1],freq_li[1])]
            ret=ret+[(isoforms_seq[2],freq_li[2])]
        else:
            freq=1/isonum
            for seq in isoforms_seq:
                ret=ret+[(seq,freq)]
    return ret


if __name__ == "__main__":
    """
    For an example of how you might call this script to run on the data provided:
    
    Usage: python proj4.py -g full_genome.txt -r shuffled_reads.txt -a DATA_PA_1100_0 -o test.out -t hw4_r_4_chr_1
    """
    parser = argparse.ArgumentParser(description='For now this starter code helps parse the files given, but leaves\n'
                                                 'the actual function that must be implemented empty')
    parser.add_argument('-g', '--genome', required=True, dest='genome_file', help='File containing the full genome')
    parser.add_argument('-r', '--reads', required=True, dest='read_file', help='File containing the shuffled reads')
    parser.add_argument('-a', '--annotation', required=True, dest='annotation_file', help='File containing gene '
                                                                                          'annotations')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file', help='Output file name')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be output on the first line of the output file so that the online\n'
                             'submission system recognizes which leaderboard this file should be submitted to. For\n'
                             'hw4, this will be hw4_r_4_chr_1')

    args = parser.parse_args()
    genome_fn = args.genome_file
    reads_fn = args.read_file
    annotation_fn = args.annotation_file
    output_fn = args.output_file

    genes = parse_annotation_file(annotation_fn)
    genome = parse_genome_file(genome_fn)
    reads = parse_reads_file(reads_fn)
    print(len(reads))
    output = quantify_isoforms(genes, genome, reads)
    with open(output_fn, 'w') as oFile:
        oFile.write('>' + args.output_header + '\n')
        oFile.write('>RNA\n')
        for isoform in output:
            out_str = '{} {}\n'.format(isoform[0], isoform[1])
            oFile.write(out_str)

    zip_fn = output_fn + '.zip'
    with zipfile.ZipFile(zip_fn, 'w') as zFile:
        zFile.write(output_fn)