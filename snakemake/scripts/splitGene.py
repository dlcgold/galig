import sys, os
import argparse
import time
import subprocess
import glob

import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gffutils
import pysam

bar_length = 50

def getTime():
    return time.strftime('[ %b %d, %Y - %l:%M:%S%p ]')

def eprint(*args, **kwargs):
    print(getTime(), *args, file=sys.stderr, **kwargs)
    
def eprint_error(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    sys.exit(1)

def print_bar(BAR, i, n):
    sys.stderr.write("[{}] {}/{}\r".format(''.join(BAR), min(i,n), n))
    sys.stderr.flush()
    
def openGTF(gtfPath, verbose=True):
    try:
        gtf = gffutils.FeatureDB("{}.db".format(gtfPath),
                                 keep_order=True)
    except ValueError:
        if verbose:
            eprint("Indexing...")
        gtf = gffutils.create_db(gtfPath,
                                 dbfn="{}.db".format(gtfPath),
                                 force=True, keep_order=True,
                                 disable_infer_genes=True,
                                 disable_infer_transcripts=True,
                                 merge_strategy='merge',
                                 sort_attribute_values=True)
        gtf = gffutils.FeatureDB("{}.db".format(gtfPath), keep_order=True)
    return gtf

def splitAnnotation(annPath, outPath, annosFold):
    eprint("Opening input annotation...")
    gtf = openGTF(annPath)

    eprint("Splitting input annotation...")
    outFold = os.path.join(outPath, annosFold)
    if not os.path.isdir(outFold):
        os.makedirs(outFold)

    genes = list()
    chr_genes_dict = {}
    tr_gene_dict = {}

    nGenes = len(list(gtf.features_of_type('gene')))
    i = 0
    count = 0
    BAR = [' ' for i in range(0, bar_length)]
    print_bar(BAR, 0, nGenes)
    for gene in gtf.features_of_type('gene'):
        chrom = gene.seqid
        geneID = gene.id.split('.')[0]
        genes.append(geneID)
        chr_genes_dict[chrom] = chr_genes_dict[chrom] + 1 if chrom in chr_genes_dict else 1
        outPath = os.path.join(outFold, "{}.gtf".format(geneID))
        with open(outPath, 'w') as out:
            out.write(str(gene) + "\n")
            for transcript in gtf.children(gene, featuretype='transcript', order_by='start'):
                transcriptID = transcript.id.split('.')[0]
                tr_gene_dict[transcriptID] = geneID
                out.write(str(transcript) + "\n")
                for exon in gtf.children(transcript, featuretype='exon', order_by='start'):
                    out.write(str(exon) + "\n")
                    
        i+=1
        count += 1
        if i >= nGenes-1:
            BAR = ['#' for i in range(0,bar_length)]
            i = nGenes - 1
        else:
            index = int(i*bar_length/nGenes)
            BAR = ['#' for i in range(0, index+1)] + [' ' for i in range(index+1,bar_length)]
        print_bar(BAR, i+1, nGenes)
    print("", file=sys.stderr)
    tmp = outPath.split('/')
    #print(tmp[0] + tmp[1] + '.txt')
    
    outcount = open(tmp[0] + '/' + tmp[1] + '.txt', 'w')
    outcount.write(str(count))
    eprint("Done.")
    return genes, chr_genes_dict, tr_gene_dict


def main():
    genes, chr_genes_dict, tr_gene_dict = splitAnnotation(snakemake.input[0], snakemake.params[1], snakemake.params[0])
if __name__ == '__main__':
    main()
