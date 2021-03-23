import sys
import os
import random
from Bio import SeqIO
import gffutils

bar_length = 50

def print_bar(BAR, i, n):
    sys.stderr.write("[{}] {}/{}\r".format(''.join(BAR), min(i,n), n))
    sys.stderr.flush()
    
def open_gtf(gtf_path, verbose=True):
    try:
        gtf = gffutils.FeatureDB("{}.db".format(gtf_path),
                                 keep_order=True)
    except ValueError:
        gtf = gffutils.create_db(gtf_path,
                                 dbfn="{}.db".format(gtf_path),
                                 force=True, keep_order=True,
                                 disable_infer_genes=True,
                                 disable_infer_transcripts=True,
                                 merge_strategy='merge',
                                 sort_attribute_values=True)
        gtf = gffutils.FeatureDB("{}.db".format(gtf_path), keep_order=True)
    return gtf

def extract_genes(gtf):
    genes = {}
    for gene in gtf.features_of_type('gene'):
        gene_idx = gene.id
        genes[gene_idx] = {}
        for transcript in gtf.children(gene, featuretype='transcript', order_by='start'):
            tr_idx = transcript.id
            exons = []
            for exon in gtf.children(transcript, featuretype='exon', order_by='start'):
                exons.append((exon.start, exon.end))
            genes[gene_idx][tr_idx] = exons
    return genes

def extract_exons(gtf):
    exons = []
    for gene in gtf.features_of_type('gene'):
        for transcript in gtf.children(gene, featuretype='transcript', order_by='start'):
            for exon in gtf.children(transcript, featuretype='exon', order_by='start'):
                exons.append((exon.start, exon.end))
    return exons

def get_exon_to_remove(exons, rangeMin, rangeMax):

    count = 0
    removable_exons = []
    exons = sorted(exons)
    for e in exons:
        if e[1] - e[0] >= rangeMin and e[1] - e[0] <= rangeMax:
            count += 1
            removable_exons.append(e)
    #print(count)
    # remove_exon_index = random.randint(0, count)
    if count == 0 or len(exons) <= 3:
        exon_to_remove = (-1, -1)
    else:
        exon_to_remove = removable_exons[int(count / 2)]
    return exon_to_remove

def new_gtf_to_file(gtf, exon_to_remove, outfile):
    if exon_to_remove != (-1, -1):
        exonlist = []
        for gene in gtf.features_of_type('gene'):
            for transcript in gtf.children(gene, featuretype='transcript', order_by='start'):
                exon_complete = []
                for exon in gtf.children(transcript, featuretype='exon', order_by='start'):
                    exon_complete.append(exon)
                for exon in gtf.children(transcript, featuretype='exon', order_by='start'):
                    if (exon.start, exon.end) == exon_to_remove and (exon.start, exon.end) != (exon_complete[0].start, exon_complete[0].end) and (exon.start, exon.end) != (exon_complete[-1].start, exon_complete[-1].end):
                        exonlist.append(exon)
        for line in gtf.all_features():
            if not line in exonlist:
                outfile.write(str(line) + "\n")
    else:
        for line in gtf.all_features():
            outfile.write(str(line) + "\n")
                
            
def main():

    #gtf_path = sys.argv[0]
    gtf_path = snakemake.input[0]
    fa_path = snakemake.input[1]
    out_path_dir = snakemake.params[0]
    out_path_name = snakemake.params[1]
    rangeMin = int(snakemake.params[2]) 
    rangeMax = int(snakemake.params[3])
    gtf_name = gtf_path.split('/')[-1]
    gtf = open_gtf(gtf_path)
    exons = extract_exons(gtf) # assuming single gene
    ref = list(SeqIO.parse(fa_path, "fasta"))[0] # Assuming single chromosome
    out_path = str.format("{}/{}", out_path_dir, out_path_name)
    outfile = open(str.format("{}/{}", out_path, gtf_name), "w")
    exon_to_remove = get_exon_to_remove(exons, rangeMin, rangeMax)
    new_gtf_to_file(gtf, exon_to_remove, outfile)
    
if __name__ == '__main__':
    main()
