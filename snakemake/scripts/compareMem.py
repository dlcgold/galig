import sys, os
import argparse
import subprocess


def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

def extractFromInfoFile(infoPath):
    lines = open(infoPath).readlines()

    line1 = lines[0].split(" ")
    ref, refLen = line1[0], int(line1[1])

    text = lines[1].strip("\n")

    exons = [(int(p[0]), int(p[1])) for p in [pos.split(",") for pos in lines[4].strip("\n").split()]]

    return ref, refLen, text, exons

def readLine(line):
    # 0: strand
    # 1: ID
    # 2: errors
    # 3 to -1: mems
    # -1: read
    line = line.strip("\n").strip(" ").split(" ")
    strand = line[0]
    readID = line[1]
    err = int(line[2])
    mems = line[3:-1]
    read = line[-1]
    return strand, readID, err, mems, read

def extractMEMs(mems, exon_length):
    MEMs = []
    for mem in mems:
        # Remove ( and ) from mem and cast to int
        mem = [int(x) for x in mem[1:-1].split(",")]
        if len(mem) == 4 and mem[3] == 1:
            mem[0] += exon_length - 1
        MEMs.append(mem)
    return MEMs

def main(mems1, mems2, gtf1, gtf2, data_path, res_name):
    out_path_dir = str.format("{}/{}", data_path, res_name)
    gtf_name = mems1.split('/')[-1].split('.')[0]
    #print(gtf_name)
    gtf1l = file_len(gtf1)
    gtf2l = file_len(gtf2)
    outfile = open(str.format("{}/{}.txt", out_path_dir, gtf_name), "w")
    #print(mems1)
    #print(mems2)
    dict1 = {}
    dict2 = {}
    ref = set()
    ref1, refLen1, textE1, exPos1 = extractFromInfoFile(gtf1 + ".sg")
    ref2, refLen2, textE2, exPos2 = extractFromInfoFile(gtf2 + ".sg")
    
    for line in open(mems1).readlines():
        val = []    
        strand1, readID1, err1, mems1, read1 = readLine(line)
        mems1list = extractMEMs(mems1, len(textE1))
        ref.add(readID1)
        #print(readID1, mems1list)
        #print(mems1list)
        for i in range(0, len(mems1list)):
            val.append((mems1list[i][1], mems1list[i][2]))
        #print(readID1, val)
        sum = 0
        for e in val:
            sum += e[1]
        if readID1 not in dict1:
            dict1[readID1] = val
        else:
            part = dict1[readID1]
            sump = 0
            for e in part:
                sump += e[1]
            if sum > sump:
                dict1[readID1] = val

    #print(dict1)
    readL = 0
    for line in open(mems2).readlines(): 
        val = []    
        strand2, readID2, err2, mems2, read2 = readLine(line)
        ref.add(readID2)
        mems2list = extractMEMs(mems2, len(textE2))
        #print(readID1, mems1list)
        #print(mems1list)
        for i in range(0, len(mems2list)):
            val.append((mems2list[i][1], mems2list[i][2]))
        #print(readID2, val)
        sum = 0
        for e in val:
            sum += e[1]
        if readID2 not in dict2:
            dict2[readID2] = val
        else:
            part = dict2[readID2]
            sump = 0
            for e in part:
                sump += e[1]
            if sum > sump:
                dict2[readID2] = val
        readL = len(read2)
    mismatch = 0
    #mismatchT = 0
    # mismatchlist = []
    semimismatch = 0
    semimatch = 0
    match = 0
    tot = 0
    err = int(readL/10)
    for r in ref:
        tot += 1
        if(r in dict1 and r in dict2):
            #print('id ', r, 'compare ', dict1[r], ' with ', dict2[r])
            if dict1[r] == dict2[r]:
                match += 1
            else:
                totl1 = 0
                totl2 = 0
                tot1 = 0
                tot2 = 0
                for l in dict1[r]:
                    totl1 += l[1]
                for l in dict2[r]:
                    totl2 += l[1]
                    
                for l in dict1[r]:
                    for e in l:
                        tot1 += e
                for l in dict2[r]:
                    for e in l:
                        tot2 += e
                        
                if tot1 == tot2 or totl1 == totl2:
                    match += 1
                elif abs(tot1 - tot2) <= err or abs(totl1 - totl2) <= err:
                    semimatch += 1
                else:
                    semimismatch += 1
                   # print('check')
        elif(r in dict1 and r not in dict2):
            #print('id ', r, 'only in mems1, with ', dict1[r])
            mismatch += 1
            #mismatchlist.append([dict1[r], 1])
        elif(r in dict2 and r not in dict1):
            #print('id ', r, 'only in mems2, with ', dict2[r])
            mismatch += 1
            #mismatchlist.append([dict2[r], 2])
        else:
            mismatch += 1
                        
    # print('match = ', match)
    # print('semimatch = ', semimatch)
    # print('semimismatch = ', semimismatch)
    # print('mismatch = ', mismatch)
    # print('tot = ', tot)
    outfile.write(str.format("len1 = {}\n", gtf1l))
    outfile.write(str.format("len2 = {}\n", gtf2l))
    outfile.write(str.format("match = {}\n", match))
    outfile.write(str.format("semimatch = {}\n", semimatch))
    outfile.write(str.format("semimismatch = {}\n", semimismatch))
    outfile.write(str.format("mismatch= {}\n", mismatch))
    #outfile.write(str.format("mismatchT = {}\n", mismatchT))
    outfile.write(str.format("tot = {}", tot))
    # for e in mismatchlist:
    #     print(e)
        
if __name__ == '__main__':
    # parser = argparse.ArgumentParser(description = "compare two mem file")
    # parser.add_argument('-mo', '--memone', required=True, help='first mem file')
    # parser.add_argument('-mt', '--memtwo', required=True, help='second mem file')
    # parser.add_argument('-ao', '--annotationone', required=True, help='first GTF input file containing the gene annotation')
    # parser.add_argument('-at', '--annotationtwo', required=True, help='second GTF input file containing the gene annotation')
    # args = parser.parse_args()
    # mems1 = args.memone
    # mems2 = args.memtwo
    # gtf1 = args.annotationone
    # gtf2 = args.annotationtwo
    main(snakemake.input[0], snakemake.input[1], snakemake.input[2], snakemake.input[3],  snakemake.params[0], snakemake.params[1])
