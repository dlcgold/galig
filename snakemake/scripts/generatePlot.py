import os
import sys
import math
import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import matplotlib.pyplot as plt

def truncate(number, decimals=0):
    """
    Returns a value truncated to a specific number of decimal places.
    """
    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer.")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more.")
    elif decimals == 0:
        return math.trunc(number)

    factor = 10.0 ** decimals
    return math.trunc(number * factor) / factor

def autolabel(rects, ax, xpos='center'):
    """
    Attach a text label above each bar in *rects*, displaying its height.

    *xpos* indicates which side to place the text w.r.t. the center of
    the bar. It can be one of the following {'center', 'right', 'left'}.
    """

    xpos = xpos.lower()  # normalize the case of the parameter
    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
                '{}'.format(height), ha=ha[xpos], va='bottom')

def main(data_path):

    filesNE = [str.format("{}/{}", data_path, f) for f in os.listdir(data_path) if f.endswith('NoExon.txt')]
    files = [str.format("{}/{}", data_path, f) for f in os.listdir(data_path) if f.endswith('.txt') and f[len(f)-10:len(f)-4] != 'NoExon']
    
    datas = []
    for file in files:
        f = open(file, 'r')
        lines = f.readlines()
        datas.append([l.split('=')[1].strip() for l in lines])
        #print(f.name,lines)
    datasNE = []
    for file in filesNE:
        f = open(file, 'r')
        lines = f.readlines()
        datasNE.append([l.split('=')[1].strip() for l in lines])
        #print(f.name,lines)
        #print(datas)
    
    noRemoved = []
    removed = []
    countzero = 0
    for data in datas:
        #print(data[6])
        if int(data[6]) != 0:
            if data[0] == data[1]:
                noRemoved.append(data)
            else:
                removed.append(data)
        else:
            countzero += 1

    #print(removed)
    # for e in noRemoved:
    #     print(e)
    # print(len(noRemoved))
    # print(len(removed))
    # print(len(removed) + len(noRemoved))
    # print(len(files))
   
    noRemovedNE = []
    removedNE = []
    countzeroNE = 0
    for data in datasNE:
        #print(data[6])
        if int(data[6]) != 0:
            if data[0] == data[1]:
                noRemovedNE.append(data)
            else:
                removedNE.append(data)
        else:
            countzeroNE += 1
    # print(noRemoved)
    # print(removed)
    
    count = 0
    noRemovedListCounter = [0] * 5
    for elem in noRemoved:
        if elem[6] != 0:
            count += 1
            for i in range(0, 5):
                noRemovedListCounter[i] += int(elem[i+2])
            
    #print(noRemovedListCounter)
    #print(count)
    for i in range(0, 5):
        noRemovedListCounter[i] = math.ceil(noRemovedListCounter[i] / count)
    #print(noRemovedListCounter)
    count = 0
    removedListCounter = [0] * 5
    for elem in removed:
        if elem[6] != 0:
            count += 1
            for i in range(0, 5):
                removedListCounter[i] += int(elem[i+2])
    # print(removedListCounter, count)
    #print(removedListCounter)
    #print(count)
    for i in range(0, 5):
        removedListCounter[i] = math.ceil(removedListCounter[i] / count)
    #print(removedListCounter)
    # print(countzero)

    count = 0
    noRemovedListCounterNE = [0] * 5
    for elem in noRemovedNE:
        if elem[6] != 0:
            count += 1
            for i in range(0, 5):
                noRemovedListCounterNE[i] += int(elem[i+2])

    for i in range(0, 5):
        noRemovedListCounterNE[i] = math.ceil(noRemovedListCounterNE[i] / count)
    count = 0
    removedListCounterNE = [0] * 5
    for elem in removedNE:       
        if elem[6] != 0:
            count += 1
            for i in range(0, 5):
                removedListCounterNE[i] += int(elem[i+2])

    # print(removedListCounterNE, count)
    for i in range(0, 5):
        removedListCounterNE[i] = math.ceil(removedListCounterNE[i] / count)
    
    removedListCounterF = [0] * 5
    for i in range(0, 5):
        if i == 0:
            removedListCounterF[0] +=  removedListCounter[i]
        elif i == 1 or i == 2:
            removedListCounterF[1] +=  removedListCounter[i]
        elif i == 3:
            removedListCounterF[2] +=  removedListCounter[i]
        else:
            removedListCounterF[3] +=  removedListCounter[i]

    #print(removedListCounter)
    noRemovedListCounterF = [0] * 5
    for i in range(0, 5):
        if i == 0:
            noRemovedListCounterF[0] +=  noRemovedListCounter[i]
        elif i == 1:
            noRemovedListCounterF[1] +=  noRemovedListCounter[i]
        elif i == 2:
            noRemovedListCounterF[2] +=  noRemovedListCounter[i]
        else:
            noRemovedListCounterF[3] +=  noRemovedListCounter[i]

    removedListCounterNEF = [0] * 5
    for i in range(0, 5):
        if i == 0:
            removedListCounterNEF[0] +=  removedListCounterNE[i]
        elif i == 1 or i == 2:
            removedListCounterNEF[1] +=  removedListCounterNE[i]
        elif i == 3:
            removedListCounterNEF[2] +=  removedListCounterNE[i]
        else:
            removedListCounterNEF[3] +=  removedListCounterNE[i]

    print(removedListCounter)
    print(noRemovedListCounter)
    print(removedListCounterNE)

    benchNE = [str.format("{}/{}", "benchmarks", f) for f in os.listdir("benchmarks") if f.endswith('NoExon.tsv')]
    bench = [str.format("{}/{}", "benchmarks", f) for f in os.listdir("benchmarks") if f.endswith('.tsv') and f[len(f)-10:len(f)-4] != 'NoExon' and f[len(f)-10:len(f)-4] != 'Origin']
    benchO = [str.format("{}/{}", "benchmarks", f) for f in os.listdir("benchmarks") if f.endswith('Origin.tsv')]

    #print(benchNE)
    benchd = []
    benchNEd = []
    benchOd = []

    for elem in bench:
        f = open(elem, 'r')
        tmp = f.readlines()[1].split()
        benchd.append(tmp)
    
    benchdCounter = [0] * 2
    for e in benchd:
        if e[2] == '-':
            e[2] = 0
        if e[3] == '-':
            e[3] = 0
        benchdCounter[0] += float(e[0])
        benchdCounter[1] += (float(e[2]) + float(e[3]))
    print(benchdCounter)
    benchdCounter[0] /= len(benchd)
    benchdCounter[1] /= len(benchd)
    #print(benchdCounter)
    for elem in benchNE:
        f = open(elem, 'r')
        tmp = f.readlines()[1].split()
        benchNEd.append(tmp)
    benchNEdCounter = [0] * 2
    for e in benchNEd:
        if e[2] == '-':
            e[2] = 0
        if e[3] == '-':
            e[3] = 0
        benchNEdCounter[0] += float(e[0])
        benchNEdCounter[1] += (float(e[2]) + float(e[3]))
    benchNEdCounter[0] /= len(benchNEd)
    benchNEdCounter[1] /= len(benchNEd)

    for elem in benchO:
        f = open(elem, 'r')
        tmp = f.readlines()[1].split()
        #print(f.name, tmp)
        benchOd.append(tmp)
    benchOdCounter = [0] * 2
    for e in benchOd:
        if e[2] == '-':
            e[2] = 0
        if e[3] == '-':
            e[3] = 0
        benchOdCounter[0] += float(e[0])
        benchOdCounter[1] += (float(e[2]) + float(e[3]))
    print(benchOdCounter)
    benchOdCounter[0] /= len(benchOd)
    benchOdCounter[1] /= len(benchOd)
    for i in range(0,2):
        benchdCounter[i] = truncate(benchdCounter[i], 3)
        benchNEdCounter[i] = truncate(benchNEdCounter[i], 3)
        benchOdCounter[i] = truncate(benchOdCounter[i],3)

    outfileD = open(str.format("{}/data.res", data_path), "w")
    outfileT = open(str.format("{}/times.res", data_path), "w")
    outfileM = open(str.format("{}/mems.res", data_path), "w")
    dataF = open(str.format("{}/data.png", data_path), "w")
    timesF = open(str.format("{}/times.png", data_path), "w")
    memsF = open(str.format("{}/mems.png", data_path), "w")
    
    outfileD.write("orig:\n")
    for e in noRemovedListCounter:
        outfileD.write(str(e) + " ")
        
    outfileD.write("\nintron:\n")
    for e in removedListCounter:
        outfileD.write(str(e) + " ")
    outfileD.write("\norig no intron:\n")
    for e in removedListCounterNE:
        outfileD.write(str(e) + " ")
        
    outfileT.write("orig:\n")
    outfileT.write(str(benchOdCounter[0]))
    outfileT.write("\nintron:\n")
    outfileT.write(str(benchdCounter[0]))
    outfileT.write("\norig no intron:\n")
    outfileT.write(str(benchNEdCounter[0]))

    outfileM.write("orig:\n")
    outfileM.write(str(benchOdCounter[1]))
    outfileM.write("\nintron:\n")
    outfileM.write(str(benchdCounter[1]))
    outfileM.write("\norig no intron:\n")
    outfileM.write(str(benchNEdCounter[1]))
    dataF.write("ok")
    timesF.write("ok")
    memsF.write("ok")
    # width = 0.3
    # labels = ['match', 'mismatch', 'tot']

    # plt.figure(1)
    # plt.xticks(np.arange(len(removedListCounterNEF)) + width, labels)
    # dataO = plt.bar(np.arange(len(noRemovedListCounterF)), noRemovedListCounterF, width=width, label="original", color="#a3be8c")
    # data = plt.bar(np.arange(len(removedListCounterF))+ width, removedListCounterF, width=width, label="without one exon", color="#ebcb8b")
    # #plt.bar(np.arange(len(noRemovedListCounterNE)) + 2*width, noRemovedListCounterNE, width=width, label= "originalNE")
    # dataNE = plt.bar(np.arange(len(removedListCounterNEF))+ 2*width, removedListCounterNEF, width=width, label="Original without one exon", color="#bf616a")
    # # plt.show()
    # ax=plt.axes()
    # ax.set_ylabel("Reads")
    # lg = plt.legend(fancybox=True, shadow=True, prop={'size':'small'}, bbox_to_anchor=(0.95,1.0,0.3,0.2), loc='upper left')
    # plt.suptitle('Average results', fontsize=16)
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # # autolabel(dataO, ax)
    # # autolabel(data, ax)
    # # autolabel(dataNE, ax)
    # plt.savefig(str.format("{}/data.png", data_path), bbox_inches='tight')
    # plt.figure(2)
    # width = 0.1
    # # labels = ['Second']
    # # plt.xticks(np.arange(1) + width, labels)
    # timeO = plt.bar(np.arange(1), benchOdCounter[0], width=width, label="original", color="#a3be8c")
    # time = plt.bar(np.arange(1)+ width, benchdCounter[0], width=width, label="without one exon", color="#ebcb8b")
    # #plt.bar(np.arange(len(noRemovedListCounterNE)) + 2*width, noRemovedListCounterNE, width=width, label= "originalNE")
    # timeNE = plt.bar(np.arange(1)+ 2*width, benchNEdCounter[0], width=width, label="Original without one exon", color="#bf616a")
    # ax=plt.axes()
    # ax.axes.xaxis.set_visible(False)
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # ax.set_ylabel("Seconds")
    # lg = plt.legend(fancybox=True, shadow=True, prop={'size':'small'}, bbox_to_anchor=(0.95,1.0,0.3,0.2), loc='upper left')
    # # autolabel(timeO, ax)
    # # autolabel(time, ax)
    # # autolabel(timeNE, ax)
    # # plt.show()
    # plt.suptitle('Average execution time', fontsize=16)
    # plt.savefig(str.format("{}/times.png", data_path), bbox_inches='tight')
    
    # plt.figure(3)
    # width = 0.1
    # # labels = ['Second']
    # # plt.xticks(np.arange(1) + width, labels)
    # memO = plt.bar(np.arange(1), benchOdCounter[1], width=width, label="original", color="#a3be8c")
    # mem = plt.bar(np.arange(1)+ width, benchdCounter[1], width=width, label="without one exon", color="#ebcb8b")
    # #plt.bar(np.arange(len(noRemovedListCounterNE)) + 2*width, noRemovedListCounterNE, width=width, label= "originalNE")
    # memNE = plt.bar(np.arange(1)+ 2*width, benchNEdCounter[1], width=width, label="Original without one exon", color="#bf616a")
    # ax=plt.axes()
    # ax.axes.xaxis.set_visible(False)
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # ax.set_ylabel("Mb")
    # lg = plt.legend(fancybox=True, shadow=True, prop={'size':'small'}, bbox_to_anchor=(0.95,1.0,0.3,0.2), loc='upper left')
    # # autolabel(memO, ax)
    # # autolabel(mem, ax)
    # # autolabel(memNE, ax)
    # # plt.show()
    # plt.suptitle('Average execution memory use', fontsize=16)
    # plt.savefig(str.format("{}/mems.png", data_path), bbox_inches='tight')
if __name__ == '__main__':
    #print(snakemake.params[0])
    main(snakemake.params[0])
    #main(sys.argv[1])
