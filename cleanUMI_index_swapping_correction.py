#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 16:57:58 2019

@author: yumingcao
"""
# adapted from alan's cleanUMI.py
import re
import numpy as np

def maxloc(num_list):
    return np.argmax(num_list)

def hammingDist(x, y):
    # find the Hamming distance between two input strings:
    if len(x)!=len(y):
        hd = len(x)
    else:
        hd = 0
        for i in range(len(x)):
            if x[i]!=y[i]:
                hd+=1    # count mismatches
    return hd

def testUmiDist(singList, multList):
    ## test each element of singList against all elements of multList to see if any
    ## singletons are one base off from any of the multis.

    badUmis = []   # output list of UMIs to be removed
    repUmis = {}   # UMIs to be incremented, with the increment count

    for bc0 in singList:
        bc1off = []   # barcodes with a Hamming distance of 1 from the singelton
        for bc1 in multList:
            if hammingDist(bc0, bc1)==1:
                bc1off.append(bc1)  # add to the list of one-offs
        # if a barcode is one off from any multi, remove it:
        if len(bc1off)>0:
            badUmis.append(bc0)
        # if bc0 is one off from only ONE multi, increment the count for that multi:
        if len(bc1off)==1:
            repUmis.setdefault(bc1off[0],0)
            repUmis[bc1off[0]]+=1

    return [badUmis, repUmis]

def sumValues(d):
    ## sums the values of the values in a key:value list:
    x=0
    for key,value in d:
        x+=value
    return x

def removeBC(bclist, the_dict):
    for key in bclist:
        if key in the_dict:
            del the_dict[key]

def mergeBarcodes(bcSampDict, the_dict):
    outDict = {}
    count = 0
    for bc, gene in the_dict.items():
        count += 1
        print(count)
        for g, umi in gene.items():
            for u, s in umi.items():
                if bc in bcSampDict.keys():
                    print(bc)
                    sample = bcSampDict[bc]
                    barcode = "_".join([str(sample),str(bc)])
                    outDict.setdefault(barcode, {})
                    outDict[barcode].setdefault(g, {})
                    outDict[barcode][g].setdefault(u, 0)
                    for n in s.values():
                        outDict[barcode][g][u] += n
                else:
                    for sample, n in s.items():
                        barcode = "_".join([str(sample),str(bc)])
                        print(barcode)
                        outDict.setdefault(barcode, {})
                        outDict[barcode].setdefault(g, {})
                        outDict[barcode][g].setdefault(u, 0)
                        outDict[barcode][g][u] += n
    return outDict

def getMergedBC(the_dict):
    mergedBC = []
    for bc in the_dict.keys():
        mergedBC.append(bc)
    return mergedBC


#def modifyDict(the_dict):
#    for bc, gene in the_dict.items():
#        for g, umi in gene.items():
#            for u, samp in umi.items():
#                for s, n in samp.items():
#                    barcode = "_".join([str(s),str(bc)])
#                    the_dict[barcode] = the_dict[bc]
#                    del the_dict[bc]

    
    
## define inputs
outFile = "flowcell1.umi.txt"
# correct index swapping:
# step 1 : combine umi.distributions.txt files together in a single donor
files = ["flu1002.umi.distributions.txt", "flu1004.umi.distributions.txt",
         "flu1201.umi.distributions.txt", "flu1202.umi.distributions.txt","flu1203.umi.distributions.txt","flu1204.umi.distributions.txt","flu1205.umi.distributions.txt",
         "flu2301.umi.distributions.txt", "flu2302.umi.distributions.txt","flu2303.umi.distributions.txt","flu2304.umi.distributions.txt","flu2305.umi.distributions.txt",
         "flu2501.umi.distributions.txt", "flu2502.umi.distributions.txt","flu2503.umi.distributions.txt","flu2504.umi.distributions.txt","flu2505.umi.distributions.txt"]
#flu12files = ["flu1201.umi.distributions.txt", "flu1202.umi.distributions.txt","flu1203.umi.distributions.txt","flu1204.umi.distributions.txt","flu1205.umi.distributions.txt"]
#flu23files = ["flu2301.umi.distributions.txt", "flu2302.umi.distributions.txt","flu2303.umi.distributions.txt","flu2304.umi.distributions.txt","flu2305.umi.distributions.txt"]
#flu25files = ["flu2501.umi.distributions.txt", "flu2502.umi.distributions.txt","flu2503.umi.distributions.txt","flu2504.umi.distributions.txt","flu2505.umi.distributions.txt"]

#umiDict = {}
umiHist = {}
#gDict = {}   # dictionary whose keys are all observed genes
bcDict = {}

for txtfile in files:
    print(txtfile)
    sample = txtfile.split(".")[0]
    sample = re.sub("..$","",sample)
    with open(txtfile, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            bc = fields[0]
            g = fields[1]
            umi = fields[2]
            n = int(fields[3])
            
#            #update the gene dictionary:
#            gDict.setdefault(g,0)
    
            # UMI counts histogram:
            umiHist.setdefault(bc,{})
            umiHist[bc].setdefault(g,{})
            umiHist[bc][g].setdefault(umi,{})
            umiHist[bc][g][umi].setdefault(sample, 0)
            umiHist[bc][g][umi][sample] += n
            
#            # keep stats on all UMIs:
#            umiDict.setdefault(umi,0)
#            umiDict[umi]+=n
            
            bcDict.setdefault(bc, {})
            bcDict[bc].setdefault(sample, 0)
            bcDict[bc][sample] += n

#with open("flowcell1_pythonvariables", "wb") as f:
#    pickle.dumps([umiHist, bcDict, gDict, umiDict] ,f)
#    
#with open("flowcell1_pythonvariables", "wb") as f:
#    pickle.dump(umiHist, f)
#
#with open("flowcell1_pythonvariables", "rb")  as f:
#    umiHist, bcDict, gDict, umiDict = pickle.load(f)
    
## correct for umis from multiple sample from bcDict
# only correct for barcode that with 95% of the reads from one donor
mergeBC = []
mergeBCwSamp = {}
conflictBC= []
for bc, samp in bcDict.items():
    if len(samp.keys()) > 1:
        total_reads = sum(samp.values())
        frac = [ x/total_reads for x in samp.values()]
        conflictBC.append(bc)
        if max(frac) >= 0.95: ## only consider merge barcode with 95% reeads from one patient
            max_idx = maxloc(frac)
            max_samp = list(samp.keys())[max_idx]
            mergeBC.append(bc)
            mergeBCwSamp.setdefault(bc, None)
            mergeBCwSamp[bc] = max_samp

discardBC = list(set(conflictBC)^set(mergeBC))

removeBC(discardBC, umiHist) ## remove BC from umiHist

## merge barcode collision to the majority sample
print("merging...")
umiHist_new = mergeBarcodes(mergeBCwSamp, umiHist) 

mergedBClist=getMergedBC(umiHist_new)

with open("flowcell1_correctedBC.txt", "w") as f:
    bc = list(umiHist_new.keys())
    for b in bc:
        f.write(b + "\n")
 ## get lists of UMIs with one count and UMIs with >nMin counts:
#meanMin = 2 # minimum number of UMI counts for use in computing UMI mean
#for bc in umiHist_out.keys():
#    for g in umiHist_out[bc].keys():
#        # for average UMI count calculations:
#        nzUmis = 0   # number of UMIs with non-zero counts
#        ntUmis = 0   # number of UMIs with >1 count
#        umiSum = 0   # sum for UMI mean calculation
#        # make lists of singlet and multis:
#        singlets = []
#        multis = []
#        for umi in umiHist_out[bc][g].keys():
#            if umiHist_out[bc][g][umi]==1:
#                singlets.append(umi)
#            elif umiHist_out[bc][g][umi]>=nMin:
#                multis.append(umi)
#                nzUmis+=1
#            else:
#                nzUmis+=1
#
#        # separate true singlets from ones that have a Hamming distance of 1 from one or more of the multis:
#        [rmUmis, incUmis] = testUmiDist(singlets, multis)
#
#        # remove bad UMIs from the dictionary
#        for i in range(len(rmUmis)):
#            del umiHist_out[bc][g][rmUmis[i]]     # delete error UMI
#        # update the counts for UMIs that were uniquely one-off from one of the singletons:
#        for u in incUmis.keys():
#            umiHist_out[bc][g][u]+=incUmis[u]
#print("merging done")


        
## For each gene, the exprssion value is len(umiHist[g].keys()) (i.e., the number of UMIs 
## remaining with at least one copy):
#fOut = open(outFile, 'w')
#
#print("UMI correcting...")
##
#### open the optional read counts file:
##if readsFile==None:
##    writeReads=False
##else:
##    writeReads=True
##    fReads = open(readsFile, 'w')
#
### remove any barcodes/cells with fewer than uMin total UMIs (if uMin>0):
#if uMin>0:
#    bcRm = []
#    for bc in umiHist_out.keys():   # loop over barcodes
#        bcSum = 0
#        for g in umiHist_out[bc].keys():   # loop over genes
#            bcSum += len(umiHist_out[bc][g].keys())  # sum total UMIs for this cell
#        ## remove this barcode if the total count is less than uMin:
#        if bcSum<uMin:
#            bcRm.append(bc)
### remove low-count cells:
#for bc in bcRm:
#    del umiHist_out[bc]
#
#print("UMI corrected")        
### file header:
#hStr = 'gene'      # header
#bcList = umiHist_out.keys()     # just to make absolutely sure they stay in the correct order
#for bc in bcList:
#    hStr = '%s\t%s' % (hStr,bc)
#fOut.write('%s\n' % hStr)
#
#### option reads count file header
##if writeReads:
##    fReads.write('%s\n' % hStr)
#print("writing table...")
## per-gene count:
#for g in gDict.keys():  # loop over all observed genes
#    oStr = '%s' % g
#    rStr = '%s' % g
#    for bc in bcList:   # UMI counts for each cell
#        if g in umiHist_out[bc]:
#            oStr = '%s\t%d' % (oStr,len(umiHist_out[bc][g].keys()))
#            rStr = '%s\t%d' % (rStr,sumValues(umiHist_out[bc][g].items()))
#        else:            
#            oStr = '%s\t0' % oStr
#            rStr = '%s\t0' % rStr
#
#    fOut.write('%s\n' % oStr)        
##    if writeReads:
##        fReads.write('%s\n' % rStr)        
#
#fOut.close()
##if writeReads:
##    fReads.close()

                  
