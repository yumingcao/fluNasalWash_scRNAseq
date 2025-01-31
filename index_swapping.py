#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 15:15:24 2019

@author: yumingcao
"""

#bcDict = {"AAAAAAAAAAAA": {"flu1002" : 50, "flu1004":20, "flu1201": 500},
#          "AAAAAAAAAACC":  {"flu1002" : 300, "flu1004":2,"flu1201": 5, "flu1202": 5},
#           "AAAAAAAGAACC":  {"flu1002" : 30, "flu1004":200}}

##### count how many reads associated with each barcode #####
f = open("conflict_BC_sample_read_table.txt", "r")
bcDict = {}
count = 0
for line in f:
    fields = line.strip().split("\t")
    count += 1
    print(count)
    bc = fields[0]
    nRead = fields[1]
    sample = fields[2]
    bcDict.setdefault(bc, {})
    bcDict[bc].setdefault(sample, {})
    bcDict[bc][sample] = nRead


fracDict = {}
for bc, samp in bcDict.items():
    fracDict.setdefault(bc, {})
    
    nReads = samp.values()
    nRead_frac = [x/sum(nReads) for x in nReads]
    idx = 0
    for samp_name in samp.keys():
        fracDict[bc].setdefault(samp_name, 0)
        fracDict[bc][samp_name] = (nRead_frac[idx], bcDict[bc][samp_name])
        idx += 1
print("fraction converted")

#df = df.reset_index()
#v = df.columns.tolist()
#v.pop(0)
#df_lng = df.melt(df, id_vars='index', value_vars = v, value_name = "fraction")

with open("fraction_conflict_barcode.txt","w") as out:
    for barcode, nested in fracDict.items():
        for sample, fraction in nested.items():
            outline = [barcode, sample, str(fraction)]
            outline_1 = "\t".join(outline)
            print(outline_1, file = out)

f.close()
        