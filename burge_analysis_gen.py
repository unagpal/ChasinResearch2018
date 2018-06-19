# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 15:24:28 2018

@author: udainagpal
"""


import numpy as np
import pandas as pd
from pandas import DataFrame
#import matplotlib.pyplot as plt

data = pd.read_csv("RBFOX2.csv", dtype=str, sep=',', index_col=0)
format = 'RBFOX2'

def excel_write (result_final):
    all_lists = spreadsheet_columns(result_final)
    df = DataFrame({'HM': all_lists[0], 'Position': all_lists[1], 'Concentration': all_lists[2], 'r': all_lists[3], 'r2': all_lists[4], 'p': all_lists[5]})
    df.to_excel('RBFOX2_correlations.xlsx', sheet_name = 'PTBP3', columns = ['HM', 'Position', 'Concentration', 'r', 'r2', 'p'], index = False)

print(len(data.values[:, 3].astype(float)))

def get7merlei (hm, posnum):
    seq = []
    pos = []
    index = []
    lei = []
    sevenmers = []
    count = 0
    wt = 0
    wt_lei = 0
    for i in range(len(hexmut)):
        if (hexmut[i]==hm):
            if count == 0:
                wt = seqlist[i]
                wt_lei = leilist[i]
            if count > 0:
                if (sevenmerposlist[i] >= posnum and sevenmerposlist[i] <= posnum + 6):
                    seq.append(seqlist[i])
                    pos.append(sevenmerposlist[i])
                    index.append(i)
            count = count + 1
    modpos = posnum + 22
    eightmers = []
    for i in seq:
        eightmers.append(i[modpos:(modpos+8)])
    sevenmers.append(wt[modpos:(modpos+7)])
    lei.append(wt_lei)
    for i in range(len(eightmers)):
        if (pos[i] < posnum + 6):
            sevenmers.append(eightmers[i][:7])
            lei.append(leilist[index[i]])
        else:
            if (eightmers[i][7] == wt[modpos+7]):
                sevenmers.append(eightmers[i][:7])
                lei.append(leilist[index[i]])
    return [sevenmers, lei];

def getposlist():
    tempposlist = []
    for i in range (-3, 47):
        if (i != 0):
            tempposlist.append(i)
    return tempposlist;

def signedRsquared (R):
    if (R > 0):
        return R**2;
    else:
        return - R**2;

def sevenmertorow (sevenmer):
    bases = ['A', 'C', 'G', 'T']
    num_str = ''
    for i in sevenmer:
        for j in range(len(bases)):
            if (i == bases[j]):
                num_str = num_str + str(j)
    num_str = str(num_str)
    ans =  int(num_str[0])* (4**6) + int(num_str[1]) * (4**5) + int(num_str[2])* (4**4) + int(num_str[3]) * (4**3) + int(num_str[4]) * (4**2) + int(num_str[5])*4 + int(num_str[6])
    #if (ans > len(allsevenmers)-1):
    #    return len(allsevenmers)-1;
    return ans;
def sevenmertorowmissing (sevenmer):
    temp = sevenmertorow(sevenmer)
    if (temp > (len(allsevenmers)-1)):
        temp = len(allsevenmers) - 1
    if (allsevenmers[temp] == sevenmer):
        return temp;
    """first = 0
    last = temp
    found = False
    while first < last and not found:
        mid = (first + last)/2
        if (allsevenmers[mid]==sevenmer):
            return mid;
        else:
            if (sevenmertorow(allsevenmers[mid]) < temp)"""
    while (allsevenmers[temp] != sevenmer and temp>=-1):
        temp = temp - 1
    return temp;
    
def spreadsheet_columns(result_final):
    poslist = [i for i in range (-3, 49)]
    hexmut = []
    pos = []
    conc = []
    r = []
    r2 = []
    p = []
    for i in range(len(result_final)):
        for j in range (len(result_final[i])):
            for k in range (len(result_final[i][j])):
                if (poslist[j] >= -2 and poslist[j]<=46):
                    hexmut.append(str(hexmutlist[i]))
                    if (poslist[j] > 0):
                        pos.append(poslist[j])
                    else:
                        pos.append(poslist[j]-1)
                    conc.append(str(concentrationlist[k]))
                    r.append(results_final[i][j][k][0])
                    r2.append(results_final[i][j][k][1])
                    p.append(results_final[i][j][k][2])
                    #print("Hexmut: " + str(hexmutlist[i]) + ', Position: ' + str(poslist[j]) + ', Concentration: ' + str(concentrationlist[k]) + ': ' + str(results_final[i][j][k])[1:-1])
    return [hexmut, pos, conc, r, r2, p];
 
def print_results (result_final):
    poslist = [i for i in range (-3, 49)]
    for i in range(len(result_final)):
        for j in range (len(result_final[i])):
            for k in range (len(result_final[i][j])):
                print("Hexmut: " + str(hexmutlist[i]) + ', Position: ' + str(poslist[j]) + ', Concentration: ' + str(concentrationlist[k]) + ': ' + str(results_final[i][j][k])[1:-1])
    return;

#print(get7merlei('A', -1))
if (format == 'missing5'):
    data_1300 = data.values[:, 3].astype(float)
    data_320 = data.values[:, 2].astype(float)
    data_80 = data.values[:, 1].astype(float)
    data_20 = data.values[:, 0].astype(float)
    concentrationlist = ['20 nM', '80 nM', '320 nM', '1300 nM']
    all_data = [data_20, data_80, data_320, data_1300]
if (format == 'missing20'):
    data_1300 = data.values[:, 3].astype(float)
    data_320 = data.values[:, 2].astype(float)
    data_80 = data.values[:, 1].astype(float)
    data_5 = data.values[:, 0].astype(float)
    concentrationlist = ['5 nM', '80 nM', '320 nM', '1300 nM']
    all_data = [data_5, data_80, data_320, data_1300]
if (format == 'missing320'):
    data_1300 = data.values[:, 3].astype(float)
    data_80 = data.values[:, 2].astype(float)
    data_20 = data.values[:, 1].astype(float)
    data_5 = data.values[:, 0].astype(float)
    concentrationlist = ['5 nM', '20 nM', '80 nM', '1300 nM']
    all_data = [data_5, data_20, data_80, data_1300]
if (format == ''):
    #data_3900 = data.values[:, 5].astype(float)
    data_1300 = data.values[:, 4].astype(float)
    data_320 = data.values[:, 3].astype(float)
    data_80 = data.values[:, 2].astype(float)
    data_20 = data.values[:, 1].astype(float)
    data_5 = data.values[:, 0].astype(float)
    concentrationlist = [ '5 nM', '20 nM', '80 nM', '320 nM', '1300 nM', '3900 nM']
    all_data = [data_5, data_20, data_80, data_320, data_1300]
if (format == '*'):
    data_0 = data.values[:, 0].astype(float)
    data_5 = data.values[:, 1].astype(float)
    data_20 = data.values[:, 2].astype(float)
    data_80 = data.values[:, 3].astype(float)
    data_320 = data.values[:, 4].astype(float)
    data_1300 = data.values[:, 5].astype(float)
    concentrationlist = ['0 nM', '5 nM', '20 nM', '80 nM', '320 nM', '1300 nM']
    all_data = [data_0, data_5, data_20, data_80, data_320, data_1300]
if (format == 'RBM47'):
    data_0 = data.values[:, 0].astype(float)
    data_5 = data.values[:, 1].astype(float)
    data_20 = data.values[:, 2].astype(float)
    data_60 = data.values[:, 3].astype(float)
    data_240 = data.values[:, 4].astype(float)
    data_960 = data.values[:, 5].astype(float)
    data_3800 = data.values[:, 6].astype(float)
    concentrationlist = ['0 nM', '5 nM', '20 nM', '60 nM', '240 nM', '960 nM', '3800 nM']
    all_data = [data_0, data_5, data_20, data_60, data_240, data_960, data_3800]
if (format == 'CELF1'):
    data_0 = data.values[:, 0].astype(float)
    data_4 = data.values[:, 1].astype(float)
    data_16 = data.values[:, 2].astype(float)
    data_64 = data.values[:, 3].astype(float)
    data_130 = data.values[:, 4].astype(float)
    data_500 = data.values[:, 5].astype(float)
    data_1000 = data.values[:, 6].astype(float)
    data_2000 = data.values[:, 7].astype(float)
    concentrationlist = ['0 nM', '4 nM', '16 nM', '64 nM', '130 nM', '500 nM', '1000 nM', '2000 nM']
    all_data = [data_0, data_4, data_16, data_64, data_130, data_500, data_1000, data_2000]
if (format == 'MBNL1'):
    data_0 = data.values[:, 0].astype(float)
    data_1 = data.values[:, 1].astype(float)
    data_4 = data.values[:, 2].astype(float)
    data_13 = data.values[:, 3].astype(float)
    data_40 = data.values[:, 4].astype(float)
    data_121 = data.values[:, 5].astype(float)
    data_365 = data.values[:, 6].astype(float)
    data_1090 = data.values[:, 7].astype(float)
    data_3280 = data.values[:, 8].astype(float)
    data_9800 = data.values[:, 9].astype(float)
    concentrationlist = ['0 nM', '1 nM', '4 nM', '13 nM', '40 nM', '121 nM', '365 nM', '1090 nM', '3280 nM', '9800 nM']
    all_data = [data_0, data_1, data_4, data_13, data_40, data_121, data_365, data_1090, data_3280, data_9800]
if (format == 'MSI1'):
    data_0 = data.values[:, 0].astype(float)
    data_05 = data.values[:, 1].astype(float)
    data_2 = data.values[:, 2].astype(float)
    data_8 = data.values[:, 3].astype(float)
    data_16 = data.values[:, 4].astype(float)
    data_64 = data.values[:, 5].astype(float)
    data_256 = data.values[:, 6].astype(float)
    data_1000 = data.values[:, 7].astype(float)
    concentrationlist = ['0 nM', '0.5 nM', '2 nM', '8 nM', '16 nM', '64 nM', '256 nM', '1000 nM']
    all_data = [data_0, data_05, data_2, data_8, data_16, data_64, data_256, data_1000]
if (format == 'RBFOX2'):
    data_0 = data.values[:, 0].astype(float)
    data_1 = data.values[:, 1].astype(float)
    data_4 = data.values[:, 2].astype(float)
    data_14 = data.values[:, 3].astype(float)
    data_40 = data.values[:, 4].astype(float)
    data_121 = data.values[:, 5].astype(float)
    data_365 = data.values[:, 6].astype(float)
    data_1100 = data.values[:, 7].astype(float)
    data_3300 = data.values[:, 8].astype(float)
    data_9800 = data.values[:, 9].astype(float)
    concentrationlist = ['0 nM', '1 nM', '4 nM', '14 nM', '40 nM', '121 nM', '365 nM', '1100 nM', '3300 nM', '9800 nM']
    all_data = [data_0, data_1, data_4, data_14, data_40, data_121, data_365, data_1100, data_3300, data_9800]
allsevenmers = data.values[:, -1].astype(str)

log_data = np.log2(all_data)

s2 = pd.read_csv("Supplemental_Table_S2.csv", dtype=str, sep=',', index_col=0)

leilist = s2.values[1:5579, 6].astype(float)
hexmut = s2.values[1:5579, 0].astype(str)
sevenmerposlist = s2.values[1:5579, 1].astype(float)
seqlist = s2.values[1:5579, 2].astype(str)

hexmutlist = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']

all_lei = []
for hm in hexmutlist:
    pos_lei = []
    for pos in range(-3, 49):
        pos_lei.append(get7merlei(hm, pos))
    all_lei.append(pos_lei)
print("all lei done")
LRLEI_full = []
for hm in all_lei:
    LRLEI_pos = []
    for pos in hm:
        LRLEI_conc = []
        sevenmers = pos[0]
        leis = pos[1]
        for concentration in log_data:
            LR = []
            LEI = []
            for s in range(len(sevenmers)):
                #if (sevenmertorowmissing(sevenmers[s]) >= 0):
                LR.append(concentration[sevenmertorow(sevenmers[s])])
                LEI.append(leis[s])
            LRLEI_init = [LR, LEI]
            LRLEI_conc.append(LRLEI_init)
        LRLEI_pos.append(LRLEI_conc)
    LRLEI_full.append(LRLEI_pos)
print("LRLEI done")
from scipy.stats.stats import pearsonr
regresults = []
for hm in LRLEI_full:
    pos = []
    for pos in hm:
        conc = []
        for conc in pos:
            x = conc[0]
            y = conc[1]
            ##print(x)
            ##print(y)
            rval = pearsonr(x, y)[0]
            signedrsqr = signedRsquared(rval)
            pval = pearsonr(x, y)[1]
            result = [rval, signedrsqr, pval]
            conc.append(result)
        pos.append(conc)
    regresults.append(pos)
print("Regresults done")
LRLEI_final = []
results_final = []
for i in range (len(all_lei)):
    LRLEI_temp = []
    results_temp = []
    for pos in all_lei[i]:
        LRLEI_conc = []
        reg_conc = []
        sevenmers = pos[0]
        leis = pos[1]
        for concentration in log_data:
            LR = []
            LEI = []
            for s in range(len(sevenmers)):
                #if (sevenmertorowmissing(sevenmers[s]) >= 0):
                LR.append(concentration[sevenmertorow(sevenmers[s])])
                LEI.append(leis[s])
            LRLEI_init = [LR, LEI]
            rval = pearsonr(LR, LEI)[0]
            signedrsqr = signedRsquared(rval)
            pval = pearsonr(LR, LEI)[1]
            results = [rval, signedrsqr, pval]
            LRLEI_conc.append(LRLEI_init)
            reg_conc.append(results)
        LRLEI_temp.append(LRLEI_conc)
        results_temp.append(reg_conc)
    LRLEI_final.append(LRLEI_temp)
    results_final.append(results_temp)
    print(i)
print("Results final done")
hexmutlist = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
#concentrationlist = [ '5 nM', '20 nM', '80 nM', '320 nM', '1300 nM', '3900 nM']

#print("Explanation: In each row, the first value is pearson's R, the second value is the signed R^2, and the third value is the p value.")
#print_results(results_final)
"""
all_counts = [0, 0, 0, 0, 0, 0]
sig_counts = [0, 0, 0, 0, 0, 0]
for hm in range(len(results_final)):
    for pos in range(len(results_final[hm])):
        for conc in range(len(results_final[hm][pos])):
            all_counts[conc] = all_counts[conc] + 1
            if (results_final[hm][pos][conc][2] < 0.01):
                sig_counts[conc] = sig_counts[conc] + 1
"""
excel_write(results_final)
#for i in range(0, len(all_counts)):
#    print(sig_counts[i]/all_counts[i]*100)