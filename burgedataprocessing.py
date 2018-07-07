# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 13:31:38 2018

@author: udainagpal
"""

import xlrd  
import numpy as np

proteins = ['A1CF', 'BOLL', 'CELF1', 'CNOT4', 'CPEB1', 'DAZ3', 'DAZAP1', 'EIF4G2', 'ELAVL4', 'FUBP1', 'FUBP3', 'HNRNPA0', 'HNRNPA2B1', 'HNRNPCL1', 'HNRNPD', 'HNRNPDL', 'HNRNPF', 'HNRNPL', 'IGF2BP1', 'IGF2BP2', 'ILF2', 'KHDRBS2', 'KHDRBS3', 'KHSRP', 'MBNL1', 'MSI1', 'NOVA1', 'PABPN1L', 'PCBP1', 'PCBP2', 'PCBP4', 'PRR3', 'PTBP3', 'PUM1', 'RALYL', 'RBFOX2', 'RBFOX3', 'RBM4', 'RBM4B', 'RBM6', 'RBM15B', 'RBM22', 'RBM24', 'RBM25', 'RBM41', 'RBM45', 'RBM47', 'RBMS2', 'RBMS3', 'RC3H1', 'SF1', 'SFPQ', 'SNRPA', 'SRSF4', 'SRSF9', 'SRSF11', 'TARDBP', 'TIA1', 'TRA2A', 'TRNAU1AP', 'UNK', 'ZCRB1', 'ZFP36', 'ZNF326']
hms = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
book = xlrd.open_workbook('Burge LEI vs LR Correlations by Protein.xlsx')  
all_p = []
all_hm = []
all_conc = []
all_positions = []
all_concentrations = []
all_r2 = []
for sheet in book.sheets():
    r2temp = sheet.col_values(4)[1:]
    hmtemp = sheet.col_values(0)[1:]
    pvaltemp = sheet.col_values(5)[1:]
    concentrationtemp = sheet.col_values(2)[1:]
    positiontemp = sheet.col_values(1)[1:]
    all_positions.append(positiontemp)
    all_r2.append(r2temp)
    all_p.append(pvaltemp)
    all_hm.append(hmtemp)
    all_conc.append(concentrationtemp)
r2= [m for sub in all_r2 for m in sub]
r2 = np.array(r2)
hm = [k for sub in all_hm for k in sub]    
p = [j for sub in all_p for j in sub]
positions = [i for sub in all_positions for i in sub]
concentration = [l for sub in all_conc for l in sub]
for i in range(len(concentration)):
    if (concentration[i] not in all_concentrations):
        all_concentrations.append(concentration[i])
def maker2positive(r2):
    for i in range(len(r2)): 
        if (r2[i]<0):
            r2[i] = -r2[i]
    return r2;
def writer2byhmtoexcel(positiver2):
    posr2 = maker2positive(r2)
    hma, hmb, hmc, hmd, hme, hmf, hmg, hmh, hmi, hmj = ([] for i in range(10))
    allr2 = [hma, hmb, hmc, hmd, hme, hmf, hmg, hmh, hmi, hmj]
    for i in range(len(posr2)):
        allr2[hms.index(hm[i])].append(posr2[i])
    return allr2;
def percentsig(pcutoff, allp):
    sigarray = []
    for protein in allp:
        tot = len(protein)
        sig = 0
        for pval in protein:
            if pval < pcutoff:
                sig = sig+1
        sigarray.append(sig/tot*100)
    return sigarray;
def benjaminihochberg(p, q):
    p = np.array(p)
    sortedp = np.sort(p)
    m = float(len(sortedp))
    i = 1
    while(sortedp[i-1] <= i/m * q):
        i = i+1
    return sortedp[i-2];
def printsigpercentages(p, all_p):
    sig_percentages = percentsig(benjaminihochberg(p, 0.05), all_p)
    for i in range(len(sig_percentages)):
        print("Protein: " + proteins[i] + ", Percent Significant: " + str(sig_percentages[i]))
    return;
def printsigbyhmuniformp ():
    pcutoff = benjaminihochberg(p, 0.05)
    sigA, sigB, sigC, sigD, sigE, sigF, sigG, sigH, sigI, sigJ = ([] for i in range(10))
    sighm = [sigA, sigB, sigC, sigD, sigE, sigF, sigG, sigH, sigI, sigJ]
    hmCount = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    for i in range(len(p)):
        hmCount[hms.index(hm[i])] = hmCount[hms.index(hm[i])] + 1
        if (p[i] < pcutoff):
            sighm[hms.index(hm[i])].append(p[i])
    for i in range(len(hms)):
        print("HM: " + hms[i] + ", Number Significant: " + str(len(sighm[i])) + ", Percent Significant: " + str(len(sighm[i])/hmCount[i]*100))
    return;
def printsigbyconcuniformp ():
    pcutoff = benjaminihochberg(p, 0.05)
    sigconc = np.zeros(len(all_concentrations))
    totconc = np.zeros(len(all_concentrations))
    for i in range(len(p)):
        totconc[all_concentrations.index(concentration[i])] += 1
        if (p[i] < pcutoff):
            sigconc[all_concentrations.index(concentration[i])] += 1
    for i in range(len(all_concentrations)):
        print("Concentration: " + all_concentrations[i] + ", number significant: " + str(sigconc[i]) + ", number total: " + str(totconc[i]) + ", percent significant: " + str(sigconc[i]/totconc[i]*100))
    return;
def getpbyhm(HM):
    pval = []
    for i in range(len(p)):
        if (hm[i] == HM):
            pval.append(p[i])
    return pval;
def printsigbyhm (q):
    #proportion_sig_by_hm = []
    for i in range(len(hms)):
        sigcount = 0
        totcount = 0
        pvals = getpbyhm(hms[i])
        sigcutoff = benjaminihochberg(pvals, q)
        for j in pvals:
            totcount = totcount + 1
            if (j < sigcutoff):
                sigcount = sigcount + 1
        print("HM: " + hms[i] + ", significant cutoff: " + str(sigcutoff) + ", number significant: " + str(int(sigcount)) + ", percent significant: " + str(sigcount/totcount*100))
    return;
def printsigbypos ():
    all_pos = [i for i in range(-3, 0)] + [i for i in range (1, 47)]
    pcutoff = benjaminihochberg(p, 0.05)
    all_by_pos = np.zeros(49)
    sig_by_pos = np.zeros(49)
    for i in range(len(p)):
        all_by_pos[all_pos.index(positions[i])] += 1
        if (p[i] < pcutoff):
            sig_by_pos[all_pos.index(positions[i])] += 1
    for j in range(len(all_by_pos)):
        print("Position: " + str(all_pos[j]) + ", Number of Regressions: " + str(all_by_pos[j]) + ", Percent Significant: " + str(sig_by_pos[j]/all_by_pos[j]*100))
    return;
sig_percentages = percentsig(benjaminihochberg(p, 0.05), all_p)
#printsigpercentages(p, all_p)
#printsigbyhmuniformp()
#printsigbyhm(0.05)
#printsigbyconcuniformp()
#printsigbypos()
sig_perc = np.array(sig_percentages)
#print("Average percent significant: " + str(np.average(sig_perc)))
