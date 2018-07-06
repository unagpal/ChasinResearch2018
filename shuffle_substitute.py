# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 11:47:05 2018

@author: udainagpal
"""
import random
from random import shuffle
from random import randint
import pandas as pd
import numpy as np
from pandas import DataFrame
from scipy.stats.stats import pearsonr

def scramble7mer (sevenmer, num_matching):
    iterations = 0
    tempsevenmer = [sevenmer[i] for i in range(len(sevenmer))]
    while (iterations < 200 and numbermatching(tempsevenmer, sevenmer)!=num_matching):        
        shuffle(tempsevenmer)        
        iterations += 1
    if (iterations == 200):
        return "ERROR";
    return ''.join(tempsevenmer);
#print(scramble7mer("ACTGCTG", 4))

def completelyscramble7mer (sevenmer):
    iterations = 0
    tempsevenmer = [sevenmer[i] for i in range(len(sevenmer))]
    while (iterations < 200 and iscompletelyscrambled(sevenmer, ''.join(tempsevenmer))==False):        
        shuffle(tempsevenmer)        
        iterations += 1
    if (iterations == 200):
        return "ERROR";
    return ''.join(tempsevenmer);
    
def numbermatchingkmerswithorder(sevenmer1, sevenmer2, k):
    num_matching = 0
    for i in range(0, 8-k):
        temp = 1
        for j in range(i, i+k):
            if (sevenmer1[j] != sevenmer2[j]):
                temp = 0
        num_matching += temp
    return num_matching;

def numbermatchingkmers(sevenmer1, sevenmer2, k):
    num_matching = 0
    for i in range(0, 8-k):
        kmer = sevenmer1[i:i+k]
        is_match = 0
        for j in range(0, 8-k):
            if (kmer == sevenmer2[j:j+k]):
                is_match = 1
        num_matching += is_match
    return num_matching;
    

def iscompletelyscrambled(sevenmer1, sevenmer2):
    return numbermatchingkmers(sevenmer1, sevenmer2, 3)==0;
    
def numbermatching(sevenmer1, sevenmer2):
    matching = 0
    for i in range(len(sevenmer1)):
        if (sevenmer1[i] == sevenmer2[i]):
            matching += 1
    return matching;

def substitutekmer (kmer, num_substitute):
    indices = [i for i in range(len(kmer))]
    shuffle(indices)
    pos_to_sub = indices[0:num_substitute]
    all_bases = ["A", "C", "G", "T"]
    all_substitutions = []
    for i in pos_to_sub:
        all_pos_sub = []
        base = kmer[i]
        for i in all_bases:
            if (i != base):
                all_pos_sub.append(i)
        rand_index = randint(0, len(all_pos_sub)-1)       
        all_substitutions.append(all_pos_sub[rand_index])
    newkmer = ""
    for i in range(len(kmer)):
        if i not in pos_to_sub:
            newkmer += kmer[i]
        else:
            newkmer += all_substitutions[pos_to_sub.index(i)]
    return newkmer;
#print(substitutekmer("ACTGCTG", 1))

def sevenmertorow (sevenmer):
    bases = ['A', 'C', 'G', 'T']
    num_str = ''
    for i in sevenmer:
        for j in range(len(bases)):
            if (i == bases[j]):
                num_str = num_str + str(j)
    num_str = str(num_str)
    ans =  int(num_str[0])* (4**6) + int(num_str[1]) * (4**5) + int(num_str[2])* (4**4) + int(num_str[3]) * (4**3) + int(num_str[4]) * (4**2) + int(num_str[5])*4 + int(num_str[6])
    return ans;

def signedRsquared(r):
    return(abs(r**2));

def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)

#Reading in data and performing analysis
def scramble_results():
    protein_count = 0
    filenames = ["A1CF.csv", "BOLL.csv", "CELF1.csv", "CNOT4.csv", "CPEB1.csv", "DAZ3.csv", "DAZAP1.csv", "EIF4G2.csv", "ELAVL4.csv", "FUBP1.csv", "FUBP3.csv", "HNRNPA0.csv", "HNRNPA2B1.csv", "HNRNPCL1.csv", "HNRNPD.csv", "HNRNPDL.csv", "HNRNPF.csv", "HNRNPL.csv", "IGF2BP1*.csv", "IGF2BP2.csv", "ILF2.csv", "KHDRBS2.csv", "KHDRBS3.csv", "KHSRP.csv", "MBNL1.csv", "MSI1.csv", "NOVA1.csv", "PABPN1L.csv", "PCBP1.csv", "PCBP2.csv", "PCBP4.csv", "PRR3.csv", "PTBP3.csv", "PUM1.csv", "RALYL.csv", "RBFOX2.csv", "RBFOX3.csv", "RBM4.csv", "RBM4B.csv", "RBM6.csv", "RBM15B.csv", "RBM22.csv", "RBM24.csv", "RBM25.csv", "RBM41.csv", "RBM45.csv", "RBM47.csv", "RBMS2.csv", "RBMS3.csv", "RC3H1.csv", "SF1.csv", "SFPQ.csv", "SNRPA.csv", "SRSF4.csv", "SRSF9.csv", "SRSF11.csv", "TARDBP.csv", "TIA1.csv", "TRA2A.csv", "TRNAU1AP.csv", "UNK.csv", "ZCRB1.csv", "ZFP36.csv", "ZNF326.csv"]
    #filenames = ["A1CF.csv", "BOLL.csv", "CELF1.csv", "CNOT4.csv", "CPEB1.csv", "DAZ3.csv", "DAZAP1.csv", "EIF4G2.csv", "ELAVL4.csv", "FUBP1.csv", "FUBP3.csv", "HNRNPA0.csv", "HNRNPA2B1.csv", "HNRNPCL1.csv", "HNRNPD.csv", "HNRNPDL.csv", "HNRNPF.csv", "HNRNPL.csv", "IGF2BP1*.csv", "IGF2BP2.csv", "ILF2.csv", "KHDRBS2.csv", "KHDRBS3.csv", "KHSRP.csv", "MBNL1.csv", "MSI1.csv", "NOVA1.csv", "PABPN1L.csv", "PCBP1.csv", "PCBP2.csv", "PCBP4.csv", "PRR3.csv", "PTBP3.csv", "PUM1.csv", "RALYL.csv", "RBFOX2.csv", "RBFOX3.csv"]
    prot, conc, num_match, r, r2, p = ([] for i in range(6))
    for f in filenames:
        data = pd.read_csv(f, dtype=str, sep=',', index_col = 0, header = None)   
        all_sevenmers = data.values[1:, -1].astype(str)
        if(len(all_sevenmers)==16384):
            protein_count += 1
            protein = f[:-4]
            print(protein)
            print(protein_count)
            header_row = data.values[0, :].astype(str)
            concentrationlist, all_data = ([] for i in range(2))
            for i in range(len(header_row)):
                if "nM" in header_row[i] or header_row[i].isdigit():
                    concentrationlist.append(header_row[i])
                    all_data.append(data.values[1:, i].astype(float))
            for i in range(len(concentrationlist)):
                for j in range(6):
                    orig_affinities = []
                    mod_affinities = []
                    for k in range(len(all_data[i])):
                        orig_sevenmer = all_sevenmers[k]
                        orig_affinity = all_data[i][k]
                        mod_sevenmer = scramble7mer(orig_sevenmer, j)
                        if (mod_sevenmer != "ERROR"):
                            orig_affinities.append(orig_affinity)
                            mod_affinity = all_data[i][sevenmertorow(mod_sevenmer)]
                            mod_affinities.append(mod_affinity)
                    results = pearsonr(orig_affinities, mod_affinities)
                    rval = results[0]
                    signedrsqr = signedRsquared(rval)
                    pval = results[1]
                    prot.append(protein)
                    conc.append(concentrationlist[i])
                    num_match.append(j)
                    r.append(rval)
                    r2.append(signedrsqr)
                    p.append(pval)
                    #print(protein + ", " + "Concentration: " + concentrationlist[i] + ", Number Matching: " + str(j) + ", R: " + str(rval) + ", R^2: " + str(signedrsqr) + ", p: " + str(pval))
    df = DataFrame({'Protein': prot, 'Concentration': conc, 'Number Matching Positions': num_match, 'r':r, 'r2':r2, 'p':p})
    df.to_excel('ScrambleResults.xlsx', sheet_name = 'final', columns = ['Protein', 'Concentration', 'Number Matching Positions', 'r', 'r2', 'p'], index = False)

def substitution_results():
    protein_count = 0
    filenames = ["A1CF.csv", "BOLL.csv", "CELF1.csv", "CNOT4.csv", "CPEB1.csv", "DAZ3.csv", "DAZAP1.csv", "EIF4G2.csv", "ELAVL4.csv", "FUBP1.csv", "FUBP3.csv", "HNRNPA0.csv", "HNRNPA2B1.csv", "HNRNPCL1.csv", "HNRNPD.csv", "HNRNPDL.csv", "HNRNPF.csv", "HNRNPL.csv", "IGF2BP1*.csv", "IGF2BP2.csv", "ILF2.csv", "KHDRBS2.csv", "KHDRBS3.csv", "KHSRP.csv", "MBNL1.csv", "MSI1.csv", "NOVA1.csv", "PABPN1L.csv", "PCBP1.csv", "PCBP2.csv", "PCBP4.csv", "PRR3.csv", "PTBP3.csv", "PUM1.csv", "RALYL.csv", "RBFOX2.csv", "RBFOX3.csv", "RBM4.csv", "RBM4B.csv", "RBM6.csv", "RBM15B.csv", "RBM22.csv", "RBM24.csv", "RBM25.csv", "RBM41.csv", "RBM45.csv", "RBM47.csv", "RBMS2.csv", "RBMS3.csv", "RC3H1.csv", "SF1.csv", "SFPQ.csv", "SNRPA.csv", "SRSF4.csv", "SRSF9.csv", "SRSF11.csv", "TARDBP.csv", "TIA1.csv", "TRA2A.csv", "TRNAU1AP.csv", "UNK.csv", "ZCRB1.csv", "ZFP36.csv", "ZNF326.csv"]
    #filenames = ["A1CF.csv", "BOLL.csv", "CELF1.csv", "CNOT4.csv", "CPEB1.csv", "DAZ3.csv", "DAZAP1.csv", "EIF4G2.csv", "ELAVL4.csv", "FUBP1.csv", "FUBP3.csv", "HNRNPA0.csv", "HNRNPA2B1.csv", "HNRNPCL1.csv", "HNRNPD.csv", "HNRNPDL.csv", "HNRNPF.csv", "HNRNPL.csv", "IGF2BP1*.csv", "IGF2BP2.csv", "ILF2.csv", "KHDRBS2.csv", "KHDRBS3.csv", "KHSRP.csv", "MBNL1.csv", "MSI1.csv", "NOVA1.csv", "PABPN1L.csv", "PCBP1.csv", "PCBP2.csv", "PCBP4.csv", "PRR3.csv", "PTBP3.csv", "PUM1.csv", "RALYL.csv", "RBFOX2.csv", "RBFOX3.csv"]
    prot, conc, num_sub, r, r2, p = ([] for i in range(6))
    for f in filenames:
        data = pd.read_csv(f, dtype=str, sep=',', index_col = 0, header = None)   
        all_sevenmers = data.values[1:, -1].astype(str)
        if(len(all_sevenmers)==16384):
            protein_count += 1
            protein = f[:-4]
            print(protein)
            print(protein_count)
            header_row = data.values[0, :].astype(str)
            concentrationlist, all_data = ([] for i in range(2))
            for i in range(len(header_row)):
                if "nM" in header_row[i] or header_row[i].isdigit():
                    concentrationlist.append(header_row[i])
                    all_data.append(data.values[1:, i].astype(float))
            for i in range(len(concentrationlist)):
                for j in range(1, 7):
                    orig_affinities = []
                    mod_affinities = []
                    for k in range(len(all_data[i])):
                        orig_sevenmer = all_sevenmers[k]
                        orig_affinity = all_data[i][k]
                        mod_sevenmer = substitutekmer(orig_sevenmer, j)
                        orig_affinities.append(orig_affinity)
                        mod_affinity = all_data[i][sevenmertorow(mod_sevenmer)]
                        mod_affinities.append(mod_affinity)
                    results = pearsonr(orig_affinities, mod_affinities)
                    rval = results[0]
                    signedrsqr = signedRsquared(rval)
                    pval = results[1]
                    prot.append(protein)
                    conc.append(concentrationlist[i])
                    num_sub.append(j)
                    r.append(rval)
                    r2.append(signedrsqr)
                    p.append(pval)
                    #print(protein + ", " + "Concentration: " + concentrationlist[i] + ", Number Matching: " + str(j) + ", R: " + str(rval) + ", R^2: " + str(signedrsqr) + ", p: " + str(pval))
    df = DataFrame({'Protein': prot, 'Concentration': conc, 'Number of Substitutions': num_sub, 'r':r, 'r2':r2, 'p':p})
    df.to_excel('SubstitutionResults.xlsx', sheet_name = 'final', columns = ['Protein', 'Concentration', 'Number of Substitutions', 'r', 'r2', 'p'], index = False)

def complete_scramble_results():
    protein_count = 0
    filenames = ["A1CF.csv", "BOLL.csv", "CELF1.csv", "CNOT4.csv", "CPEB1.csv", "DAZ3.csv", "DAZAP1.csv", "EIF4G2.csv", "ELAVL4.csv", "FUBP1.csv", "FUBP3.csv", "HNRNPA0.csv", "HNRNPA2B1.csv", "HNRNPCL1.csv", "HNRNPD.csv", "HNRNPDL.csv", "HNRNPF.csv", "HNRNPL.csv", "IGF2BP1*.csv", "IGF2BP2.csv", "ILF2.csv", "KHDRBS2.csv", "KHDRBS3.csv", "KHSRP.csv", "MBNL1.csv", "MSI1.csv", "NOVA1.csv", "PABPN1L.csv", "PCBP1.csv", "PCBP2.csv", "PCBP4.csv", "PRR3.csv", "PTBP3.csv", "PUM1.csv", "RALYL.csv", "RBFOX2.csv", "RBFOX3.csv", "RBM4.csv", "RBM4B.csv", "RBM6.csv", "RBM15B.csv", "RBM22.csv", "RBM24.csv", "RBM25.csv", "RBM41.csv", "RBM45.csv", "RBM47.csv", "RBMS2.csv", "RBMS3.csv", "RC3H1.csv", "SF1.csv", "SFPQ.csv", "SNRPA.csv", "SRSF4.csv", "SRSF9.csv", "SRSF11.csv", "TARDBP.csv", "TIA1.csv", "TRA2A.csv", "TRNAU1AP.csv", "UNK.csv", "ZCRB1.csv", "ZFP36.csv", "ZNF326.csv"]
    #filenames = ["A1CF.csv", "BOLL.csv", "CELF1.csv", "CNOT4.csv", "CPEB1.csv", "DAZ3.csv", "DAZAP1.csv", "EIF4G2.csv", "ELAVL4.csv", "FUBP1.csv", "FUBP3.csv", "HNRNPA0.csv", "HNRNPA2B1.csv", "HNRNPCL1.csv", "HNRNPD.csv", "HNRNPDL.csv", "HNRNPF.csv", "HNRNPL.csv", "IGF2BP1*.csv", "IGF2BP2.csv", "ILF2.csv", "KHDRBS2.csv", "KHDRBS3.csv", "KHSRP.csv", "MBNL1.csv", "MSI1.csv", "NOVA1.csv", "PABPN1L.csv", "PCBP1.csv", "PCBP2.csv", "PCBP4.csv", "PRR3.csv", "PTBP3.csv", "PUM1.csv", "RALYL.csv", "RBFOX2.csv", "RBFOX3.csv"]
    prot, conc, r, r2, p = ([] for i in range(5))
    for f in filenames:
        data = pd.read_csv(f, dtype=str, sep=',', index_col = 0, header = None)   
        all_sevenmers = data.values[1:, -1].astype(str)
        if(len(all_sevenmers)==16384):
            protein_count += 1
            protein = f[:-4]
            print(protein)
            print(protein_count)
            header_row = data.values[0, :].astype(str)
            concentrationlist, all_data = ([] for i in range(2))
            for i in range(len(header_row)):
                if "nM" in header_row[i] or header_row[i].isdigit():
                    concentrationlist.append(header_row[i])
                    all_data.append(data.values[1:, i].astype(float))  
            for i in range(len(concentrationlist)):
                orig_affinities = []
                mod_affinities = []
                for k in random.sample(range(len(all_data[i])), 2500):
                    orig_sevenmer = all_sevenmers[k]
                    orig_affinity = all_data[i][k]
                    mod_sevenmer = completelyscramble7mer(orig_sevenmer)
                    if (mod_sevenmer != "ERROR"):
                        orig_affinities.append(orig_affinity)
                        mod_affinity = all_data[i][sevenmertorow(mod_sevenmer)]
                        mod_affinities.append(mod_affinity)
                results = pearsonr(orig_affinities, mod_affinities)
                rval = results[0]
                signedrsqr = signedRsquared(rval)
                pval = results[1]
                prot.append(protein)
                conc.append(concentrationlist[i])
                r.append(rval)
                r2.append(signedrsqr)
                p.append(pval)
                #print(protein + ", " + "Concentration: " + concentrationlist[i] + ", Number Matching: " + str(j) + ", R: " + str(rval) + ", R^2: " + str(signedrsqr) + ", p: " + str(pval))
    df = DataFrame({'Protein': prot, 'Concentration': conc, 'r':r, 'r2':r2, 'p':p})
    df.to_excel('ScrambleResults.xlsx', sheet_name = 'final', columns = ['Protein', 'Concentration', 'r', 'r2', 'p'], index = False)

def mult_reg(X, y):
    X = np.array(X)
    X = X.transpose()
    y = np.array(y)
    from sklearn.model_selection import train_test_split  
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)
    from sklearn.linear_model import LinearRegression  
    regressor = LinearRegression()  
    regressor.fit(X_train, y_train)
    y_pred = regressor.predict(X_test)
    return [regressor.intercept_, regressor.coef_, regressor.score(X_test, y_test)]
#print(mult_reg([[1, 2, 3, 4, 5], [1, 3, 5, 7, 9]], [0, -1, -2, -3, -4])[1])

def multi_regression_results():
    protein_count = 0
    filenames = ["A1CF.csv", "BOLL.csv", "CELF1.csv", "CNOT4.csv", "CPEB1.csv", "DAZ3.csv", "DAZAP1.csv", "EIF4G2.csv", "ELAVL4.csv", "FUBP1.csv", "FUBP3.csv", "HNRNPA0.csv", "HNRNPA2B1.csv", "HNRNPCL1.csv", "HNRNPD.csv", "HNRNPDL.csv", "HNRNPF.csv", "HNRNPL.csv", "IGF2BP1*.csv", "IGF2BP2.csv", "ILF2.csv", "KHDRBS2.csv", "KHDRBS3.csv", "KHSRP.csv", "MBNL1.csv", "MSI1.csv", "NOVA1.csv", "PABPN1L.csv", "PCBP1.csv", "PCBP2.csv", "PCBP4.csv", "PRR3.csv", "PTBP3.csv", "PUM1.csv", "RALYL.csv", "RBFOX2.csv", "RBFOX3.csv", "RBM4.csv", "RBM4B.csv", "RBM6.csv", "RBM15B.csv", "RBM22.csv", "RBM24.csv", "RBM25.csv", "RBM41.csv", "RBM45.csv", "RBM47.csv", "RBMS2.csv", "RBMS3.csv", "RC3H1.csv", "SF1.csv", "SFPQ.csv", "SNRPA.csv", "SRSF4.csv", "SRSF9.csv", "SRSF11.csv", "TARDBP.csv", "TIA1.csv", "TRA2A.csv", "TRNAU1AP.csv", "UNK.csv", "ZCRB1.csv", "ZFP36.csv", "ZNF326.csv"]
    prot, conc, r2, intercept, coeff1, coeff2, coeff3, coeff4, coeff5, coeff6 = ([] for i in range(10))
    coeffs = [coeff1, coeff2, coeff3, coeff4, coeff5, coeff6]
    for f in filenames:
        data = pd.read_csv(f, dtype=str, sep=',', index_col = 0, header = None)   
        all_sevenmers = data.values[1:, -1].astype(str)
        if(len(all_sevenmers)==16384):
            protein_count += 1
            protein = f[:-4]
            print(protein)
            print(protein_count)
            header_row = data.values[0, :].astype(str)
            concentrationlist, all_data = ([] for i in range(2))
            for i in range(len(header_row)):
                if "nM" in header_row[i] or header_row[i].isdigit():
                    concentrationlist.append(header_row[i])
                    all_data.append(data.values[1:, i].astype(float))  
            for i in range(len(concentrationlist)):
                matching_1mer, matching_2mer, matching_3mer, matching_4mer, matching_5mer, matching_6mer, absdiff = ([] for counter in range(7))
                all_matching = [matching_1mer, matching_2mer, matching_3mer, matching_4mer, matching_5mer, matching_6mer]
                for k in range(100000):
                    sevenmer_indices = random.sample(range(len(all_data[i])), 2)                    
                    sevenmer_one = all_sevenmers[sevenmer_indices[0]]
                    sevenmer_two = all_sevenmers[sevenmer_indices[1]]
                    one_affinity = all_data[i][sevenmer_indices[0]]
                    two_affinity = all_data[i][sevenmer_indices[1]]
                    abs_diff = abs(one_affinity - two_affinity)
                    for kmer in range(1, 7):
                        all_matching[kmer-1].append(numbermatchingkmers(sevenmer_one, sevenmer_two, kmer))
                    absdiff.append(abs_diff)
                reg_results = mult_reg(all_matching, absdiff)
                prot.append(protein)
                conc.append(concentrationlist[i])
                r2.append(reg_results[2])
                intercept.append(reg_results[0])
                for i in range(len(coeffs)):
                    coeffs[i].append(reg_results[1][i])
    df = DataFrame({'Protein': prot, 'Concentration': conc, 'r2':r2, 'Intercept':intercept, 'coeff#1mers':coeff1, 'coeff#2mers':coeff2, 'coeff#3mers':coeff3, 'coeff#4mers':coeff4, 'coeff#5mers':coeff5, 'coeff#6mers':coeff6})
    df.to_excel('MultipleRegressionResults100000.xlsx', sheet_name = 'final', columns = ['Protein', 'Concentration', 'r2', 'Intercept', 'coeff#1mers', 'coeff#2mers', 'coeff#3mers', 'coeff#4mers', 'coeff#5mers', 'coeff#6mers'], index = False)
    return;

#multi_regression_results()
#systematic_sampling_results()
#complete_scramble_results()
#substitution_results()