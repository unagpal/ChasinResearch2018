#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 09:51:39 2018

@author: rishabbhatt
"""
import os, sys, string 
import pandas as pd
from pandas import ExcelFile
from pandas import ExcelWriter
import scipy
from scipy.stats.stats import pearsonr

dirname = "/Users/rishabbhatt/Desktop/DifferentProteins"
outfile = "StatsGraphsConcenComparisonByProtein1.xlsx"
tabName = "Statistics" 


def run_Combos(dataframe, protein):
    SumData = pd.DataFrame() #sumarized data with r, r2 and p values
    rCorr = []
    r2Corr = []
    comp = [] # to ensure which values are being checked
    pVal = []
    headerList  = len(dataframe.columns) #run the correlations against all possible combinations
    i = 1; 
    while (i  < headerList - 1): # calculates all values for every combination 
        j = i + 1
        while(j < headerList):
            comp.append("(" + str(protein) + ")"+ str(dataframe.columns[i]) + "vs." + str( dataframe.columns[j]))
            rCorr.append(pearsonr(dataframe[dataframe.columns[i]], dataframe[dataframe.columns[j]])[0])
            r2Corr.append(calcr2((pearsonr(dataframe[dataframe.columns[i]], dataframe[dataframe.columns[j]])[0])))
            pVal.append(pearsonr(dataframe[dataframe.columns[i]], dataframe[dataframe.columns[j]])[1])
            j += 1
        i+= 1
    #put all these lists as columns into the dataframe   
    
    s1 = pd.Series(comp)
    #print(s1)
    s2 = pd.Series(rCorr)
    #print (s2)
    s3 = pd.Series(r2Corr)
    s4 = pd.Series(pVal)
    SumData = pd.concat([s1, s2, s3, s4],  axis=1)
    #print(SumData)
    return SumData
    #write_excel(SumData)
    
def JoinExcelFiles(filelist, outfile, tabname):
    #print("Into Join Excel Files")
    AllDF = pd.DataFrame()
    
    for file in filelist:
        protein = os.path.split(file)[1].split(".")[0]
        #print("Reading and Appending file content for " + file)
        df = pd.read_excel(file)
        run_Combos(df, protein)
        AllDF = AllDF.append(run_Combos(df, protein))
    #return AllDF
    write_excel(AllDF, outfile)


def getExcelFiles(dirname): #stores all the directories of the desired excel files
    print("Into Excel File " + dirname)
    xlFileList = []
    files = os.listdir(dirname)
    for file in files:
        # print(file)
        ext = os.path.splitext(file)[1]
        if ( (ext == ".xls") or (ext == ".xlsx") ):
            fullfile = os.path.join(dirname, file)
            xlFileList.append(fullfile)
            print("Excel File Name is " + fullfile)
    return xlFileList

def calcr2(r):
    if (r > 0):
        return r**2
    elif (r < 0):
        return - (r**2)
    else:
        return 0;

def write_excel(dataframe, outfile): 
    writer = ExcelWriter(outfile)
    dataframe.to_excel(writer, tabName, index = False)
    writer.save()
    writer.close()


def main(): 
    JoinExcelFiles(getExcelFiles(dirname), outfile, tabName)


if __name__ == '__main__': 
    if(len(sys.argv) == 2):
        dirname = sys.argv[1]
        if(os.path.isdir(dirname)):
            print("Valid Directory Name " + dirname)
        else:
            print("Invalid Input Directory Name " + dirname)
    else:
        dirname = "/Users/rishabbhatt/Desktop/DifferentProteins"
    main()