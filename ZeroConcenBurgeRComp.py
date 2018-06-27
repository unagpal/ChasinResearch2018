#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 15:25:05 2018

@author: rishabbhatt
"""
import os, sys, string
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile

def read_excel(infile, tabname):
    df = pd.read_excel(infile, sheetname=tabname)  
    print("Column headings:")
    print(df.columns)
    return df
    
def write_excel(outfile, tabname, df):
    print("Writting the Content from Input file")
    writer = ExcelWriter(outfile)
    df.to_excel(writer, tabname, index=False)
    writer.save()
    writer.close()

def append_excel(outfile, tabname, df):
    print("Appending the Content from Input file")
    writer = ExcelWriter(outfile)
    df.to_excel(writer,tabname, index=False)
    writer.save()
    writer.close()
    
''' 
List File Names 
'''
def getFiles(dirname):
    files = os.listdir(dirname)
    #for file in files:
        #print(file)

'''
Join all excel file data into single output file
'''
def JoinExcelFiles(filelist, outfile, tabname):
    #print("Into Join Excel Files")
    AllDF = pd.DataFrame()
    
    for file in filelist:
        protein = os.path.split(file)[1].split(".")[0]
        #print("Reading and Appending file content for " + file)
        data = pd.read_excel(file)
        ldata = data.values.tolist()
        newlist = []
        for delement in ldata:
            delement.insert(0, protein)
            newlist.append(delement)
        AllDF = AllDF.append(pd.DataFrame(newlist))
    return AllDF

''' 
List File Names 
'''
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
        
'''
Main Method
Input - Directory Name
'''
if __name__ == '__main__':
    if(len(sys.argv) == 2):
        dirname = sys.argv[1]
        if(os.path.isdir(dirname)):
            print("Valid Directory Name " + dirname)
        else:
            print("Invalid Input Directory Name " + dirname)
    else:
        dirname = "/Users/rishabbhatt/Desktop/ChasinResearch"
        # dirname = os.path.curdir
       # print("Assuming current directory is input directory" + dirname)
        
    filelist = getExcelFiles(dirname)
 
    ofile = "Test.xlsx"
    outfile = os.path.join(dirname, ofile)
    outtabname = "Comparison"
    
    allData = JoinExcelFiles(filelist, outfile, outtabname)
    print(allData)
    
    write_excel(outfile, outtabname, allData)
    
