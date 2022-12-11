"""
    This is the main part of the BLAST heuristic algorithm.
    This is where all the functions are executed and files are imported.
"""

import pandas as pd
import numpy as np

import blastFunctions
import fastaProcessing
import matrixProcessing
import regex

# Defining file and library variables for passing in
# read files.
scoringMatrixFile = ''
queryFasta = ''
databaseFasta = ''
databaseDictionary = {}

try:
    # Opening matrix file for reading.
    fMatrix = open("matrices/BLOSUM62.matrix", mode='r', encoding='utf-8')

    # Opening query FASTA file for reading.
    fQuery = open("query.fasta")

    # Opening FASTA protein databse file for reading.
    fDatabase = open("proteinDatabaseSmaller.fasta")

    # Reading each file and storing it in a variable.
    scoringMatrixFile = fMatrix.read()
    queryFasta = fQuery.read()
    databaseFasta = fDatabase.read()
except Exception as e:
    print(e)
    print("Something went wrong when opening the file to read.")
finally:
    #Closing each file after reading was done.
    fMatrix.close()
    fQuery.close()
    fDatabase.close()

# Passing matrix file to the matrixProcessing function that is described in
# matrixProcessing.py. This function creates a pandas module dataFrame.
scoringMatrixDataFrame = matrixProcessing.readScoringMatrix(scoringMatrixFile)

# Passing FASTA database file to the fastaProcessing function that is described in
# fastaProcessing.py. This function converts FASTA file into a single FASTA string
# that is going to be used as a database. The second parameter can be set to TRUE
# to also return a dictionary that defines where each protein molecule sequence
# starts and ends.
proteinDatabase, databaseDictionary, databaseEndIndices = fastaProcessing.readFasta(databaseFasta, True)

# Passing query FASTA file to the fastaProcessing function that is described in
# fastaProcessing.py. This funciton converts FASTA file into a single FASTA string
# that is going to be used as a database.
querySequence = fastaProcessing.readFasta(queryFasta)

# Generating FASTA database sequence k-mers with their indexes.
kmerArray, kmerDictionary = blastFunctions.kmerGeneration(proteinDatabase, True)

# Generating high scoring sequence pairs between the database and the query.
HSSPdict = blastFunctions.generatingHSSP(scoringMatrixDataFrame, kmerDictionary, querySequence)

# This function returns a dictionary with the results after running the actual pBLAST algorithm.
resultDictionary = blastFunctions.searchAlgorithm(databaseEndIndices, HSSPdict, databaseDictionary, scoringMatrixDataFrame, proteinDatabase, querySequence)

# This function sorts the resultDictionary in ascending manner by the ['Score'] value.
dict = dict(sorted(resultDictionary.items(), key=lambda item: item[1]['Score'], reverse=True))



##########################################################################
alreadyPrinted = []
for result in dict:
    existBoolean = False
    object = dict[result]
    if (object['Score'] == 0):
        continue
    for header in alreadyPrinted:
        if (object['Header'] == header):
            existBoolean = True
    if (existBoolean == True):
        continue
    alreadyPrinted.append(object['Header'])

    print(object['Header'] + ' with score of:', object['Score'])
    print('Database hit: ' + object['AlignedDatabase'])
    print('              ' + object['Seperator'])
    print('Query hit:    ' + object['AlignedQuery'], '\n')


# This part is required for the test.py file.
#def getMatrix():
    #return scoringMatrixDataFrame
