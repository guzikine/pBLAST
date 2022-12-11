# pBLAST
Python based pBLAST algorithm.

This program written in python executes the main steps in the BLAST algorithm.
To run this program, open the **main.py** file. The file contains 3 parameter variables 
that are crucial for determining how sensitive the BLAST algorithm will be.

The files contained in this project, that have the **.fasta** extension are used for
database and query sequence extraction.

Folder **matrices** contains a BLOSUM62 matrix with the **.matrix** extension. The program
uses this matrix for calculating match values between amino acids.
