# STP_sprint3

## **"THE HOMOLOGENATOR"**

Introduction to Programming - SPRINT3 Paired Coding Project\
Thomas Scott-Adams :  9627185\
Jay Miles          : 10806682


### **What does this app do?**

Given a query transcript (nucleotide) sequence in a fasta file:

1. Get user-defined query name, create new directory for output.
2. Get user-defined nucleotide sequence as string or fasta file.
3. Confirm that any .fasta file is of the correct format.
4. Submit sequence to BLASTN using Biopython NCBIWWW.
5. Extract relevant information from BLASTN XML output.
6. Identify closest sequence match and corresponding gene symbol.
7. Identify homologous protein sequences using Homologene API.
8. Format homologous protein sequences as a .fasta file.
9. Construct a multiple sequence alignment using Clustal Omega.


### **What are typical use cases for this app?**
Given a nucleotide sequence, generating a multiple sequence alignment (MSA) can be time-consuming and requires multiple steps (identifying the gene, finding orthologues, sequence alignment and clustering). Following multiple individual steps also introduces the potential for errors.

User story:\
-As a bioinformatics trainee, I want to automate the process of getting an MSA from a single sequence so I can reduce the potential for errors and focus on using/interpreting the MSA rather than generating it

User requirements:\
-Take a transcript (nucleotide) sequence as input (in FASTA format or as direct sequence entry)\ 
-Provide an MSA and phylogenetic tree from Clustal Omega as output


### **What data are required for this app to run?**
This app requires an internet connection and is written in Python 3. Python modules required for this app are listed in 'requirements.txt', and can be installed using the following console command:\
python -m pip install -r requirements.txt\
('python' may need to be replaced with the name of your specific Python version, e.g. 'python3')


### **How should this app be used?**
-If you are in the working directory, the app can be run from the command line with the command 'python homologenator.py'. Depending on your version of Python, you may need to replace 'python' with 'python3'.
-A valid email address must be provided for access to the web API services.\
-A unique user-defined query identifier must be provided in order to create a new folder for output files.\
-A transcript nucleotide sequence in a .fasta file. Sequences shorter than 1 line can be pasted in directly as a string.\
-Example sequences have been provided in the example_1.fasta and example_2.fasta files.


### **What does this app output?**
-All output files will be created in a new folder in the working directory. This folder will be named based on the query name.\
-Raw BLASTN output in .xml format (generated by BLASTN query using input sequence)\
-Raw .fasta file of homologous protein sequences from Homologene\
-Formatted .fasta file of homologous protein sequences from Homologene\
-Text file containing multiple sequence alignment (generated from homologous protein sequences by Clustal Omega)\
-Text file containing information for phylogenetic tree construction in Newick format\
-Image file showing phylogenetic tree


### **Limitations of the app**
-BLASTN search can take a considerable amount of time depending on server load\
-The Homologene database identifies a limited number of homologous protein sequences, which can also vary depending on the gene\
-The MSA output by Clustal Omega is not in .fasta format and so may not be appropriate input for downstream tools


*This app was made by Jay Miles and Tom Scott-Adams*
