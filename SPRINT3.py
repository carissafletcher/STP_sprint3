#!/c/Users/Jay/anaconda3/python

"""
Given a query nucleotide sequence in a fasta file:

1. Get input: nucleotide sequence as string or fasta
2. Query BLASTN with the sequence
3. Extract relevant information from BLASTN XML output
4. Identify best sequence match and corresponding gene
5. Identify gene orthologues and obtain their sequences
6. Construct a multiple sequence alignment and return in FASTA format

"""

#Import necessary modules
import re
from Bio.Blast import NCBIWWW
import xml.etree.ElementTree as ET
from Bio import Entrez, SeqIO

sequence_data = "GCTGTTCAGCGTTCTGCTGGAGCAGGGCCCCGGACGGCCAGGCGACGCCCCGCACACCGG"

#Get query sequence into .fasta format
def acquire_input():
    request_text = 'Sequence format is (1) string (2) fasta file: '
    input_type = input(request_text)
    while not ((input_type == 1) or (input_type == 2)):
        print('Please enter either 1 or 2.')
        input_type = input(request_text)

    if input_type == 1:
        query_sequence = input('Paste sequence as string here: ')
        ((query_sequence.replace(' ', '')).strip()).upper()
        ok_input = ['A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', '-']
        if [character not in ok_input for character in query_sequence]:
            print('Warning, unsupported character: ' + character)

    elif input_type == 2:
        filename = input('Enter full filename: ')
        while not (filename[-6:] == '.fasta'):
            print('File must be in .fasta format.')
            filename = input('Enter full filename: ')
        query_sequence = open(filename).read()

    return query_sequence

#Perform blastn search on query sequence, get XML output
def blast_search(sequence_data):
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data)
    print('BLAST query complete.')
    blast_results = result_handle.read()
    with open('results.xml', 'w') as save_file:
        save_file.write(blast_results)
    return blast_results

#Extract relevant info for each BLAST hit
def get_hit_list(raw_blast_output):
    list_of_hits = []
    root = raw_blast_output.getroot()

    for hit in root.findall('Hit'):
        hit_list = []
        gene_description = hit.findtext('Hit_def')
        gene_symbol = re.search('[(](.*)[)]', gene_description).group(1)
        transcript = hit.findtext('Hit_accession')
        percent_identity = hit/Hit_hsps/Hsp.findtext('Hsp_identity')
        hit_list.append([gene_symbol, transcript, percent_identity])
        list_of_hits.append(hit_list)

    return list_of_hits



blast_search(sequence_data)
