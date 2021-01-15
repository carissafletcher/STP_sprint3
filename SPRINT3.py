#!/c/Users/Jay/anaconda3/python

"""
Given a query nucleotide sequence in a fasta file:

1. Get input: nucleotide sequence as string or fasta
2. Query BLASTN with the sequence
3. Extract relevant information from BLASTN XML output
4. Identify best sequence match and corresponding gene symbol
5. Identify gene orthologues and obtain their protein sequences
6. Construct a multiple sequence alignment and return in FASTA format

"""

#Import necessary modules
import re
from Bio.Blast import NCBIWWW
import xml.etree.ElementTree as ET
from Bio import Entrez, SeqIO #might not need this?

#Function 1: Get query sequence from string or .fasta file
def acquire_input():
    input_type = input('Sequence format is (1) string (2) fasta file? ')
    while not ((input_type == '1') or (input_type == '2')):
        print('Please enter either 1 or 2.')
        input_type = input('Sequence format is (1) string (2) fasta file? ')

    if input_type == '1':
        query_sequence = input('Paste sequence here: ')
        ((query_sequence.replace(' ', '')).strip()).upper() #NEED TO TEST
        ok_input = ['A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', '-']
        for character in query_sequence:
            if character not in ok_input:
                print('Warning, unsupported character: ' + character) #NEED TO TEST

    elif input_type == '2':
        filename = input('Enter filename: ')
        while not (filename[-6:] == '.fasta'):
            print('File must be .fasta type.')
            filename = input('Enter filename: ')
        query_sequence = open(filename).read()
        print(query_sequence)

    query_name = input('Please enter a name for this query: ')
    return query_sequence, query_name

#Function 2: Perform blastn search on query sequence, get XML output
def blast_search(query_sequence, query_name):
    print('Querying BLAST - this may take some time...')
    result_handle = NCBIWWW.qblast("blastn", "nt", query_sequence)
    print('BLAST query complete.')
    blast_results = result_handle.read()
    output_file = 'results_' + query_name + '.xml'
    with open(output_file, 'w') as save_file:
        save_file.write(blast_results)
    return output_file

#Function 3: Extract relevant info for each BLAST hit
def get_blast_hits(filename):
    list_of_hits = []
    tree = ET.parse(filename)
    root = tree.getroot()

    for hit in root.iter('Hit'):
        gene_description = hit.find('Hit_def').text
        transcript = hit.find('Hit_accession').text

        length = int(hit.find('Hit_len').text)
        for hsp in hit.iter('Hsp'):
            identity = int(hsp.find('Hsp_identity').text)
        percent_identity = 100*(identity/length)

        list_of_hits.append([gene_description, transcript, percent_identity])

    print('There are ' + str(len(list_of_hits)) + ' sequence matches.')
    return list_of_hits

#Function 4: Find sequence match with highest % identity
def find_best_match(list_of_hits):
    current_best = 0
    best_matches = []
    for hit in list_of_hits:
        if hit[2] > current_best:
            current_best = hit[2]

    for hit in list_of_hits:
        if hit[2] == current_best:
            best_matches.append(hit)
    
    gene_symbols = []
    for hit in best_matches:
        gene_desc = hit[0]
        split1 = ((gene_desc.split('('))[1:])[0]
        split2 = (split1.split(')'))[:1]
        gene_symbols.append(split2)

    unique_genes = []
    for gene in gene_symbols:
        if gene not in unique_genes:
            unique_genes.append(gene)
    if len(unique_genes) == 1:
        gene_output = unique_genes[0][0]
    else: #NEEDS TESTING - we need a query which outputs 2+ potential genes
        print('Multiple possible gene matches: ')
        print([gene for gene in unique_genes])
        selection = input('Please select one gene symbol: ')
        while selection not in unique_genes:
            print('Selection must be in list above.')
            selection = input('Please select one gene symbol: ')
        gene_output = selection

    print(gene_output)
    return gene_output

#Function 5: Identify gene orthologues
#Function 6: Get protein sequences into a single fasta
#Function 7: Construct an MSA


example_sequence = "GCTGTTCAGCGTTCTGCTGGAGCAGGGCCCCGGACGGCCAGGCGACGCCCCGCACACCGG" #from CACNA1F

#Main function to call other functions
def main():
    sequence, name = acquire_input()
    blast_output = blast_search(sequence, name)
    output_list = get_blast_hits(blast_output)
    find_best_match(output_list)

#Call to main function
if __name__ == '__main__':
    main()
