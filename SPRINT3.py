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
import os
import requests
import time
import xml.etree.ElementTree as ET

from Bio import Entrez, Phylo, SeqIO
from Bio.Blast import NCBIWWW

#Function 1: Function to check whether a file is in fasta format
def fasta_check(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)

#
def get_query_name():
    query_name = input('Please enter a name for this query: ')
    while not query_name:
        query_name = input('Please enter a name for this query: ')
    
    current_dir = os.getcwd()
    path_name = current_dir + '/' + query_name

    sentry = False
    while sentry == False:
        if not os.path.exists(path_name):
            os.mkdir(path_name)
            sentry = True
        else:
            print('Folder already exists, please choose another name.')
            query_name = input('Please enter a name for this query: ')
            path_name = current_dir + '/' + query_name
        
    return query_name, path_name

#Function 2: Get query sequence from string or .fasta file
def acquire_input():
    input_type = input('DNA sequence format is (1) string (2) fasta file? ')
    while not ((input_type == '1') or (input_type == '2')):
        print('Please enter either 1 or 2.')
        input_type = input('DNA sequence format is (1) string (2) fasta file? ')

    if input_type == '1':
        string_input = (input('Paste sequence here: ')).upper()
        query_sequence = (string_input.replace(' ', '')).strip()
        print('Sequence to query: ' + query_sequence)
        ok_input = ['A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', '-']
        for character in query_sequence:
            if character not in ok_input:
                print('Warning, unsupported character: ' + character) #warning only, stills runs query

    elif input_type == '2':
        filename = input('Enter filename: ')
        while not (filename.endswith('.fasta')):
            print('File must be .fasta type.')
            filename = input('Enter filename: ')
        while fasta_check(filename) == False:
            print('Data must be in fasta format.')
            filename = input('Enter filename: ')
        query_sequence = open(filename).read()
        print(query_sequence)

    return query_sequence

#Function 3: Perform blastn search on query sequence, get XML output
def blast_search(query_sequence, query_name, path_name):
    print('Querying BLAST - this may take some time...')
    result_handle = NCBIWWW.qblast("blastn", "nt", query_sequence)
    print('BLAST query complete.')
    blast_results = result_handle.read()
    result_handle.close()
    output_file = query_name + '_blast_output.xml'
    output_path = os.path.join(path_name, output_file)
    with open(output_path, 'w') as save_file:
        save_file.write(blast_results)
    return output_path

#Function 4: Extract relevant info from BLAST output
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

#Function 5: Find sequence match with highest % identity and get gene symbol
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
        gene_symbols.append(split2[0])

    unique_genes = []
    for gene in gene_symbols:
        if gene not in unique_genes:
            unique_genes.append(gene)
    if len(unique_genes) == 1:
        gene_output = unique_genes[0]
        print('Closest gene match: ' + gene_output)
    else: #NEEDS TESTING - not sure about this section
        print('Multiple possible gene matches: ')
        print([gene for gene in unique_genes])
        gene_output = input('Please select one gene symbol: ')
        while gene_output not in unique_genes:
            print('Selection must be in list above.')
            gene_output = input('Please select one gene symbol: ')

    return gene_output

#Function 6: Identify gene homologues and return protein sequences as .fasta
def find_homologues(gene_output):
    Entrez.email = 'jomiles26@googlemail.com'
    search_handle = Entrez.esearch(db='homologene', term=(gene_output + '[Gene Name] AND Homo sapiens[Organism]'))
    search_record = Entrez.read(search_handle)
    search_handle.close()

    hg_id = search_record['IdList'][0]
    print('Homologene ID: ' + hg_id)

    hg_handle = Entrez.efetch(db='homologene', id = hg_id, rettype = 'homologene', retmode='text')
    hg_record = hg_handle.readlines()
    hg_handle.close()
    accession_list = []
    for line in hg_record:
        line_strip = line.strip()
        if line_strip[-3:] == ' aa':
            accession = (line_strip.split(' '))[-2]
            accession_list.append(accession)
    print('Transcript accession numbers: ', [acs for acs in accession_list])

    fasta_handle = Entrez.efetch(db='homologene', id = hg_id, rettype = 'fasta', retmode='text')
    fasta_record = fasta_handle.read()
    fasta_handle.close()

    return fasta_record

#Function 7: Format fasta file with species names
def format_fasta(gene_output, fasta_record, path_name):
    first_seq = int(fasta_record.find('>'))
    clipped_fasta_record = (fasta_record[first_seq:]).strip()
    fasta_file = gene_output + '_raw_homologues.fasta'

    og_fasta_path = os.path.join(path_name, fasta_file)
    with open(og_fasta_path, 'w') as output_object:
        output_object.write(clipped_fasta_record)
    with open(og_fasta_path, 'r') as input_object:
        fasta_homologues = input_object.readlines()

    i = 0
    for line in fasta_homologues:
        if line[0] == '>':
            split_line = line.split('|')
            handle = Entrez.efetch(db="protein", id=split_line[3], rettype="gb", retmode="text")
            record = SeqIO.read(handle, "gb")
            handle.close()
            species = (record.annotations['organism']).replace(' ', '_')
            line_end = line.split('ref|')
            if '[' in line_end[-1]:
                remove_dup_species = line_end[-1].split('[')
                new_line = '>' + gene_output + '|' +  str(species) + '|' + remove_dup_species[0] + '\n'
            else:
                new_line = '>' + gene_output + '|' +  str(species) + '|' + line_end[-1]
            fasta_homologues[i] = new_line
        i += 1

    homologues_file = gene_output + '_formatted_homologues.fasta'
    homologues_path = os.path.join(path_name, homologues_file)
    with open(homologues_path, 'w') as output_object:
        [output_object.write(line) for line in fasta_homologues]

    return homologues_path

#Function 8: Construct an MSA using Clustal Omega
def construct_MSA(homologues_file, gene_symbol, path_name):
    
    email_address = input("Please enter a valid email address: ")
    while "@" not in email_address:
        print("Web API requires a valid email address.")
        email_address = input("Please enter a valid email address: ")

    with open(homologues_file,'r') as homologues_object:
        fasta_transcripts = homologues_object.read()

    run = ''
    while not str(run) == '<Response [200]>':
        run = requests.post('https://www.ebi.ac.uk/Tools/services/rest/clustalo/run', data = {'email':email_address, 'sequence':fasta_transcripts})
        #print(run)
        try:
            run.raise_for_status()
        except requests.exceptions.HTTPError:
            print('Error whilst submitting query. Please ensure email address is valid.')
            email_address = input("Please enter a valid email address: ")            
    job_id_bytes = run.content
    job_id = job_id_bytes.decode('utf-8')

    status = requests.get('https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/'+ str(job_id))
    status_bytes = status.content
    status_decoded = status_bytes.decode('utf-8')

    while status_decoded != "FINISHED":
        status = requests.get('https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/'+ str(job_id))
        status_bytes = status.content
        status_decoded = status_bytes.decode('utf-8')
        print(status_decoded)
        if status_decoded == "RUNNING":
            time.sleep(10)

    msa_request = requests.get('https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/'+ str(job_id) +'/aln-clustal_num')
    msa_bytes = msa_request.content
    msa = msa_bytes.decode('utf-8')
    msa_file = gene_symbol + '_msa.fasta'
    msa_path = os.path.join(path_name, msa_file)
    with open(msa_path, 'w') as msa_file_object:
        msa_file_object.write(msa)

    phylo_request = requests.get('https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/'+ str(job_id) +'/phylotree')
    phylo_bytes = phylo_request.content
    phylo = phylo_bytes.decode('utf-8')
    phylo_file = gene_symbol + '_phylo.txt'
    phylo_path = os.path.join(path_name, phylo_file)
    with open(phylo_path, 'w') as phylo_file_object:
        phylo_file_object.write(phylo)

    return msa, phylo

example_sequence = "GCTGTTCAGCGTTCTGCTGGAGCAGGGCCCCGGACGGCCAGGCGACGCCCCGCACACCGG" #from CACNA1F
test_string = "  gCTGTTCAGCGTTCTGCtggAGCa GGGCEFCGGACGGCCAGGCGAC  GCCCCICACACCgg " #some gaps, lowercase, and unsupported characters

#Main function to call other functions
def main():
    query_name, path_name = get_query_name()
    sequence = acquire_input()
    blast_output = blast_search(sequence, query_name, path_name)
    list_of_hits = get_blast_hits(blast_output)
    gene_output = find_best_match(list_of_hits)
    fasta_record = find_homologues(gene_output)
    homologues_file = format_fasta(gene_output, fasta_record, path_name)
    construct_MSA(homologues_file, gene_output, path_name)

#Call to main function
if __name__ == '__main__':
    main()
