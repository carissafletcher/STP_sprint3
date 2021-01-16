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
from Bio import AlignIO, Entrez, Phylo
from Bio.Blast import NCBIWWW
from Bio.Align.Applications import ClustalwCommandline
import xml.etree.ElementTree as ET

#Function 1: Get query sequence from string or .fasta file
def acquire_input():
    input_type = input('Sequence format is (1) string (2) fasta file? ')
    while not ((input_type == '1') or (input_type == '2')):
        print('Please enter either 1 or 2.')
        input_type = input('Sequence format is (1) string (2) fasta file? ')

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
        while not (filename[-6:] == '.fasta'):
            print('File must be .fasta type.')
            filename = input('Enter filename: ')
        query_sequence = open(filename).read()
        print(query_sequence)

    query_name = input('Please enter a name for this query: ')
    while not query_name:
        query_name = input('Please enter a name for this query: ')
    return query_sequence, query_name

#Function 2: Perform blastn search on query sequence, get XML output
def blast_search(query_sequence, query_name):
    print('Querying BLAST - this may take some time...')
    result_handle = NCBIWWW.qblast("blastn", "nt", query_sequence)
    print('BLAST query complete.')
    blast_results = result_handle.read()
    result_handle.close()
    output_file = 'results_' + query_name + '.xml'
    with open(output_file, 'w') as save_file:
        save_file.write(blast_results)
    return output_file

#Function 3: Extract relevant info from BLAST output
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

#Function 4: Find sequence match with highest % identity and get gene symbol
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

#Function 5: Identify gene homologues and return protein sequences as .fasta
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
    first_seq = int(fasta_record.find('>'))
    clipped_fasta_record = (fasta_record[first_seq:]).strip()
    fasta_file = gene_output + '_hg.fasta'
    with open(fasta_file, 'w') as output_object:
        output_object.write(clipped_fasta_record)

    return fasta_file

#Function 7: Construct an MSA
#def construct_MSA(fasta_transcripts):
    #cline = ClustalwCommandline("clustalw2", infile = fasta_transcripts) #clustalW needs local installation or api...


example_sequence = "GCTGTTCAGCGTTCTGCTGGAGCAGGGCCCCGGACGGCCAGGCGACGCCCCGCACACCGG" #from CACNA1F
test_string = "  gCTGTTCAGCGTTCTGCtggAGCa GGGCEFCGGACGGCCAGGCGAC  GCCCCICACACCgg " #some gaps, lowercase, and unsupported characters

#Main function to call other functions
def main():
    sequence, name = acquire_input()
    blast_output = blast_search(sequence, name)
    output_list = get_blast_hits(blast_output)
    gene_symbol = find_best_match(output_list)
    find_homologues(gene_symbol)

#Call to main function
if __name__ == '__main__':
    main()
