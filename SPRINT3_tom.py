from Bio.Blast import NCBIWWW, NCBIXML
import xml.etree.ElementTree as ET

query_name = input("Please enter query name: ")

sequence_data = "ACATGGAAGCGGAGGAGGACCCATCGCCCTGTGCGTCTTCGGGATCAGGGCGTGCGTGCA"

def blast_search(sequence_data):
    print("Sending request to BLAST, this may take several minutes during peak times!")
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data)
    print(result_handle)

    with open('results_' + query_name + '.xml', 'w') as save_file: 
        blast_results = result_handle.read() 
        save_file.write(blast_results)
    return blast_results

def parse_blast_xml(blast_results):
    #tree = ET.parse('results.xml')
    root = blast_results.getroot()



    for hit_def in root.findall('hit_def'):
        
        print(hit_def.findtext())
    




output = blast_search(sequence_data)

parse_blast_xml(output)

