from Bio.Blast import NCBIWWW

sequence_data = "GCTGTTCAGCGTTCTGCTGGAGCAGGGCCCCGGACGGCCAGGCGACGCCCCGCACACCGG"

def blast_search(sequence_data):
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data)
    print(result_handle)

    with open('results.xml', 'w') as save_file: 
        blast_results = result_handle.read() 
        save_file.write(blast_results)
    return blast_results

blast_search(sequence_data)