import argparse
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

#define a dictionary to store base counts for seq sorting
base_counts = {
    'A': 0, 'G': 0, 'C': 0, 'T': 0, 'U': 0, 'N': 0,
    'R': 0, 'Y': 0, 'S': 0, 'W': 0, 'K': 0, 'M': 0,
    'D': 0, 'H': 0, 'B': 0, 'V': 0, 'L': 0, 'P': 0,
    'Q': 0, 'F': 0, 'E': 0, 'I': 0
}

#count instances of each base and store in dictionary
def base_count(sequence):
    counts = {base: sequence.count(base) for base in base_counts}
    return counts

#sort seq into their own category (protein, dna, rna)
def sequence_sort(sequences):
    sorted_sequences = []
    for record in sequences:
        sequence = str(record.seq) 
        counts = base_count(sequence) #perform base counts
        a_percentage = counts['A'] / sum(counts.values())
        if a_percentage > 0.10 and counts['U'] == 0: #using percentage of a and instances of U to differentiate btwn dna and mrna and dna and protein
            category = "DNA"
        elif counts['U'] > 0.1: #mRNA has lots of U, more than protein
            category = "mRNA"
        else:
            category = "Protein" #everything else goes into the protein category
        sorted_sequences.append((record.id, category, sequence))
    return sorted_sequences

#search for sequences on blast
def perform_blast(sorted_sequences):
    blast_results = [] #list of blast results
    for header, category, sequence in sorted_sequences:
        print(f"Performing BLAST for sequence {header} in category {category}...") #error check
        if category == "DNA":
            result_handle = NCBIWWW.qblast("blastn", "nt", sequence) #search in nucelotide sequences
        elif category == "mRNA":
            result_handle = NCBIWWW.qblast("blastn", "nr", sequence)
        else:
            result_handle = NCBIWWW.qblast("blastp", "nr", sequence)

        blast_records = NCBIXML.parse(result_handle)
        blast_record = next(blast_records, None)
        if blast_record and blast_record.alignments: #if there is a record for the sequence, find the best alignment
            best_alignment = blast_record.alignments[0]
            blast_results.append((header, best_alignment.title, category, best_alignment.length, best_alignment.hsps[0].expect))
        else:
            blast_results.append((header, "No hits found", category, 0, "N/A")) #otherwise, if no hits found, print that in output file
    return blast_results
#main
if __name__ == "__main__":
    #command-line arguments with argparse
    parser = argparse.ArgumentParser(description='Sort sequences into DNA, mRNA, and protein categories and perform BLAST search.')
    parser.add_argument('input_file', help='Input FASTA file') #input file of sequences
    parser.add_argument('output_file', help='Output FASTA file') #output file
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file

    print("Starting sequence sorting...") #print for error checking

    #read FASTA file using SeqIO
    sequences = list(SeqIO.parse(input_file, "fasta"))

    #sort seqs
    sorted_sequences = sequence_sort(sequences)

    #blast function
    blast_results = perform_blast(sorted_sequences)

    print("BLAST search completed.") #print for error checking

    #write BLAST results to file
    with open(output_file, 'w') as file:
        for header, hit_title, category, hit_length, hit_evalue in blast_results:
            file.write(f">{header} ({category})\nBest Hit: {hit_title}\nLength: {hit_length}\nE-value: {hit_evalue}\n\n")

    print(f"BLAST results written to {output_file}.") #so i know it's done
