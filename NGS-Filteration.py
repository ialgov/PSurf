# This script allows choosing a specific NGS results raw data file (.fastq) and filtration of the results in several steps.
# The output files would be saved in the same folder as the input fastq file.
#Each sample would have a folder containing sequences from the different filtration steps and a folder named "Final" would contain the filtered data that can be used for further analyses.
# Make sure to install all required packages before running this script.

import tkinter as tk
from tkinter import filedialog
import os
from Bio import SeqIO

# Define barcode sequences and sample names.
# For barcodes, it is recommended to include the sequence downstream to the barcode for better specificity.
barcodes = {
    "Sample 1": "ATCACGACCGCCTCCACTAGCATATG",
    "Sample 2": "ACAGTGGTCGCCTCCACTAGCATATG",
    "Sample 3": "CAGATCCACGCCTCCACTAGCATATG",
    "Sample 4": "ACAAACGGCGCCTCCACTAGCATATG",
    "Sample 5": "ACCCAGCACGCCTCCACTAGCATATG",
    "Sample 6": "AACCCCTCCGCCTCCACTAGCATATG",
    "Sample 7": "CCCAACCTCGCCTCCACTAGCATATG",
    "Sample 8": "CACCACACCGCCTCCACTAGCATATG",
    "Sample 9": "GAAACCCACGCCTCCACTAGCATATG",
    "Sample 10": "TGTGACCACGCCTCCACTAGCATATG",
    "Sample 11": "AGGGTCAACGCCTCCACTAGCATATG",
    "Sample 12": "AGGAGTGGCGCCTCCACTAGCATATG"
    
}

# Create a file dialog box to choose the input file
root = tk.Tk()
root.withdraw()
input_file = filedialog.askopenfilename()

# Get the directory and filename from the chosen input file
input_dir, input_filename = os.path.split(input_file)

# Create the output directory path
output_dir = os.path.join(input_dir, "DeMultiplexed_" + os.path.splitext(input_filename)[0])

# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# Initialize a dictionary to store output file handles
output_files = {bc: open(os.path.join(output_dir, bc + ".fastq"), "w") for bc in barcodes}

# Loop through input file and write to appropriate output file based on barcode
for record in SeqIO.parse(input_file, "fastq"):
    seq = str(record.seq)
    for bc, bc_seq in barcodes.items():
        if seq.startswith(bc_seq):
            SeqIO.write(record, output_files[bc], "fastq")
            break

# Close all output file handles
for f in output_files.values():
    f.close()

from Bio import SeqIO, Seq
from collections import defaultdict
import re
import os
import tkinter as tk
from tkinter import filedialog

# Prompt the user to choose a folder
root = tk.Tk()
root.withdraw()
folder_path = output_dir

# Find all the fastq files in the selected folder
fastq_files = [f for f in os.listdir(folder_path) if f.endswith('.fastq')]

# Check if there are any fastq files in the folder
if not fastq_files:
    print("No fastq files found in the selected folder.")
    exit(1)

# Process each fastq file in the folder
for fastq_file in fastq_files:
    # Get the base name of the fastq file
    base_name, ext = os.path.splitext(fastq_file)

    # Create a new directory with the same name as the fastq file
    output_dir = os.path.join(folder_path, base_name)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # Create a new directory for the final files
    output_dir_final = os.path.join(folder_path, 'Final')
    if not os.path.exists(output_dir_final):
        os.mkdir(output_dir_final)

    # Define input and output file paths - rename the files as desired.
    input_file_1 = os.path.join(folder_path, fastq_file)
    output_file_1 = os.path.join(output_dir, base_name + '-Filtered.fastq')
    input_file_2 = os.path.join(output_dir, base_name + '-Filtered.fastq')
    output_file_2 = os.path.join(output_dir, base_name + '-21.fastq')
    input_file_3 = os.path.join(output_dir, base_name + '-21.fastq')
    output_file_3 = os.path.join(output_dir, base_name + '-collapsed.fasta')
    input_file_4 = os.path.join(output_dir, base_name + '-collapsed.fasta')
    output_file_4 = os.path.join(output_dir, base_name + '-noStop.fasta')
    input_file_5 = os.path.join(output_dir, base_name + '-noStop.fasta')
    output_file_5 = os.path.join(output_dir, base_name + '-noTEV.fasta')
    input_file_6 = os.path.join(output_dir, base_name + '-noTEV.fasta')
    output_file_6 = os.path.join(output_dir, base_name + '-noMotif2.fasta')
    input_file_7 = os.path.join(output_dir, base_name + '-noMotif2.fasta')
    output_file_7 = os.path.join(output_dir, base_name + '-noMotif2-multiply.fasta')
    input_file_8 = os.path.join(output_dir, base_name + '-noMotif2-multiply.fasta')
    output_file_8 = os.path.join(output_dir_final, base_name + '.fasta')

    from Bio import SeqIO, Seq
    from collections import defaultdict
    import re
    import os
    import tkinter as tk
    from tkinter import filedialog

    # Set the motif to search for.
    # motif_1 would be a sequence upstream to the substrate region, to include only sequences containing that sequence of interest.
    # motif_2 can be defined for removal of specific unwanted motifs or sequences from the data.
    motif_1 = "GTAGCTAGCCCCACTACCGCC"
    motif_2 = re.compile(r"GGGGGGG")


    # Function to reverse complement a sequence - should be used in case the data used is from the reverse sequencing.
    def reverse_complement(seq):
        return str(Seq.Seq(seq).reverse_complement())


    # Step 1: Filter reads by motif - this step make sure to include in the data only sequences that contains a specific sequence
    num_reads_processed = 0
    num_reads_with_motif_1 = 0
    with open(output_file_1, "w") as out_handle:
        for record in SeqIO.parse(input_file_1, "fastq"):
            num_reads_processed += 1
            if motif_1 in str(record.seq):
                SeqIO.write(record, out_handle, "fastq")
                num_reads_with_motif_1 += 1
    print(f"Processed {num_reads_processed} reads, found {num_reads_with_motif_1} reads with the motif_1")

    # Step 2: Extract 21 bp sequences after motif and create new reads - this step extracts 21bp (which corresponds to 7 amino acids substrate used) to get a list of dna sequences of the substrates only.
    input_reads = SeqIO.parse(output_file_1, "fastq")
    output_reads = []
    for read in input_reads:
        if motif_1 in read.seq:
            start = read.seq.find(motif_1) + len(motif_1)
            end = start + 21
            seq = str(read.seq[start:end])
            new_read = read[start:end]
            new_read.id = read.id + "_21bp"
            new_read.description = ""
            output_reads.append(new_read)
    SeqIO.write(output_reads, output_file_2, "fastq")

    # Step 3: Reverse complement the sequences and collapse identical sequences - in this step the 21bp are processed to result in a fasta file with ranked dna sequences based on their count in the data
    seq_counts = defaultdict(int)
    with open(output_file_2, "r") as in_file:
        for record in SeqIO.parse(in_file, "fastq"):
            rev_comp_seq = reverse_complement(str(record.seq))
            seq_counts[rev_comp_seq] += 1
        sequences = [(seq, count) for seq, count in seq_counts.items()]
    sequences.sort(key=lambda x: x[1], reverse=True)
    with open(output_file_3, "w") as out_file:
        for i, (seq, count) in enumerate(sequences, 1):
            SeqIO.write(SeqIO.SeqRecord(Seq.Seq(seq), id=f"{i}-{count}", description=""), out_file, "fasta")

    # Step 4: Translate DNA sequences to protein sequences and remove stop codons.
    with open(output_file_3, "r") as in_file, open(output_file_4, "w") as out_file:
        for record in SeqIO.parse(in_file, "fasta"):
            protein_seq = record.seq.translate()
            if "*" not in protein_seq:
                SeqIO.write(SeqIO.SeqRecord(protein_seq, id=record.id, description=""), out_file, "fasta")

    # Step 5: Remove TEV substrates - this step removes any sequences that shares >70% similarity the substrate used as the template for library construction.

    import re
    from Bio.Seq import Seq
    from Bio.Align import PairwiseAligner

    # Open input and output files
    with open(output_file_4, "r") as input_file, open(output_file_5, "w") as output_file:
        # Initialize variables
        sequence = ""
        header = ""
        aligner = PairwiseAligner()
        # Loop through input file
        for line in input_file:
            line = line.strip()
            # If line is a header line
            if line.startswith(">"):
                # Write the previous sequence to output file if it doesn't contain the motif
                if not any(aligner.score(Seq("ENLYFQM"), Seq(sequence[start:end])) / 7 >= 0.7
                           for start in range(len(sequence) - 6) for end in range(start + 7, len(sequence) + 1)):
                    output_file.write(header + "\n" + sequence + "\n")
                # Update header and reset sequence
                header = line
                sequence = ""
            # If line is a sequence line
            else:
                sequence += line
        # Write the last sequence to output file if it doesn't contain the motif
        if not any(aligner.score(Seq("ENLYFQM"), Seq(sequence[start:end])) / 7 >= 0.7
                   for start in range(len(sequence) - 6) for end in range(start + 7, len(sequence) + 1)):
            output_file.write(header + "\n" + sequence + "\n")

    # Step 6: Remove sequences with unwanted motif
    # Define the motif_2 to search for
    motif_2 = re.compile(r"GGGGGGG") #GGGGGGG can be replaced with an unwanted motif such as LxxR in yeast systems

    # Define the name of the input fasta file
    input_file = input_file_6

    # Define the name of the output fasta file
    output_file = output_file_6

    # Open the output file for writing
    with open(output_file, "w") as f_out:

        # Iterate over the sequences in the input fasta file
        for seq_record in SeqIO.parse(input_file, "fasta"):

            # Flag to keep track of whether the motif_2_2 is present
            motif_2_present = False

            # Iterate over all possible 4 amino acid peptide sequences in the protein sequence
            for i in range(len(seq_record.seq) - 3):
                peptide = str(seq_record.seq[i:i + 4])

                # Check if the motif_2 is present in the peptide
                if motif_2.search(peptide):
                    motif_2_present = True
                    break

            # If the motif_2 is present, skip this sequence and move on to the next one
            if motif_2_present:
                continue

            # If the motif_2 is not present, write the sequence to the output file
            SeqIO.write(seq_record, f_out, "fasta")

    # Step 7: Multiply (deCollapse) - this step is done to find similar substrate sequences that were counted seperately as  a result of synonymous mutations
    # Define the name of the input fasta file
    input_file = input_file_7
    # Define the name of the output fasta file
    output_file = output_file_7

    with open(input_file, "r") as f:
        data = f.readlines()

    sequences = {}
    for line in data:
        if line.startswith(">"):
            seq_id = line.strip()[1:]
            rank = int(seq_id.split("-")[1].split("_")[0])
            sequences[seq_id] = {"rank": rank, "seq": ""}
        else:
            sequences[seq_id]["seq"] += line.strip()

    output_sequences = {}
    for seq_id in sequences:
        rank = sequences[seq_id]["rank"]
        seq = sequences[seq_id]["seq"]
        for i in range(1, rank + 1):
            new_seq_id = f"{seq_id.split('-')[0]}-s_{i}"
            output_sequences[new_seq_id] = seq

    with open(output_file, "w") as f:
        for seq_id in output_sequences:
            seq = output_sequences[seq_id]
            f.write(f">{seq_id}\n{seq}\n")

    # CollapseAA sequences - in this step all similar amino acids substrates are being collapsed and ranked according to their abundance in the data.
    from collections import Counter
    import operator

    input_file = input_file_8
    output_file = output_file_8

    # Read input file and count sequences
    seq_counts = Counter()
    with open(input_file, "r") as infile:
        seq = ""
        for line in infile:
            if line.startswith(">"):
                if seq:
                    seq_counts[seq] += 1
                    seq = ""
            else:
                seq += line.strip()
        if seq:
            seq_counts[seq] += 1

    # Sort sequences by count, most abundant first
    sorted_seqs = sorted(seq_counts.items(), key=operator.itemgetter(1), reverse=True)

    # Write sorted sequences to output file, excluding sequences that appear only once
    with open(output_file, "w") as outfile:
        for i, (seq, count) in enumerate(sorted_seqs):
            if count > 1:  # Skip sequences that appear only once
                outfile.write(">seq{0}_{1}\n{2}\n".format(i + 1, count, seq))

    print("Collapsing and filtering sequences complete.")


print("File organization complete.")

