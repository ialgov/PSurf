import os
import csv
import tkinter as tk
from tkinter import filedialog
from Bio import SeqIO
import numpy as np

SEQUENCE_LENGTH = 7  # Specify the fixed sequence length


# def choose_reference_file():
#     root = tk.Tk()
#     root.withdraw()
#     reference_path = filedialog.askopenfilename(title="Choose a reference fasta file")
#     return reference_path
#
#
# def choose_sample_files():
#     root = tk.Tk()
#     root.withdraw()
#     sample_paths = filedialog.askopenfilenames(title="Choose sample fasta files")
#     return sample_paths
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Your script description")
    parser.add_argument("--reference", required=True, help="Path to the reference fasta file")
    parser.add_argument("--samples", nargs="+", required=True, help="Paths to sample fasta files")
    return parser.parse_args()

args = parse_arguments()

def calculate_total_sequences(file_path):
    total_sequences = sum(1 for record in SeqIO.parse(file_path, "fasta"))
    return total_sequences


import numpy as np  # Ensure numpy is imported for log2

def calculate_frequencies_per_position(file_path):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    frequencies = {aa: [0] * SEQUENCE_LENGTH for aa in amino_acids}
    total_weighted_sequences = 0  # Total weighted sequences

    for record in SeqIO.parse(file_path, "fasta"):
        # Extract count from sequence name (e.g., '>seq1_1079' -> 1079)
        count = 1  # Default count if not specified
        if "_" in record.id:
            try:
                count = int(record.id.split("_")[-1])  # Extract the last part as count
            except ValueError:
                pass  # If the last part is not a number, use count = 1

        # Skip sequences with count < 2
        if count < 2:
            continue

        total_weighted_sequences += count  # Add to total weighted sequences

        # Count amino acids at each position, weighted by the sequence count
        for i, aa in enumerate(record.seq):
            frequencies[aa][i] += count

    # Apply log2 transformation to weighted counts (add pseudocount to avoid log2(0))
    pseudocount = 1  # Add 1 to avoid log2(0)
    log2_counts = {aa: [np.log2(count + pseudocount) for count in frequencies[aa]] for aa in amino_acids}

    # Calculate frequencies using log2-transformed counts
    total_log2_counts = sum(sum(log2_counts[aa]) for aa in amino_acids)
    log2_frequencies = {aa: [log2_count / total_log2_counts for log2_count in log2_counts[aa]] for aa in amino_acids}

    print("Total Weighted Sequences:", total_weighted_sequences)
    print("Log2-Transformed Frequencies per Position:")
    for aa in log2_frequencies:
        print(f"{aa}: {log2_frequencies[aa]}")

    return log2_frequencies, total_weighted_sequences

def calculate_stderror_per_position(reference_path, sample_total_sequences):
    reference_frequencies, total_sequences = calculate_frequencies_per_position(reference_path)

    stderrors = {
        aa: [np.sqrt(reference_frequencies[aa][i] * (1 - reference_frequencies[aa][i]) * (1 / sample_total_sequences))
             for i in range(SEQUENCE_LENGTH)] for aa in reference_frequencies}

    print("Standard Errors:")
    for aa, stderror_list in stderrors.items():
        print(f"{aa}: {stderror_list}")

    return reference_frequencies, stderrors, total_sequences


def process_single_sample(reference_path, sample_path):
    sample_total_sequences = calculate_total_sequences(sample_path)
    control_frequencies, stderror, _ = calculate_stderror_per_position(reference_path, sample_total_sequences)

    frequencies, _ = calculate_frequencies_per_position(sample_path)

    output_file_path = os.path.join(os.path.dirname(sample_path), f"aa-freq-{os.path.basename(sample_path)}.csv")
    with open(output_file_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        header = ["Amino Acid", "Position", "Frequency"]
        writer.writerow(header)

        for aa in frequencies:
            for i, freq in enumerate(frequencies[aa]):
                writer.writerow([aa, i + 1, freq])

    # Calculate z-scores for amino acids and create a new CSV file
    z_scores = calculate_z_scores_positions(frequencies, control_frequencies, stderror)
    z_output_file_path = os.path.join(os.path.dirname(sample_path), f"z-score-aa-{os.path.basename(sample_path)}.csv")
    with open(z_output_file_path, "w", newline="") as z_csvfile:
        z_writer = csv.writer(z_csvfile)
        header = ["Amino Acid"] + [f"P{i}" for i in range(1, SEQUENCE_LENGTH + 1)]
        z_writer.writerow(header)

        for aa in frequencies:
            row_data = [aa]
            for i in range(SEQUENCE_LENGTH):
                row_data.append(z_scores[aa][i])
            z_writer.writerow(row_data)

    print(f"Output files created for {os.path.basename(sample_path)}: {output_file_path}, {z_output_file_path}")

def calculate_z_scores_positions(sample_frequencies, control_frequencies, stderrors):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    z_scores = {aa: [0] * SEQUENCE_LENGTH for aa in amino_acids}

    # Calculate z-scores using the formula: (samplefreq - controlfreq) / stderror
    for aa in sample_frequencies:
        for i in range(SEQUENCE_LENGTH):
            z_scores[aa][i] = (sample_frequencies[aa][i] - control_frequencies[aa][i]) / stderrors[aa][i]

    return z_scores

def main():
    # # Step 1: Choose a reference fasta file
    # reference_path = choose_reference_file()
    #
    # # Step 2: Choose multiple sample fasta files
    # sample_paths = choose_sample_files()

    # Get input file paths from command-line arguments
    reference_path = args.reference
    sample_paths = args.samples

    # Step 3: Process each sample fasta file
    for sample_path in sample_paths:
        process_single_sample(reference_path, sample_path)



if __name__ == "__main__":
    main()
