import os
import csv
import tkinter as tk
from tkinter import filedialog
from Bio import SeqIO
import numpy as np
from itertools import product

import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Your script description")
    parser.add_argument("--reference", required=True, help="Path to the reference fasta file")
    parser.add_argument("--samples", nargs="+", required=True, help="Paths to sample fasta files")
    return parser.parse_args()

args = parse_arguments()


def calculate_total_sequences(file_path):
    # Count the total number of sequences in the file
    total_sequences = sum(7 for record in SeqIO.parse(file_path, "fasta"))
    return total_sequences

import numpy as np  # Ensure numpy is imported for log2

def calculate_frequencies_single_aa(file_path):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    frequencies = {aa: 0 for aa in amino_acids}
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

        # Count amino acids, weighted by the sequence count
        for aa in record.seq:
            frequencies[aa] += count

    # Apply log2 transformation to weighted counts (add pseudocount to avoid log2(0))
    pseudocount = 1  # Add 1 to avoid log2(0)
    log2_counts = {aa: np.log2(count + pseudocount) for aa, count in frequencies.items()}

    # Calculate frequencies using log2-transformed counts
    total_log2_counts = sum(log2_counts.values())
    log2_frequencies = {aa: log2_count / total_log2_counts for aa, log2_count in log2_counts.items()}

    print("Total Weighted Sequences:", total_weighted_sequences)
    print("Log2-Transformed Frequencies:")
    for aa, freq in log2_frequencies.items():
        print(f"{aa}: {freq}")

    return log2_frequencies, total_weighted_sequences

def calculate_stderror(reference_frequencies, total_sequences):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    # Calculate stderror for all amino acids
    stderrors = {aa: np.sqrt(reference_frequencies[aa] * (1 - reference_frequencies[aa]) * (1 / total_sequences)) for aa in amino_acids}

    print("Standard Errors:")
    for aa, stderror in stderrors.items():
        print(f"{aa}: {stderror}")

    return stderrors

def calculate_stderror_single_aa(reference_path, sample_total_sequences):
    reference_frequencies, total_sequences = calculate_frequencies_single_aa(reference_path)

    # Print the actual calculation for amino acids "A" and "T"
    for aa in ["A", "F"]:
        actual_calculation = f"sqrt({reference_frequencies[aa]} * (1 - {reference_frequencies[aa]}) * (1 / {sample_total_sequences}))"
        print(f"Actual Calculation for {aa}: {actual_calculation}")

    stderrors = calculate_stderror(reference_frequencies, sample_total_sequences)

    return reference_frequencies, stderrors, total_sequences

# ... (rest of the code remains unchanged)

def process_single_sample(reference_path, sample_path):
    sample_total_sequences = calculate_total_sequences(sample_path)
    control_frequencies, stderror, _ = calculate_stderror_single_aa(reference_path, sample_total_sequences)

    frequencies, _ = calculate_frequencies_single_aa(sample_path)

    output_file_path = os.path.join(os.path.dirname(sample_path), f"aa-freq-{os.path.basename(sample_path)}.csv")
    with open(output_file_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        header = ["Amino Acid", "Frequency"]
        writer.writerow(header)

        for aa, freq in frequencies.items():
            writer.writerow([aa, freq])

    # Calculate z-scores for amino acids and create a new CSV file
    z_scores = calculate_z_scores_single_aa(frequencies, control_frequencies, stderror)
    z_output_file_path = os.path.join(os.path.dirname(sample_path), f"z-score-aa-{os.path.basename(sample_path)}.csv")
    with open(z_output_file_path, "w", newline="") as z_csvfile:
        z_writer = csv.writer(z_csvfile)
        header = ["Amino Acid", "Z-Score"]
        z_writer.writerow(header)

        for aa, z_score in z_scores.items():
            z_writer.writerow([aa, z_score])

    print(f"Output files created for {os.path.basename(sample_path)}: {output_file_path}, {z_output_file_path}")

def calculate_z_scores_single_aa(sample_frequencies, control_frequencies, stderrors):
    z_scores = {aa: 0 for aa in sample_frequencies}

    # Calculate z-scores using the formula: (samplefreq - controlfreq) / stderror
    for aa in sample_frequencies:
        z_scores[aa] = (sample_frequencies[aa] - control_frequencies[aa]) / stderrors[aa]

    # Print the z-score calculation for amino acids "A" and "F"
    for aa in ["A", "F"]:
        z_calculation = f"({sample_frequencies[aa]} - {control_frequencies[aa]}) / {stderrors[aa]}"
        print(f"Z-Score Calculation for {aa}: {z_calculation} = {z_scores[aa]}")

    return z_scores

def combine_z_scores():
    current_directory = os.path.dirname(os.path.abspath(__file__))
    samples_folder = os.path.join(current_directory)

    all_z_scores = []
    motifs = list("ACDEFGHIKLMNPQRSTVWY")

    for file_name in os.listdir(samples_folder):
        if file_name.startswith("z-score-aa-") and file_name.endswith(".csv"):
            file_path = os.path.join(samples_folder, file_name)
            with open(file_path, newline="") as csvfile:
                reader = csv.DictReader(csvfile)
                z_scores = {row["Amino Acid"]: float(row["Z-Score"]) for row in reader}

            all_z_scores.append({"Sample": file_name.replace("z-score-aa-", "").replace(".csv", ""), **z_scores})

    output_file_path = os.path.join(current_directory, "combined_z_scores_aa.csv")
    with open(output_file_path, "w", newline="") as csvfile:
        fieldnames = ["Amino Acid"] + [z_scores["Sample"] for z_scores in all_z_scores]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for motif in motifs:
            row_data = {"Amino Acid": motif}
            for z_scores in all_z_scores:
                row_data[z_scores["Sample"]] = z_scores.get(motif, 0)
            writer.writerow(row_data)

    print(f"Combined Z-scores file created: {output_file_path}")

def main():
    # Get input file paths from command-line arguments
    reference_path = args.reference
    sample_paths = args.samples

    # Step 3: Process each sample fasta file
    for sample_path in sample_paths:
        process_single_sample(reference_path, sample_path)

    # Step 4: Combine Z-scores for all samples
    combine_z_scores()

if __name__ == "__main__":
    main()
