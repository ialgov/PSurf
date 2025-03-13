import os
import csv
import tkinter as tk
from tkinter import filedialog
from Bio import SeqIO
import numpy as np
from itertools import product

# def choose_reference_file():
#     root = tk.Tk()
#     root.withdraw()
#     reference_path = filedialog.askopenfilename(title="Choose a reference fasta file")
#     return reference_path
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
    # Count the total number of sequences in the file
    total_sequences = sum(1 for record in SeqIO.parse(file_path, "fasta"))
    return total_sequences

import numpy as np  # Ensure numpy is imported for log2

def calculate_frequencies_dipeptides(file_path):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    dipeptides = [f"{aa1}{aa2}" for aa1, aa2 in product(amino_acids, repeat=2)]
    frequencies = {dipeptide: 0 for dipeptide in dipeptides}
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

        # Count dipeptides, weighted by the sequence count
        sequence = str(record.seq)
        for i in range(len(sequence) - 1):
            dipeptide = sequence[i:i+2]
            frequencies[dipeptide] += count

    # Apply log2 transformation to weighted counts (add pseudocount to avoid log2(0))
    pseudocount = 1  # Add 1 to avoid log2(0)
    log2_counts = {dipeptide: np.log2(count + pseudocount) for dipeptide, count in frequencies.items()}

    # Calculate frequencies using log2-transformed counts
    total_log2_counts = sum(log2_counts.values())
    log2_frequencies = {dipeptide: log2_count / total_log2_counts for dipeptide, log2_count in log2_counts.items()}

    print("Total Weighted Sequences:", total_weighted_sequences)
    print("Log2-Transformed Dipeptide Frequencies:")
    for dipeptide, freq in log2_frequencies.items():
        print(f"{dipeptide}: {freq}")

    return log2_frequencies, total_weighted_sequences

def calculate_stderror(reference_frequencies, total_sequences):
    dipeptides = ["".join(aa) for aa in product("ACDEFGHIKLMNPQRSTVWY", repeat=2)]

    # Calculate stderror for all dipeptides
    stderrors = {dipeptide: np.sqrt(reference_frequencies[dipeptide] * (1 - reference_frequencies[dipeptide]) * (1 / total_sequences)) for dipeptide in dipeptides}

    print("Standard Errors:")
    for dipeptide, stderror in stderrors.items():
        print(f"{dipeptide}: {stderror}")

    return stderrors

def calculate_stderror_dipeptides(reference_path, sample_total_sequences):
    reference_frequencies, total_sequences = calculate_frequencies_dipeptides(reference_path)

    # Print the actual calculation for dipeptides "AA" and "AT"
    for dipeptide in ["AA", "AT"]:
        actual_calculation = f"sqrt({reference_frequencies[dipeptide]} * (1 - {reference_frequencies[dipeptide]}) * (1 / {sample_total_sequences}))"
        print(f"Actual Calculation for {dipeptide}: {actual_calculation}")

    stderrors = calculate_stderror(reference_frequencies, sample_total_sequences)

    return reference_frequencies, stderrors, total_sequences

# ... (rest of the code remains unchanged)

def process_single_sample(reference_path, sample_path):
    sample_total_sequences = calculate_total_sequences(sample_path)
    control_frequencies, stderror, _ = calculate_stderror_dipeptides(reference_path, sample_total_sequences)

    frequencies, _ = calculate_frequencies_dipeptides(sample_path)

    output_file_path = os.path.join(os.path.dirname(sample_path), f"dipeptide-freq-{os.path.basename(sample_path)}.csv")
    with open(output_file_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        header = ["Dipeptide", "Frequency"]
        writer.writerow(header)

        for dipeptide, freq in frequencies.items():
            writer.writerow([dipeptide, freq])

        # Calculate z-scores for dipeptides and create a new CSV file
        z_scores = calculate_z_scores_dipeptides(frequencies, control_frequencies, stderror)
        z_output_file_path = os.path.join(os.path.dirname(sample_path),
                                          f"z-score-dipeptide-{os.path.basename(sample_path)}.csv")
        with open(z_output_file_path, "w", newline="") as z_csvfile:
            z_writer = csv.writer(z_csvfile)
            header = ["Dipeptide", "Z-Score"]
            z_writer.writerow(header)

            for dipeptide, z_score in z_scores.items():
                z_writer.writerow([dipeptide, z_score])

        print(f"Output files created for {os.path.basename(sample_path)}: {output_file_path}, {z_output_file_path}")

        # Create "uniques" file for each sample
        unique_dipeptides = set(dipeptide for dipeptide in frequencies if
                                np.isnan(z_scores[dipeptide]) or z_scores[dipeptide] == float('inf'))
        unique_file_path = os.path.join(os.path.dirname(sample_path), f"uniques_{os.path.basename(sample_path)}.csv")
        with open(unique_file_path, "w", newline="") as unique_csvfile:
            unique_writer = csv.writer(unique_csvfile)
            header = ["Dipeptide", "Frequency"]
            unique_writer.writerow(header)

            for dipeptide in unique_dipeptides:
                unique_writer.writerow([dipeptide, frequencies[dipeptide]])

        print(f"Unique dipeptides file created: {unique_file_path}")

def calculate_z_scores_dipeptides(sample_frequencies, control_frequencies, stderrors):
    z_scores = {dipeptide: 0 for dipeptide in sample_frequencies}

    # Calculate z-scores using the formula: (samplefreq - controlfreq) / stderror
    for dipeptide in sample_frequencies:
        z_scores[dipeptide] = (sample_frequencies[dipeptide] - control_frequencies[dipeptide]) / stderrors[dipeptide]

    # Print the z-score calculation for dipeptides "AA" and "AT"
    for dipeptide in ["AA", "AT"]:
        z_calculation = f"({sample_frequencies[dipeptide]} - {control_frequencies[dipeptide]}) / {stderrors[dipeptide]}"
        print(f"Z-Score Calculation for {dipeptide}: {z_calculation} = {z_scores[dipeptide]}")

    return z_scores

def combine_z_scores():
    current_directory = os.path.dirname(os.path.abspath(__file__))
    samples_folder = os.path.join(current_directory)

    all_z_scores = []
    dipeptides = ["".join(aa) for aa in product("ACDEFGHIKLMNPQRSTVWY", repeat=2)]

    for file_name in os.listdir(samples_folder):
        if file_name.startswith("z-score-dipeptide-") and file_name.endswith(".csv"):
            file_path = os.path.join(samples_folder, file_name)
            with open(file_path, newline="") as csvfile:
                reader = csv.DictReader(csvfile)
                z_scores = {row["Dipeptide"]: float(row["Z-Score"]) for row in reader}

            all_z_scores.append({"Sample": file_name.replace("z-score-dipeptide-", "").replace(".csv", ""), **z_scores})

    output_file_path = os.path.join(current_directory, "combined_z_scores_dipeptide.csv")
    with open(output_file_path, "w", newline="") as csvfile:
        fieldnames = ["Dipeptide"] + [z_scores["Sample"] for z_scores in all_z_scores]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for dipeptide in dipeptides:
            row_data = {"Dipeptide": dipeptide}
            for z_scores in all_z_scores:
                row_data[z_scores["Sample"]] = z_scores.get(dipeptide, 0)

                # Replace 'inf' and 'nan' with 0
                if row_data[z_scores["Sample"]] == float('inf') or np.isnan(row_data[z_scores["Sample"]]):
                    row_data[z_scores["Sample"]] = 0
            writer.writerow(row_data)

    print(f"Combined Z-scores file created: {output_file_path}")

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

    # Step 4: Combine Z-scores for all samples
    combine_z_scores()

if __name__ == "__main__":
    main()
