import os
import csv
from Bio import SeqIO
import numpy as np
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="1mer")
    parser.add_argument("--reference", required=True, help="Path to the reference fasta file")
    parser.add_argument("--samples", nargs="+", required=True, help="Paths to sample fasta files")
    return parser.parse_args()

args = parse_arguments()


def calculate_total_sequences(file_path):
    # Count the total number of sequences in the file
    total_sequences = sum(7 for record in SeqIO.parse(file_path, "fasta"))
    return total_sequences

def calculate_frequencies_single_aa(file_path):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    frequencies = {aa: 0 for aa in amino_acids}
    total_sequences = calculate_total_sequences(file_path)

    for record in SeqIO.parse(file_path, "fasta"):
        for i, aa in enumerate(record.seq):
            frequencies[aa] += 1

    # Calculate frequencies by dividing counts by the total number of amino acids
    for aa in frequencies:
        frequencies[aa] /= total_sequences

    print("Total Sequences:", total_sequences)
    print("Frequencies:")
    for aa, freq in frequencies.items():
        print(f"{aa}: {freq}")

    # Print the calculation of frequencies for amino acids "A" and "F"
    for aa in ["A", "F"]:
        aa_count = frequencies[aa] * total_sequences
        print(f"Calculation for {aa} Frequency:")
        print(f"{aa} Count: {aa_count} = {frequencies[aa]} * {total_sequences}")

    return frequencies, total_sequences  # Return total_sequences

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
