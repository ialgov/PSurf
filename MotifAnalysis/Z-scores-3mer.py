import os
import csv
from Bio import SeqIO
import numpy as np
from itertools import product
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="3mer")
    parser.add_argument("--reference", required=True, help="Path to the reference fasta file")
    parser.add_argument("--samples", nargs="+", required=True, help="Paths to sample fasta files")
    return parser.parse_args()

args = parse_arguments()

def calculate_total_sequences(file_path):
    # Count the total number of sequences in the file
    total_sequences = sum(1 for record in SeqIO.parse(file_path, "fasta"))
    return total_sequences

def calculate_frequencies_triipeptides(file_path):
    triipeptides = ["".join(aa) for aa in product("ACDEFGHIKLMNPQRSTVWY", repeat=3)]
    frequencies = {triipeptide: 0 for triipeptide in triipeptides}
    total_sequences = calculate_total_sequences(file_path)

    for record in SeqIO.parse(file_path, "fasta"):
        sequence = record.seq
        for i in range(len(sequence) - 2):
            triipeptide = sequence[i:i+3]
            frequencies[triipeptide] += 1

    # Calculate frequencies by dividing counts by the total number of triipeptides
    for triipeptide in frequencies:
        frequencies[triipeptide] /= total_sequences

    print("Total Sequences:", total_sequences)
    print("triipeptide Frequencies:")
    for triipeptide, freq in frequencies.items():
        print(f"{triipeptide}: {freq}")

    return frequencies, total_sequences

def calculate_stderror(reference_frequencies, total_sequences):
    triipeptides = ["".join(aa) for aa in product("ACDEFGHIKLMNPQRSTVWY", repeat=3)]

    # Calculate stderror for all triipeptides
    stderrors = {triipeptide: np.sqrt(reference_frequencies[triipeptide] * (1 - reference_frequencies[triipeptide]) * (1 / total_sequences)) for triipeptide in triipeptides}

    print("Standard Errors:")
    for triipeptide, stderror in stderrors.items():
        print(f"{triipeptide}: {stderror}")

    return stderrors

def calculate_stderror_triipeptides(reference_path, sample_total_sequences):
    reference_frequencies, total_sequences = calculate_frequencies_triipeptides(reference_path)

    # Print the actual calculation for triipeptides "AAA" and "AAT"
    for triipeptide in ["AAA", "AAT"]:
        actual_calculation = f"sqrt({reference_frequencies[triipeptide]} * (1 - {reference_frequencies[triipeptide]}) * (1 / {sample_total_sequences}))"
        print(f"Actual Calculation for {triipeptide}: {actual_calculation}")

    stderrors = calculate_stderror(reference_frequencies, sample_total_sequences)

    return reference_frequencies, stderrors, total_sequences

def process_single_sample(reference_path, sample_path):
    sample_total_sequences = calculate_total_sequences(sample_path)
    control_frequencies, stderror, _ = calculate_stderror_triipeptides(reference_path, sample_total_sequences)

    frequencies, _ = calculate_frequencies_triipeptides(sample_path)

    output_file_path = os.path.join(os.path.dirname(sample_path), f"triipeptide-freq-{os.path.basename(sample_path)}.csv")
    with open(output_file_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        header = ["triipeptide", "Frequency"]
        writer.writerow(header)

        for triipeptide, freq in frequencies.items():
            writer.writerow([triipeptide, freq])

    # Calculate z-scores for triipeptides and create a new CSV file
    z_scores = calculate_z_scores_triipeptides(frequencies, control_frequencies, stderror)
    z_output_file_path = os.path.join(os.path.dirname(sample_path), f"z-score-triipeptide-{os.path.basename(sample_path)}.csv")
    with open(z_output_file_path, "w", newline="") as z_csvfile:
        z_writer = csv.writer(z_csvfile)
        header = ["triipeptide", "Z-Score"]
        z_writer.writerow(header)

        for triipeptide, z_score in z_scores.items():
            z_writer.writerow([triipeptide, z_score])

    print(f"Output files created for {os.path.basename(sample_path)}: {output_file_path}, {z_output_file_path}")

    # Create "uniques" file for each sample
    unique_triipeptides = set(triipeptide for triipeptide in frequencies if
                              np.isnan(z_scores[triipeptide]) or z_scores[triipeptide] == float('inf'))
    unique_file_path = os.path.join(os.path.dirname(sample_path), f"uniques_{os.path.basename(sample_path)}.csv")
    with open(unique_file_path, "w", newline="") as unique_csvfile:
        unique_writer = csv.writer(unique_csvfile)
        header = ["triipeptide", "Frequency"]
        unique_writer.writerow(header)

        for triipeptide in unique_triipeptides:
            unique_writer.writerow([triipeptide, frequencies[triipeptide]])

    print(f"Unique triipeptides file created: {unique_file_path}")

def calculate_z_scores_triipeptides(sample_frequencies, control_frequencies, stderrors):
    z_scores = {triipeptide: 0 for triipeptide in sample_frequencies}

    # Calculate z-scores using the formula: (samplefreq - controlfreq) / stderror
    for triipeptide in sample_frequencies:
        z_scores[triipeptide] = (sample_frequencies[triipeptide] - control_frequencies[triipeptide]) / stderrors[triipeptide]

    # Print the z-score calculation for triipeptides "AA" and "AT"
    for triipeptide in ["AAA", "AAT"]:
        z_calculation = f"({sample_frequencies[triipeptide]} - {control_frequencies[triipeptide]}) / {stderrors[triipeptide]}"
        print(f"Z-Score Calculation for {triipeptide}: {z_calculation} = {z_scores[triipeptide]}")

    return z_scores

def combine_z_scores():
    current_directory = os.path.dirname(os.path.abspath(__file__))
    samples_folder = os.path.join(current_directory)

    all_z_scores = []
    triipeptides = ["".join(aa) for aa in product("ACDEFGHIKLMNPQRSTVWY", repeat=3)]

    for file_name in os.listdir(samples_folder):
        if file_name.startswith("z-score-triipeptide-") and file_name.endswith(".csv"):
            file_path = os.path.join(samples_folder, file_name)
            with open(file_path, newline="") as csvfile:
                reader = csv.DictReader(csvfile)
                z_scores = {row["triipeptide"]: float(row["Z-Score"]) for row in reader}

            all_z_scores.append({"Sample": file_name.replace("z-score-triipeptide-", "").replace(".csv", ""), **z_scores})

    output_file_path = os.path.join(current_directory, "combined_z_scores_triipeptide.csv")
    with open(output_file_path, "w", newline="") as csvfile:
        fieldnames = ["triipeptide"] + [z_scores["Sample"] for z_scores in all_z_scores]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for triipeptide in triipeptides:
            row_data = {"triipeptide": triipeptide}
            for z_scores in all_z_scores:
                row_data[z_scores["Sample"]] = z_scores.get(triipeptide, 0)

                # Replace 'inf' and 'nan' with 0
                if row_data[z_scores["Sample"]] == float('inf') or np.isnan(row_data[z_scores["Sample"]]):
                    row_data[z_scores["Sample"]] = 0
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
