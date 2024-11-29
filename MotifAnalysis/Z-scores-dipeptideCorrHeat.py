import os
import csv
from Bio import SeqIO
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="dipeptideCorrHeat")
    parser.add_argument("--reference", required=True, help="Path to the reference fasta file")
    parser.add_argument("--samples", nargs="+", required=True, help="Paths to sample fasta files")
    return parser.parse_args()

args = parse_arguments()

def calculate_total_sequences(file_path):
    # Count the total number of sequences in the file
    total_sequences = sum(1 for record in SeqIO.parse(file_path, "fasta"))
    return total_sequences

def calculate_frequencies_dipeptides(file_path):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    dipeptides = [f"{aa1}{aa2}" for aa1, aa2 in product(amino_acids, repeat=2)]
    frequencies = {dipeptide: 0 for dipeptide in dipeptides}
    total_sequences = calculate_total_sequences(file_path)

    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq)
        for i in range(len(sequence) - 1):
            dipeptide = sequence[i:i+2]
            frequencies[dipeptide] += 1

    # Calculate frequencies by dividing counts by the total number of dipeptides
    for dipeptide in frequencies:
        frequencies[dipeptide] /= total_sequences

    print("Total Sequences:", total_sequences)
    print("Dipeptide Frequencies:")
    for dipeptide, freq in frequencies.items():
        print(f"{dipeptide}: {freq}")

    return frequencies, total_sequences

def calculate_stderror_dipeptides(reference_frequencies, total_sequences):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    dipeptides = [f"{aa1}{aa2}" for aa1, aa2 in product(amino_acids, repeat=2)]
    # Calculate stderror for all dipeptides
    stderrors = {dipeptide: np.sqrt(reference_frequencies[dipeptide] * (1 - reference_frequencies[dipeptide]) * (1 / total_sequences)) for dipeptide in dipeptides}

    print("Standard Errors for Dipeptides:")
    for dipeptide, stderror in stderrors.items():
        print(f"{dipeptide}: {stderror}")

    return stderrors

def calculate_stderror_dipeptides(reference_path, sample_total_sequences):
    reference_frequencies, total_sequences = calculate_frequencies_dipeptides(reference_path)

    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    dipeptides = [f"{aa1}{aa2}" for aa1, aa2 in product(amino_acids, repeat=2)]

    # Calculate stderror for all dipeptides
    stderrors = {dipeptide: np.sqrt(reference_frequencies[dipeptide] * (1 - reference_frequencies[dipeptide]) * (1 / sample_total_sequences)) for dipeptide in dipeptides}

    print("Standard Errors for Dipeptides:")
    for dipeptide, stderror in stderrors.items():
        print(f"{dipeptide}: {stderror}")

    return reference_frequencies, stderrors, total_sequences

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
    z_output_file_path = os.path.join(os.path.dirname(sample_path), f"z-score-dipeptide-{os.path.basename(sample_path)}.csv")
    with open(z_output_file_path, "w", newline="") as z_csvfile:
        z_writer = csv.writer(z_csvfile)
        header = ["Dipeptide", "Z-Score"]
        z_writer.writerow(header)

        for dipeptide, z_score in z_scores.items():
            z_writer.writerow([dipeptide, z_score])

    # Save the z-scores heatmap as a PNG file
    save_heatmap_as_png(z_scores, sample_path)

    print(f"Output files created for {os.path.basename(sample_path)}: {output_file_path}, {z_output_file_path}")

def save_heatmap_as_png(z_scores, sample_path):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    dipeptides = [f"{aa1}{aa2}" for aa1, aa2 in product(amino_acids, repeat=2)]

    # Initialize a 2D array to store z-scores
    z_score_array = np.zeros((20, 20))

    for dipeptide in dipeptides:
        # Extract amino acids at positions P1 and P2 from the dipeptide
        p1, p2 = dipeptide
        # Find the index of the amino acids in the list
        idx_p1 = amino_acids.index(p1)
        idx_p2 = amino_acids.index(p2)

        # Populate the z_score_array with z-scores
        z_score_value = z_scores.get(dipeptide, 0)
        z_score_array[idx_p1, idx_p2] = z_score_value

    # Save the heatmap as a PNG file
    output_png_path = os.path.join(os.path.dirname(sample_path), f"DipepHM-{os.path.basename(sample_path)}.png")
    plt.imshow(z_score_array, cmap="coolwarm", interpolation="nearest")
    plt.colorbar(label="Z-Score")
    plt.xticks(np.arange(20), list(amino_acids))
    plt.yticks(np.arange(20), list(amino_acids))
    plt.xlabel("Amino Acid at P2")
    plt.ylabel("Amino Acid at P1")
    plt.title("Z-Scores for Dipeptides")
    plt.savefig(output_png_path, format="png")
    plt.close()

    print(f"Heatmap saved as PNG: {output_png_path}")

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
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    dipeptides = [f"{aa1}{aa2}" for aa1, aa2 in product(amino_acids, repeat=2)]

    for file_name in os.listdir(samples_folder):
        if file_name.startswith("z-score-dipeptide-") and file_name.endswith(".csv"):
            file_path = os.path.join(samples_folder, file_name)
            with open(file_path, newline="") as csvfile:
                reader = csv.DictReader(csvfile)
                z_scores = {row["Dipeptide"]: float(row["Z-Score"]) for row in reader}

            all_z_scores.append({"Sample": file_name.replace("z-score-dipeptide-", "").replace(".csv", ""), **z_scores})

    # Initialize a 2D array to store z-scores
    z_score_array = np.zeros((20, 20))

    for dipeptide in dipeptides:
        # Extract amino acids at positions P1 and P2 from the dipeptide
        p1, p2 = dipeptide
        # Find the index of the amino acids in the list
        idx_p1 = amino_acids.index(p1)
        idx_p2 = amino_acids.index(p2)

        # Populate the z_score_array with z-scores
        for z_scores in all_z_scores:
            sample_name = z_scores["Sample"]
            z_score_value = z_scores.get(dipeptide, 0)
            z_score_array[idx_p1, idx_p2] = z_score_value

    # Save the z_score_array as a CSV file with amino acid labels for both rows and columns
    output_csv_path = os.path.join(current_directory, "z_score_array.csv")
    np.savetxt(output_csv_path, z_score_array, delimiter=",", fmt="%1.3f", header=",".join(list(amino_acids)))

    print(f"Z-score heatmap and array saved to {output_csv_path}.")


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
