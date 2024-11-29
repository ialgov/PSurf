import os
import csv
from Bio import SeqIO
import numpy as np
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="1-spaced")
    parser.add_argument("--reference", required=True, help="Path to the reference fasta file")
    parser.add_argument("--samples", nargs="+", required=True, help="Paths to sample fasta files")
    return parser.parse_args()

args = parse_arguments()

def calculate_total_sequences(file_path):
    # Count the total number of sequences in the file
    total_sequences = sum(1 for record in SeqIO.parse(file_path, "fasta"))
    return total_sequences

def calculate_frequencies_spaced_2mer(file_path):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    spaced_2mers = [f"{aa1}.{aa2}" for aa1 in amino_acids for aa2 in amino_acids]
    frequencies = {spaced_2mer: 0 for spaced_2mer in spaced_2mers}
    total_sequences = calculate_total_sequences(file_path)

    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq)
        for i in range(len(sequence) - 2):
            spaced_2mer = f"{sequence[i]}.{sequence[i + 2]}"
            frequencies[spaced_2mer] += 1

        # Calculate frequencies by dividing counts by the total number of spaced 2-mers
    for spaced_2mer in frequencies:
        frequencies[spaced_2mer] /= total_sequences

    print("Total Sequences:", total_sequences)
    print("spaced_2mer Frequencies:")
    for spaced_2mer, freq in frequencies.items():
        print(f"{spaced_2mer}: {freq}")

    return frequencies, total_sequences

def calculate_stderror(reference_frequencies, total_sequences):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    spaced_2mer = [f"{aa1}.{aa2}" for aa1 in amino_acids for aa2 in amino_acids]
    # Calculate stderror for all spaced_2mer
    stderrors = {spaced_2mer: np.sqrt(reference_frequencies[spaced_2mer] * (1 - reference_frequencies[spaced_2mer]) * (1 / total_sequences)) for spaced_2mer in spaced_2mer}

    print("Standard Errors:")
    for spaced_2mer, stderror in stderrors.items():
        print(f"{spaced_2mer}: {stderror}")

    return stderrors

def calculate_stderror_spaced_2mer(reference_path, sample_total_sequences):
    reference_frequencies, total_sequences = calculate_frequencies_spaced_2mer(reference_path)

    # Print the actual calculation for spaced_2mer "AA" and "AT"
    for spaced_2mer in ["A.A", "A.T"]:
        actual_calculation = f"sqrt({reference_frequencies[spaced_2mer]} * (1 - {reference_frequencies[spaced_2mer]}) * (1 / {sample_total_sequences}))"
        print(f"Actual Calculation for {spaced_2mer}: {actual_calculation}")

    stderrors = calculate_stderror(reference_frequencies, sample_total_sequences)

    return reference_frequencies, stderrors, total_sequences

def process_single_sample(reference_path, sample_path):
    sample_total_sequences = calculate_total_sequences(sample_path)
    control_frequencies, stderror, _ = calculate_stderror_spaced_2mer(reference_path, sample_total_sequences)

    frequencies, _ = calculate_frequencies_spaced_2mer(sample_path)

    output_file_path = os.path.join(os.path.dirname(sample_path), f"spaced_2mer-freq-{os.path.basename(sample_path)}.csv")
    with open(output_file_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        header = ["spaced_2mer", "Frequency"]
        writer.writerow(header)

        for spaced_2mer, freq in frequencies.items():
            writer.writerow([spaced_2mer, freq])

        # Calculate z-scores for spaced_2mer and create a new CSV file
        z_scores = calculate_z_scores_spaced_2mer(frequencies, control_frequencies, stderror)
        z_output_file_path = os.path.join(os.path.dirname(sample_path),
                                          f"z-score-spaced_2mer-{os.path.basename(sample_path)}.csv")
        with open(z_output_file_path, "w", newline="") as z_csvfile:
            z_writer = csv.writer(z_csvfile)
            header = ["spaced_2mer", "Z-Score"]
            z_writer.writerow(header)

            for spaced_2mer, z_score in z_scores.items():
                z_writer.writerow([spaced_2mer, z_score])

        print(f"Output files created for {os.path.basename(sample_path)}: {output_file_path}, {z_output_file_path}")

        # Create "uniques" file for each sample
        unique_spaced_2mers = set(spaced_2mer for spaced_2mer in frequencies if
                                  np.isnan(z_scores[spaced_2mer]) or z_scores[spaced_2mer] == float('inf'))
        unique_file_path = os.path.join(os.path.dirname(sample_path), f"uniques_{os.path.basename(sample_path)}.csv")
        with open(unique_file_path, "w", newline="") as unique_csvfile:
            unique_writer = csv.writer(unique_csvfile)
            header = ["spaced_2mer", "Frequency"]
            unique_writer.writerow(header)

            for spaced_2mer in unique_spaced_2mers:
                unique_writer.writerow([spaced_2mer, frequencies[spaced_2mer]])

        print(f"Unique spaced_2mers file created: {unique_file_path}")

def calculate_z_scores_spaced_2mer(sample_frequencies, control_frequencies, stderrors):
    z_scores = {spaced_2mer: 0 for spaced_2mer in sample_frequencies}

    # Calculate z-scores using the formula: (samplefreq - controlfreq) / stderror
    for spaced_2mer in sample_frequencies:
        z_scores[spaced_2mer] = (sample_frequencies[spaced_2mer] - control_frequencies[spaced_2mer]) / stderrors[spaced_2mer]

    # Print the z-score calculation for spaced_2mer "A.A" and "A.T"
    for spaced_2mer in ["A.A", "A.T"]:
        z_calculation = f"({sample_frequencies[spaced_2mer]} - {control_frequencies[spaced_2mer]}) / {stderrors[spaced_2mer]}"
        print(f"Z-Score Calculation for {spaced_2mer}: {z_calculation} = {z_scores[spaced_2mer]}")

    return z_scores

def combine_z_scores():
    current_directory = os.path.dirname(os.path.abspath(__file__))
    samples_folder = os.path.join(current_directory)

    all_z_scores = []
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    spaced_2mer = [f"{aa1}.{aa2}" for aa1 in amino_acids for aa2 in amino_acids]

    for file_name in os.listdir(samples_folder):
        if file_name.startswith("z-score-spaced_2mer-") and file_name.endswith(".csv"):
            file_path = os.path.join(samples_folder, file_name)
            with open(file_path, newline="") as csvfile:
                reader = csv.DictReader(csvfile)
                z_scores = {row["spaced_2mer"]: float(row["Z-Score"]) for row in reader}

            all_z_scores.append({"Sample": file_name.replace("z-score-spaced_2mer-", "").replace(".csv", ""), **z_scores})

    output_file_path = os.path.join(current_directory, "combined_z_scores_spaced_2mer.csv")
    with open(output_file_path, "w", newline="") as csvfile:
        fieldnames = ["spaced_2mer"] + [z_scores["Sample"] for z_scores in all_z_scores]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for spaced_2mer in spaced_2mer:
            row_data = {"spaced_2mer": spaced_2mer}
            for z_scores in all_z_scores:
                row_data[z_scores["Sample"]] = z_scores.get(spaced_2mer, 0)

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
