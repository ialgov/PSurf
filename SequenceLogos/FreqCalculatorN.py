import os
import re
from collections import defaultdict


def parse_fasta(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        sequence = ''
        header = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    weight = int(header.split('_')[-1])
                    if weight >= 2:
                        sequences.append((sequence, weight))
                    sequence = ''
                header = line[1:]
            else:
                sequence += line
        if sequence:
            weight = int(header.split('_')[-1])
            if weight >= 2:
                sequences.append((sequence, weight))
    return sequences


def calculate_frequencies(sequences):
    position_counts = defaultdict(lambda: defaultdict(int))
    total_counts = defaultdict(int)

    for sequence, weight in sequences:
        for i, aa in enumerate(sequence, start=1):
            position_counts[i][aa] += weight
            total_counts[i] += weight

    frequencies = []
    for position in sorted(position_counts.keys()):
        for aa in position_counts[position]:
            frequency = position_counts[position][aa] / total_counts[position]
            frequencies.append((aa, position, frequency))

    return frequencies


def write_csv(output_path, frequencies):
    with open(output_path, 'w') as file:
        file.write("Amino Acid,Position,Frequency\n")
        for aa, position, frequency in frequencies:
            file.write(f"{aa},{position},{frequency:.4f}\n")


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    fasta_files = [f for f in os.listdir(script_dir) if f.endswith('.fasta')]

    for fasta_file in fasta_files:
        file_path = os.path.join(script_dir, fasta_file)
        sequences = parse_fasta(file_path)
        frequencies = calculate_frequencies(sequences)

        output_file = fasta_file.replace('.fasta', '_frequencies.csv')
        output_path = os.path.join(script_dir, output_file)
        write_csv(output_path, frequencies)
        print(f"Processed {fasta_file} and saved results to {output_file}")


if __name__ == "__main__":
    main()