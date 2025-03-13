import subprocess
import argparse
import os
import tkinter as tk
from tkinter import filedialog

def choose_reference_file():
    root = tk.Tk()
    root.withdraw()
    reference_path = filedialog.askopenfilename(title="Choose a reference fasta file")
    return reference_path

def choose_sample_files():
    root = tk.Tk()
    root.withdraw()
    sample_paths = filedialog.askopenfilenames(title="Choose sample fasta files")
    return sample_paths

def main():
    # Step 1: Choose a reference fasta file
    reference_path = choose_reference_file()

    # Step 2: Choose multiple sample fasta files
    sample_paths = choose_sample_files()


    # Command for Z-scores-1mer.py
    command_1mer = [
        "python",
        "Z-scores-1mer.py",
        "--reference", reference_path,
        "--samples", *sample_paths
    ]


    # Command for Z-scores-2mer.py
    command_2mer = [
        "python",
        "Z-scores-2mer.py",
        "--reference", reference_path,
        "--samples", *sample_paths
    ]


    # Command for Z-scores-aaPerPosition.py
    command_aa_per_position = [
        "python",
        "Z-scores-aaPerPosition.py",
        "--reference", reference_path,
        "--samples", *sample_paths
    ]

    # Command for Z-scores-dipeptideCorrHeat.py
    command_dipeptide_corr_heat = [
        "python",
        "Z-scores-dipeptideCorrHeat.py",
        "--reference", reference_path,
        "--samples", *sample_paths
    ]

    # Execute the commands
    subprocess.run(command_1mer)
    subprocess.run(command_2mer)
    subprocess.run(command_aa_per_position)
    subprocess.run(command_dipeptide_corr_heat)

if __name__ == "__main__":
    main()
