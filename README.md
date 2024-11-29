NGS data was processed using this Python script with multiple filtration steps.
The analysis begins with de-multiplexing sequencing reads into individual samples based on their DNA barcodes. 
Next, the script filters the DNA reads to retain only those containing the linker sequence downstream of the substrate region.
21 bp substrate sequences located upstream of the linker are extracted. 
The DNA sequences are being translated into amino acid sequences. 
Any substrate exhibiting 70% or greater similarity to the template sequence used for library construction (TEVs in this case) is excluded. 
Identical sequences are collapsed and counted, with the resulting data saved in a FASTA file. 
Sequence names indicated their frequency, such as seq1_249 representing the most abundant sequence, found 249 times in the dataset.

![image](https://github.com/user-attachments/assets/db9e97ed-eed1-4575-b500-819d623f4a33)
