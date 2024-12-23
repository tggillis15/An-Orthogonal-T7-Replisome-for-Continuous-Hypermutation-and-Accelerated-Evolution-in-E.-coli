import pandas as pd
from collections import Counter, defaultdict
from multiprocessing import Pool
import os

# Define constants
input_file = "/gpfs/home/tgillis/mutspectra_pilot/merge_and_filter_starterpack_v7/Fastq2/Tem1-48/processed_data/filtered_sequences_only.csv"
output_dir = "/gpfs/home/tgillis/mutspectra_pilot/merge_and_filter_starterpack_v7/Fastq2/Tem1-48/processed_data"
os.makedirs(output_dir, exist_ok=True)

mutation_spectra_output_file = os.path.join(output_dir, "mutation_spectra.csv")
heatmap_output_file = os.path.join(output_dir, "mutation_heatmap.csv")

# WT reference sequence
reference_sequence = "TTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAA"

# Define the downstream sequence for trimming
downstream_sequence = "GCTAATCAGC"
barcode_length = 7

# Function to trim 3' barcodes based on downstream sequence
def trim_sequence(sequence, downstream_seq, barcode_len):
    downstream_pos = sequence.find(downstream_seq)
    if downstream_pos != -1:
        return sequence[:downstream_pos - barcode_len]
    return sequence  # Return the sequence as is if the downstream sequence is not found

# Function to identify mutations compared to the WT reference
def identify_mutations(reference, sequence):
    mutations = []
    for i, (ref_base, seq_base) in enumerate(zip(reference, sequence)):
        if ref_base != seq_base:
            mutation = f"{ref_base}{i + 1}{seq_base}"  # Format: "A1T"
            mutations.append(mutation)
    return mutations

# Load the filtered sequences
filtered_sequences = pd.read_csv(input_file, header=None)
print(f"Loaded {len(filtered_sequences)} sequences from {input_file}")

# Trim the sequences
trimmed_sequences = filtered_sequences[0].apply(
    lambda seq: trim_sequence(seq, downstream_sequence, barcode_length)
)

# Analyze mutations
mutation_counts = Counter()
def process_sequence(sequence):
    if len(sequence) == len(reference_sequence):  # Only compare sequences of equal length
        return identify_mutations(reference_sequence, sequence)
    return []

# Use multiprocessing to speed up mutation analysis
with Pool() as pool:
    all_mutations = pool.map(process_sequence, trimmed_sequences)

# Aggregate mutations
for mutations in all_mutations:
    mutation_counts.update(mutations)

# Save mutation spectra
mutation_spectra_df = pd.DataFrame.from_dict(mutation_counts, orient="index", columns=["Count"]).reset_index()
mutation_spectra_df.columns = ["Mutation", "Count"]
mutation_spectra_df.to_csv(mutation_spectra_output_file, index=False)
print(f"Mutation spectra saved to {mutation_spectra_output_file}")

# Generate mutation heatmap data (4x4 grid for A, T, G, C mutations)
bases = ["A", "T", "G", "C"]
heatmap = pd.DataFrame(0, index=bases, columns=bases)

for mutation, count in mutation_counts.items():
    ref_base = mutation[0]
    mut_base = mutation[-1]
    if ref_base in bases and mut_base in bases:
        heatmap.at[ref_base, mut_base] += count

# Save heatmap data
heatmap.to_csv(heatmap_output_file)
print(f"Mutation heatmap saved to {heatmap_output_file}")
