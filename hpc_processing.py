import pandas as pd
from collections import Counter, defaultdict
import os
from multiprocessing import Pool

# Define constants
input_file = os.path.expanduser("~/mutspectra_pilot/merge_and_filter_starterpack_v7/Fastq2/Tem1-48/good_reads_Tem1-48.csv")
output_dir = os.path.expanduser("~/mutspectra_pilot/merge_and_filter_starterpack_v7/Fastq2/Tem1-48/processed_data")
os.makedirs(output_dir, exist_ok=True)

normalized_output_file = os.path.join(output_dir, "normalized_unique_sequences_Tem49-96.csv")
mutation_output_file = os.path.join(output_dir, "mutation_spectra_Tem49-96.csv")
barcode_sequences_output_file = os.path.join(output_dir, "barcode_sequences_Tem49-96.csv")
filtered_output_file = os.path.join(output_dir, "filtered_barcode_sequences_Tem49-96.csv")
mutated_sequences_output_file = os.path.join(output_dir, "mutated_sequences_Tem49-96.csv")
filtered_sequences_only_output_file = os.path.join(output_dir, "filtered_sequences_only.csv")

# Define the upstream and downstream constant sequences
upstream_sequence = "GGCATTCTGCTG"  # Replace with your upstream sequence
downstream_sequence = "GCTAATCAGC"  # Replace with your downstream sequence
barcode_length = 14  # Length of each barcode (7 bp each)

# Load all rows from the CSV file (no headers)
df = pd.read_csv(input_file, header=None)  # Full 450 million reads
print(f"Loaded {len(df)} rows from {input_file}")

# Assume sequences are in the first column
sequences = df[0]

# Function to extract the sequence between the barcodes, excluding both barcodes
def extract_barcodes_and_sequence(sequence, upstream_seq, downstream_seq, barcode_len):
    # Find the upstream barcode position
    upstream_pos = sequence.find(upstream_seq)
    # Find the downstream barcode position
    downstream_pos = sequence.find(downstream_seq)

    if upstream_pos != -1 and downstream_pos != -1:
        # Extract the sequence between the two barcodes, excluding both barcodes
        sequence_between_barcodes = sequence[upstream_pos + len(upstream_seq) + barcode_len : downstream_pos]
        return sequence_between_barcodes, upstream_pos, downstream_pos  # Return also the positions for reference
    return None, None, None  # Return None if either barcode is not found

# Parallelize the sequence processing across CPUs
def process_sequences(sequences):
    barcode_to_sequences = defaultdict(list)
    for seq in sequences:
        between_sequence, upstream_pos, downstream_pos = extract_barcodes_and_sequence(seq, upstream_sequence, downstream_sequence, barcode_length)
        if between_sequence:
            # Combine both barcodes to create a 14-bp barcode (upstream + downstream)
            barcode1 = seq[upstream_pos + len(upstream_sequence):upstream_pos + len(upstream_sequence) + barcode_length]
            barcode2 = seq[downstream_pos - barcode_length:downstream_pos]  # Fixed downstream barcode extraction
            combined_barcode = barcode1 + barcode2
            barcode_to_sequences[combined_barcode].append(between_sequence)
    return barcode_to_sequences

# Split the sequences into chunks for parallel processing
num_chunks = 50  # Number of chunks based on available CPUs
chunk_size = len(sequences) // num_chunks
chunks = [sequences[i:i + chunk_size] for i in range(0, len(sequences), chunk_size)]

# Use multiprocessing to process the sequences in parallel
with Pool(processes=num_chunks) as pool:
    results = pool.map(process_sequences, chunks)

# Combine the results from all chunks
barcode_to_sequences = defaultdict(list)
for result in results:
    for barcode, seqs in result.items():
        barcode_to_sequences[barcode].extend(seqs)

# Count how many times each barcode appears
barcode_counts = {barcode: len(seqs) for barcode, seqs in barcode_to_sequences.items()}
barcode_df = pd.DataFrame(list(barcode_counts.items()), columns=["Barcode", "Count"])
barcode_df.to_csv(barcode_sequences_output_file, index=False)

# Create a list of (barcode, sequence, sequence_count) for the detailed DataFrame
barcode_sequences_data = []
for barcode, seqs in barcode_to_sequences.items():
    seq_counts = Counter(seqs)  # Count occurrences of each sequence
    for sequence, count in seq_counts.items():
        barcode_sequences_data.append((barcode, sequence, count))

# Create the detailed DataFrame
barcode_sequences_df = pd.DataFrame(barcode_sequences_data, columns=["Barcode", "Sequence", "Sequence_Count"])

# Filter the sequences by barcode count (minimum count of 7 or more)

filtered_barcode_sequences_df = barcode_sequences_df[barcode_sequences_df['Sequence_Count'] >= 7]
filtered_barcode_sequences_df.to_csv(filtered_output_file, index=False)

# Remove "Barcode" and "Sequence_Count" columns and recount unique sequences
filtered_sequences_only_df = filtered_barcode_sequences_df.drop(columns=["Barcode", "Sequence_Count"])
unique_sequences = filtered_sequences_only_df["Sequence"].nunique()
print(f"Number of unique sequences after filtering: {unique_sequences}")

# Save filtered sequences only
filtered_sequences_only_df.to_csv(filtered_sequences_only_output_file, index=False)

# Normalize by unique barcode and sequence combinations
normalized_df = filtered_barcode_sequences_df.drop_duplicates(subset=["Barcode", "Sequence"])
normalized_df.to_csv(normalized_output_file, index=False)

# Find the most abundant sequence (WT reference sequence) based on unique sequences
sequence_counts = normalized_df["Sequence"].value_counts()
reference_sequence = sequence_counts.idxmax()
print(f"Reference Sequence (WT): {reference_sequence}")

# Function to identify mutations compared to the reference sequence
def identify_mutations(reference, sequence):
    mutations = []
    for i, (ref_base, seq_base) in enumerate(zip(reference, sequence)):
        if ref_base != seq_base:
            mutation = f"{ref_base}{i+1}{seq_base}"  # Format: "A1T"
            mutations.append(mutation)
    return mutations

# Analyze mutations
mutation_counts = Counter()
mutated_sequences = []

for _, row in normalized_df.iterrows():
    sequence = row['Sequence']
    if len(sequence) == len(reference_sequence):  # Only compare sequences of equal length
        mutations = identify_mutations(reference_sequence, sequence)
        if mutations:  # Sequence has mutations
            mutated_sequences.append((row['Barcode'], sequence, mutations))
            mutation_counts.update(mutations)

# Create a DataFrame for mutated sequences
mutated_sequences_df = pd.DataFrame(mutated_sequences, columns=["Barcode", "Sequence", "Mutations"])
mutated_sequences_df.to_csv(mutated_sequences_output_file, index=False)

# Generate mutation spectra
mutation_spectra = pd.DataFrame.from_dict(mutation_counts, orient='index', columns=["Count"]).reset_index()
mutation_spectra.columns = ["Mutation", "Count"]
mutation_spectra.to_csv(mutation_output_file, index=False)

print(f"Mutation spectra saved to {mutation_output_file}")
print(f"Mutated sequences saved to {mutated_sequences_output_file}")
print("We are done!")
