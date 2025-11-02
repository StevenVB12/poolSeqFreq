import sys

def summarize_positions(file_path):
    # Initialize a dictionary to hold the count of unique positions for each sample
    unique_positions_count = {}

    # Open and process the file
    with open(file_path, 'r') as file:
        for line in file:
            sample, scaffold, position, genotype, max_allele, min_allele = line.strip().split('\t')
            # Concatenate scaffold and position to create a unique identifier for each position
            unique_position = f"{scaffold}_{position}"
            if sample not in unique_positions_count:
                unique_positions_count[sample] = set()
            unique_positions_count[sample].add(unique_position)

    # Output the summary to stdout
    for sample, positions in unique_positions_count.items():
        print(f"{sample}\t{len(positions)}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python summarize_positions.py <input_file>")
        sys.exit(1)
    input_file = sys.argv[1]
    summarize_positions(input_file)
