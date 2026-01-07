import os
import argparse

def split_and_merge(input_file, output_dir):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Ensure the output directory exists.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Split the file into chunks of 22 lines.
    chunks = [lines[i:i+22] for i in range(0, len(lines), 22)]

    # Create a dictionary to store the merged counts.
    merged_data = {}

    # Process each chunk.
    for chunk in chunks:
        if len(chunk) == 22:
            # The last line of each chunk is used as the key.
            key = chunk[-1].strip()
            if key not in merged_data:
                merged_data[key] = chunk[:21]
            else:
                # Merge the counts for the same key.
                merged_data[key] = merge_counts(merged_data[key], chunk[:21])

    # Write the merged data to files.
    for key, value in merged_data.items():
        filename = key + '.txt'
        with open(os.path.join(output_dir, filename), 'w') as new_file:
            new_file.writelines(value)

def merge_counts(chunk1, chunk2):
    merged_chunk = []
    for line1, line2 in zip(chunk1, chunk2):
        parts1 = line1.split()
        parts2 = line2.split()
        merged_line = []
        for i in range(0, len(parts1), 2):
            amino_acid = parts1[i]
            count1 = int(parts1[i+1])
            count2 = int(parts2[i+1])
            merged_line.append(f"{amino_acid} {count1 + count2}")
        merged_chunk.append(' '.join(merged_line) + '\n')
    return merged_chunk

def format_conversion(input_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Format conversion from 8 columns to 2 columns.
    formatted_lines = []
    for line in lines:
        parts = line.split()
        for i in range(0, len(parts), 2):
            formatted_lines.append(f"{parts[i]} {parts[i+1]}\n")

    # Write the formatted content back to the same file.
    with open(input_file, 'w') as file:
        file.writelines(formatted_lines)

def calculate_percentages(input_file):
    amino_acid_counts = {}
    total_counts = 0

    with open(input_file, 'r') as file:
        lines = file.readlines()

    for line in lines:
        parts = line.split()
        amino_acid = parts[0]
        count = int(parts[1])
        amino_acid_counts[amino_acid] = amino_acid_counts.get(amino_acid, 0) + count
        total_counts += count

    percentages = {aa: count / total_counts for aa, count in amino_acid_counts.items()}
    return percentages

def main(input_file, output_dir):
    split_and_merge(input_file, output_dir)
    
    # Process each file in the output directory.
    for filename in os.listdir(output_dir):
        if filename.endswith('.txt'):
            input_path = os.path.join(output_dir, filename)
            format_conversion(input_path)
            percentages = calculate_percentages(input_path)
            
            # Write the percentages to a new file.
            output_path = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}_percentages.txt")
            sorted_percentages = sorted(percentages.items())
            with open(output_path, 'w') as file:
                for amino_acid, percentage in sorted_percentages:
                    file.write(f"{amino_acid}\t {percentage:.15f}\n")  # Display all decimal places

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a text file by splitting it into smaller files, merging duplicates, formatting the data, and calculating percentages.")
    parser.add_argument("input_file", help="The path to the input text file.")
    parser.add_argument("output_dir", help="The directory where the output files will be stored.")
    args = parser.parse_args()

    main(args.input_file, args.output_dir)
