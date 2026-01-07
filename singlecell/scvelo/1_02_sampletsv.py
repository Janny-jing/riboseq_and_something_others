import re
import sys

def process_samples(input_file):
    current_sample = None
    current_r1 = None
    current_r2 = None
    output = []
    
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            sample_match = re.match(r'^Processing sample (\d+):$', line)
            if sample_match:
                current_sample = sample_match.group(1)
                continue
                
            r1_match = re.match(r'^R1:\s+(.+)$', line)
            if r1_match and current_sample:
                current_r1 = r1_match.group(1)
                continue
                
            r2_match = re.match(r'^R2:\s+(.+)$', line)
            if r2_match and current_sample and current_r1:
                current_r2 = r2_match.group(1)
                output_line = f"sample {current_sample}\t{current_r1};{current_r2}"
                output.append(output_line)
                current_sample = current_r1 = current_r2 = None
                
    return '\n'.join(output)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input_file>")
        sys.exit(1)
    
    input_filename = sys.argv[1]
    try:
        result = process_samples(input_filename)
        print(result)
    except FileNotFoundError:
        print(f"Error: File '{input_filename}' not found")
        sys.exit(1)

