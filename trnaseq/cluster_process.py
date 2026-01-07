import sys

def process_file(input_file, output_file):
    seen = {}
    unique_lines = []

    with open(input_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')  # 假设列之间用制表符分隔
            key = parts[0]
            value = '\t'.join(parts[1:])
            
            if key not in seen:
                seen[key] = True
                unique_lines.append((key, value))

    with open(output_file, 'w') as file:
        for key, value in unique_lines:
            file.write(f"{key}\t{value}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_file.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    process_file(input_file, output_file)
