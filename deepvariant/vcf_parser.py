import sys
import gzip

# Define the allowed chromosomes for V2P (1-22, X, Y).
ALLOWED_CHROMS = set([str(i) for i in range(1, 23)] + ['X', 'Y'])

def process_vcf(file_path):
    """
    Reads a VCF (plain text or .gz).
    Filters:
      1. Allowed chromosomes (1-22, X, Y)
      2. FILTER column must be 'PASS' (Removes RefCall, LowQual, etc.)
    Outputs:
      Valid VCF lines to stdout with 'chr' prefixes stripped.
    """
    # 1. Write the Header required by V2P
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    try:
        # Automatically handle gzip or plain text
        if file_path.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'

        with open_func(file_path, mode) as vcf:
            for line in vcf:
                # Skip existing header lines
                if line.startswith('#'):
                    continue
                
                cols = line.strip().split('\t')
                
                # Ensure we have enough columns (FILTER is index 6)
                if len(cols) < 7:
                    continue
                
                # --- NEW FILTERING LOGIC ---
                # 2. Check the FILTER column (Index 6)
                # If it is not 'PASS', skip it. This removes 'RefCall'.
                filter_status = cols[6]
                if filter_status != 'PASS':
                    continue
                # ---------------------------

                # 3. Normalize Chromosome Name
                raw_chrom = cols[0]
                norm_chrom = raw_chrom.lower().replace('chrom', '').replace('chr', '').upper()
                
                # 4. Filter "Poison" Chromosomes (M, Un, etc.)
                if norm_chrom not in ALLOWED_CHROMS:
                    continue

                # 5. Set Output to ONLY the number/letter (e.g. '1' not 'chr1')
                cols[0] = norm_chrom
                
                # Extract first 8 fields and print
                selected_fields = cols[:8]
                print("\t".join(selected_fields))

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python vcf_parser.py <input_vcf_file>")
    else:
        input_file = sys.argv[1]
        process_vcf(input_file)

