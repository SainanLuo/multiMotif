import argparse
from Bio import SeqIO
from itertools import product
import multiprocessing

def generate_motif_variants(motif):
    variants = []
    ambiguous_bases = {
        'W': 'AT',
        'S': 'GC',
        'M': 'AC',
        'K': 'GT',
        'R': 'AG',
        'Y': 'CT',
        'B': 'CGT',
        'D': 'AGT',
        'H': 'ACT',
        'V': 'ACG',
        'N': 'ACGT'
    }
    variant_bases = []
    motif_variants = []

    for base in motif:
        if base in ambiguous_bases:
            variant_bases.append(list(ambiguous_bases[base]))
        else:
            variant_bases.append([base])

    motif_variants = product(*variant_bases)

    for variant in motif_variants:
        variants.append(''.join(variant))

    return variants

def calculate_mismatches(fragment, variants):
    mismatches_list = []

    for variant in variants:
        mismatches = sum([1 for base, variant_base in zip(fragment, variant) if variant_base != base])
        mismatches_list.append(mismatches)

    return min(mismatches_list)

def is_matching_position(base, motif_base):
    degenerate_bases = {
        'W': 'AT',
        'S': 'GC',
        'M': 'AC',
        'K': 'GT',
        'R': 'AG',
        'Y': 'CT',
        'B': 'CGT',
        'D': 'AGT',
        'H': 'ACT',
        'V': 'ACG',
        'N': 'ACGT'
    }

    if motif_base in degenerate_bases:
        allowed_bases = degenerate_bases[motif_base]
        return base in allowed_bases
    else:
        return base == motif_base

def search_fasta_record(record, motif, max_mismatches, direction):
    motif_variants = generate_motif_variants(motif)
    motif_length = len(motif)

    results = []
    sequence = str(record.seq)
    sequence_length = len(sequence)

    for i in range(sequence_length - motif_length + 1):
        if direction == '+' or direction == '+,-':
            fragment = sequence[i:i+motif_length]
            mismatches = calculate_mismatches(fragment, motif_variants)
            mismatch_positions = [i for i, (base, motif_base) in enumerate(zip(fragment, motif)) if not is_matching_position(base, motif_base)]

            if mismatches <= max_mismatches:
                result = {
                    'motif': motif,
                    'sequence_id': record.id,
                    'sequence_length':sequence_length,
                    'start': i,
                    'end': i + motif_length - 1,
                    'sequence': fragment,
                    'strand': '+',
                    'mismatch': mismatches,
                    'mismatch_positions': mismatch_positions
                }
                results.append(result)

        if direction == '-' or direction == '+,-':
            reverse_complement = str(record.seq.reverse_complement())
            reverse_start = sequence_length - (i + motif_length)
            reverse_end = sequence_length - i - 1
            start = sequence_length - reverse_end - 1
            end = sequence_length - reverse_start - 1
            reverse_fragment = reverse_complement[reverse_start:reverse_end+1]
            reverse_mismatches = calculate_mismatches(reverse_fragment, motif_variants)
            mismatch_positions = [i for i, (base, motif_base) in enumerate(zip(reverse_fragment, motif)) if not is_matching_position(base, motif_base)]

            if reverse_mismatches <= max_mismatches:
                result = {
                    'motif': motif,
                    'sequence_id': record.id,
                    'sequence_length':sequence_length,
                    'start': start,
                    'end': end,
                    'sequence': reverse_fragment,
                    'strand': '-',
                    'mismatch': reverse_mismatches,
                    'mismatch_positions': mismatch_positions
                }
                results.append(result)

    return results

def parse_arguments():
    parser = argparse.ArgumentParser(description='Search for motif in a fasta file with allowed mismatches.')
    parser.add_argument('-f', '--fasta', type=str, help='Path to the fasta file')
    parser.add_argument('-motif', type=str, help='Target motif')
    parser.add_argument('-m', '--mismatches', default=0, type=int, help='Maximum allowed mismatches')
    parser.add_argument('-d', '--direction', default='+,-', type=str, choices=['+,-', '+', '-'], help='Search direction: both(default), forward, or reverse')

    return parser.parse_args()

def main():
    args = parse_arguments()

    if not args.fasta or not args.motif:
        print("Please input fasta or motif")
    else:
        fasta_file_path = args.fasta
        motif = args.motif
        max_mismatches = args.mismatches
        direction = args.direction

        records = list(SeqIO.parse(fasta_file_path, "fasta"))
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        results = pool.starmap(search_fasta_record, [(record, motif, max_mismatches, direction) for record in records])
        results = [item for sublist in results for item in sublist]  # Flatten the list

        if not results:
            print("#No result")
        else:
            print("Sequence ID\tSequence Length\tmotif\tstart\tend\tstrand\tmismatch\tmismatch_positions\tprediction sequence") 
            for result in results:
                mismatch_positions_str = ",".join(str(pos) for pos in result['mismatch_positions'])
                print(f"{result['sequence_id']}\t{result['sequence_length']}\t{result['motif']}\t{result['start']}\t{result['end']}\t{result['strand']}\t{result['mismatch']}\t{mismatch_positions_str}\t{result['sequence']}")

if __name__ == "__main__":
    main()