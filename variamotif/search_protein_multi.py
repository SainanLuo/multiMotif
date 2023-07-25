import argparse
from Bio import SeqIO
import multiprocessing

def calculate_mismatches(fragment, motif):
    return sum([1 for aa, motif_aa in zip(fragment, motif) if motif_aa != aa])

def search_fasta_record(record, motif, max_mismatches):
    results = []
    sequence = str(record.seq)
    motif_length = len(motif)
    sequence_length = len(sequence)

    for i in range(sequence_length - motif_length + 1):
        fragment = sequence[i:i+motif_length]
        mismatches = calculate_mismatches(fragment, motif)

        if mismatches <= max_mismatches:
            result = {
                'motif': motif,
                'sequence_id': record.id,
                'start': i,
                'end': i + motif_length - 1,
                'strand': "+",
                'sequence': fragment,
                'mismatch': mismatches
            }
            results.append(result)

    return results


# Command-line argument parsing
def parse_arguments():
    parser = argparse.ArgumentParser(description='Search for motif in a fasta file with allowed mismatches.')
    parser.add_argument('-f', '--fasta', type=str, help='Path to the fasta file')
    parser.add_argument('-motif', type=str, help='Target motif')
    parser.add_argument('-m', '--mismatches', default=0, type=int, help='Maximum allowed mismatches')

    return parser.parse_args()

def main():
    args = parse_arguments()

    if not args.fasta or not args.motif:
        print("Please input fasta or motif")
    else:
        fasta_file_path = args.fasta
        motif = args.motif
        max_mismatches = args.mismatches

        records = list(SeqIO.parse(fasta_file_path, "fasta"))
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        results = pool.starmap(search_fasta_record, [(record, motif, max_mismatches) for record in records])
        results = [item for sublist in results for item in sublist]  # Flatten the list

        if not results:
            print("#No result")
        else:
            print("Sequence ID\tmotif\tstart\tend\tstrand\tmismatch\tprediction sequence") 
            for result in results:
                print(f"{result['sequence_id']}\t{result['motif']}\t{result['start']}\t{result['end']}\t{result['strand']}\t{result['mismatch']}\t{result['sequence']}")

if __name__ == "__main__":
    main()