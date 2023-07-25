import argparse
from Bio import SeqIO
from itertools import product
import multiprocessing
import time

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
    variant_bases = [ambiguous_bases.get(base, base) for base in motif]
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

def insert_gap(sequence, gap_length):
    return sequence + 'N' * gap_length
    
def forward_search(record, sequence, motif1, motif2, motif1_variants, motif2_variants, motif_length1, motif_length2, gap_length, max_mismatches, i):
    fragment1 = sequence[i:i+motif_length1]
    fragment2 = sequence[i+motif_length1+gap_length:i+motif_length1+gap_length+motif_length2]
    mismatches1 = calculate_mismatches(fragment1, motif1_variants)
    mismatches2 = calculate_mismatches(fragment2, motif2_variants)
    if mismatches1 + mismatches2 <= max_mismatches:
        result = {
            'motif':motif1 + ',' + motif2,
            'sequence_id': record.id,
            'start': i,
            'end': i + motif_length1 + gap_length + motif_length2 - 1,
            'mismatches1':mismatches1,
            'fragment1':fragment1,
            'mismatches2':mismatches2,
            'fragment2':fragment2,
            'mismatches':mismatches1 + mismatches2,
            'sequence': fragment1 + sequence[i+motif_length1:i+motif_length1+gap_length] + fragment2,
            'strand': '+',
            'gap': gap_length
        }
        return result

def reverse_result(result, sequence_length):
    result['strand'] = '-'
    end = result['end']
    start = result['start']
    result['start'] = sequence_length - 1 - end
    result['end'] = sequence_length - 1 - start
    return result
    
def search_fasta_record(record, motif1, motif2, min_gap_length, max_gap_length, max_mismatches, direction):
    motif1_variants = generate_motif_variants(motif1)
    motif2_variants = generate_motif_variants(motif2)
    motif_length1 = len(motif1)
    motif_length2 = len(motif2)

    results = []
    
    sequence = str(record.seq)
    sequence_length = len(sequence)

    for i in range(sequence_length - motif_length1 - max_gap_length - motif_length2 + 2):
        for gap_length in range(min_gap_length, max_gap_length+1):
            if direction == '+' or direction == '+,-':
                result = forward_search(record, sequence, motif1, motif2, motif1_variants, motif2_variants, motif_length1, motif_length2, gap_length, max_mismatches, i)
                if result is not None:
                    results.append(result)

            if direction == '-' or direction == '+,-':
                reverse_seq = str(record.seq.reverse_complement())
                result = forward_search(record, reverse_seq, motif1, motif2, motif1_variants, motif2_variants, motif_length1, motif_length2, gap_length, max_mismatches, i)
                if result is not None:
                    results.append(reverse_result(result, sequence_length))

    return sequence_length, results

def parse_arguments():
    parser = argparse.ArgumentParser(description='Search for motif in a fasta file with allowed mismatches.')
    parser.add_argument('-f', '--fasta', type=str, help='Path to the fasta file')
    parser.add_argument('-motif1',required=True, type=str, help='First target motif')
    parser.add_argument('-motif2',default="None", type=str, help='Second target motif')
    parser.add_argument('-min_g', '--min_gap', default=0, type=int, help='Minimum number of gaps between motifs')
    parser.add_argument('-max_g', '--max_gap', default=0, type=int, help='Maximum number of gaps between motifs')
    parser.add_argument('-m', '--mismatches', default=0, type=int, help='Maximum allowed mismatches')
    parser.add_argument('-d', '--direction', type=str, default='+,-', choices=['+,-', '+', '-'], help='Search direction: both (default), forward, or reverse')
    return parser.parse_args()

def main():
    args = parse_arguments()

    if not args.fasta or not args.motif1 or not args.motif2:
        print("请输入fasta文件路径、motif1和motif2")
    else:
        fasta_file_path = args.fasta
        motif1 = args.motif1
        motif2 = args.motif2
        min_gap_length = args.min_gap
        max_gap_length = args.max_gap
        max_mismatches = args.mismatches
        direction = args.direction

        records = list(SeqIO.parse(fasta_file_path, "fasta"))
        pool = multiprocessing.Pool(multiprocessing.cpu_count())

        # 使用starmap同时获取sequence_length和results
        results_with_lengths = pool.starmap(search_fasta_record, [(record, motif1, motif2, min_gap_length, max_gap_length, max_mismatches,direction) for record in records])

        # 解包sequence_length和results
        sequence_lengths, results_list = zip(*results_with_lengths)

        results = [item for sublist in results_list for item in sublist]  # 展平列表

        # 打印包含sequence_length的输出
        if not results:
            print("#没有结果")
        else:
            print("Sequence ID\tSequence Length\tmotif\tstart\tend\tstrand\tmotif1_mismatch\tfragment1\tmotif2_mismatch\tfragment2\tmismatch\tgap_length\tprediction sequence")
            for sequence_length, result in zip(sequence_lengths, results):
                print(f"{result['sequence_id']}\t{sequence_length}\t{result['motif']}\t{result['start']}\t{result['end']}\t{result['strand']}\t{result['mismatches1']}\t{result['fragment1']}\t{result['mismatches2']}\t{result['fragment2']}\t{result['mismatches']}\t{result['gap']}\t{result['sequence']}")

if __name__ == "__main__":
    main()