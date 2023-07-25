# variamotif_cli.py
import argparse
from Bio import SeqIO
from itertools import product
import multiprocessing
import search_motif_multi
import search_protein_multi
import search_motif_gap_multi
import motif_location_many

def main():
    # Your code to handle the command-line arguments goes here

    parser.add_argument('-T', action='store_true', help='Search in traditional mode')
    parser.add_argument('-V', action='store_true', help='Search in variable mode')
    parser.add_argument('-P', action='store_true', help='Search protein sequence')
    parser.add_argument('-D', action='store_true', help='Display your table')

    parser = argparse.ArgumentParser(description='Description of your tool')
    parser.add_argument('-f', '--fasta', type=str, help='Path to the fasta file')
    parser.add_argument('-motif', type=str, help='Target motif')
    parser.add_argument('-m', '--mismatches', default=0, type=int, help='Maximum allowed mismatches')
    parser.add_argument('-d', '--direction', default='+,-', type=str, choices=['+,-', '+', '-'], help='Search direction: both(default), forward, or reverse')
    
    #search_protein_multi.py #argument: -f -motif -m
 
    #-V (search variable length motif,only for DNA)
    #DNA motif variable
    parser.add_argument('-motif1',required=True, type=str, help='First target motif')
    parser.add_argument('-motif2',default="None", type=str, help='Second target motif')
    parser.add_argument('-min_g', '--min_gap', default=0, type=int, help='Minimum number of gaps between motifs')
    parser.add_argument('-max_g', '--max_gap', default=0, type=int, help='Maximum number of gaps between motifs')
 
    #V
    parser.add_argument('-t', '--table', dest='table_file', required=True, help='Input table file.')
    parser.add_argument('-r', '--display_both_directions', action='store_true', help='Display motifs from both forward and reverse strands.')


    args = parser.parse_args()

    # Call the appropriate function from your package based on the user inputs
    # For example: variamotif.some_function(args.input_file, args.output)
    if args.T:
        search_in_traditional_mode()
    elif args.V:
        search_in_variable_mode()
    elif args.P:
        search_protein_sequence()
    elif args.D:
        display_your_table()
    else:
        print("Please select a valid option.")

if __name__ == '__main__':
    main()
