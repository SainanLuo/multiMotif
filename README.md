# VariaMotif
 Scanning and visualization for variable gap motifs.

## Installation

To install VariaMotif, you can use [pip](https://pip.pypa.io/en/stable/installation/):

```
pip install variamotif
```


Or clone the repository and then install it local:

```
git clone https://github.com/Luo-Sainan/VariaMotif.git
cd VariaMotif
pip install .
```

If you are a Windows user, make sure you installed [Python3](https://www.python.org/downloads/), then run this command in your `cmd` or `PowerShell`:

```
C:\Users\lsn> set PATH=C:\Users\lsn\AppData\Local\Programs\Python\Python311 #set your Python path in cmd window
C:\Users\lsn> set PATH=C:\Users\lsn\AppData\Local\Programs\Python\Python311\Scripts #set your pip path in cmd window
C:\Users\lsn> pip install variamotif
C:\Users\lsn> variamotif -h
```

## Dependencies

VariaMotif requires the following packages:

`Biopython (>= 1.78)`

`Matplotlib (>= 3.3)`

`NumPy (>= 1.19)`

`Pandas (>= 1.1)`

These dependencies will be installed automatically when you install VariaMotif using pip.

## Usage

After installation, you can use the `variamotif` command to run VariaMotif. Here is the basic usage:

```
variamotif -h

usage: VariaMotif.py [-h] [-extract_sequences] [-fna FNA] [-gff GFF]
                     [-up UPSTREAM] [-down DOWNSTREAM] [--promoter] [--orf]
                     [-VariaMotif] [-f FASTA] [-motif1 MOTIF1]
                     [-motif2 MOTIF2] [-min_g MIN_GAP] [-max_g MAX_GAP]
                     [-m MISMATCHES] [-d {+,-,+,-}] [-fix] [-variable] [-DNA]
                     [-RNA] [-protein] [-o OUTPUT] [-VisualMotif] [-i] [-r]
                     [-t TABLE_FILE]

VariaMotif for motif scanning

optional arguments:
  -h, --help            show this help message and exit
  -extract_sequences, --extract_sequences
                        extract promoter or orf sequences
  -fna FNA, --fna FNA   Input FNA file
  -gff GFF, --gff GFF   Input GFF file
  -up UPSTREAM, --upstream UPSTREAM
                        Gene start location upstream length (optional, default
                        is 400)
  -down DOWNSTREAM, --downstream DOWNSTREAM
                        Gene start location downstream length (optional,
                        default is 0)
  --promoter            Extract promoters
  --orf                 Extract ORFs
  -VariaMotif, --VariaMotif
                        motif scanning
  -f FASTA, --fasta FASTA
                        FASTA file path
  -motif1 MOTIF1        motif1,required=True
  -motif2 MOTIF2        motif2,default="None"
  -min_g MIN_GAP, --min_gap MIN_GAP
                        mix gap length between motif1 and motif2
  -max_g MAX_GAP, --max_gap MAX_GAP
                        max gap length between motif1 and motif2
  -m MISMATCHES, --mismatches MISMATCHES
                        max mismatches
  -d {+,-,+,-}, --direction {+,-,+,-}
                        Search direction: both, forward (default), or reverse
  -fix                  For fixed length motif
  -variable             For variable length motif
  -DNA                  For DNA variable motif
  -RNA                  For RNA variable motif
  -protein              For protein variable motif
  -o OUTPUT, --output OUTPUT
                        Output file for motif scanning result and Output file
                        prefix for display
  -VisualMotif, --VisualMotif
                        Display motif in sequence
  -i, --image           Display motif in sequence
  -r, --display_both_directions
                        Display motifs from both + and - strands.
  -t TABLE_FILE, --table TABLE_FILE
                        Input table file.
```

## Examples
(1) Fixed length DNA motif (CodY, AATTTTCWGAAAATT)

```
variamotif -VariaMotif -f GCA_000009045.1.promoter.fa -fix -DNA -motif1 AATTTTCWGAAAATT -m 2 -d +,- -o CodY.fix.out -i
```
(2) Variable length DNA motif (CcpA, TGTAAA-N(0-40)-TTTACA)

```
variamotif -VariaMotif -f GCA_000008765.1.promoter.fa -variable -DNA -motif1 TGTAAA -motif2 TTTACA -min_g 0 -max_g 40 -m 0 -d + -o CcpA.variable.out -i
```

(3) Fixed length protein motif

```
variamotif -VariaMotif -f 105.HTH.fa -fix -protein -motif1 GXTRSVIVN -m 0 -o one_motif.out -i
```

(4) Variable length protein motif

```
variamotif -VariaMotif -f 105.HTH.fa -variable -protein -motif1 GXTRSVIVN -motif2 LGMKGT -min_g 15 -max_g 15 -m 0 -o motif_gap.out -i
```
(5) RNA

```
variamotif -VariaMotif -f RNA_random.fa -fix -RNA -motif1 ACCGUUUUGAAAGGCG -m 0 -o motif.out -i
```

## Visualization Tool: VisualMotif

For visualization, if there are multiple motifs (3 or more), it is recommended to search separately according to the fixed length motif, merge the result files, and use the VisualMotif option for visualization.

```
variamotif -VisualMotif -t VariaMotif.result.txt -o single -r
```

To merge multiple motif scanning results:

```
variamotif -VisualMotif -t many_files.txt -o multi -r
```

## Utility for Extracting Gene Sequences

This utility can be used to extract either promoter sequences or open reading frames (ORFs) from a given genomic sequence. It takes as input a FASTA file containing the genomic sequence and a GFF file containing the gene annotations. Users can specify the desired upstream and downstream lengths relative to the start of each gene.


This example demonstrates how to extract the promoter sequences of genes from a given genomic sequence. The upstream length is set to 400 base pairs, and the downstream length is set to 0. The output is saved in a file named `GCA_000009045.1.promoter.fa`

```
variamotif -extract_sequences -fna GCA_000009045.1_ASM904v1_genomic.fna -gff genomic.gff --promoter -up 400 -down 0 -o GCA_000009045.1.promoter.fa
```


## Runtime Records

CodY, genome size 4.1M, 15bp motif, both strands, 2 mismatches, runtime 12 seconds.
CcpA, genome size 3.9M, two motifs, gap 0~40, both strands, 0 mismatches, runtime 13.851 seconds.

## Data Sources

CcpA flexible binding site (PMID: 28119470), Genome used: GCA_000008765.1 (Clostridium acetobutylicum ATCC 824)
CodY, fixed length motif (PMID: 18083814), GCA_000009045.1, AATTTTCWGAAAATT, Bacillus subtilis subsp. subtilis str. 168
