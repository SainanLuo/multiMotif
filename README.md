# VariaMotif
 A Computational Tool for Variable Motif scanning and Sequence-based Relative Position Visualization of Search Results in Sequences.

## Installation

To install VariaMotif, you can clone the repository and then install it using [pip](https://pip.pypa.io/en/stable/installation/):

```
git clone https://github.com/Luo-Sainan/VariaMotif.git
cd VariaMotif
pip install .
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

usage: variamotif [-h] [-VariaMotif] [-f FASTA] [-motif1 MOTIF1] [-motif2 MOTIF2] [-min_g MIN_GAP] [-max_g MAX_GAP]
                  [-m MISMATCHES] [-d {+,-,+,-}] [-fix] [-variable] [-DNA] [-RNA] [-protein] [-o OUTPUT] [-VisualMotif]
                  [-i] [-r] [-t TABLE_FILE]

VariaMotif for motif scanning

optional arguments:
  -h, --help            show this help message and exit
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
                        Output file for motif scanning result and Output file prefix for display
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
variamotif -f GCA_000009045.1.promoter.fa -fix -DNA -motif1 AATTTTCWGAAAATT -m 2 -d +,- -o CodY.fix.out -i
```
(2) Variable length DNA motif (CcpA, TGTAAA-N(0-40)-TTTACA)

```
variamotif -f GCA_000008765.1.promoter.fa -variable -DNA -motif1 TGTAAA -motif2 TTTACA -min_g 0 -max_g 40 -m 0 -d + -o CcpA.variable.out -i
```

(3) Fixed length protein motif

```
variamotif -f 105.HTH.fa -fix -protein -motif1 GXTRSVIVN -m 0 -o one_motif.out -i
```

(4) Variable length protein motif

```
variamotif -f 105.HTH.fa -variable -protein -motif1 GXTRSVIVN -motif2 LGMKGT -min_g 15 -max_g 15 -m 0 -o motif_gap.out -i
```
(5) RNA

```
variamotif -f RNA_random.fa -fix -RNA -motif1 ACCGUUUUGAAAGGCG -m 0 -o motif.out -i
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

## Runtime Records

CodY, genome size 4.1M, 15bp motif, both strands, 2 mismatches, runtime 12 seconds.
CcpA, genome size 3.9M, two motifs, gap 0~40, both strands, 0 mismatches, runtime 19 seconds.

## Data Sources

CcpA flexible binding site (PMID: 28119470), Genome used: GCA_000008765.1 (Clostridium acetobutylicum ATCC 824)
CodY, fixed length motif (PMID: 18083814), GCA_000009045.1, AATTTTCWGAAAATT, Bacillus subtilis subsp. subtilis str. 168