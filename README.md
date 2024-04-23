# multiMotif
 A tool for Scanning and Visualization of Diverse and Distant Multiple Motifs.

## Installation

To install multiMotif, you can use [pip](https://pip.pypa.io/en/stable/installation/):

```
pip install multiMotif
```


Or clone the repository and then install it local:

```
git clone https://github.com/SainanLuo/multiMotif.git
cd multiMotif
pip install .
```

If you are a Windows user, make sure you installed [Python3](https://www.python.org/downloads/), then run this command in your `cmd` or `PowerShell`:

```
C:\Users\lsn> set PATH=C:\Users\lsn\AppData\Local\Programs\Python\Python311 #set your Python path in cmd window
C:\Users\lsn> set PATH=C:\Users\lsn\AppData\Local\Programs\Python\Python311\Scripts #set your pip path in cmd window
C:\Users\lsn> pip install multiMotif
C:\Users\lsn> multiMotif -h
```

## Dependencies

multiMotif requires the following packages:

`Biopython (>= 1.78)`

`Matplotlib (>= 3.3)`

`NumPy (>= 1.19)`

`Pandas (>= 1.1)`

`plotly (>=5.18.0)`

These dependencies will be installed automatically when you install multiMotif using pip.

## Usage

After installation, you can use the `multiMotif` command to run multiMotif. Here is the basic usage:

```
usage: multiMotif [-h]
                  {singleMotif,multiMotif,visualMotif,extract_sequences} ...

VariaMotif for motif scanning

positional arguments:
  {singleMotif,multiMotif,visualMotif,extract_sequences}
                        sub-command help
    singleMotif         Scanning for single motif
    multiMotif          Scanning more than two ordered motifs
    visualMotif         Visualization:display motif in sequence
    extract_sequences   Extract Sequences

optional arguments:
  -h, --help            show this help message and exit
```

multiMotif tool features modules for identifying single motifs(‘singleMotif’), multiple motifs(‘multiMotif’), a visualization module (‘visualMotif’) to enhance results interpretation, and a sequence extraction module (‘extract_sequences’) for obtaining promoters or open reading frames, or user-specified locations. 

We provide comprehensive user manuals in both Chinese and English. You can download the English version of the manual at [English Language User Manual Download Link](https://github.com/SainanLuo/multiMotif/releases/download/v1.2.1/multiMotif_Installation_UsageInstructions_EnglishLanguage.pdf). For the Chinese version, please visit this link:[Chinese Language User Manual Download Link](https://github.com/SainanLuo/multiMotif/releases/download/v1.2.1/multiMotif_Installation_UsageInstructions_ChineseLanguage.pdf).


