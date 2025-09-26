<p align="center">
<img src="logo.jpg" />
</p>

[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=CryptoGenotyper)

## CryptoGenotyper

The `CryptoGenotyper` is a fast and reproducible tool that can be used to classify the genotype of *Cryptosporidium* samples by directly analyzing the DNA electropherogram files (.ab1) and DNA sequences in FASTA format that correspond to two of its characteristic gene markers: SSU rRNA and *gp60*. 

*Cryptosporidium* is a protozoan parasite that causes the enteric disease, cryptosporidiosis. It is transmitted to both humans and animals through zoonotic or anthroponotic means, and these dynamics can be studied through the analysis of SSU rRNA or *gp60* gene locus. Although due to the nature of these gene targets, manual analysis can be repetitive and difficult, allowing for the potential of inaccurate or incomplete results to be reported. 

CryptoGenotyper is able to analyze both well-defined and poorly-resolved peaks from SANGER sequencing files to ultimately output the corresponding sequence along with the *Cryptosporidium* genotype in standard nomenclature.

Since v1.5.0 CryptoGenotyper introduces an ability to process popular and widely available FASTA format both for *18S* and *gp60* markers. Importantly the gp60 marker module now incorporates the most recent typing nomenclature defined by the gp60 typing table [here](https://cryptodb.org/cryptodb/app/static-content/gp60.html) of the [Robinson, Guy, et al. "Deciphering a cryptic minefield: A guide to Cryptosporidium gp60 subtyping." Current Research in Parasitology & Vector-Borne Diseases (2025): 100257](https://www.sciencedirect.com/science/article/pii/S2667114X25000172) publication.

## Requirements
* `biopython >= 1.70,<1.78`
* `numpy >= 1.1`
* `python > 3.6`
* `blast == 2.9.0`
* `clustalw >= 2.1`

If running Ubuntu these requirements could be installed via one-liner command (might need preceed with `sudo`). Due to version restriction, biopython is installed via `pip3`

```
apt install clustalw ncbi-blast+ python3 python3-numpy python3-pip && pip3 install "biopython>=>=1.7,<1.78"
```

## Installation
The `CryptoGenotyper` can be installed by pulling source code from this repository or via `conda` package management system from the `bioconcda` channel.

* `git clone https://github.com/phac-nml/CryptoGenotyper.git && cd CryptoGenotyper && python3 setup.py install`
* `conda install -c bioconda cryptogenotyper`

## Usage
Only few parameters required to run `cryptogenotyper` including input `*.ab1` file, marker (`18S` or `gp60`), sequence type (`forward`, `reverse`, `contig`), output prefix. 


```
usage: cryptogenotyper [-h] [--verbose] -i INPUT -m MARKER -t SEQTYPE
                       [-f FORWARDPRIMERNAME] [-r REVERSEPRIMERNAME]
                       [-o OUTPUTPREFIX] [-d DATABASEFILE] [-v]
                       [--noheaderline]

In silico type cryptosporidium from sanger reads in AB1 format

optional arguments:
  -h, --help            show this help message and exit
  --verbose             Turn on verbose logging [False].
  -i INPUT, --input INPUT
                        Path to directory with AB1/FASTA forward and reverse files
                        OR path to a single AB1/FASTA file
  -m MARKER, --marker MARKER
                        Name of the marker. Currently gp60 and 18S markers are
                        supported
  -t SEQTYPE, --seqtype SEQTYPE
                        Input sequences type. Select one option out of these
                        three: contig - both F and R sequences provided
                        forward - forward only sequence provided reverse -
                        reverse only sequence provided
  -f FORWARDPRIMERNAME, --forwardprimername FORWARDPRIMERNAME
                        Name of the forward primer to identify forward read
                        (e.g. gp60F, SSUF)
  -r REVERSEPRIMERNAME, --reverseprimername REVERSEPRIMERNAME
                        Name of the reverse primer to identify forward read
                        (e.g. gp60R, SSUR)
  -o OUTPUTPREFIX, --outputprefix OUTPUTPREFIX
                        Output name prefix for the results (e.g. test results
                        in test_report.fa)
  -d DATABASEFILE, --databasefile DATABASEFILE
                        Custom database reference file
  -v, --version         show program's version number and exit
  --noheaderline        Display header on tab-delimited file [False]
```

## Examples
The `example` folder contains a couple of sequences to try out by executing the following commands. The `-i` accepts either a single file input or directory. The `-f` and `-r` parameters allow to filter input folder files in case multiple files need to be processed, they are also used to create file pairs used for `contig` mode (i.e. `-t contig`).

```
#Use a single Sanger file in a forward-only mode for the 18S marker
cryptogenotyper -i example/P17705_Crypto16-2F-20170927_SSUF_G12_084.ab1 -m 18S -t forward -f SSUF -o test
cryptogenotyper -i example/

#Use a single Sanger file in a forward-only mode for the gp60 marker
P17705_gp60-Crypt14-1F-20170927_gp60F_G07_051.ab1 -m gp60 -t forward -f gp60F -o test

#use a pair of Sanger files in contig mode for the GP60 marker, using forward and reverse file filters specified by -f and -r
cryptogenotyper -i example/ -m 18S -t contig -f SSUF -r SSUR -o test
cryptogenotyper -i example/ -m gp60 -t contig -f gp60F -r gp60R -o test

#Use a pair of FASTA files in contig mode for the GP60 marker found in the current directory
cryptogenotyper -i . -f "forward.fasta" -r "reverse.fasta" -o tmp -t contig -m gp60

#Use a single FASTA file in a forward-only mode for the GP60 marker
cryptogenotyper -i ./tests/Data/Cmortiferum_gp60_f2.fasta  -o tmp -t forward -m gp60

#to filter multiple files found in a directory, such as tests/Data/, use the filters specified by -f and/or -r
cryptogenotyper -i ./tests/Data/  -f "_gp60_f" -o tmp -t forward -m gp60

#to pair multiple files in contig mode, use the -f and/or -r to select forward and reverse files from a directory, such as tests/Data/ 
cryptogenotyper -i ./tests/Data/  -f "Cmortiferum_gp60_f" -r "Cmortiferum_gp60_r" -o tmp -t contig -m gp60
```

## Galaxy workflows
Galaxy workflows for running multiple files and instructions are found in the [`GalaxyWorkflows`](./Cryptogenotyper/GalaxyWorkflows/README.md) folder, along with installation instructions. 

## Citation
Please cite the following publication if you find this subtyping tool useful in your work.


Yanta, C.A., Bessonov, K., Robinson, G., Troell, K., Guy, R.A. (2021) CryptoGenotyper: A new bioinformatics tool for rapid _Cryptosporidium_ identification. _Food and Waterborne Parasitology_, 23:e00115. doi: [10.1016/j.fawpar.2021.e00115](https://doi.org/10.1016/j.fawpar.2021.e00115)

