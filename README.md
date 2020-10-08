<div style="text-align:center"><img src="logo.jpg" /></div>

## CryptoGenotyper

The `CryptoGenotyper` is a fast and reproducible tool that can be used to classify the genotype of *Cryptosporidium* samples by directly analyzing the DNA electropherogram files (.ab1) that correspond to two of its characteristic gene markers: SSU rRNA and *gp60*. 

*Cryptosporidium* is a protozoan parasite that causes the enteric disease, cryptosporidiosis. It is transmitted to both humans and animals through zoonotic or anthroponotic means, and these dynamics can be studied through the analysis of SSU rRNA or *gp60* gene locus. Although due to the nature of these gene targets, manual analysis can be repetitive and difficult, allowing for the potential of inaccurate or incomplete results to be reported. 

CryptoGenotyper is able to analyze both well-defined and poorly-resolved peaks to ultimately output the corresponding sequence along with the *Cryptosporidium* genotype in standard nomenclature.

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
* `conda install -c bioconda cryptotyper`

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
                        Path to directory with AB1 forward and reverse files
                        OR path to a single AB1 file
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
The `example` folder contains a couple of sequences to try out by executing the following commands.

```
cryptogenotyper -i example/P17705_Crypto16-2F-20170927_SSUF_G12_084.ab1 -m 18S -t forward -f SSUF -o test
cryptogenotyper -i example/P17705_gp60-Crypt14-1F-20170927_gp60F_G07_051.ab1 -m gp60 -t forward -f gp60F -o test
cryptogenotyper -i example/ -m 18S -t contig -f SSUF -r SSUR -o test
cryptogenotyper -i example/ -m gp60 -t contig -f gp60F -r gp60R -o test

```

## Citation
Please cite the following publication if you find useful this subtyping tool in your work.


Christine A. Yanta, Kyrylo Bessonov, Guy Robinson, Karin Troell, Rebecca A. Guy. *CryptoGenotyper: a new bioinformatics tool for rapid Cryptosporidium identification.* (Submitted to [FAWPAR journal](https://www.journals.elsevier.com/food-and-waterborne-parasitology))

