#!/usr/bin/env python

import os
import textwrap

from CryptoGenotyper.msr import msr_main
from CryptoGenotyper.gp60 import gp60_main
import argparse
from CryptoGenotyper.version import  __version__

def make_custom_database(input_fasta):
    print("Making custom database from fasta file")
    print(os.getcwd(), print(input_fasta))
    os.system("makeblastdb  -dbtype nucl -in "+input_fasta+" -out custom_db")


def parse_cli_arguments():
    parser = argparse.ArgumentParser(prog="cryptogenotyper",
                                     description='In silico type cryptosporidium from sanger reads in AB1 format\n',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog = textwrap.dedent('''Example usage using example ab1 files included in /example folder:
    cryptogenotyper -i example/P17705_Crypto16-2F-20170927_SSUF_G12_084.ab1 -m 18S -t forward -f SSUF -o test
    cryptogenotyper -i example/P17705_gp60-Crypt14-1F-20170927_gp60F_G07_051.ab1 -m gp60 -t forward -f gp60F -o test
    cryptogenotyper -i example/ -m 18S -t contig -f SSUF -r SSUR -o test
    cryptogenotyper -i example/ -m gp60 -t contig -f gp60F -r gp60R -o test
                                                              '''))

    parser.add_argument('--verbose', action='store_true', dest='verbose',
                        help='Turn on verbose logging [False].', required=False)
    parser.add_argument('-i','--input',  nargs=1,  required=True,
                        help="Path to directory with AB1 forward and reverse files OR path to a single AB1 file")
    parser.add_argument('-m', '--marker', type=str, required=True,
                        help="Name of the marker. Currently gp60 and 18S markers are supported")
    parser.add_argument('-t', '--seqtype', type=str, required=True,
                        help="Input sequences type. Select one option out of these three:\n"
                             "contig - both F and R sequences provided\n "
                             "forward - forward only sequence provided\n"
                             "reverse - reverse only sequence provided\n")
    parser.add_argument('-f', '--forwardprimername', type=str, required=False,
                        help="Name of the forward primer to identify forward read (e.g. gp60F, SSUF)")
    parser.add_argument('-r', '--reverseprimername', type=str, required=False,
                        help="Name of the reverse primer to identify forward read (e.g. gp60R, SSUR)")
    parser.add_argument('-o', '--outputprefix', type=str, required=False, default="cryptorun",
                        help="Output name prefix for the results (e.g. test results in test_report.fa)")
    parser.add_argument('-d', '--databasefile',type=str, required=False,
                        help="Path to custom database reference FASTA file")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__))
    parser.add_argument('--noheaderline',action='store_true',dest='header',help='Display header on tab-delimited file [False]', required=False)

    return parser.parse_args()

#in command line: sequences, marker, contig/f/r, fname, rname, expName
def main():
    args= parse_cli_arguments()

    print("Running cryptogenotyper v{}".format(__version__))

    if args.databasefile:
        make_custom_database(input_fasta=args.databasefile)

    seq_dir = args.input[0]

    marker = args.marker
    typeSeq = args.seqtype

    expName = args.outputprefix

    header = args.header

    pathlist=list()
    if os.path.isdir(seq_dir):
        seq_dir = os.path.abspath(seq_dir) #get directory name

        for file in os.listdir(seq_dir):
            if file.endswith("ab1"):
                pathlist.append(os.path.abspath(seq_dir + "/" + file))
    else:
        pathlist.append(os.path.abspath(seq_dir)) #get absolute path of a single file


    if typeSeq == "contig":
        fPrimer = args.forwardprimername
        rPrimer = args.reverseprimername

        if marker == "18S":
            return msr_main(pathlist, fPrimer, rPrimer, typeSeq, expName, args.databasefile, header)

        elif marker == "gp60":
            return gp60_main(pathlist, fPrimer, rPrimer, typeSeq, expName,args.databasefile, header)

    elif typeSeq == "forward":
        fPrimer = args.forwardprimername

        if marker == "18S":
            return msr_main(pathlist, fPrimer, "", typeSeq, expName, args.databasefile, header)

        elif marker == "gp60":
            return gp60_main(pathlist, fPrimer, "", typeSeq, expName, args.databasefile,header)


    elif typeSeq == "reverse":
        rPrimer = args.reverseprimername


        if marker == "18S":
            return msr_main(pathlist, '', rPrimer, typeSeq, expName, args.databasefile, header)

        elif marker == "gp60":
            return gp60_main(pathlist, '', rPrimer, typeSeq, expName, args.databasefile, header)




if __name__ == "__main__":
    main()
