#!/usr/bin/env python

import os, logging
import textwrap

from CryptoGenotyper.msr import msr_main
from CryptoGenotyper.gp60 import gp60_main
import argparse
from CryptoGenotyper.version import  __version__
from CryptoGenotyper.logging import create_logger
import CryptoGenotyper.definitions as definitions
from CryptoGenotyper import utilities

LOG = create_logger() #create a root logger with settings that will be used by other module loggers


def make_custom_database(input_fasta):
    LOG.info("Making custom database from FASTA file")
    LOG.info(f"Working in {os.getcwd()} on {input_fasta}")
    os.system("makeblastdb  -dbtype nucl -in "+input_fasta+" -out custom_db -parse_seqids")


def parse_cli_arguments():
    parser = argparse.ArgumentParser(prog="cryptogenotyper",
                                     description='In silico type cryptosporidium from sanger reads in AB1 format\n',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog = textwrap.dedent('''Example usage using example ab1 files included in /example folder:
    cryptogenotyper -i example/P17705_Crypto16-2F-20170927_SSUF_G12_084.ab1 -m 18S -t forward -f SSUF -o test
    cryptogenotyper -i example/P17705_gp60-Crypt14-1F-20170927_gp60F_G07_051.ab1 -m gp60 -t forward -f gp60F -o test
    cryptogenotyper -i example/ -m 18S -t contig -f SSUF -r SSUR -o test
    cryptogenotyper -i example/ -m gp60 -t contig -f gp60F -r gp60R -o test
    cryptogenotyper -i example/input.fasta -m gp60 -o test
    cryptogenotyper -i example/ -m gp60 -o test -f input                                                          
                                                              '''))

    parser.add_argument('--verbose', action='store_true', dest='verbose',
                        help='Turn on verbose logging [False].', required=False)
    parser.add_argument('-i','--input',  nargs=1,  required=True,
                        help='''Path to SINGLE directory with AB1/FASTA forward and reverse files OR path to a SINGLE AB1/FASTA file. 
                        Use -f and/or -r to filter inputs''')
    parser.add_argument('-m', '--marker', type=str, required=True,
                        help="Name of the marker. Currently gp60 and 18S markers are supported")
    parser.add_argument('-t', '--seqtype', type=str, required=False, default="forward",
                        help="Input sequences type. Select one option out of these three:\n"
                             "contig - both F and R sequences provided\n "
                             "forward - forward only sequence provided\n"
                             "reverse - reverse only sequence provided\n")
    parser.add_argument('-f', '--forwardprimername', type=str,  required=False,
                        help="Name of the forward primer to identify forward read (e.g. gp60F, SSUF)")
    parser.add_argument('-r', '--reverseprimername', type=str, required=False,
                        help="Name of the reverse primer to identify forward read (e.g. gp60R, SSUR)")
    parser.add_argument('-s', '--suffix', type=str, required=False, default=None,
                        help="Optional suffix to filter filenames (e.g. only include files ending with a specific pattern)")
    parser.add_argument('-o', '--outputprefix', type=str, required=False, default="cryptorun",
                        help="Output name prefix for the results (e.g. test results in test_report.fa)")
    parser.add_argument('-d', '--databasefile',type=str, required=False, default=None,
                        help="Path to custom database reference FASTA file")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__))
    parser.add_argument('--noheaderline',action='store_true',dest='header',help='Display header on tab-delimited file [False]', required=False)
    
    return parser.parse_args()

#in command line: sequences, marker, contig/f/r, fname, rname, expName
def main():
    args= parse_cli_arguments()
    LOG.debug(args)

    if args.verbose:
        for handler in logging.getLogger().handlers:
            if type(handler) is logging.FileHandler:
                handler.setLevel(logging.DEBUG)  # Ensure file logs DEBUG too
            if type(handler) is logging.StreamHandler:
                handler.setLevel(logging.DEBUG)  # Optional: console to DEBUG 
    else:
        for handler in logging.getLogger().handlers:
            if type(handler) is logging.FileHandler:
                handler.setLevel(logging.DEBUG)  # Set file handler to DEBUG
            elif type(handler) is logging.StreamHandler:
                handler.setLevel(logging.INFO)   # Set console handler to INFO

    LOG.info("Running cryptogenotyper v{}".format(__version__))
    LOG.info(args)

    if utilities.is_databases_initialized() == False:
        utilities.init_blast_databases()
    else:
        LOG.info("Databases already initialized and in good health. Continue on")    

    if args.databasefile:
        make_custom_database(input_fasta=args.databasefile)  

    seq_dir = args.input[0]

    marker = args.marker
    if marker not in definitions.MARKERS:
        msg=f"marker value provided {marker} is not supported (supported markers {definitions.MARKERS})"
        LOG.error(msg)
        raise ValueError(msg)
    
    typeSeq = args.seqtype
    expName = args.outputprefix
    header = args.header

    pathlist=list()
    if os.path.isdir(seq_dir):
        seq_dir = os.path.abspath(seq_dir) #get directory name
        for file in os.listdir(seq_dir): #filter by extension type
            full_path = os.path.join(seq_dir, file)
            if os.path.isfile(full_path) and any(file.endswith(file_ext) for file_ext in definitions.FILETYPES):
                pathlist.append(full_path)
        pathlist = utilities.filter_files_by_suffix(pathlist,args.suffix)        
        #try to pair any ab1 files
        if typeSeq == "contig":
            file_pairs, other_files = utilities.pair_files(pathlist, args.forwardprimername, args.reverseprimername)
            flat_file_list = [file for pair in file_pairs for file in pair] #Now `flat_list` is a single list of all file paths
            pathlist = flat_file_list     
    else:
        #get absolute path of a single file provided
        if any([seq_dir.endswith(file_ext) for file_ext in definitions.FILETYPES]):
            pathlist.append(os.path.abspath(seq_dir)) 
        pathlist = utilities.filter_files_by_suffix(pathlist,args.suffix)    
    
    
    if not pathlist:
        raise ValueError("No input files found to process. Please check your input path or suffix filters.")


    #check if files exists actually
    for path in pathlist:
        if os.path.exists(path) == False:
            msg=f"File does not exist {path}. Check input file path"
            raise Exception(msg)    

    fPrimer=""; rPrimer=""
    if args.forwardprimername:
        fPrimer = args.forwardprimername
    if args.reverseprimername:
        rPrimer = args.reverseprimername

    if typeSeq == "contig":
        if marker == "18S":
            return msr_main(pathlist, fPrimer, rPrimer, typeSeq, expName, args.databasefile, header, args.verbose)

        elif marker == "gp60":
            return gp60_main(pathlist, fPrimer, rPrimer, typeSeq, expName, args.databasefile, header, args.verbose)

    elif typeSeq == "forward":
        #fPrimer = args.forwardprimername
        if marker == "18S":
            return msr_main(pathlist, fPrimer, rPrimer, typeSeq, expName, args.databasefile, header, args.verbose)

        elif marker == "gp60":
            return gp60_main(pathlist, fPrimer, rPrimer, typeSeq, expName, args.databasefile, header, args.verbose)


    elif typeSeq == "reverse":
        #rPrimer = args.reverseprimername
        if marker == "18S":
            return msr_main(pathlist, '', rPrimer, typeSeq, expName, args.databasefile, header, args.verbose)

        elif marker == "gp60":
            return gp60_main(pathlist, '', rPrimer, typeSeq, expName, args.databasefile, header, args.verbose)




if __name__ == "__main__":
    main()
