import os,sys
from CryptoGenotyper.typer import main as cryptogenotyper_main

def read_report_file(report_file_path):
     with open(report_file_path) as outfp:
         lines = [line.rstrip() for line in outfp.readlines()]
         return lines

#test the main function
def test_default_database_gp60_contiginput(input_dir=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","example"))):
    print(os.getcwd(),os.listdir(input_dir))

    args = [
        "-i", input_dir ,
        "-m", "gp60",
        "-t", "contig",
        "-f", "gp60F",
        "-r", "gp60R",
        "-o", "test"

    ]

    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("test_cryptogenotyper_report.txt")
    secondrow = lines[1].split("\t")
    assert 'Could not classify repeat region. Check manually.' in secondrow, f'Could not find the expected QC message'
    assert 'AY262034' in secondrow
    assert 'C.parvum' in secondrow
    assert 'IIa' in secondrow

def test_custom_database_gp60_contig_input(input_dir=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","example"))):
    args = [
        "-i", input_dir,
        "-m", "gp60",
        "-t", "contig",
        "-f", "gp60F",
        "-r", "gp60R",
        "-d", os.path.abspath(os.path.join(os.path.dirname(__file__),"..","reference_database","gp60_ref.fa")),
        "-o", "test"

    ]

    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("test_cryptogenotyper_report.txt")
    secondrow = lines[1].split("\t")
    assert 'Could not classify repeat region. Check manually.' in secondrow, f'Could not find the expected QC message'
    assert 'AY262034' in secondrow
    assert 'C.parvum' in secondrow
    assert 'IIa' in secondrow

def test_default_database_18S(input_dir=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "example"))):
    args = [
        "-i", input_dir,
        "-m", "18S",
        "-t", "contig",
        "-f", "P17705_Crypto16-2F",
        "-r", "P17705_Crypto16-2R",
        "-o", "test"

    ]

    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("test_cryptogenotyper_report.txt")
    
    secondrow = lines[1].split("\t")
    thirdrow = lines[2].split("\t")
    forthrow = lines[3].split("\t")

    assert 'KT948751.1' in secondrow 
    assert 'C.parvum' in secondrow
    assert 'KM012040.1' in thirdrow
    assert 'C.parvum' in secondrow
    assert 'L16997'  in forthrow
    assert 'C.hominis' in forthrow



def test_custom_database_18S(input_dir=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","example"))):

    args = [
        "-i", input_dir,
        "-m", "18S",
        "-t", "contig",
        "-f", "P17705_Crypto16-2F",
        "-r", "P17705_Crypto16-2R",
        "-d", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "reference_database", "msr_ref.fa")),
        "-o", "test"

    ]

    sys.argv[1:] = args
    cryptogenotyper_main()
    lines=read_report_file("test_cryptogenotyper_report.txt")
    
    secondrow = lines[1].split("\t")
    thirdrow = lines[2].split("\t")
    forthrow = lines[3].split("\t")

    assert 'KT948751.1' in secondrow 
    assert 'C.parvum' in secondrow
    assert 'KM012040.1' in thirdrow
    assert 'C.parvum' in secondrow
    assert 'L16997'  in forthrow
    assert 'C.hominis' in forthrow



def test_default_singlefile(input_dir=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "example"))):
    os.chdir(input_dir)

    args = [
        "-i", "P17705_Crypto16-2F-20170927_SSUF_G12_084.ab1",
        "-m", "18S",
        "-t", "forward",
        "-f", "SSUF",
        "-r", "SSUR",
        "-o", "test"

    ]

    sys.argv[1:] = args

    cryptogenotyper_main()

    lines=read_report_file("test_cryptogenotyper_report.txt")
    
    secondrow = lines[1].split("\t")
    thirdrow = lines[2].split("\t")
    forthrow = lines[3].split("\t")

    assert 'KT948751.1' in secondrow 
    assert 'C.parvum' in secondrow
    assert 'KM012040.1' in thirdrow
    assert 'C.parvum' in secondrow
    assert 'L16997'  in forthrow
    assert 'C.hominis' in forthrow

