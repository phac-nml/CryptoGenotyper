import os,sys
from CryptoGenotyper.typer import main as cryptogenotyper_main



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

def test_custom_database_gp60_contiginput(input_dir=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","example"))):
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

def test_custom_database_18S(input_dir=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","example"))):

    args = [
        "-i", input_dir,
        "-m", "18S",
        "-t", "contig",
        "-f", "SSUF",
        "-r", "SSUR",
        "-d", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "reference_database", "msr_ref.fa")),
        "-o", "test"

    ]

    sys.argv[1:] = args
    cryptogenotyper_main()

def test_default_database_18S(input_dir=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "example"))):
    args = [
        "-i", input_dir,
        "-m", "18S",
        "-t", "contig",
        "-f", "SSUF",
        "-r", "SSUR",
        "-o", "test"

    ]

    sys.argv[1:] = args
    cryptogenotyper_main()

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

