import os,sys, logging, pytest
from CryptoGenotyper.typer import main as cryptogenotyper_main


TEST_DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),"Data"))
EXAMPLE_DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),"..","example"))

def read_report_file(report_file_path):
     with open(report_file_path) as outfp:
         lines = [line.rstrip() for line in outfp.readlines()]
         return lines

#test the contig gp60 typing with a pair of sanger files
def test_default_database_gp60_contig_input(input_dir=EXAMPLE_DATA_DIR):

    args = [
        "-i", input_dir ,
        "-m", "gp60",
        "-t", "contig",
        "-f", 'P17705_gp60-Crypt14-1F',
        "-r", "P17705_gp60-Crypt14-1R",
        "-s", "P17705",
        "-o", "test"

    ]

    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("test_cryptogenotyper_report.txt")
 
    secondrow = lines[1].split("\t")
    assert 'AY262034' in secondrow, secondrow[::-1]
    assert 'C.parvum' in secondrow
    assert 'IIaA15G2R1' in secondrow
    
#test the contig gp60 typing with a pair of sanger files with a custom database
def test_custom_database_gp60_contig_input(input_dir=os.path.abspath(EXAMPLE_DATA_DIR)):
    args = [
        "-i", input_dir,
        "-m", "gp60",
        "-t", "contig",
        "-f", "P17705_gp60-Crypt14-1F",
        "-r", "P17705_gp60-Crypt14-1R",
        "-d", os.path.abspath(os.path.join(os.path.dirname(__file__),"..","reference_database","gp60_ref.fa")),
        "-o", "test"

    ]

    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("test_cryptogenotyper_report.txt")
    secondrow = lines[1].split("\t")
    assert any('Found double peaks in repeat region' in item for item in secondrow), f'Could not find the expected QC message'
    assert any('AY262034' in item for item in secondrow), "AY262034 not found in any field."
    assert any('AY262034' in item for item in secondrow), "AY262034 not found in any field."
    assert any('C.parvum' in item for item in secondrow), "C.parvum not found in any field."
    assert any('IIa' in item for item in secondrow), "IIa not found in any field."

# test contig mode with sanger file with default 18S database
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
    print(lines)
    secondrow = lines[1].split("\t")
    thirdrow = lines[2].split("\t")
    


    assert 'C.parvum' in secondrow
    assert 'KT948751.1' in secondrow, secondrow[::-1]
    assert 'KM012040.1'  in thirdrow, thirdrow[::-1]
    assert 'C.parvum' in thirdrow


# test contig mode for sanger files if user supplies a custom database 
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

    assert 'C.parvum' in secondrow
    assert 'KT948751.1' in secondrow, secondrow[::-1]
    assert 'KM012040.1'  in thirdrow, thirdrow[::-1]
    assert 'C.parvum' in thirdrow

# test sanger 18S file in forward orientation
def test_default_singlefile(input_dir=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "example"))):
   

    args = [
        "-i", os.path.join(input_dir,"P17705_Crypto16-2F-20170927_SSUF_G12_084.ab1"),
        "-m", "18S",
        "-t", "forward",
        "-f", "SSUF",
        "-o", "test"

    ]

    sys.argv[1:] = args

    cryptogenotyper_main()

    lines=read_report_file("test_cryptogenotyper_report.txt")
    assert len(lines) == 3, f"Expected 3 lines in output but got {len(lines)}"
    
    secondrow = lines[1].split("\t")
    assert 'C.parvum' in secondrow
    assert 'KT948751.1' in secondrow
    thirdrow = lines[2].split("\t")
    assert 'KM012040.1' in thirdrow
   
            
 
# test if the same fasta file can be correctly 18S typed provided in the forward and reverse orientation
def test_gp60_fasta_single_sequence(input_dir = os.path.join(TEST_DATA_DIR) ):
    args = [
        ["-i", os.path.join(input_dir,"P17705_gp60-Crypt16-1F-20170927_gp60F_G08_052.fasta"),
        "-m", "gp60",
        "-t", "forward",
        "-o", "P17705_gp60_fasta"],
        [
        "-i", os.path.join(input_dir,"P17705_gp60-Crypt16-1F-20170927_gp60F_G08_052.fasta"),
        "-m", "gp60",
        "-t", "reverse",
        "-o", "P17705_gp60_fasta"
        ]
    ]

    for arg in args:
        sys.argv[1:] = arg
        cryptogenotyper_main()

        lines=read_report_file("P17705_gp60_fasta_cryptogenotyper_report.txt")  
        secondrow = lines[1].split("\t")

        assert "C.parvum" in secondrow
        assert "IIaA15G2R1" in secondrow



# test if the same fasta file can be correctly 18S typed provided in the forward and reverse orientation and as a contig
def test_18S_fasta_single_sequence(input_dir = os.path.join(TEST_DATA_DIR) ):
    args = [
        ["-i", os.path.join(input_dir,"P17705_Crypto16-20170927_SSUR.fasta"),
        "-t", "forward",
        "-m", "18S",
        "-o", "results_18S_fasta"],
        ["-i", os.path.join(input_dir,"P17705_Crypto16-20170927_SSUR.fasta"),
        "-t", "reverse",
        "-m", "18S",
        "-o", "results_18S_fasta"],
        ["-i", input_dir,
        "-t", "contig",
        "-f", "P17705_Crypto16-20170927_SSUF.fasta",
        "-r", "P17705_Crypto16-20170927_SSUR.fasta",
        "-m", "18S",
        "-o", "results_18S_fasta"], #contig for Crypto16 successfully formed
        ["-i", input_dir,
        "-t", "contig",
        "-f", "P16182_MD123e_SSUF.fasta",
        "-r", "P16182_MD123e_SSUR.fasta",
        "-m", "18S",
        "-o", "results_18S_fasta" #contig for MD123e successfully formed
        ]

    ]        

    for arg in args:
        sys.argv[1:] = arg
        print(f"*** Running {arg} ***")
        cryptogenotyper_main()
        lines=read_report_file("results_18S_fasta_cryptogenotyper_report.txt")  
        secondrow = lines[1].split("\t")

        assert "C.parvum" in secondrow or "C.hominis" in secondrow
        assert "KT948751.1" in secondrow or "KM012040.1" in secondrow or "GQ983348.1" in secondrow

# test if single read fasta file in reverse orientation correctly
def test_18S_fasta_single_sequence_сustom_database(input_dir = os.path.join(TEST_DATA_DIR) ):
    args = [
        "-i", os.path.join(input_dir,"P17705_Crypto16-20170927_SSUR.fasta"),
        "-m", "18S",
        "-d", os.path.abspath(os.path.join(os.path.dirname(__file__),"..","reference_database","msr_ref.fa")),
        "-o", "P17705_Crypto16_18S_fasta"
    ]        

    sys.argv[1:] = args
    print(os.getcwd())
    cryptogenotyper_main()
    lines=read_report_file("P17705_Crypto16_18S_fasta_cryptogenotyper_report.txt")  
    secondrow = lines[1].split("\t")

    assert "C.parvum" in secondrow

#test if multi-fasta file with several sequences can by typed well in forward mode
def test_18S_fasta_multi_sequence(input_dir = os.path.join(TEST_DATA_DIR) ):
    args = [
        "-i", os.path.join(input_dir,"18S_multifasta.fasta"),
        "-m", "18S",
        "-d", os.path.abspath(os.path.join(os.path.dirname(__file__),"..","reference_database","msr_ref.fa")),
        "-o", "18S_multifasta"
    ]        

    sys.argv[1:] = args
    cryptogenotyper_main()
    lines=read_report_file("18S_multifasta_cryptogenotyper_report.txt")   

    assert "C.parvum" in lines[1].split("\t")
    assert "C.hominis" in lines[2].split("\t") 

# Test sanger two files in a contig mode when contig is formed
def test_sanger_contig_mode_gp60_two_reverse(caplog, input_dir = TEST_DATA_DIR):
    caplog.set_level(logging.DEBUG)
    
    expected_message = "Detected 411bp overlap region between sequences"
    args = [
        "-i", input_dir,
        "-f", "P27277_Crypt242_Cmortiferum_R2.ab1",
        "-r", "P27277_Crypt242_Cmortiferum_R3.ab1",
        "-m", "gp60",
        "-t", "contig",
        "-o", "P17705_Crypto242_gp60_sanger_contig",
        "--verbose"
    ]  
    
    sys.argv[1:] = args   
    cryptogenotyper_main()  
  
    lines=read_report_file("P17705_Crypto242_gp60_sanger_contig_cryptogenotyper_report.txt")   
    assert "C.mortiferum" in lines[1].split("\t")
    assert "XIVaA18G2T2a" in lines[1].split("\t")
    assert expected_message in caplog.text , f"Expected message not found in stdout"

# Test if in contig mode user provides valid input files but of different types (Sanger and FASTA). 
# Raise error in that case 
def test_different_valid_file_types_contig(caplog, input_dir = TEST_DATA_DIR):
    caplog.set_level(logging.DEBUG)
    # 1. Create a symbolic link to a file
    source_file = f"{EXAMPLE_DATA_DIR}/P17705_Crypto16-2R-20170927_SSUR_H12_082.ab1"
    link_file = "P17705_Crypto16-2R-20170927_SSUR_H12_082.ab1"

    # os.symlink(source, destination)
    if os.path.exists(source_file) == False:
        os.symlink(source_file, link_file)
 
    args = [
        "-i", input_dir,
        "-f", "P17705_Crypto16-20170927_SSUF.fasta",
        "-r", "P17705_Crypto16-2R-20170927_SSUR_H12_082.ab1",
        "-m", "18S",
        "-t", "contig",
        "-o", "P17705_Crypto16_18S_sanger_contig",
        "--verbose"
    ]  
    sys.argv[1:] = args   
    with pytest.raises(ValueError) as excinfo:
        cryptogenotyper_main()
    expected_message = "No input files found to process. Please check your input path or suffix filters."
    assert expected_message in str(excinfo.value)     

# test using a distant actin sequence how 18S and gp60 will handle this scenario
def test_using_distant_actin_seq(caplog, input_dir = TEST_DATA_DIR):
    caplog.set_level(logging.DEBUG)
    args = [
        ["-i", os.path.join(input_dir,"XM_001388245.1_Cryptosporidium_parvum_actin.fa"),
        "-t", "forward",
        "-m", "gp60",
        "-o", "actin_fasta"],
        ["-i", os.path.join(input_dir,"XM_001388245.1_Cryptosporidium_parvum_actin.fa"),
        "-t", "forward",
        "-m", "18S",
        "-o", "actin_fasta"]
        
    ]    

    for arg in args:
        sys.argv[1:] = arg
        print(f"*** Running {arg} ***")
        cryptogenotyper_main()
        lines=read_report_file("actin_fasta_cryptogenotyper_report.txt")  
        secondrow = lines[1].split("\t")
        assert "Poor Sequence Quality. Check manually." in secondrow or \
            "Could not analyze. No species detected (potential reasons: not an 18S sequence, poor sequence quality, ref. database limitations, BLAST failure). Check manually." in secondrow
        assert "Maybe not be an 18S Crypto sequence or outdated database" in caplog.text  or \
            "Maybe not be a gp60 Crypto sequence or outdated database" in caplog.text
        break
    