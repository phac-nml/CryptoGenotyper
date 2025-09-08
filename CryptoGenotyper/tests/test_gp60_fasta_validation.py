import os,sys, json
from CryptoGenotyper.typer import main as cryptogenotyper_main

TEST_DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),"Data"))

def read_report_file(report_file_path):
     with open(report_file_path) as outfp:
         lines = [line.rstrip() for line in outfp.readlines()]
         return lines

#test the main function
def test_validation_gp60_dataset(input_fasta_file=os.path.join(TEST_DATA_DIR,"dataset_validation_gp60.fasta")):
    results_dict={}
    mismatch_sampleid_list = [] #stores all mismatched samples for debug (if any)
    
    args = [
        "-i", input_fasta_file ,
        "-m", "gp60",
        "-o", "test_validation_gp60"
    ]

    sys.argv[1:] = args
    cryptogenotyper_main()


    lines=read_report_file(f"{TEST_DATA_DIR}/dataset_validation_gp60_metadata.txt")
    for ln, line in enumerate(lines[1:]):
        line_list = line.split("\t")
        if len(line_list) != 3:
            print(f"Line {ln} did not have 3 fields but {len(line_list)}")
        results_dict[line_list[0]]={"predicted_species":"","predicted_subtype":"",
                                    "expected_species":line_list[1],"expected_subtype":line_list[2]}    
        

    lines=read_report_file("test_validation_gp60_cryptogenotyper_report.txt")
    header_fields_list = lines[0].split("\t")
    sampleid_idx = [idx for idx, field in enumerate(header_fields_list) if field=="Sample Name"][0]
    species_idx = [idx for idx, field in enumerate(header_fields_list) if field=="Species"][0]
    subtype_idx = [idx for idx, field in enumerate(header_fields_list) if field=="Subtype"][0]

    #add CryptoGentoyper results to dictionary
    for line in lines[1:]:
        line_list = line.split("\t")
        results_dict[line_list[sampleid_idx]]["predicted_species"] = line_list[species_idx]
        results_dict[line_list[sampleid_idx]]["predicted_subtype"] = line_list[subtype_idx]
        
    
    
    with open('test_validation_gp60_dataset_results_dict.json', 'w') as fp:
        json.dump(results_dict, fp, indent=4)

    for sample_id in results_dict.keys():
        if results_dict[sample_id]["predicted_subtype"] != results_dict[sample_id]["expected_subtype"]:
            mismatch_sampleid_list.append(sample_id)
            
    if len(mismatch_sampleid_list) != 0:
        print(f"Found {len(mismatch_sampleid_list)} samples that did not match expected gp60 subtype values")
        for sample_id in mismatch_sampleid_list:
            print(sample_id, results_dict[sample_id])
    
    print(f"Matched {len(results_dict) - len(mismatch_sampleid_list)} out of {len(results_dict) } gp60 marker samples in FASTA inputs mode")             
    
    assert len(mismatch_sampleid_list) == 0
        
#test short truncated gp60 sequence of only 86bp
def test_short_gp60_sequence(input_fasta_file=os.path.join(TEST_DATA_DIR,"OP925781.1_gp60_short.fasta")):
    args = [
        "-i", input_fasta_file ,
        "-m", "gp60",
        "-o", "short_sequence_validation_gp60"

    ]

    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("short_sequence_validation_gp60_cryptogenotyper_report.txt")
  
    
    secondrow = lines[1].split("\t")

    assert 'C.ryanae' in secondrow 
    assert 'XXIj' in secondrow

#test long 2863bp gp60 sequence 3x longer than normal size to make sure results are generated and QC warnings do not prevent results generation   
def test_long_gp60_sequence(input_fasta_file=os.path.join(TEST_DATA_DIR,"CM098023.1_gp60_long.fasta")):
    args = [
        "-i", input_fasta_file ,
        "-m", "gp60",
        "-o", "long_sequence_validation_gp60"

    ]

    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("long_sequence_validation_gp60_cryptogenotyper_report.txt")
    secondrow = lines[1].split("\t")

    assert 'C.hominis' in secondrow 
    assert 'IbA12G3' in secondrow

def test_very_shot_gp60_sequence_blast_identical_hits(input_fasta_file=os.path.join(TEST_DATA_DIR,"very_short_gp60_sequence_identical_blast_hits.fasta")):
    args = [
        "-i", input_fasta_file ,
        "-m", "gp60",
        "-o", "very_short_identical_blast_hits_validation_gp60"

    ]

    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("very_short_identical_blast_hits_validation_gp60_cryptogenotyper_report.txt")
  
 
    secondrow = lines[1].split("\t")
    assert 'C.hominis' in secondrow 
    assert any('Check manually' in s for s in secondrow), secondrow

#provide two sequences in reverse orientation in contig mode and test sequence orientation detection algorithm
def test_both_sequences_in_reverse_orientation_contig_mode(input_fasta_file=TEST_DATA_DIR):
    args = [
        "-i", input_fasta_file ,
        "-m", "gp60",
        "-f", "Cmortiferum_gp60_r2",
        "-r", "Cmortiferum_gp60_r3",
        "-t", "contig",
        "-o", "Cmortiferum_gp60_both_reverse"
    ]
    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("Cmortiferum_gp60_both_reverse_cryptogenotyper_report.txt")
    secondrow = lines[1].split("\t")
    assert 'C.mortiferum' in secondrow 
    assert "XIVaA16G2T2a" in secondrow

def test_both_sequences_in_forward_orientation_contig_mode(input_fasta_file=TEST_DATA_DIR):
    args = [
        "-i", input_fasta_file ,
        "-m", "gp60",
        "-f", "Cmortiferum_gp60_f2",
        "-r", "Cmortiferum_gp60_f3",
        "-t", "contig",
        "-o", "Cmortiferum_gp60_both_forward"
    ]
    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("Cmortiferum_gp60_both_forward_cryptogenotyper_report.txt")
    secondrow = lines[1].split("\t")
    assert 'C.mortiferum' in secondrow 
    assert "XIVaA16G2T2a" in secondrow

# test contig mode given a single read fasta files in correct orientation (forward and reverse)
def test_typical_fasta_single_read_contig_mode(input_fasta_file=TEST_DATA_DIR):
    args = [
        "-i", input_fasta_file ,
        "-m", "gp60",
        "-f", "Cmortiferum_gp60_f2",
        "-r", "Cmortiferum_gp60_r3",
        "-t", "contig",
        "-o", "Cmortiferum_gp60_typical_contig"

    ]
    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("Cmortiferum_gp60_typical_contig_cryptogenotyper_report.txt")
    secondrow = lines[1].split("\t")
    assert 'C.mortiferum' in secondrow 
    assert "XIVaA16G2T2a" in secondrow

# test contig mode for multifasta file with 3 reads each from Illumina MiSeq
def test_multifasta_contig_mode_illumina(input_fasta_file=TEST_DATA_DIR):
    args = [
        "-i", input_fasta_file ,
        "-m", "gp60",
        "-f", "test_illumina_gp60_F",
        "-r", "test_illumina_gp60_R",
        "-t", "contig",
        "-o", "test_illumina_gp60_contig"
     
    ]
    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("test_illumina_gp60_contig_cryptogenotyper_report.txt")
    secondrow = lines[1].split("\t")
    thirdrow = lines[2].split("\t")
    assert 'C.parvum' in secondrow 
    assert "IIaA16G3R1" in secondrow
    assert 'C.parvum' in thirdrow 
    assert "IIaA15G2R2" in thirdrow

def test_short_gp60_with_bad_ref_coverage(input_fasta_file=TEST_DATA_DIR):
    args = [
        "-i", input_fasta_file ,
        "-m", "gp60",
        "-f", "Cparvum_IIa_gp60_AY262034.fa",
        "-t", "forward",
        "-o", "Cparvum_IIa_gp60_AY262034_short_gp60_forward"
     
    ]
    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("Cparvum_IIa_gp60_AY262034_short_gp60_forward_cryptogenotyper_report.txt")
    for line in lines[1:]:
        assert 'C.parvum' in line 
        assert "IIaA15G2R1" in line
        assert "Reference allele coverage is < 60%" in line
