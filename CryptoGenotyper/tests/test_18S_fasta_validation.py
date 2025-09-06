import os,sys, json
from CryptoGenotyper.typer import main as cryptogenotyper_main

TEST_DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),"Data"))

def read_report_file(report_file_path):
     with open(report_file_path) as outfp:
         lines = [line.rstrip() for line in outfp.readlines()]
         return lines

#test the main function
def test_validation_18S_dataset(input_fasta_file=os.path.join(TEST_DATA_DIR,"dataset_18S_validation.fasta")):
    results_dict={}
    mismatch_sampleid_list = [] #stores all mismatched samples for debug (if any)
    
    args = [
        "-i", input_fasta_file ,
        "-m", "18S",
        "-o", "test_validation_18S"

    ]

    sys.argv[1:] = args
    cryptogenotyper_main()


    lines=read_report_file(f"{TEST_DATA_DIR}/dataset_validation_18S_metadata.txt")
    for ln, line in enumerate(lines[1:]):
        line_list = line.split("\t")
        if len(line_list) != 2:
            raise Exception(f"Line {ln} did not have 2 fields but {len(line_list)}")
        results_dict[line_list[0]]={"predicted_species":"",
                                    "expected_species":line_list[1]}    

    lines=read_report_file("test_validation_18S_cryptogenotyper_report.txt")
    header_fields_list = lines[0].split("\t")
    sampleid_idx = [idx for idx, field in enumerate(header_fields_list) if field=="Sample Name"][0]
    species_idx = [idx for idx, field in enumerate(header_fields_list) if field=="Species"][0]

    
    #add CryptoGentoyper 18S results to dictionary and find how many matched and how many did not (if any)
    mismatch_sampleid_list = []
    for line in lines[1:]:
        line_list = line.split("\t")
        sample_id = line_list[sampleid_idx]
        #print(f"CGT predicted: {line_list[sampleid_idx]}\t{line_list[species_idx]}\t{line_list[subtype_idx]}")
        results_dict[sample_id]["predicted_species"] = line_list[species_idx]
        #print(results_dict[line_list[sampleid_idx]])
        if results_dict[sample_id]["predicted_species"] != results_dict[sample_id]["expected_species"]:
            mismatch_sampleid_list.append(sample_id)
    
    print(f"Overall results of 18S marker testing on {len(results_dict)} samples:\n{json.dumps(results_dict, indent=4, sort_keys=True)}\n")

    if len(mismatch_sampleid_list) != 0:
        print(f"Found {len(mismatch_sampleid_list)} samples that did not match expected 18S subtype values")
        for sample_id in mismatch_sampleid_list:
            print(sample_id, results_dict[sample_id])
    
    print(f"Matched {len(results_dict) - len(mismatch_sampleid_list)} out of {len(results_dict) } 18S marker samples in FASTA inputs mode") 
    assert len(mismatch_sampleid_list) == 0, f"Found {len(mismatch_sampleid_list)} mismatch(es) ({mismatch_sampleid_list}) between expected and predicted species values"

def test_both_sequences_in_forward_orientation_contig_mode(input_fasta_file=TEST_DATA_DIR):
    output_name="test_hominis_18S_both_forward"
    args = [
        "-i", input_fasta_file ,
        "-m", "18S",
        "-f", "Chominis_18S_f2",
        "-r", "Chominis_18S_f3",
        "-t", "contig",
        "-o", output_name
    ]
    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file(f"{output_name}_cryptogenotyper_report.txt")
    assert len(lines) == 2
    header_fields_list = lines[0].split("\t")
    sampleid_idx = [idx for idx, field in enumerate(header_fields_list) if field=="Sample Name"][0]
    line_list = lines[1].split("\t")
    sample_id = line_list[sampleid_idx]
    assert "Chominis_18S_f2" == sample_id
    assert "C.parvum" in lines[1].split("\t")
    assert "1746" in lines[1].split("\t")

def test_both_sequences_in_reverse_orientation_contig_mode(input_fasta_file=TEST_DATA_DIR):
    output_name="test_hominis_18S_both_reverse"
    args = [
        "-i", input_fasta_file ,
        "-m", "18S",
        "-f", "Chominis_18S_r2",
        "-r", "Chominis_18S_r3",
        "-t", "contig",
        "-o", output_name
    ]
    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file(f"{output_name}_cryptogenotyper_report.txt")
    assert len(lines) == 2
    header_fields_list = lines[0].split("\t")
    sampleid_idx = [idx for idx, field in enumerate(header_fields_list) if field=="Sample Name"][0]
    line_list = lines[1].split("\t")
    sample_id = line_list[sampleid_idx]
    assert "Chominis_18S_r2" == sample_id
    assert "C.parvum" in lines[1].split("\t")
    assert "1746" in lines[1].split("\t")

def test_typical_contig_mode(input_fasta_file=TEST_DATA_DIR):
    args = [
        "-i", input_fasta_file ,
        "-m", "18S",
        "-f", "Chominis_18S_f2.fasta",
        "-r", "Chominis_18S_r3.fasta",
        "-t", "contig",
        "-o", "Chominis_18S_typical_contig"
    ]
    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("Chominis_18S_typical_contig_cryptogenotyper_report.txt")
    secondrow = lines[1].split("\t")
    assert "C.parvum" in secondrow
    assert "1746" in secondrow


# test contig mode for multifasta file with several reads from Illumina MiSeq sequencer
def test_multifasta_contig_mode_illumina(input_fasta_file=TEST_DATA_DIR):
    args = [
        "-i", input_fasta_file,
        "-m", "18S",
        "-f", "test_illumina_18S_F",
        "-r", "test_illumina_18S_R",
        "-t", "contig",
        "-o", "test_illumina_18S_contig"
     
    ]
    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("test_illumina_18S_contig_cryptogenotyper_report.txt")
    for line in lines[1:]:
        assert 'C.hominis' in line 
        assert "L16997" in line or "GQ983348.1" in line

# test contig mode for multifasta file with several reads each from MinIon sequencer
def test_multifasta_contig_mode_nanopore(input_fasta_file=TEST_DATA_DIR):
    args = [
        "-i", input_fasta_file ,
        "-m", "18S",
        "-f", "test_nanopore_18S_forward.fasta",
        "-t", "forward",
        "-o", "test_nanopore_18S_forward"
     
    ]
    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("test_nanopore_18S_forward_cryptogenotyper_report.txt")
    for line in lines[1:]:
        assert 'C.parvum' in line 
        assert "KM012040.1" in line
    
def test_short_18S_with_bad_ref_coverage(input_fasta_file=TEST_DATA_DIR):

    args = [
        "-i", input_fasta_file ,
        "-m", "18S",
        "-f", "C.parvum_AF093492_short.fa",
        "-t", "forward",
        "-o", "Cparvum_AF093492_short_18S_forward"
     
    ]
    sys.argv[1:] = args
    cryptogenotyper_main()

    lines=read_report_file("Cparvum_AF093492_short_18S_forward_cryptogenotyper_report.txt")
    for line in lines[1:]:
        assert 'C.parvum' in line 
        assert "Reference allele coverage is < 60%" in line

