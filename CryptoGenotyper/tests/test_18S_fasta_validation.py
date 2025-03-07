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
            print(f"Line {ln} did not have 2 fields but {len(line_list)}")
        results_dict[line_list[0]]={"predicted_species":"",
                                    "expected_species":line_list[1]}    

    lines=read_report_file("test_validation_18S_cryptogenotyper_report.txt")
    header_fields_list = lines[0].split("\t")
    sampleid_idx = [idx for idx, field in enumerate(header_fields_list) if field=="Sample Name"][0]
    species_idx = [idx for idx, field in enumerate(header_fields_list) if field=="Species"][0]

    
    #add CryptoGentoyper 18S results to dictionary
    for line in lines[1:]:
        line_list = line.split("\t")
        #print(f"CGT predicted: {line_list[sampleid_idx]}\t{line_list[species_idx]}\t{line_list[subtype_idx]}")
        results_dict[line_list[sampleid_idx]]["predicted_species"] = line_list[species_idx]
        #print(results_dict[line_list[sampleid_idx]])
    print(results_dict)