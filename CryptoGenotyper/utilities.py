import os, shutil, glob, itertools, subprocess, datetime, logging, statistics
from CryptoGenotyper import definitions
from CryptoGenotyper.logging import create_logger

# setup the application logging
LOG = create_logger(__name__)
DATABASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),"reference_database"))

def createTempFastaFiles(experiment_prefix="", record=None):
    LOG.info("Creating temporary FASTA files from multi-FASTA file ...")
    temp_dir = os.path.join("./",f"tmp_fasta_files_{experiment_prefix}")
    if os.path.exists(temp_dir == False):
        os.makedirs(temp_dir, exist_ok=True)
    tempFastaFilePath=os.path.join(temp_dir,record.name+".fasta") 
    with open(tempFastaFilePath, "w") as fp:   
        fp.write(record.format("fasta"))    
    return tempFastaFilePath

def getFileType(path):
    filetype = [filetype for filetype in definitions.FILETYPES if path.endswith(filetype) ]
    if filetype and len(filetype) == 1:
        filetype = filetype[0]
        if 'ab1' == filetype:
            filetype = "abi"
        return filetype 
    else:
        return None
    
def cleanTempFastaFilesDir(temp_dir="tmp_fasta_files"):
    LOG.info("Cleaning temporary files after the run")
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir, ignore_errors=True)
    tmpfiles2remove = list(itertools.chain.from_iterable([glob.glob(e) for e in 
                                                          ["align.*","custom_db.*", "query*.txt", "result*.txt", "SSUresult*",
                                                          "result*.xml", "refseq.fa","blast*.xml", "*_tmp.fasta", "*_tmp.xml"]]))
    for file in tmpfiles2remove:
        try:
            os.remove(file)
        except FileNotFoundError:
            pass       

def init_blast_databases():
    LOG.info("Initializing databases ...") 
    files_list = os.listdir(DATABASE_DIR)
    fasta_files = [file for file in files_list if file.endswith(tuple(definitions.FASTA_FILETYPES))]
    tmpfiles2remove = [file for file in files_list if file.endswith((".ndb",".nhr",".nin",".nog",".nos",".not",".nsq",".ntf",".nto"))]

    for file in tmpfiles2remove:
        try:
            os.remove(file)
        except FileNotFoundError:
            pass       

    for fasta in fasta_files:
        subprocess.run(f"cd {DATABASE_DIR} && makeblastdb -dbtype nucl -in {fasta}  -out {fasta}  -parse_seqids", shell=True)

    with open(os.path.join(DATABASE_DIR,"db_status.txt"),"w") as fp:
        data_and_time_str = datetime.datetime.today().strftime('%Y-%m-%d')
        fp.write(f"Databases initialized on {data_and_time_str}")
        LOG.info("Databases initialized successfully")


def is_databases_initialized():
     db_status_filepath = os.path.join(DATABASE_DIR,"db_status.txt")
     if os.path.exists(db_status_filepath) and os.stat(db_status_filepath).st_size != 0:
        return True
     else:
        return False 
#sort BLAST hits based on %identity first and if a tie, 
#then by the bitscore (default BLAST only sorts by the bitscore) and reference coverage     
def sort_blast_hits_by_id_and_bitscore(blast_record):
    if len(blast_record.alignments) > 1:
        # Calculate coverage for each alignment
        coverages = {a: a.hsps[0].align_length / a.length for a in blast_record.alignments}

        # Compute mean and standard deviation
        coverage_values = list(coverages.values())
        mean_cov = statistics.mean(coverage_values)
        std_dev_cov = statistics.stdev(coverage_values)
        threshold = mean_cov - 2 * std_dev_cov
        coverages = {a: a.hsps[0].align_length / a.length for a in blast_record.alignments}
        average_query_coverage = statistics.mean(coverages.values())

        LOG.debug(f"Average query coverage for BLAST hits is {average_query_coverage:.4f} and coverage threshold is {threshold}")

        # Filter alignments (only alignment objects)
        kept_alignments = [a for a in blast_record.alignments if coverages[a] >= threshold]
        filtered_alignments = [a for a in blast_record.alignments if coverages[a] < threshold]

        # Log filtered hits info for debug
        if filtered_alignments:
            LOG.debug("Filtered out alignments:")
            for a in filtered_alignments:
                LOG.debug(f"  {a.hit_id}: coverage={coverages[a]:.4f}")
        
         # Sort kept alignments using composite sort key
        return sorted(
            kept_alignments,
            key=lambda a: (
                a.hsps[0].identities / a.hsps[0].align_length,
                a.hsps[0].score,
                a.hsps[0].align_length / a.length
            ),
            reverse=True
        )
    else:
        return blast_record.alignments

            

