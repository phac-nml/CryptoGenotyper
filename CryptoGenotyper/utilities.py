import os, shutil, glob, itertools, subprocess, datetime, logging, re
from CryptoGenotyper import definitions
from collections import defaultdict

# setup the application logging
LOG = logging.getLogger(__name__)
DATABASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),"reference_database"))

def createTempFastaFiles(experiment_prefix="", record=None):
    LOG.debug("Creating temporary FASTA files from multi-FASTA file ...")
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
    LOG.debug("Cleaning temporary files after the run")
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



def quantile(data, q):
    """
    Calculate the q-th quantile of a list of numbers.
    
    Parameters:
    - data: list of numeric values
    - q: quantile to compute (e.g., 0.25 for Q1, 0.5 for median, 0.75 for Q3)

    Returns:
    - The quantile value
    """
    if not 0 <= q <= 1:
        raise ValueError("q must be between 0 and 1")

    sorted_data = sorted(data)
    n = len(sorted_data)
    pos = q * (n - 1)
    lower = int(pos)
    upper = min(lower + 1, n - 1)
    fraction = pos - lower

    return sorted_data[lower] * (1 - fraction) + sorted_data[upper] * fraction

#sort BLAST hits based on %identity first and if a tie, 
#then by the bitscore (default BLAST only sorts by the bitscore) and reference coverage     
def sort_blast_hits_by_id_and_bitscore(blast_record):
    if len(blast_record.alignments) > 1:
        # Calculate coverage for each alignment
        coverages = {a: a.hsps[0].align_length / blast_record.query_length for a in blast_record.alignments}
        identities = {a: a.hsps[0].identities / a.hsps[0].align_length for a in blast_record.alignments}
        scores = {a: a.hsps[0].score for a in blast_record.alignments}
        
        # Compute %identity and %coverage (reference allele)
        coverage_values = list(coverages.values())
        identity_values = list(identities.values())
        threshold_coverage = quantile(coverage_values,0.75)
        threshold_identity = quantile(identity_values,0.25)
        threshold_score = quantile(scores.values(),0.95)

        LOG.debug(f"HSP score threshold 95th quantile {threshold_score}")
        # Filter alignments (only alignment objects)
        kept_alignments = [a for a in blast_record.alignments if scores[a] >= threshold_score]
        # Fallback: if no alignments pass the threshold, keep the top-scoring alignment
        if not kept_alignments:
            top_alignment = max(blast_record.alignments, key=lambda a: a.hsps[0].score)
            kept_alignments = [top_alignment]
            LOG.debug(f"No alignments passed threshold — using the top scoring alignment: {top_alignment.hit_id} with score {top_alignment.hsps[0].score}")
        filtered_alignments = [a for a in blast_record.alignments if coverages[a] < threshold_coverage and identities[a] < threshold_identity ]

        # Log filtered hits info for debug
        if filtered_alignments:
            LOG.debug(f"Filtered out {len(filtered_alignments)} alignments:{[a.hit_id for a in filtered_alignments]}")
        LOG.debug(f"Kept {len(kept_alignments)} alignments")    
        
        # Sort kept alignments using composite sort key by %identity, score, reference allele
        sorted_blast_hits = sorted(
            kept_alignments,
            key=lambda a: (
                a.hsps[0].identities / a.hsps[0].align_length,
                a.hsps[0].score,
                a.hsps[0].align_length / a.length
            ),
            reverse=True
            )
        # Log first 10 hits
        table_lines = [
            f"{'Idx':<3} │ {'ID':<50} │ {'Score':>6} │ {'Ident %':>8} │ {'Ref Cov %':>8} │ {'Gaps':>4} │ {'QLen':>6} │ {'ALen':>6}",
            "-" * 115
        ]
        for idx, alignment in enumerate(sorted_blast_hits[:100]):
                hsp = alignment.hsps[0]
                identity = hsp.identities / hsp.align_length * 100
                coverage = min(hsp.align_length / blast_record.query_length * 100, 100)
                table_lines.append(
                f"{idx + 1:<3} │ "
                f"{alignment.hit_id:<50} │ "
                f"{hsp.score:6.1f} │ "
                f"{identity:8.2f} │ "
                f"{coverage:9.2f} │ "
                f"{hsp.gaps:4d} │ "
                f"{blast_record.query_length:6d} │ "
                f"{hsp.align_length:6d}"
                )
        LOG.debug("Top 10 BLAST hits ranked based on perecent identity, then bitscore and then reference allele coverage:\n" + "\n".join(table_lines))        
        return sorted_blast_hits 
    else:
        return blast_record.alignments

            

def pair_files(file_paths, forward_suffix, reverse_suffix):
    """
    Pairs sequencing files into forward and reverse read tuples based on specified suffix patterns.

    This function processes a list of AB1 file paths and groups them into forward and reverse read pairs
    based on the provided suffix strings. Each file name is examined for the presence of either the forward
    or reverse suffix. If a match is found, the suffix is used to identify and group files by sample ID 
    (by stripping characters starting from the suffix match). Files with matching sample IDs are then paired.

    :param file_paths: List[str]
        A list of file paths to files to be paired.
    
    :param forward_suffix: str
        A string pattern indicating that a file is a forward read (e.g., 'SSUF_' or '_F1_') specified by --forwardprimername.

    :param reverse_suffix: str
        A string pattern indicating that a file is a reverse read (e.g., 'SSUR_' or '_R2_') specified by --reverseprimername.

    :return: List[Tuple[str, str]]
        A list of tuples each consisting of two values
        Each tuple contains the file paths for a forward and reverse read pair (paired_files list) and occasionally
        Files that cannot be paired are returned as a separate unpaired_files lies .
    """

    def strip_suffix(filename, suffix):
        best_match = None
        best_span = (0, 0)

        #for suffix in suffixes:
        match = re.search(re.escape(suffix), filename)
        if match:
            if (match.end() - match.start()) > (best_span[1] - best_span[0]):
                best_match = match
                best_span = match.span()

            if best_match:
                return filename[:best_span[0]]
        return filename
    
    grouped_files = defaultdict(dict)
    unpaired_files = []
    
    for path in sorted(file_paths):
        filename = os.path.basename(path) #use only filenames but not paths 
        name_no_ext = os.path.splitext(filename)[0] #remove file extension from the file name

        if forward_suffix in name_no_ext:
            sample_id = strip_suffix(name_no_ext, forward_suffix)
            grouped_files[sample_id]['F'] = path
        elif reverse_suffix in name_no_ext:
            sample_id = strip_suffix(name_no_ext, reverse_suffix)
            grouped_files[sample_id]['R'] = path
        else:
            unpaired_files.append(path)
            LOG.warning(f"Could not identify direction for file: {filename}")

    file_pairs = []
    for sample_id, files in sorted(grouped_files.items()):
        f = files.get('F')
        r = files.get('R')
        if f and r:
            file_pairs.append((f, r))
        else:
            LOG.warning(f"Incomplete pair for '{sample_id}': {files}")
    
    
    if file_pairs:
        table_lines = [
        f"\n{'Idx':<3} │ {'Forward Read':<50} │ {'Reverse Read':<50}", "─" * 115]
        for idx, (fwd, rev) in enumerate(file_pairs, 1):
            table_lines.append(f"{idx:<3} │ {os.path.basename(fwd):<50} │ {os.path.basename(rev):<50}")  
        # Print the table
        LOG.info("\n".join(table_lines))   

    # Print unpaired files as a simple list
    if unpaired_files:
        LOG.info(
        f"\nUnpaired {len(unpaired_files)} files:\n" +
        "\n".join(f"{i+1:>2}. {file}" for i, file in enumerate(sorted(unpaired_files))) +
        "\n"
    )
    return  file_pairs, unpaired_files


def filter_files_by_suffix(pathlist, suffix):
    """
    Filters a list of file paths based on a provided suffix.

    :param pathlist: List of file paths
    :param suffix: Suffix string to filter files by (if None or empty, no filtering is applied)
    :return: Filtered list of file paths
    """
    if suffix:
        original_count = len(pathlist)
        pathlist = [path for path in pathlist if suffix in os.path.basename(path)]
        LOG.info(f"Filtered {original_count} files using suffix '{suffix}' and {len(pathlist)} files remaining after filtering.")
    
    return pathlist