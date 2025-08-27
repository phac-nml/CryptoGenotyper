import os, shutil, glob, itertools, subprocess, datetime, logging, re, time, itertools
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Blast import NCBIXML
from CryptoGenotyper import definitions
from collections import defaultdict

# ANSI escape codes for text formatting
RED = "\033[31m"
BOLD = "\033[1m"
RESET = "\033[0m"


# setup the application logging
LOG = logging.getLogger(__name__)
DATABASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),"reference_database"))

def createTempFastaFiles(experiment_prefix="", record=None):
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
def sort_blast_hits_by_id_and_bitscore(blast_record, mode="default"):
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

        LOG.debug(f"HSP score threshold 95th quantile {int(threshold_score)}")
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
        
        # Sort kept alignments using composite sort key by %identity, score, reference allele coverage
        if mode=="default":
            sorted_blast_hits = sorted(
                kept_alignments,
                key=lambda a: (
                    a.hsps[0].identities / a.hsps[0].align_length,
                    a.hsps[0].score,
                    a.hsps[0].align_length / a.length #coverage of the reference allele
                ),
                reverse=True
                )
        elif mode == "bitscore":
            sorted_blast_hits = sorted(
                kept_alignments,
                key=lambda a: (
                    a.hsps[0].score,
                    a.hsps[0].align_length / a.length #length of the subject (reference allele)
                ),
                reverse=True
                )
        else:
            sorted_blast_hits = kept_alignments    

        # Log first 10 hits
        table_lines = [
            f"{'Idx':<3} │ {'ID':<50} │ {'Score':>6} │ {'Ident %':>8} │ {'Ref Cov %':>9} │ {'Gaps':>4} │ {'QLen':>6} │ {'ALen':>6} │ {'Qstart':>7} │ {'Qend':>5}",
            "-" * 132
        ]
        for idx, alignment in enumerate(sorted_blast_hits[:100]):
                hsp = alignment.hsps[0]
                identity = hsp.identities / hsp.align_length * 100
                coverage = min(alignment.hsps[0].align_length / alignment.length * 100, 100)
                table_lines.append(
                f"{idx + 1:<3} │ "
                f"{alignment.hit_id:<50} │ "
                f"{hsp.score:6.1f} │ "
                f"{identity:8.2f} │ "
                f"{coverage:9.2f} │ "
                f"{hsp.gaps:4d} │ "
                f"{blast_record.query_length:6d} │ "
                f"{hsp.align_length:6d} │ "
                f"{hsp.query_start:7d} │ "
                f"{hsp.query_end:5d}"
                )
        LOG.debug(f"Top 10 BLAST hits for {blast_record.query} ranked based on percent identity, then bitscore and then reference allele coverage:\n" + "\n".join(table_lines))        
        return sorted_blast_hits 
    else:
        return blast_record.alignments

            
def pair_files(file_paths, forward_suffix: str = "", reverse_suffix: str = ""):
    """
    Pairs sequencing files into forward/reverse read tuples.

    Operates in two modes:
    1. Direct File Mode: If `forward_suffix` and `reverse_suffix` are complete
       filenames of existing files in `file_paths`, they are paired directly.
    2. Suffix Mode: Otherwise, the function treats them as patterns and searches
       for them within the filenames to group and pair files by a common sample ID.

    :param file_paths: A list of file paths to be paired.
    :param forward_suffix: A full filename or a string pattern for a forward read.
    :param reverse_suffix: A full filename or a string pattern for a reverse read.
    :return: A tuple containing a list of (forward, reverse) pairs and a list of unpaired files.
    """

    def strip_suffix(filename: str, suffix: str) -> str:
        """Removes the given suffix and preceding text from a filename."""
        match = re.search(re.escape(suffix), filename)
        if match:
            common_name_before = filename[:match.start()]
            common_name_after = filename[match.end():]
            if common_name_before:
                return common_name_before
            if common_name_after:
                return common_name_after
            else:
                filename
        return filename

    # --- Step 0: Check for Direct File Mode ---
    # This handles the case where the user provides full filenames instead of suffixes.
    f_basename=""; r_basename=""
    if forward_suffix:
        f_basename = os.path.basename(forward_suffix)
    if reverse_suffix:  
        r_basename = os.path.basename(reverse_suffix)
    # Create a map of basenames to full paths for efficient lookup
    path_map = {os.path.basename(p): p for p in file_paths}

    if f_basename in path_map and r_basename in path_map:
        LOG.info("Direct file mode activated: Provided forward and reverse names are existing files.")
        
        f_path = path_map[f_basename]
        r_path = path_map[r_basename]
        
        file_pairs = [(f_path, r_path)]
        unpaired_files = [p for p in file_paths if p != f_path and p != r_path]
        
        return file_pairs, unpaired_files

    # --- If not in Direct File Mode, proceed with Suffix Mode ---
    LOG.info("Suffix mode activated: Searching for files containing specified suffix patterns.")

    # --- Step 1: Categorize all files in a single pass ---
    forward_reads = []
    reverse_reads = []
    other_files = []
    for path in file_paths:
        name_no_ext = os.path.splitext(os.path.basename(path))[0]
        if forward_suffix in name_no_ext:
            forward_reads.append(path)
        elif reverse_suffix in name_no_ext:
            reverse_reads.append(path)
        else:
            other_files.append(path)

    # --- Step 2: Check for the special case (single F/R pair found via suffix) ---
    if len(forward_reads) == 1 and len(reverse_reads) == 1:
        LOG.info(
            "Special suffix case: Found exactly one forward and one reverse file via suffix. "
            "Pairing them directly."
        )
        file_pairs = [(forward_reads[0], reverse_reads[0])]
    
        return file_pairs, other_files

    # --- Step 3: Proceed with standard sample_id based pairing ---
    grouped_files = defaultdict(dict)
    unpaired_files = list(other_files)

    for path in sorted(forward_reads):
        name_no_ext = os.path.splitext(os.path.basename(path))[0]
        
        sample_id = strip_suffix(name_no_ext, forward_suffix)
        
        if 'F' in grouped_files[sample_id]:
            LOG.warning(f"Duplicate forward read for sample ID '{sample_id}'. Keeping first, discarding {os.path.basename(path)}")
            unpaired_files.append(path)
        else:
            grouped_files[sample_id]['F'] = path
    
    for path in sorted(reverse_reads):
        name_no_ext = os.path.splitext(os.path.basename(path))[0]
        sample_id = strip_suffix(name_no_ext, reverse_suffix)
        if 'R' in grouped_files[sample_id]:
            LOG.warning(f"Duplicate reverse read for sample ID '{sample_id}'. Keeping first, discarding {os.path.basename(path)}")
            unpaired_files.append(path)
        else:
            grouped_files[sample_id]['R'] = path
          

    # --- Step 4: Create pairs and identify remaining unpaired files ---
    file_pairs = []
    for sample_id, files in sorted(grouped_files.items()):
        f_path = files.get('F')
        r_path = files.get('R')

        if f_path and r_path:
            file_pairs.append((f_path, r_path))
        else:
            if f_path: unpaired_files.append(f_path)
            if r_path: unpaired_files.append(r_path)
            LOG.warning(f"Incomplete pair for sample ID '{sample_id}': Found files {list(files.keys())} but not a complete F/R set.")

    # --- Step 5: Log results and return ---
    if file_pairs:
        table_lines = [
        f"\n{'Idx':<3} │ {'Forward Read':<50} │ {'Reverse Read':<50}", "─" * 115]
        for idx, (fwd, rev) in enumerate(file_pairs, 1):
            table_lines.append(f"{idx:<3} │ {os.path.basename(fwd):<50} │ {os.path.basename(rev):<50}")  
        # Print the table
        LOG.info("\n".join(table_lines))   

    # Print unpaired files as a simple list
    # Print unpaired files as a simple list with highlighting
    # Print unpaired files as a simple list with highlighting
    if unpaired_files:
        highlighted_unpaired_lines = []
        sorted_unpaired = sorted(unpaired_files) 
        
        for i, file_path in enumerate(sorted_unpaired):
            display_name = os.path.basename(file_path)
            
            highlighted_display_name = display_name # Start with the original name

            # Only highlight if the filename contains either the forward or reverse suffix.
            # This covers duplicates and files with missing partners.
            if forward_suffix in display_name:
                highlighted_display_name = re.sub(re.escape(forward_suffix), f"{RED}{BOLD}{forward_suffix}{RESET}", display_name, count=1)
            elif reverse_suffix in display_name:
                highlighted_display_name = re.sub(re.escape(reverse_suffix), f"{RED}{BOLD}{reverse_suffix}{RESET}", display_name, count=1)
            # Files in 'other_files' will not contain these suffixes and thus won't be highlighted,
            # which is the desired behavior for "pattern fails the pattern matches"
            
            highlighted_unpaired_lines.append(f"{i+1:>2}. {highlighted_display_name}")

        LOG.info(
            f"\nUnpaired {len(unpaired_files)} files:\n" +
            "\n".join(highlighted_unpaired_lines) +
            "\n"
        )
    LOG.info("Will start analyzing files in 10 seconds ...")    
    time.sleep(10)
    return sorted(file_pairs), sorted(unpaired_files)

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


def checkInputOrientation(sequence,  blastdbpath):
    LOG.info(f"Checking input sequence orientation")
    blast_output_filename = "blastn_orientation_results_tmp.xml"
    query_filename = "query.txt"

        # Open the file with writing permission
    myfile = open(query_filename, 'w')
    myfile.write('>query\n')    
    myfile.write(''.join(sequence))
    myfile.close()
 
    blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn',query="query.txt", dust='yes',
                                                    db=blastdbpath, reward=1, 
                                                    penalty=-2, gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, 
                                                    out=blast_output_filename)

    LOG.info(f"Running BLASTN to find if input path is Forward or Reverse orientation")
    stdout, stderr = blastn_cline()
    if stdout:
        LOG.debug(f"BLASTN stdout: {stdout}")
    if stderr:
        LOG.warning(f"BLASTN stderr: {stderr}")

    if (os.stat(blast_output_filename).st_size == 0):
        LOG.warning(f"No Significant Hits Found")
        return "No Significant Hits Found"

    result_handle = open(blast_output_filename , 'r')
    blast_records = NCBIXML.parse(result_handle)

    total_alignments = sum(len(blast_record.alignments) for blast_record in NCBIXML.parse(result_handle))

    LOG.info(f"Found {total_alignments} alignments")
    if total_alignments == 0:
        LOG.warning("No BLAST records found in the XML output.")
        return ""
    

    
    

    # 4. Parse the BLAST XML results
    forward_frame_matches = 0
    reverse_frame_matches = 0

    with open(blast_output_filename, 'r') as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        
        try:
            blast_record = next(blast_records) # Get the first (and usually only) query record
        except StopIteration:
            LOG.warning("No BLAST records found in the XML output. No significant hits")
            return ""
        
       

        # Iterate through alignments and HSPs to check frames
        # Iterate through alignments and HSPs to check strands
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                    # hsp.strand is a tuple (query_strand, subject_strand) e.g. ('Plus', 'Minus')
                    # For BLASTN, Plus indicates forward, Minus indicates reverse.
                    if hsp.strand[0] == "Plus" and hsp.strand[1] == "Plus":
                        forward_frame_matches += 1
                    elif hsp.strand[0] == "Plus" and hsp.strand[1] == "Minus":
                        reverse_frame_matches += 1
                    # Other strand combinations (e.g., query_strand -1) might occur
                    # if the query itself is being treated as reverse complemented by BLAST,
                    # but for typical query orientation checks, we assume the input query
                    # is "forward" and check subject alignment relative to it.

        # Clean up temporary files
        if os.path.exists(query_filename):
            os.remove(query_filename)
            LOG.debug(f"Cleaned up {query_filename}")
        if os.path.exists(blast_output_filename):
            os.remove(blast_output_filename)
            LOG.debug(f"Cleaned up {blast_output_filename}")
        
        # 5. Determine the overall orientation
        if forward_frame_matches > reverse_frame_matches:
            LOG.info(f"Determined orientation: Forward (Forward matches: {forward_frame_matches}, Reverse matches: {reverse_frame_matches})")
            return "Forward"
        elif reverse_frame_matches > forward_frame_matches:
            LOG.info(f"Determined orientation: Reverse (Forward matches: {forward_frame_matches}, Reverse matches: {reverse_frame_matches})")
            return "Reverse"
        elif forward_frame_matches > 0 and forward_frame_matches == reverse_frame_matches:
            LOG.info(f"Determined orientation: Ambiguous (Equal Forward and Reverse matches: {forward_frame_matches})")
            return ""
        else:
            LOG.info("No clear orientation could be determined from significant BLASTN hits.")
            return ""

def add_check_manually_str(string_list):
    msg_check = "Check manually."
   
    if string_list == []:
        return "-"
    if msg_check not in string_list:
        string_list.append(msg_check)
    return " ".join(string_list)

def add_to_string_if_missing(self, str2add):
    if str2add not in self.qcMsgsStr:
        self.qcMsgsStr += f"{str2add}"

def SSU_final_result_qc_checker(self, species, percent_identity,query_coverage,query_length,subject_length): 
    LOG.debug(f"Analyzing {self.name} file by final result QC checker species={species}, %id={percent_identity}, %query_cov={query_coverage} {query_length} {subject_length} ...")
    is_good_result = True; #qc_failure_reasons = []
    #need to have extra QC for C.hominis as it is the most represented species and similar to C.parvum
    if "C.hominis" in species:
        if percent_identity < 95:
            msg=f'Identity < 95% ({percent_identity}%).'
            LOG.debug(msg)
            add_to_string_if_missing(self, msg)
            is_good_result=False
            #qc_failure_reasons.append("Percent identity is less than 95.")
        elif query_coverage < 90:
            msg=f'Coverage < 90% ({query_coverage}%)'
            LOG.debug(msg)
            add_to_string_if_missing(self, msg)
            is_good_result=False
            #qc_failure_reasons.append("Percent query coverage is less than 60%.")
        else:
            is_good_result=True
    else:
       if query_coverage < 60 and query_length <= subject_length*1.05: #add 5% to the subject length for gaps and alignments 
            msg=f'Coverage < 60% ({query_coverage}%)'
            LOG.debug(msg)
            add_to_string_if_missing(self, msg)
            is_good_result=False
       elif percent_identity < 90:
            msg=f'Identity < 90% ({percent_identity}%)'
            LOG.debug(msg)
            add_to_string_if_missing(self, msg)
            is_good_result=False   
                    
 
    return is_good_result

def slice_multifasta(typeSeq, fasta_paths, pathlist, expName):
    # fasta file slicing
    if typeSeq == "contig":
        step_size = 2
    else:
        step_size = 1
            
    for i in range(0, len(fasta_paths), step_size):
        if step_size == 1:
            fasta_path = fasta_paths[i]
            with open(fasta_path,"r") as handle:
                records=list(SeqIO.parse(handle, "fasta"))
                LOG.info(f"File {handle.name} has {len(records)} sequences")
                #this is multifasta format as multiple records in fasta file
                if len(records) > 1:
                    for fasta_record in records:
                        tempFastaFilePath = createTempFastaFiles(expName,fasta_record)
                        pathlist.append(tempFastaFilePath)
                    pathlist.remove(fasta_path) #remove original fasta path that was split       
        elif step_size == 2:
            forward_path = fasta_paths[i]
            reverse_path = fasta_paths[i + 1]
           
            with open(forward_path,"r") as handle_fwd, open(reverse_path,"r") as handle_rev:
                records_fwd=list(SeqIO.parse(handle_fwd, "fasta"))
                records_rev=list(SeqIO.parse(handle_rev, "fasta"))
                LOG.info(f"Forward file {handle_fwd.name} has {len(records_fwd)} sequences and reverse file {handle_rev.name} has {len(records_rev)} sequences")
                
                if len(records_fwd) > 1 and len(records_rev) > 1 :
                    for fasta_record_fwd, fasta_record_rev in zip(records_fwd, records_rev):
                        tempFastaFilePath = createTempFastaFiles(expName+"_1",fasta_record_fwd)
                        pathlist.append(tempFastaFilePath)
                        tempFastaFilePath = createTempFastaFiles(expName+"_2",fasta_record_rev)
                        pathlist.append(tempFastaFilePath)
                    pathlist.remove(forward_path) #remove original fasta path that was split
                    pathlist.remove(reverse_path) #remove original fasta path that was split
            if pathlist:
                width_col = 60
                LOG.info(f"{os.path.basename(handle_fwd.name)} and {os.path.basename(handle_rev.name)} ..." )
                table_lines = [
                    f"\n{'Idx':<3} │ {'Forward':<{width_col}} │ {'Reverse':<{width_col}}",
                    "─" * (width_col*2+3)
                ]
                for i in range(0, len(pathlist), 2):
                    if i + 1 < len(pathlist):
                        fwd = pathlist[i]
                        rev = pathlist[i+1]
                        # Corrected line: use {width_col} instead of width_col
                        table_lines.append(f"{int(i/2) + 1:<3} │ {os.path.basename(fwd):<{width_col}} │ {os.path.basename(rev):<{width_col}}")
                    else:
                        fwd = pathlist[i]
                        table_lines.append(f"{int(i/2) + 1:<3} │ {os.path.basename(fwd):<{width_col}} │ {'':<{width_col}}")
                LOG.info("\n".join(table_lines))     