#******************************************************************************
#******************************************************************************
#                 Cryptosporidium Gp60 Analyzer (Version 1.1)
#             Written by: Christine Yanta (christine.yanta@canada.ca)
#                           Date: May 2, 2017
#   Public Health Agency of Canada - National Microbiology Laboratory Guelph
#******************************************************************************
#******************************************************************************


import io
import os, sys
import re
import logging
from CryptoGenotyper.logging import create_logger
from CryptoGenotyper import definitions, utilities

from Bio import SeqIO
from Bio.Seq import Seq

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio.Blast import NCBIXML

from Bio.Alphabet import IUPAC

import numpy as np

import copy
import math
import shutil, glob, itertools, collections



global TESTING
TESTING = False


# setup the application logging
LOG = logging.getLogger(__name__)

#analyzingGp60 Class Description:
#   This class keeps all the information of each Cryptosporidium sample's
#   gp60 Sanger sequence and analyzes it to classify.
#   The data is obtained from reading the given ab1 file and is then aligned
#   with the reference sequences to determine family and subtype.  It then
#   counts the number of repeats within the sequence, depending on the
#   family determined.
class analyzingGp60(object):

    def __init__(self):
        self.g = []
        self.a = []
        self.t = []
        self.c = []

        self.name=""  #name of the file
        self.species=""
        self.repeats=""

        self.oldseq = []  #before trimming the ends of the sequence
        self.seq = []    #after trimming the ends of the sequence
        self.seqLength = 0

        self.peakLoc = [] #peak locations
        self.phred_qual = [] #phred quality: 10 = 90% base call accuracy, 20=99% (log scale)

        self.beginSeq = 0
        self.endSeq = 0

        self.repeatStarts=0
        self.repeatEnds=0

        self.forwardSeq = True
        self.file = None
        self.tabfile = None

        self.checkRepeatManually = False
        self.averagePhredQuality = 0
        self.doublePeaksinRepeat = False
        self.ambigSpeciesBlast = [] #multiple BLAST hits with identical ranking making typing ambiguous
        self.fileType = "" # filetype (sanger or fasta) that determines QC messages displayed and final report





    def readFiles(self, dataFile, forw, output, tabFile, filetype = "abi", customdatabasename = None):
        #setting whether the sequence is the forward or reverse strand
        #(used later on know whether the reverse complement is needed)
        if forw:
            self.forwardSeq = True #forward input file
        else:
            self.forwardSeq = False #reverse input file

        self.file = output
        self.tabfile = tabFile
        

        #**********Beginning to work with the ab1 file given:**********#

        #retrieves the sample file name (removes directory pathway)
        self.name = dataFile.split("/")[len(dataFile.split("/"))-1]

        LOG.debug(f"Reading file {self.name} ...") #Lets user know which sequence the program is on

        #opens the ab1 file
        if  filetype in definitions.SANGER_FILETYPES:
            handle=open(dataFile,"rb")
            record=SeqIO.read(handle, "abi")
            #retrieving base amplitude data
            self.g=np.array(record.annotations['abif_raw']['DATA9'])
            self.a=np.array(record.annotations['abif_raw']['DATA10'])
            self.t=np.array(record.annotations['abif_raw']['DATA11'])
            self.c=np.array(record.annotations['abif_raw']['DATA12'])

            #obtaining primary bases, includes N's
            raw_seq = list(record.annotations['abif_raw']['PBAS2'].decode('UTF-8'))

            #peak locations
            self.peakLoc=np.array(record.annotations['abif_raw']['PLOC2'])

            #phred quality of the bases
            self.phred_qual=record.letter_annotations['phred_quality']
            self.seq = raw_seq
        elif filetype in definitions.FASTA_FILETYPES:
            handle=open(dataFile,"r")   
            record=SeqIO.read(handle, "fasta")
            raw_seq = list(record.seq.upper()) #making sure all bases converted to upper case so matching works
            self.seq = raw_seq
            self.phred_qual = [60] * len(raw_seq)
        else:
            sys.exit(f"Unsupported file type: '{filetype}'. Please provide a valid file. Supported extensions are: {', '.join(definitions.FILETYPES)}")

        # check input sequence orientation
        sequence_orientation_check_result  = utilities.checkInputOrientation(self.seq, os.path.dirname(__file__)+"/reference_database/gp60_ref.fa")
     
        if sequence_orientation_check_result:
            LOG.info(f"Input {self.name} sequence orientation is {sequence_orientation_check_result}")
            if sequence_orientation_check_result  == "Reverse":
                self.forwardSeq = False 
            elif sequence_orientation_check_result  == "Forward":
                self.forwardSeq = True
        else:
            LOG.info(f"Could not verify orientation from BLAST searches for {self.name}")  
        
        self.seqLength = len(self.seq)


        self.repeatStarts = 0
        self.repeatEnds = 0

        self.b = 0
        self.e = 0

        #calculate inital average quality taking into account entire sequence and later update value with selected region (if not failing condition is hit)
        self.averagePhredQuality = round(sum(self.phred_qual) / len(self.phred_qual),2) 

        

        


        #if any([isinstance(element,str) == False for element in self.seq]):
        #    self.seq = list(record.annotations['abif_raw']['PBAS2'].decode('UTF-8'))
            #self.oldseq = list(record.annotations['abif_raw']['PBAS2'].decode('UTF-8'))
       

        b = 0
        e = 0

        filename = f"query_vs_blast_gp60.txt"
        myfile = open(filename, 'w')

        LOG.debug(f"Initial {len(self.seq)}bp long input sequence from {self.name}:"+''.join(self.seq)) #debug
        myfile.write(f">{self.name}\n")
        myfile.write(''.join(self.seq))

        # Close the file
        myfile.close()
        
        
        dbpath = os.path.join(os.path.dirname(__file__),"reference_database/blast_gp60.fasta")
        blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn',query=filename, dust='yes',
                                             db= dbpath,
                                             reward=1, penalty=-2,gapopen=5, gapextend=2, 
                                             evalue=0.00001, outfmt=5, out="gp60result2.xml")
        input_seq_length = len(self.seq)
        LOG.info(f"Querying {self.name} seq of {input_seq_length}bp against global {os.path.basename(blastn_cline.db)} gp60 database ...")
        stdout, stderr = blastn_cline()
        LOG.debug(f"BLASTN stdout={stdout} and stderr={stderr}")
      
        

        if (os.stat("gp60result2.xml").st_size == 0):
            LOG.error(f"Generated an empty gp60 BLAST result for {self.name}. Maybe not be a gp60 Crypto sequence or outdated database at {dbpath}?")
            return False
        
        result_handle = open("gp60result2.xml", 'r')
        blast_records = NCBIXML.parse(result_handle)
        blast_record = next(blast_records)
        if len(blast_record.alignments) == 0:
            LOG.warning(f"Found 0 alignments for {self.name}. Maybe not be a gp60 Crypto sequence or outdated database at {dbpath}?")
            return False
        blast_record.alignments = utilities.sort_blast_hits_by_id_and_bitscore(blast_record, "bitscore") 
        br_alignment = blast_record.alignments[0]
        hsp = br_alignment.hsps[0]
        self.species = br_alignment.hit_id
       


        # use full sequence length in FASTA format (no sequence trimming) except if sequence is larger then 2000 bp 
        # E.g. longest gp60 sequences such as C.fayeri|IVc(FJ490069), Length: 1505, Sequence ID: C.suis|XXVa(MH187874), Length: 1338)
        # trim input seqeunce to the top reference in case user supplies entire genome or large sequence that could negatively impact repeat region determination
        if filetype in definitions.FASTA_FILETYPES:
            if self.seqLength <= 2000:
                b=0
                e=self.seqLength
            else: #sequence larger then 2000 bp then it can be trimmed if a top reference is > 700 bp, otherwise do nothing
                if abs(e-b) > 2000:
                    b = hsp.query_start-1 #trim index
                    e = hsp.query_end     #trim index
                else:
                    b=0
                    e = self.seqLength-1

        else:
            #works even if sequence is in reverse orientation as results are always given with respect to the query sequence so
            #trimming coordinates are always valid
            dbpath = os.path.join(os.path.dirname(__file__),"reference_database/gp60_ref.fa")
            LOG.info(f"Querying {self.name} against the {dbpath} gp60 database ...")

            filename = f"query.txt"

            # Open the file with writing permission
            myfile = open(filename, 'w')

            t = round(len(self.seq)*0.5)
            myfile.write(f">{self.name}\n")
            myfile.write(''.join(self.seq[0:t]))

            # Close the file
            myfile.close()

            if customdatabasename:
                    dbpath = "custom_db"
            else:
                    dbpath = os.path.dirname(__file__)+"/reference_database/gp60_ref.fa"
            blastn_cline = NcbiblastnCommandline(cmd ='blastn', task='blastn',query=filename, dust='yes',
                                                     db = dbpath, 
                                                     reward=1, penalty=-2,gapopen=5, gapextend=2,evalue=1e-5, outfmt=5, out="gp60result.xml")
            LOG.info(f"Querying {self.name} against the {blastn_cline.db} gp60 database ...")
            stdout, stderr = blastn_cline()
            
            b = hsp.query_start-1 #trim index
            e = hsp.query_end #trim index
 
            if b > (len(self.seq)*0.25):
                LOG.debug(f"Position b={b} and is GREATER than 25% of the original sequence length ({len(self.seq)}). Threshold: {len(self.seq)*0.25:.2f}")
                
                if (os.stat("gp60result.xml").st_size !=0):
                    result_handle = open("gp60result.xml", 'r')
                    blast_records = NCBIXML.parse(result_handle)
                    blast_record = next(blast_records)

                    if len(blast_record.alignments) == 0:
                        return False

                    br_alignment = blast_record.alignments[0]
                    hsp = br_alignment.hsps[0]

                    self.species = br_alignment.hit_id
                    b = hsp.query_start-1
                    e = hsp.query_end
            else:
                LOG.debug(f"Position b={b} and is LESS than 25% of the original sequence length ({len(self.seq)}). Threshold: {len(self.seq)*0.25:.2f}")        
                if (os.stat("gp60result.xml").st_size !=0):
                    result_handle = open("gp60result.xml", 'r')
                    blast_records = NCBIXML.parse(result_handle)
                    blast_record = next(blast_records)
                    if len(blast_record.alignments) == 0:
                        return False

                    br_alignment = blast_record.alignments[0]
                    self.species = br_alignment.hit_id
                  
        #Checking quality of the sequence.  If the avg. phred quality < 20 (99% base
        #   calling certainty), then the sequence cannot be analyzed.
        
        quality=0
        full_quality = 0
  
        for i in range(b ,e):
            quality += self.phred_qual[i]

        for i in range(0, self.seqLength):
            full_quality += self.phred_qual[i]

        full_quality = full_quality/self.seqLength
        

        quality = quality / len(self.phred_qual[b:e])
        self.averagePhredQuality = round(quality,2)
        LOG.info(f"Average PHRED {self.name} sequence quality in selected range ({b} - {e}) is {self.averagePhredQuality}")

        if quality < 20 or full_quality < 10:
            #self.file.write("\n")
            #self.file.write("\n;>Sequence: " + self.name.split(".ab1")[0])
            #self.file.write("\n;>\tPoor sequence data.  Could not analyze.  Refer to sequence manually.")

            #self.tabfile.write(self.name.split(".ab1")[0] + "\t-\t-\t-\t" + "Poor sequence quality. Check manually.\t" + str(self.averagePhredQuality) + "\n")

            LOG.warning(f"\tPoor sequence data for {self.name}. Could not analyze. Refer to sequence manually.")

            return False
        else:
            self.peakLoc = copy.copy(self.peakLoc[b:e])
            self.seq = copy.copy(self.seq[b:e])
            self.phred_qual = copy.copy(self.phred_qual[b:e])
            self.seqLength = len(self.seq)
            

            self.b = b
            self.e = e



            LOG.info(f"Final sequence length {self.seqLength}bp (coordinates {b}:{e}) selected from raw {len(raw_seq)} bp sequence ({round(self.seqLength/len(raw_seq)*100,1)}% left)")
            #calling trimSeq to trim the ends of the Sanger sequence
            #goodTrim = self.trimSeq()
            
            return True







    #trimSeq() trims both ends of the Sanger sequence since the sequencing is
    #   never clear (usually messy) at the beginning and the end.  Trims based
    #   on phred_quality (10 bases in a row with 99% base calling accuracy)
    def trimSeq(self, filetype="ab1"):
        if filetype in definitions.FASTA_FILETYPES:
            self.seq = list(self.seq)
            self.findRepeatRegion()
            return True
        
        LOG.info(f"Trimming sequence {self.name} ({len(self.seq)}bp) of filetype {filetype}")
        qual_cutOff = 20 #99% certainty
        foundBegin = False
        foundEnd = False

        #obtain non-ambiguous sequence by fixing the N's
        self.fixN()
        #finding the repeat region
        self.findRepeatRegion()
     

        goodRepeatFind = True

        for i in range(self.repeatEnds, len(self.seq)):
            if i == len(self.seq) -1:
                goodRepeatFind = True
            elif abs(self.peakLoc[i] - self.peakLoc[i+1]) <10:
                goodRepeatFind = False
                break

        #Starting to trim the start of the sequence
        #   Looking for where the phred_quality is greater the cutoff (99% accuracy) for
        #   10 bases in a row.
        for x in range (0, self.seqLength):
            #if the start of the repeat region, stop trimming the sequence
            #  update where the repeat region starts
            if x == self.repeatStarts:
                if x == 0:
                    self.beginSeq = 0
                    break

                for i in range (self.repeatStarts-1, 0, -1):
                    if self.phred_qual[i] < qual_cutOff:
                        self.beginSeq = i
                        break

                break

            #if phred_quality is greater than 99% base accuracy and not an N
            if self.phred_qual[x] >= qual_cutOff:
                for i in range(0,5):
                    if self.phred_qual[x+i] < qual_cutOff:
                        foundBegin=False
                        break
                    else:
                        foundBegin = True

            #if it passed 10+ bases at 99%, start trim there
            if foundBegin:
                self.beginSeq = x
                break

        #reverse input logic        
        if self.forwardSeq == False:
            averageA = 0
            averageC = 0
            averageG = 0
            averageT = 0
            for x in range(self.repeatStarts, self.repeatEnds):
                averageA += self.a[self.peakLoc[x]]
                averageC += self.c[self.peakLoc[x]]
                averageT += self.t[self.peakLoc[x]]
                averageG += self.g[self.peakLoc[x]]
            
            denom = self.repeatEnds - self.repeatStarts
            averageA = (averageA / denom) *0.75
            averageC = (averageC / denom) *0.75
            averageG = (averageG / denom) *0.75
            averageT = (averageT / denom) * 0.75

            if self.repeatEnds == len(self.seq):
                self.endSeq = self.repeatEnds

            for x in range(self.repeatEnds, self.seqLength, 1):
                #print(x, self.a[self.peakLoc[x+i]], self.c[self.peakLoc[x+i]], self.g[self.peakLoc[x+i]],self.t[self.peakLoc[x+i]])
                if x == self.seqLength-1 or x+5 >= self.seqLength-1:
                    self.endSeq = x
                    break
                if self.phred_qual[x] < qual_cutOff:
                    for i in range(0, 5, 1):
                        if self.phred_qual[x+i] >= qual_cutOff:

                            if self.a[self.peakLoc[x+i]] >= averageA or  self.g[self.peakLoc[x+i]] >= averageG or self.c[self.peakLoc[x+i]] >= averageC or self.t[self.peakLoc[x+i]] >= averageT:
                                foundEnd=False
                                break
                            else:
                                foundEnd=True
                        else:
                            foundEnd=True
                    if foundEnd:
                        self.endSeq = x
                        break

        #forward input logic               
        else:
            if self.repeatEnds == len(self.seq):
                self.endSeq = self.repeatEnds
            #Starting to trim the end of the sequence
            #  Same algorithm as the beginning, but starting from the end
            for x in range (self.seqLength-1, self.beginSeq, -1):
                #if at the end of the repeat region, stop cutting there
                if x == self.repeatEnds:
                    self.endSeq = x
                    break

                #if phred quality is higher than the quality cut off, must be when 10+bases
                if self.phred_qual[x] >= qual_cutOff:
                    for i in range (0, -10, -1):
                        if self.phred_qual[x+i] < qual_cutOff:
                            foundEnd=False
                            break

                        else:
                            foundEnd = True

                #if found where to cut off
                if foundEnd:
                    self.endSeq = x
                    break

        if abs(self.endSeq-self.repeatEnds) <= 2 and goodRepeatFind==False:
            goodRepeatFind = False
        else:
            goodRepeatFind = True


        cutSeqLength = self.endSeq - self.beginSeq
        self.repeatStarts = self.repeatStarts - self.beginSeq
        self.repeatEnds = self.repeatEnds - self.beginSeq

        


        tempSeq = ['']*(cutSeqLength)
        tempPeakLoc = ['']*(cutSeqLength)
        tempPhredQual = ['']*(cutSeqLength)

        #shift everything in sequence over to trim the sequence
        for t in range (self.beginSeq, self.endSeq):
            tempSeq[t-self.beginSeq] = self.seq[t]
            tempPeakLoc[t-self.beginSeq] = self.peakLoc[t]
            tempPhredQual[t-self.beginSeq] = self.phred_qual[t]

        length = len(tempSeq)


        tempA = ['']*length
        tempG = ['']*length
        tempC = ['']*length
        tempT = ['']*length

        #copy over the amplitude values
        for i in range(0, length):
            #print(self.peakLoc[self.beginSeq+i])
            tempA[i] = copy.copy(self.a[self.peakLoc[self.beginSeq+i]])
            tempG[i] = copy.copy(self.g[self.peakLoc[self.beginSeq+i]])
            tempC[i] = copy.copy(self.c[self.peakLoc[self.beginSeq+i]])
            tempT[i] = copy.copy(self.t[self.peakLoc[self.beginSeq+i]])



        LOG.debug(f"Cut new sequence length {cutSeqLength}bp  cut from position {self.beginSeq} and ending at {self.endSeq} removing {len(self.seq)-cutSeqLength} bp")
      
        #copy all temporary variables for analysis and store for that class
        self.seq = copy.copy(tempSeq)
        self.peakLoc = copy.copy(tempPeakLoc)
        self.phred_qual = copy.copy(tempPhredQual)
        self.a = copy.copy(tempA)
        self.g = copy.copy(tempG)
        self.c = copy.copy(tempC)
        self.t = copy.copy(tempT)
        self.seqLength = cutSeqLength


        #if the reverse sequence, make the reverse complement
        if not self.forwardSeq:
            #self.peakLoc = self.peakLoc[::-1] #reverse the array
            self.phred_qual = self.phred_qual[::-1]

            a = self.a
            t = self.t
            g = self.g
            c = self.c

            self.a = t[::-1]
            self.t = a[::-1]
            self.g = c[::-1]
            self.c = g[::-1]

            revCseq=''.join(self.seq)

            temp = Seq(revCseq, IUPAC.ambiguous_dna)
            rev_complement = temp.reverse_complement()

            self.seq = copy.copy(rev_complement)

            #update where the repeat region starts as the sequence is being
            #  reversed

            start = len(self.seq) - self.repeatEnds
            end = len(self.seq) - self.repeatStarts

            
            self.repeatStarts = start
            self.repeatEnds = end
            self.forwardSeq = True





        #if forward sequence, just make the sequence in string form instead of array
        else:
            forwSeq = ''.join(self.seq)
            temp = Seq(forwSeq, IUPAC.ambiguous_dna)
            self.seq = copy.copy(temp)
        

        return goodRepeatFind




    #fixN() finds the best base when an 'N' is encountered in the sequence by finding
    #   the maximum amplitude of all four bases at that position
    def fixN(self):
        l = len(self.seq)
        for i in range (0, l):
            base=self.seq[i]


            #if the base is not one of the 4 DNA base codes
            if (base!='A' and base!='G' and base!='C' and base!='T'):
                
                index = self.peakLoc[i] #finding peak location

                #get amplitudes for all four bases and find the max
                amplitudes=[self.a[index], self.g[index], self.c[index], self.t[index]]
                maxAmp = max(amplitudes)

                if self.a[index] == maxAmp:
                    self.seq[i] = 'A'

                elif self.g[index] == maxAmp:
                    self.seq[i] = 'G'

                elif self.c[index] == maxAmp:
                    self.seq[i] = 'C'

                elif self.t[index] == maxAmp:
                    self.seq[i] = 'T'

            #print(base, self.seq[i])



    #*************FINDING POSSIBLE REPEAT REGION***************************
    #finding TCA region algorithm:
    #   Looks through entire sequence for the repeat region, always keeping
    #   track of where the longest TC_ (or AG_ for reverse) repeat is
    def findRepeatRegion(self):
        LOG.info(f"Trying to find repeat region in the sequence {len(self.seq)}bp ...")
      
        seq = list(self.seq)
    


        repeatStart = 0
        maxCount = 0 #length of the repeat region

        tempCounter = 0
        inRepeat = False
        
        repeatForw = ['T','C']
        repeatRev = ['G', 'A']
        ACA_case_inCparvum = ['A','C','A']+ ['T','C','A']*2 #see section 2.3 The r repeat designation in IIa of "Deciphering a cryptic minefield: a guide to Cryptosporidium gp60 subtyping"

        index=0
        
        #matched_positions = [match for match in re.finditer(r"TC\w{1,1}","".join(self.seq))]
        
        #print([match for match in re.finditer(r"TC\w{1,1}",self.seq)])
        #repeats_dict = {}; inRepeat = False; repeatsCounter=0; 
        #for match in matched_positions:
        #    if match.span()[0]+3 in [m.span()[0] for m in matched_positions]:
        #        if inRepeat == False:
        #            repeatsCounter += 1
        #            repeats_dict[repeatsCounter]=[]
        #        inRepeat = True
        #        repeats_dict[repeatsCounter].append(match.span()[0])
        #    else:
        #        #not in repeat region or repeat region has just finished
        #        if inRepeat == True:
        #            repeats_dict[repeatsCounter].append(match.span()[0]+3)    #add the end position of the repeat region running in triplets
        #        inRepeat=False
        #print(repeats_dict)        
        #longest_repeat_codons_pos = repeats_dict[sorted(repeats_dict, key=lambda key: len(repeats_dict[key]),reverse=True)[0]]     
        #repeatStart = min(longest_repeat_codons_pos)
        #repeatEnd = max(longest_repeat_codons_pos)
        #repeatSeq = self.seq[repeatStart:repeatEnd]
        
        #LOG.debug(f"Found longest repeat {repeatSeq} of {len(repeatSeq)}bp (positions {repeatStart}-{repeatEnd}) composed of {len(longest_repeat_codons_pos)} triplets")        
        #LOG.debug(f"Matched positions binned into repeats dictionary {repeats_dict}")
        #self.repeatStarts = repeatStart
        #self.repeatEnds = repeatEnd
        #return  True
        

        # OLD CODE from v1.0.0
        #loops through entire sequence once
        while True:
            #once entire sequence has been looked at, exit loop
            if index >= self.seqLength:
                #if the last part of the sequence was a repeat region, update stats
                #if appropriate
                if inRepeat and tempCounter > maxCount:
                    maxCount = tempCounter
                    repeatStart = index - (maxCount*3)
                break


            #if its a TC_ or AG_ part of the sequence (in a potential repeat region)
            if (self.forwardSeq and seq[index:index+2]==repeatForw) or \
                (not self.forwardSeq and seq[index:index+2]==repeatRev) or ("C.parvum" in self.species and seq[index:index+9] == ACA_case_inCparvum):
                #print(f"if {index+1} and inRepeat {inRepeat} trying to find repeat and next 2 codons {seq[index+3:index+9]}")
                #if found previous consecutive TC_ and AG_, keep adding onto how
                #many there is and go to the next 3 bases (next maybe repeat)
                if inRepeat:
                    tempCounter += 1
                    index += 3

                #start of a new potential repeat region
                else:
                    tempCounter = 1
                    index+=3
                    inRepeat=True

                
            #if it's the end of the repeat region        
            elif inRepeat:
                inRepeat = False
               
                candidate_repeat = "".join(seq[index - (tempCounter*3):index - (tempCounter*3) + (tempCounter*3)])
                before_repeat_2codons = "".join(seq[index - (tempCounter*3)-6:index - (tempCounter*3)])
                 
                #C.viatorum family XVa has GAAAGC sequence before repeat region
                if re.match(r"C.viatorum\|XV[a-z]{1,1}",self.species) and before_repeat_2codons == "GAAAGC":
                    maxCount = tempCounter
                    repeatStart  = index - (maxCount*3)
                    break  

                #if that was the longest region found, keep track of where it was and report it
                if tempCounter > maxCount:
                    maxCount = tempCounter
                    repeatStart  = index - (maxCount*3)
                    
                   
                

                index+=1

            #if not in a repeat region
            else:
                index+=1
        
        self.repeatStarts = repeatStart

        #can not exceed the length of the entire sequence
        if repeatStart + (maxCount*3) > self.seqLength:
            self.repeatEnds = self.seqLength
        else:    
            self.repeatEnds = repeatStart + (maxCount*3)


        #if it's the reverse sequence, want to subtract one since the reverse sequence
        #  reads _GA and looking for the 'GA' portions (meaning you want to start at the
        #  _ part)
        if not self.forwardSeq:
            self.repeatStarts -= 1
            self.repeatEnds -= 1
        LOG.debug(f"Repeat starts at position {self.repeatStarts} and ends at {self.repeatEnds} {''.join(seq[self.repeatStarts:self.repeatEnds])}")    







    #determineRepeats() finds the number of repeats found in the gp60 sequence
    #   Looks for TCA, TCG, TCT, and R (defined below).  Will need to modify
    #   this function if different repeat regions are required to be found.
    def determineRepeats(self):
        LOG.info(f"Trying to determine the repeat region ({self.repeatStarts}-{self.repeatEnds}) and encode it in the gp60 standard nomenclature for {self.name}")
        self.seq = "".join(self.seq)

        self.repeatLen = len(self.seq[self.repeatStarts:self.repeatEnds])
        seq = copy.copy(self.seq)
        

        #Seeing if good enough quality of sequence to determine repeats
        quality_cutOff=15
        goodQuality=False
        self.checkRepeatManually = False
  
        for i in range(self.repeatStarts, self.repeatEnds):
            if i < len(self.a):
                a= self.a[i]
                c= self.c[i]
                g= self.g[i]
                t= self.t[i]
            else:
                a=0
                c=0
                g=0
                t=0
            
            maxAmp = max(a,c,g,t)
            lri = 3.1

            if a != 0 and a != maxAmp:
                    lri = math.log((maxAmp/a),2)
            if g != 0 and g != maxAmp:
                    lri = math.log((maxAmp/g),2)
            if c != 0 and c!=maxAmp:
                    lri = math.log((maxAmp/c),2)
            if t != 0 and t!=maxAmp:
                    lri = math.log((maxAmp/t),2)
            
            if lri <= 2.0:
                self.doublePeaksinRepeat = True
 
            if self.phred_qual[i] < quality_cutOff:  
                self.checkRepeatManually = True
                for x in range(0, 3):
                    if i+x >= self.repeatEnds:
                        break
                    elif i+x >= len(self.phred_qual):
                        break
                    elif self.phred_qual[i+x] < quality_cutOff and (i+x+1) < len(self.peakLoc):
                        if abs(self.peakLoc[i+x] - self.peakLoc[i+x+1]) < 10:
                            goodQuality = True
                        else:
                            goodQuality = False


                    else:
                        goodQuality = True
                        break


            else:
                goodQuality=True

            if not goodQuality:
                break
            
        if goodQuality:
            LOG.info(f"Sequence {self.name} of good quality to detect repeats")
            numAcatca = 0   #C. parvum R
            numACA = 0 # for C. cuniculus Vb that has ACA repeats
            numHomR = 0     #C. hominis R= A(A/G)(A/G)ACGGTGGTAAGG (15bp minisatellite)
            numHomR2 = 0    #C. hominis R=A(A/G)(A/G)ACGGTGAAGG (13bp minisatellite)
            numHomfR = 0    #C. hominis If R=AAGAAGGCAAAGAAG

            LOG.info(f"Determining the R repeat after position {self.repeatEnds} ...")
            #after repeat region sequence, looking for R
            afterRepeat = "".join(copy.copy(self.seq[self.repeatEnds:]))
            
            if (afterRepeat.find("ACATCA")!= -1) : #returns -1 if substring is not found
                numAcatca = afterRepeat.count("ACATCA")
                self.repeatLen = self.repeatLen+numAcatca*len("ACATCA")

            if re.search(r"(ACA){1,}",afterRepeat):
                numACA = re.search(r"(ACA){1,}",afterRepeat).group(0).count('ACA')    

            if len(self.species.split("|")) >= 2:
                numHomR = 0; numHomR2 = 0
                if "Ia" in self.species.split("|")[1]:
                    numHomR = afterRepeat.count("AAAACGGTGGTAAGG") + afterRepeat.count("AGAACGGTGGTAAGG") + afterRepeat.count("AAGACGGTGGTAAGG") + afterRepeat.count("AGGACGGTGGTAAGG")

                    #If 2 numHomR are found, sometimes the third is different (still counted as same R)
                    #   Reference:
                    if numHomR >= 2:
                        numHomR2 = afterRepeat.count("AAAACGGTGAAGG") + afterRepeat.count("AAGACGGTGAAGG") + afterRepeat.count("AGAACGGTGAAGG") + afterRepeat.count("AGGACGGTGAAGG")
                        numHomR += numHomR2
                    
                elif "If" in self.species.split("|")[1]:
                    #numHomfR = afterRepeat.count("AAGAAGGCAAAGAAG") + afterRepeat.count("CAGAAGGCAAAGAAG") + afterRepeat.count("AAGAAGGCAAGAGAAG") + afterRepeat.count("AAGAGGGCAAAGAAG") + afterRepeat.count("AAGAGGGCAGTGAAG")
                    #numHomfR = afterRepeat.count("AAGAAGGCAAGAGAAG") + afterRepeat.count("CAGAAGGCAAAGAAG")+ afterRepeat.count("AAGAGGGCAAAGAAG") + afterRepeat.count("AAGAGGGCAGTGAAG")
                    numHomfR = len(re.findall(r"AAGAAGGCAAGAGAAG",afterRepeat))
                    numHomfR += len(re.findall(r"\w{1,1}AGA\w{1,1}GGCA\w{1,1}\w{1,1}GAAG",afterRepeat))

                self.repeatLen = self.repeatLen+numHomR*len("AAAACGGTGGTAAGG")+numHomR2*len("AAAACGGTGAAGG")    
            
            LOG.info(f"Found {numAcatca} ACATCA repeats")
            #repeat region sequence
            seq = copy.copy(self.seq[self.repeatStarts:self.repeatEnds])
            LOG.info(f"Found repeat {self.repeatLen} bp sequence in {len(self.seq)} bp {self.name} with coordinates {self.repeatStarts+1}-{self.repeatStarts+self.repeatLen}: {self.seq[self.repeatStarts:self.repeatStarts+self.repeatLen]}")
            
            numTCA = seq.count("TCA")
            numTCG = seq.count("TCG")
            numTCT = seq.count("TCT")
            
            LOG.info(f"In repeat region found {numTCA} TCA  {numTCG} TCG  {numTCT} TCT ...")
            
            repeat=""
            family=""
            goodRepeat = True
        
            if len(self.species.split("|")) >= 2:
                family=self.species.split("|")[1].split("(")[0]
            
            #piecing together the number of repeats and writing in standard format
            #  Format: A#G#T#R#
            if (numTCA!=0):
                repeat += "A" + str(numTCA)
         
            if (numTCG!=0):
                repeat += "G" + str(numTCG)
          
            if (numTCT!=0):
                if family == "Ib" or family == "Ie" or re.match(r"IV[a-z]{1,1}",family) or family == "VIIa" \
                or family == "XIa" or "XIV" in family:
                    repeat += "T" + str(numTCT)
                elif family == "XVIIa" and re.search(r"GGTGTTACCACTGCTCCTGTGGCA", afterRepeat):
                    repeat += "R1"    
                else:
                    goodRepeat = False
            
           
            if re.match(r"XXV[a-z]{1,1}|XXIV[a-z]{1,1}",family):
                repeat = ""        
            
            if 'Apodemus' in self.species:
                numRepeatApodemus = len(re.findall(r"GGTGTTACCACTGCTCCTGTGGCA",afterRepeat))
                if numRepeatApodemus:
                    repeat = "R"+str(numRepeatApodemus)
                else:
                    repeat = ""           
            
            if re.match(r"XXV[a-z]{1,1}",family):
                RepeatsSuis = re.findall(r"GGTG\w{1,1}TCAAG\w{1,1}GAATGC\w{1,1}CAG",self.seq)
                numRepeatsSuis = len(RepeatsSuis)
                if RepeatsSuis:
                    repeat = "R"+ str(numRepeatsSuis)

            if family == "Vb":
                repeat += "R"+ str(numACA)

            if (numAcatca!=0) and re.match(r"II[l,a,t]{1,1}|XIII[a-z]{1,1}",family):
                repeat += "R" + str(numAcatca)

            if (numHomR != 0) and family =="Ia":
                repeat += "R" + str(numHomR)

            
            if (numHomfR != 0) and family == "If":
                repeat += "R" + str(numHomfR)

                
            if not goodRepeat or re.match(r"XX[a-z]{1,1}|XXI[a-z]{1,1}|XXII[a-z]{1,1}|XXIII[a-z]{1,1}|XXVI[a-z]{1,1}",family):
                repeat = ""
            
            self.repeats = repeat
            if repeat:
                LOG.info(f"Final REPEAT REGION standard nomenclature encoded value is '{self.repeats}'")
        
            return True

        else:
            LOG.warning("Sequence of BAD quality to determine repeats. Repeat nomenclature is skipped ...")
            self.repeats = ""
            return False


    def determineFamily(self,customdatabasename=None):
        LOG.info(f"Determine family and species of {self.name} of {len(self.seq)}bp sequence ...")
        # Filename to write
        filename = "query.txt"

        # Open the file with writing permission
        myfile = open(filename, 'w')
        
        # Write a line to the file after repeat region or full length
        #if self.repeatEnds != 0:
        #    myfile.write(''.join(self.seq[self.repeatEnds:]))
        #else:
        myfile.write(f">{self.name}\n")
        myfile.write(''.join(self.seq))


        # Close the file
        myfile.close()
        
        if customdatabasename:
            blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn', query="query.txt", dust='no',
                                                 db="custom_db",
                                                 reward=1, penalty=-2, gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="gp60result.xml")
        else:
            blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn',query="query.txt", dust='yes',
                                             db=os.path.dirname(__file__)+"/reference_database/gp60_ref.fa", 
                                             reward=1, penalty=-2, gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="gp60result.xml")
        
        LOG.info(f"Querying {self.name} against the {os.path.basename(blastn_cline.db)} gp60 database ...")
        stdout, stderr = blastn_cline()
        LOG.debug(f"BLAST stdout {stdout} and stderr {stderr}...")


        if (os.stat("gp60result.xml").st_size == 0):
            #self.tabfile.write("Could not determine sequence. Check manually")
            self.species="No blast hits."

        else:
            result_handle = open("gp60result.xml", 'r')
            blast_records = list(NCBIXML.parse(result_handle)) #blast records are sorted by bitscore by default
            
            blast_record = blast_records[0] #take the first record only 
            #print([(a.hsps[0].identities/a.hsps[0].align_length, a.hsps[0].bits) for a in blast_record.alignments])
            
            if len(blast_record.alignments) > 0:
                #sort BLAST hists based on IDENTITY and if ties then by BITSCORE as this works better for reference alleles of diff length
                blast_record.alignments = utilities.sort_blast_hits_by_id_and_bitscore(blast_record)
                top10hits=[f"\n{idx+1} - {a.hit_id}: accession:{a.accession}, ref_length:{a.length}, bitscore:{a.hsps[0].bits}, score:{a.hsps[0].score}, identity:{round((a.hsps[0].identities/a.hsps[0].align_length)*100,1)}%, query_coverage: {int((a.hsps[0].align_length/blast_record.query_length)*100)}%, gaps:{a.hsps[0].gaps}, strand:{a.hsps[0].strand}, query match coordinates (start-end):{a.hsps[0].query_start}-{a.hsps[0].query_end}" for r in blast_records for idx,a in enumerate(r.alignments) if idx < 10]
                LOG.debug("Top 10 hits in species and determine family identification:"+"".join(top10hits))
                top10bitscores = [a.hsps[0].bits for r in blast_records for idx,a in enumerate(r.alignments) if idx <= 10]
                #find hits with identical bitscore and warn user
             
                if len(top10bitscores) != len(set(top10bitscores)):
                    duplicated_bitscores = [item for item, count in collections.Counter(top10bitscores).items() if count > 1]
                    identicalHits = [a.hit_id for r in blast_records for idx,a in enumerate(r.alignments) if idx <= 10 and a.hsps[0].bits in duplicated_bitscores]
                    identicalHits_str = ", ".join(identicalHits)
                    LOG.warning("Hits with identical BLAST bitscore are found. The list of duplicated bitscore(s): " \
                        f"{duplicated_bitscores} ({identicalHits_str})")
                    
                    
                    self.ambigSpeciesBlast = identicalHits
                else:
                    self.ambigSpeciesBlast = []   
                    
            
                if len(blast_record.alignments) >= 2:
                    br_alignment1 = blast_record.alignments[0]
                    hsp1 = br_alignment1.hsps[0]

                    br_alignment2 = blast_record.alignments[1]
                    hsp2 = br_alignment2.hsps[0]

                    if hsp1.score == hsp2.score or hsp2.align_length-hsp1.align_length >=25:
                        query_from1 = hsp1.query_start - 1
                        query_to1 = hsp1.query_end
                        sbjctseq1 = hsp1.sbjct
                        sbjct_from1 = hsp1.sbjct_start - 1
                        sbjct_to1 = hsp1.sbjct_end
                        num_identities1 = hsp1.identities
                        align_len1 = hsp1.align_length
                        identity1 = num_identities1/align_len1

                        query_from2 = hsp2.query_start - 1
                        query_to2 = hsp2.query_end - 1
                        sbjctseq2 = hsp2.sbjct
                        sbjct_from2 = hsp2.sbjct_start - 1
                        sbjct_to2 = hsp2.sbjct_end - 1
                        num_identities2 = hsp2.identities
                        align_len2 = hsp2.align_length
                        identity2 = num_identities2/align_len2




                        if identity1 >= identity2:
                            self.species = br_alignment1.hit_id
                            data = br_alignment1.hit_id.split(" ")[0]


                        else:
                            self.species = br_alignment2.hit_id
                            data = br_alignment2.hit_id.split(" ")[0]

                    else:
                        hit = blast_record.alignments[0].hit_id
                        data = hit.split(" ")[0]
                        self.species = blast_record.alignments[0].hit_id

                else:
                    hit = blast_record.alignments[0].hit_id
                    data = hit.split(" ")[0]
                    self.species = blast_record.alignments[0].hit_id
                
                LOG.info(f"Species identified {self.species} from {os.path.basename(blastn_cline.db)} gp60 database")    
            else:
                LOG.warning("No blast hits.")
                self.species="No blast hits."
            
    def blast(self, sequence, contig, filetype, customdatabasename=None):
        if customdatabasename:
            blastdbpath="custom_db"
        else:    
            blastdbpath=os.path.dirname(__file__)+"/reference_database/gp60_ref.fa"

        
        if not contig:
            # Filename to write
            filename = "query.txt"

            # Open the file with writing permission
            myfile = open(filename, 'w')

            # Write a line to the file
            myfile.write(f">{self.name}\n")
            myfile.write(sequence)
            query = sequence

            # Close the file
            myfile.close()

        
            blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn', query="query.txt", dust='yes',
                                             db=os.path.dirname(__file__)+"/reference_database/blast_gp60.fasta", 
                                             reward=1, penalty=-2,gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="gp60result.xml")
        
            LOG.info(f"Running BLAST on {len(sequence)}bp sequence from {self.name} on {os.path.basename(blastn_cline.db)} and contig value is {contig}")
            blastn_cline()
            if (os.stat("gp60result.xml").st_size == 0):
                LOG.warning("No BLAST hits found! Check database global blast_gp60.fasta and if input sequence belongs to Cryptosporidium species")
                return "","","","","","","","",""

            result_handle = open("gp60result.xml", 'r')
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)
            
       
            if len(blast_record.alignments) == 0:
                LOG.warning("No BLAST hits found! Check database global blast_gp60.fasta and if input sequence belongs to Cryptosporidium species")
                return "","","","","","","","",""
            else:
                blast_record.alignments = utilities.sort_blast_hits_by_id_and_bitscore(blast_record)    
                br_alignment = blast_record.alignments[0]
                hsp = br_alignment.hsps[0]
            
                subjectSeq = hsp.sbjct

                file = open("align.fa", 'w')
                file.write("\n>Seq1\n")
                file.write(sequence)
                file.write("\n>Seq2\n")
                file.write(subjectSeq)
                file.close()
                
                if len(sequence) < 10000: #do not run more than 10000bp on ClustalW as it will be slow
                    clustalw_cline = ClustalwCommandline("clustalw", infile="align.fa")
                    clustalw_cline()

                    align = AlignIO.read("align.aln", "clustal")

                    for record in align:
                        if record.id == "Seq1":
                            query = record.seq
                        else:
                            subject = record.seq
                else:
                    LOG.warning(f"Submitted sequence is larger that 10000 bp ({len(sequence)}bp) so ClustalW protein alignment will be skipped")
                    query = sequence
                    subject = subjectSeq    


            

            if filetype in definitions.SANGER_FILETYPES:
                newseq=""
                tempA = []
                tempC = []
                tempT = []
                tempG = []
                tempPL = []

                for i in range(0, len(query)):
                    if abs(self.peakLoc[i]-self.peakLoc[i+1]) < 10:
                        pass

                    if subject[i] == '-':
                        newseq += query[i]
                    else:
                        query = query[i:]
                        subject = subject[i:]

                        tempA = self.a
                        tempG = self.g
                        tempC = self.c
                        tempT = self.t
                        tempPL = self.peakLoc

                        self.a = self.a[i:]
                        self.g = self.g[i:]
                        self.t = self.t[i:]
                        self.c = self.c[i:]
                        self.peakLoc = self.peakLoc[i:]
                        break

                l = len(query)
            
                queryOffset=0
                subjectOffset=0

                for i in range(0, l):
                    #print(l, query[i+queryOffset], subject[i+subjectOffset], self.a[i], self.g[i], self.c[i], self.t[i])

                    if (i+queryOffset) >= len(self.peakLoc) or (i+subjectOffset)>=len(subject):
                        break
                    if query[i+queryOffset] != subject[i+subjectOffset]:
                        if queryOffset == -1:
                            if (i+queryOffset+1) >= len(self.peakLoc):
                                pass
                            elif abs(self.peakLoc[i+queryOffset]-self.peakLoc[i+queryOffset+1]) <= 10:
                                queryOffset = 0
                                subjectOffset = -1

                        elif (i+queryOffset)>= self.repeatStarts and (i+queryOffset)<=self.repeatEnds:
                            #print("IN")
                            if (i+queryOffset+1) >= len(self.peakLoc):
                                pass
                            elif abs(self.peakLoc[i+queryOffset]-self.peakLoc[i+queryOffset+1]) <= 10 or abs(self.peakLoc[i+queryOffset]-self.peakLoc[i+queryOffset-1]) <= 10:
                                if subject[i+subjectOffset] != '-' and query[i+queryOffset]!='-' and (query[i+queryOffset] == query[i+queryOffset-1] or query[i+queryOffset] == query[i+queryOffset+1]):
                                    newseq+= subject[i+subjectOffset]
                                else:
                                    if query[i+queryOffset] != '-':
                                        newseq += query[i+queryOffset]
                            else:
                                if query[i+queryOffset] != '-':

                                    a= self.a[i+queryOffset]
                                    c= self.c[i+queryOffset]
                                    g= self.g[i+queryOffset]
                                    t= self.t[i+queryOffset]

                                    maxAmp = max(a,c,g,t)
                                    lri = 1.6 #2.1


                                    if subject[i+subjectOffset] == 'A':
                                        if a != 0 and a != maxAmp:
                                            lri = math.log((maxAmp/a),2)
                                    elif subject[i+subjectOffset] == "G":
                                        if g != 0 and g != maxAmp:
                                            lri = math.log((maxAmp/g),2)
                                    elif subject[i+subjectOffset] == "C":
                                        if c != 0 and c!=maxAmp:
                                            lri = math.log((maxAmp/c),2)
                                    elif subject[i+subjectOffset] == "T":
                                        if t != 0 and t!=maxAmp:
                                            lri = math.log((maxAmp/t),2)

                                    if lri <= 1.5: #2.0
                                        newseq += subject[i+subjectOffset]
                                    else:
                                        newseq += query[i+queryOffset]

                        elif query[i+queryOffset] == "-":
                            pass
                        elif subject[i+subjectOffset] == "-":
                            if (i+queryOffset+1) >= len(self.peakLoc):
                                pass
                            elif abs(self.peakLoc[i+queryOffset]-self.peakLoc[i+queryOffset+1]) < 10 : #>9
                                queryOffset = -1

                            else:
                                pass
                        else:
                            if i+queryOffset < len(self.a):
                                a= self.a[i+queryOffset]
                                c= self.c[i+queryOffset]
                                g= self.g[i+queryOffset]
                                t= self.t[i+queryOffset]

                                #print(a,c,g,t, query[i+queryOffset], subject[i+subjectOffset])

                                maxAmp = max(a,c,g,t)
                                lri = 3.1

                                if subject[i+subjectOffset] == 'A':
                                    if a != 0 and a != maxAmp:
                                        lri = math.log((maxAmp/a),2)
                                elif subject[i+subjectOffset] == "G":
                                    if g != 0 and g != maxAmp:
                                        lri = math.log((maxAmp/g),2)
                                elif subject[i+subjectOffset] == "C":
                                    if c != 0 and c!=maxAmp:
                                        lri = math.log((maxAmp/c),2)
                                elif subject[i+subjectOffset] == "T":
                                    if t != 0 and t!=maxAmp:
                                        lri = math.log((maxAmp/t),2)

                                if lri <= 3.0:
                                    newseq += subject[i+subjectOffset]
                                else:
                                    rest = 0
                                    if maxAmp == a:
                                        rest = max(c,g,t)
                                    elif maxAmp == g:
                                        rest = max(c,a,t)
                                    elif maxAmp == c:
                                        rest = max(g,a,t)
                                    elif maxAmp == t:
                                        rest = max(c,a,g)

                                    beforeMax = max(self.a[i+queryOffset-1], self.g[i+queryOffset-1], self.c[i+queryOffset-1],self.t[i+queryOffset-1])
                                    beforeBase = query[i+queryOffset-1]

                                    tocontinue = False

                                    if beforeBase == 'A' and beforeMax == self.a[i+queryOffset-1] or beforeBase == 'G' and beforeMax == self.g[i+queryOffset-1] or beforeBase == 'C' and beforeMax == self.c[i+queryOffset-1] or beforeBase == 'T' and beforeMax == self.t[i+queryOffset-1]:
                                        if (i+queryOffset+1) < len(self.a):
                                            afterMax = max(self.a[i+queryOffset+1], self.g[i+queryOffset+1], self.c[i+queryOffset+1],self.t[i+queryOffset+1])
                                            afterBase = query[i+queryOffset+1]

                                            if afterBase == 'A' and afterMax == self.a[i+queryOffset+1] or afterBase == 'G' and afterMax == self.g[i+queryOffset+1] or afterBase == 'C' and afterMax == self.c[i+queryOffset+1] or afterBase == 'T' and afterMax == self.t[i+queryOffset+1]:
                                                #newseq += query[i+queryOffset]
                                                tocontinue = False

                                            else:
                                                tocontinue = True
                                        else:
                                            tocontinue=True
                                    else:
                                        tocontinue = True

                                    if tocontinue:
                                        if subject[i+subjectOffset] == 'A':
                                            if a == rest:
                                                newseq += subject[i+subjectOffset]
                                            else:
                                                newseq += query[i+queryOffset]

                                        elif subject[i+subjectOffset] == 'C':
                                            if c == rest:
                                                newseq += subject[i+subjectOffset]
                                            else:
                                                newseq += query[i+queryOffset]

                                        elif subject[i+subjectOffset] == 'G':
                                            if g == rest:
                                                newseq += subject[i+subjectOffset]
                                            else:
                                                newseq += query[i+queryOffset]
                                        elif subject[i+subjectOffset] == 'T':
                                            if t == rest:
                                                newseq += subject[i+subjectOffset]
                                            else:
                                                newseq += query[i+queryOffset]
                                    else:
                                        newseq += query[i+queryOffset]
                            else:
                                newseq += query[i+queryOffset]

                    else:
                        newseq += query[i+queryOffset]
            
                if tempA != [] and tempC != [] and tempG != [] and tempT != []:
                    self.a = tempA
                    self.c = tempC
                    self.g = tempG
                    self.t = tempT
                    self.peakLoc = tempPL

                self.seq = newseq

            # Filename to write
            filename = "query.txt"

            # Open the file with writing permission
            myfile = open(filename, 'w')
            self.findRepeatRegion()
            # Write a line to the file
            #myfile.write(''.join(self.seq[self.repeatEnds:])) #original only align from repeat onwards
            myfile.write(''.join(self.seq))
       
                         
            # Close the file
            myfile.close()
            LOG.debug(f"Running BLAST on {len(self.seq)}bp sequence ...")
            blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn', query="query.txt", dust='yes',
                                                 db=blastdbpath, reward=1, penalty=-2,gapopen=5, gapextend=2,
                                                 evalue=0.00001, outfmt=5, out="gp60result.xml")

            blastn_cline()

            if (os.stat("gp60result.xml").st_size == 0):
                LOG.warning(f"BLAST result file is empty (zero size)! Check database {blastdbpath} or inputs!")
                return "","","","","","","","",""

            else:
                result_handle = open("gp60result.xml", 'r')
                blast_records = list(NCBIXML.parse(result_handle))
                #blast_record = next(blast_records)
                
                if not blast_records:
                    LOG.warning(f"No BLAST hits found! Check database {blastdbpath} or inputs")
                    return "","","","","","","","",""
                
                if len(blast_records[0].alignments) == 0:
                    LOG.warning(f"No BLAST hits found! Check database {blastdbpath} or inputs")
                    return "","","","","","","","",""

                blast_record = blast_records[0]
                blast_record.alignments = utilities.sort_blast_hits_by_id_and_bitscore(blast_record)


                br_alignment = blast_record.alignments[0]
                hsp = br_alignment.hsps[0]
                
                top10hits=[f"\n{idx+1} - {a.hit_id}: {a.accession}, length:{a.length}, bitscore:{a.hsps[0].bits}, score:{a.hsps[0].score},"\
                            f"identity: {round((a.hsps[0].identities/a.hsps[0].align_length)*100,1)}%, query_cov: {round((a.hsps[0].align_length/blast_record.query_length)*100)}%," \
                            f"ref_cov: {round(a.hsps[0].align_length/a.length*100,2)}% ,"\
                            f"gaps:{a.hsps[0].gaps}, strand:{a.hsps[0].strand}, query coord(start-end):{a.hsps[0].query_start}-{a.hsps[0].query_end}" 
                            for r in blast_records for idx,a in enumerate(r.alignments) if idx < 10]            
                LOG.debug("Top 10 BLAST hits:"+"".join(top10hits))
            

                LOG.info(f"Top hit strand orientation (query={self.name}, target={br_alignment.accession}) is {hsp.strand}")

                if hsp.strand == ("Plus", "Minus"):
                    LOG.warning(f"The reverse complement query {self.name} is being submitted in 3'->5' orientation. Performing reverse complement of the sequence (5'->3')")
                    sequence = str(Seq(sequence).reverse_complement())
                    self.seq=sequence
                
                percent_identity = round(hsp.identities/hsp.align_length*100,2)
                evalue = hsp.expect
                bitscore = hsp.score

                if hsp.align_length <= blast_record.query_length :
                    query_coverage = round(hsp.align_length/blast_record.query_length*100, 1)
                    # query_coverage = round( ((hsp.query_end - hsp.query_start) + 1)/blast_record.query_length * 100,1)
                else:
                    query_coverage = 100
                query_length = blast_record.query_length
                subject_coverage = round( abs(hsp.sbjct_end - hsp.sbjct_start) / br_alignment.length, 3) * 100
                accession = blast_record.alignments[0].hit_id
                species = br_alignment.hit_id.split("|")[0]

                if "|" in accession:
                    accession = accession.split("|")[-1].split("(")[1].split(")")[0]

                return bitscore,evalue,query_coverage,query_length,percent_identity, accession, self.seq, species, subject_coverage

    
        else:
            # Filename to write
            filename = "query.txt"

            # Open the file with writing permission
            myfile = open(filename, 'w')
       
            self.seq = sequence
            self.findRepeatRegion()
          
            # Write a line to the file
            #myfile.write(''.join(sequence[self.repeatEnds:]))
            myfile.write("".join(sequence))

            # Close the file
            myfile.close()

            blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn',query="query.txt", dust='yes',
                                                 db=blastdbpath, reward=1, 
                                                 penalty=-2,gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="gp60result.xml")

            LOG.info(f"Running BLASTN on filtered sequence {len(self.seq)}bp from {self.name} on {blastn_cline.db}")
            stdout, stderr = blastn_cline()

            if (os.stat("gp60result.xml").st_size == 0):
                return "","","","","","","","",""

            result_handle = open("gp60result.xml", 'r')
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)

            if len(blast_record.alignments) == 0:
                return "","","","","","","","",""

            br_alignment = blast_record.alignments[0]
            hsp = br_alignment.hsps[0]

            percent_identity = round((hsp.identities/hsp.align_length)*100,2)
            evalue = hsp.expect
            bitscore = hsp.score
            if hsp.align_length <= blast_record.query_length :
                query_coverage = round(hsp.align_length/blast_record.query_length*100)
            else:
                query_coverage = 100
            subject_coverage = round( abs(hsp.sbjct_end - hsp.sbjct_start) / br_alignment.length, 3) * 100    
            query_length = blast_record.query_length
            accession = blast_record.alignments[0].hit_id
            species = br_alignment.hit_id.split("|")[0]

            if "|" in accession:
                accession = accession.split("|")[-1].split("(")[1].split(")")[0]
            LOG.info(f"{self.name} BLAST results on {query_length} bp {query_coverage}% coverage, {percent_identity}% identity, top hit accession {accession}")    

            return bitscore,evalue,query_coverage,query_length,percent_identity, accession, sequence, species, subject_coverage



    #printFasta() prints the results of the alignment and the repeat
    #   region for each sample. It then prints the fasta sequence for
    #   each sample that can successfully be analyzed
    def printFasta(self, contig, mode, sampleName, filetype="ab1", customdatabasename = None):
        LOG.info(f"Running printFasta() on {self.name} with species '{self.species}' and repeat encoding '{self.repeats}' ...")
        bitscore = 0.0
        evalue = 0.0
        query_coverage = 0.0
        subject_coverage = 0.0 #ref allele coverage
        query_length = 0
        percent_identity = 0.0
        accession = ""
        seq = ""
        species = ""
        self.file.write("\n>" + sampleName)
        

        #Output sample name
        #self.tabfile.write(self.name.split(".ab1")[0] + "\t")
        self.tabfile.write(sampleName + "\t")

        #Output type of sequences
        self.tabfile.write(mode + "\t")
        
        #for forward or reverse inputs

        if contig == "":
            sequence=str("".join(self.seq))
            if sequence != "Poor Sequence Quality":
                LOG.info(f"Running BLAST on {len(sequence)}bp sequence of acceptable quality and getting top hit accession ...")
                bitscore,evalue,query_coverage,query_length,percent_identity, accession, seq, species, subject_coverage = self.blast(sequence,False, filetype, customdatabasename)
                if self.species != species:
                    LOG.warning(f"Species mismatch w.r.t to accession number {self.species} vs {species}. Fixing it")
                    self.species = species #making sure the accession number matches the species
            else:
                LOG.info("Poor Sequence Quality. BLAST would not be run")

        else:
            sequence = contig

            if sequence != "Poor Sequence Quality":
                LOG.info(f"Running BLAST on acceptable quality {len(sequence)}bp contig and getting top hit accession ...")
                bitscore,evalue,query_coverage,query_length,percent_identity, accession, seq, species, subject_coverage = self.blast(sequence, True, customdatabasename) 
                
                if self.species != species:
                    LOG.warning(f"Species mismatch w.r.t to accession number {self.species} vs {species}. Fixing it")
                    self.species = species #making sure the accession number matches the species
            else:
                LOG.info("Poor Sequence Quality. BLAST would not be run")
                
            
        
        #1.Sample Name
        #2.Type of Sequences
        #3.Species
        #4.Subtype
        #5.Sequence
        #6.Comments
        #7.Avg. Phred Quality
        #8.Bit Score
        #9.Query Length (bp)
        #10.Query Coverage
        #11.E-value
        #12.Percent Identity
        #13. Accession Number
        if self.seq == "Poor Sequence Quality":
            avg_phred_msg = ""
            if self.fileType == "sanger":
                avg_phred_msg = f"Average PHRED Quality = {self.averagePhredQuality}" 
            self.tabfile.write("\t\t\tPoor Sequence Quality. Check manually." + avg_phred_msg + "\t\t\t\t\t\t\t\n")
            self.file.write(" | Poor Sequence Quality (" + avg_phred_msg + "). Check manually.")
        
        #elif evalue > 1e-75 or query_coverage < 50:
        #    self.tabfile.write("\t\t\tCould not determine species from chromatogram. Check manually.\t" + str(self.averagePhredQuality)+ "\t\t\t\t\t\t\n")
        #    self.file.write("\n;Could not determine species from chromatogram. Check manually. Sequence Phred Quality is "+ str(self.averagePhredQuality))    
        
        elif seq == "" or self.species == "No blast hits.":
            avg_phred_msg = ""
            if self.fileType == "sanger":
                avg_phred_msg = f"Average PHRED Quality = {self.averagePhredQuality}"
            self.tabfile.write("\t\t\tNo blast hits. Check manually." + avg_phred_msg + "\t\t\t\t\t\t\t\n")
            self.file.write("\n;No blast hits. Check manually. Sequence PHRED Quality is " + avg_phred_msg)

        else:
            self.seq = seq
            self.determineFamily(customdatabasename)

            #Output Species and Subfamily(ex. C.parvum\tIIa)
            if len(self.species.split("(")) > 0:
                speciesName = self.species.split("(")[0]
                
                if len(speciesName.split("|")) >= 2:
                    self.tabfile.write(speciesName.split("|")[0] + "\t" + speciesName.split("|")[1])
                    self.file.write(" | " + speciesName.split("|")[0] + " " + speciesName.split("|")[1])
                else:
                    self.tabfile.write(self.species + "\t")
                    self.file.write(" | " + self.species)
            else:
                self.tabfile.write(self.species + "\t")
                self.file.write(" | " + self.species)


            foundRepeat=True
            #Find repeats
            if "ubiquitum" not in speciesName and "felis" not in speciesName:
                self.forwardSeq = True
                self.seq = seq
                #again find repeat region
                self.findRepeatRegion()
                
            
                #if self.checkRepeatManually and self.doublePeaksinRepeat: #self.repeats == "Could not classify repeat region. Check manually." or (
                #    foundRepeat = False
                #elif self.repeatStarts < 3: #need some sequence before the repeat region to trust it
                #    foundRepeat = False
                #else:
                foundRepeat = self.determineRepeats()
                #print(foundRepeat, self.repeats, self.tabfile.getvalue());exit()
                if foundRepeat:
                    self.tabfile.write(self.repeats)
                    self.file.write(self.repeats)
                


                #Still outputting repeats if there's a subfamily
                if len(speciesName.split("|")) == 3 and foundRepeat == True:
                    if len(self.repeats):
                        subfamily = speciesName.split("|")[2]
                    else:
                        subfamily = "|"+speciesName.split("|")[2]    
                    
                    LOG.debug(f'Appending family subtype {subfamily} to {speciesName.split("|")[1]}{self.repeats}')    
                    self.file.write(subfamily)
                    self.tabfile.write(subfamily)
                

            #Output sequence
            self.tabfile.write("\t" + str(seq))
            self.file.write("\n" +str(seq))

            #Quality check in comments section of the report
            qc_messages = [] # Initialize a list to collect all QC messages
            if self.ambigSpeciesBlast:
                qc_messages.append(f"Check manually.")
            if not foundRepeat:
                qc_messages.append("Could not classify repeat region.")

            elif len(self.species.split("|")) > 1 and "If" in self.species.split("|")[1]:
                    begin = seq.find("CACCCCAACTC")
                    end = seq.find("AAGCTACTCCAAAGG")
                    if begin != -1 and end != -1:
                        begin += 10

                    insert = end - begin
                    if insert > 1 and self.checkRepeatManually:
                        qc_messages.append(f"Not all bases in repeat region had phred quality >= 20. A {insert}bp insert was found.")
                    elif insert > 1:
                        qc_messages.append(f"A {insert}bp insert was found.")
                    elif percent_identity < 99.1:
                        qc_messages.append("BLAST percent identity less than 99.1%. Check manually in case of new gp60 family.")
    
            if self.doublePeaksinRepeat:
                qc_messages.append("Found double peaks in repeat region.")    
            if self.checkRepeatManually:
                qc_messages.append("Check repeat region manually. Not all bases in repeat region had phred quality >= 20.")
            if percent_identity < 99.1:
                qc_messages.append("BLAST percent identity less than 99.1%. Check manually in case of new gp60 family.")
            if subject_coverage < 60:
                qc_messages.append(f"Reference allele coverage is < 60% ({subject_coverage :.1f}%).")    
            if foundRepeat == False:
                qc_messages.append("Could not classify repeat region.")   
            
            if not qc_messages:
                final_comment = "-"
            else:
                if "Check manually." not in qc_messages:
                    qc_messages.append("Check manually.")
                # Join all collected messages with a space
                final_comment = " ".join(qc_messages)
            
            # Write the final consolidated comment to the tabfile, followed by a single tab
            self.tabfile.write(f"\t{final_comment}\t")

            #write the average phred quality of the sequence
            #if its a contig, will output as #(F), #(R) to show phred of both
            self.tabfile.write(str(self.averagePhredQuality) + "\t")
            self.tabfile.write(str(bitscore)+"\t")
            self.tabfile.write(str(query_length)+"\t")
            self.tabfile.write(str(round(query_coverage,2)) + "%\t")
            self.tabfile.write(str(evalue) +"\t")
            self.tabfile.write(str(percent_identity) + "%\t")
            self.tabfile.write(str(accession)+"\n")



#Function to build the contig of forward and reverse sequences
#assumes that both s1 and s2 are in forward 5'->3' orientation
#falls back to BLAST if CLUSTALW is giving aligment with identity lower than 90%
def buildContig(s1, s2):
    contig=""
    LOG.info("Building a gp60 contig from sequences ...")
    if s1 == "":
        return s2
    elif s2 == "":
        return s1
    else:
        file = open("align.fa", 'w')
        file.write("\n>Seq1\n")
        file.write(str(s1))
        file.write("\n>Seq2\n")
        file.write(str(s2))
        file.close()

        clustalw_cline = ClustalwCommandline("clustalw", infile="align.fa")
        stdout, stderr = clustalw_cline()
        LOG.debug(f"stdout:\n{stdout}\nstderr:\n{stderr}")

        align = AlignIO.read("align.aln", "clustal")
        alignment_length = align.get_alignment_length()
        seq1_aligned = str(align[0].seq)
        seq2_aligned = str(align[1].seq)
        gaps_in_seq1 = 0
        gaps_in_seq2 = 0
        total_bases = 0
        exact_matches = 0
        # Iterate through each column of the alignment
        for i in range(alignment_length):
            char1 = seq1_aligned[i]
            char2 = seq2_aligned[i]
            
            # Count gaps
            if char1 == '-':
                gaps_in_seq1 += 1
            if char2 == '-':
                gaps_in_seq2 += 1
                
            # Count positions where both sequences are not gaps (counts only regions were BOTH sequences have a base at a given position)
            if char1 != '-' and char2 != '-':
                total_bases += 1
                # Count exact matches
                if char1 == char2:
                    exact_matches += 1
        percentage_identity = (exact_matches / total_bases) * 100 if total_bases > 0 else 0
        LOG.info(f"ClustalW: identity {percentage_identity :.2f} %, exact matches: {exact_matches}, gaps in seq1: {gaps_in_seq1}, gaps in seq2: {gaps_in_seq2}")

        if percentage_identity > 90:
            for record in align:
                if record.id == "Seq1":
                    seq1 = record.seq
                else:
                    seq2 = record.seq

            # Calculate the length of the overlap based on the alignment
            overlap_length = 0
            for i in range(len(seq1)):
                if seq1[i] != "-" and seq2[i] != "-":
                    overlap_length += 1
            
            # Log the overlap length and other info
            

            l = len(seq1)

            for i in range(0, l):
                if seq1[i] == "-":
                    contig += seq2[i]
                elif seq2[i] == "-":
                    contig+=seq1[i]
                elif seq1[i] == seq2[i]:
                    contig += seq1[i]
                elif i < l/2:
                    contig += seq1[i]
                else:
                    contig += seq2[i]

            LOG.info(f"A {len(contig)}bp contig formed ({len(s1)}bp of forward + {len(s2)}bp reverse) using ClustalW")
            LOG.info(f"Detected {overlap_length}bp overlap region between sequences ({overlap_length/len(contig)*100 :.1f}% of contig length).")
        else:
            #build contig with BLASTN as a 2nd fallback option. Both sequences are in forward orientation
            with open("query_tmp.fasta", "w") as f:
                f.write(f">query_seq\n{s1}\n")

            with open("subject_tmp.fasta", "w") as f:
                f.write(f">subject_seq\n{s2}\n")
    

            blastn_cline = NcbiblastnCommandline(
                query="query_tmp.fasta",
                subject="subject_tmp.fasta",
                task="dc-megablast", 
                evalue=10, word_size=11,
                outfmt=5,
                out="blastn_contig_results.xml"
            )
            blastn_stdout, blastn_stderr = blastn_cline()

            QUERY_LEN = len(s1)
            SBJCT_LEN = len(s2) # Length of the original reverse sequence
            LOG.info(f"Building contig from forward {QUERY_LEN}bp and reverse {SBJCT_LEN} sequences using BLASTN instead of ClustalW as alignment was poor (%identity < 90)");

            # Read and parse the BLAST XML output
            with open("blastn_contig_results.xml") as blast_results_file:
                blast_records = NCBIXML.parse(blast_results_file)
                    
                # Get the first record, which contains our search results
                first_record = next(blast_records)
                    
                # Check if there are any alignments found
                if first_record.alignments:
                    # Get the first (best) alignment
                    best_alignment = first_record.alignments[0]
                    hsp = best_alignment.hsps[0]

                    query_start_0based = hsp.query_start - 1  #forward seq
                    query_end_0based = hsp.query_end  #forward seq

                    sbjct_start_0based = hsp.sbjct_start - 1
                    sbjct_end_0based = hsp.sbjct_end


                    contig = s1[0 : query_end_0based] + s2[sbjct_end_0based : ]
                    overlap_length = hsp.align_length      
            
                        
                    LOG.info(f"A {len(contig)}bp contig formed (used {len(s1[0 : query_end_0based])}bp of forward (incl. {overlap_length}bp overlap region) + non-overlap {len(s2[sbjct_end_0based : ])}bp of reverse)")            
                else:
                    print("No significant alignments found by BLAST.")



  
        return contig



def writeResults(expName, file, tabfile):
    LOG.info(f"Writing output files ...")
    experimentName = expName + "_"
    output_report_file_name = f"{experimentName}cryptogenotyper_report.fa"
    filename = os.path.join('.', output_report_file_name)
    with open(filename, 'w') as resultFile:
        resultFile.write(file.getvalue())

    output_tabreport_file_name = f"{experimentName}cryptogenotyper_report.txt"
    tab_filename = os.path.join('.', output_tabreport_file_name)

    
    reportdata2write = ""
    with open(tab_filename, 'w') as tabFile:
        reportdata2write = tabfile.getvalue()
        tabFile.write(reportdata2write)
    
    print("\n>>>GP60 RESULTS  REPORT (only first 10 lines are printed)")
    for idx, line in enumerate(reportdata2write.split("\n")):
        print(line)
        if idx == 10:
            LOG.warning(f"Please check the {output_tabreport_file_name} for the completed output ...")
            break

    print(">>> FASTA report written to " + os.getcwd()+"/"+output_report_file_name + "")
    print(">>> Tab-delimited report written to " + os.getcwd() + "/" + output_tabreport_file_name + "")
    LOG.info("FASTA report written to " + os.getcwd()+"/"+output_report_file_name)
    LOG.info("Tab-delimited report written to " + os.getcwd() + "/" + output_tabreport_file_name)



def cleanTempFastaFilesDir(temp_dir="tmp_fasta_files"):
    LOG.info("Cleaning the temporary FASTA and BLAST database files (if any)")
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir, ignore_errors=True)
    tmpfiles2remove = list(itertools.chain.from_iterable([glob.glob(e) for e in ["align.*","custom_db.*"]]))
    for file in tmpfiles2remove:
        try:
            os.remove(file)
        except FileNotFoundError:
            pass       

def gp60_main(pathlist_unfiltered, fPrimer, rPrimer, typeSeq, expName, customdatabasename, noheader, verbose):

    fPrimer = fPrimer.replace(' ', '')
    rPrimer = rPrimer.replace(' ', '')

    
    pathlist = [path for path in pathlist_unfiltered if any([True if path.endswith(type) else False for type in definitions.FILETYPES ])]
    
    
    LOG.info(f"Filtering paths by the primer names (fPrimer:'{fPrimer}' rPrimer:'{rPrimer})'")
    if fPrimer and rPrimer:
        pathlist = [path for path in pathlist if re.search(fPrimer, path) or re.search(rPrimer, path)]  # select only files matching the primers
    elif fPrimer:
        pathlist = [path for path in pathlist if re.search(fPrimer, path)]
        pathlist.sort()
    elif rPrimer:
        pathlist = [path for path in pathlist if re.search(rPrimer, path)]
        pathlist.sort()

    
    #if multi-FASTA file is present in the list slice it up into individual files https://www.metagenomics.wiki/tools/fastq/multi-fasta-format 
    fasta_paths = [path for path in pathlist for fasta_extension in definitions.FASTA_FILETYPES if path.endswith(fasta_extension)]
    utilities.slice_multifasta(typeSeq,fasta_paths,pathlist,expName)

    
    if pathlist == []:
        msg=f"No supported input file(s) found in pathlist {pathlist_unfiltered}. Supported input filetypes are {definitions.FILETYPES}"
        LOG.error(msg)
        sys.exit(msg)
  
    
    
    LOG.info(f"Processing {len(pathlist)} file(s) in {typeSeq} mode.")

    contig = False
    onlyForwards = False
    onlyReverse = False

    if typeSeq == 'forward':
        onlyForwards = True

    elif typeSeq == 'reverse':
        onlyReverse = True

    elif typeSeq == 'contig':
        contig = True



    tabfile = io.StringIO()

    if not noheader:
        tabfile.write("Sample Name\tType of Sequences\tSpecies\tSubtype\tSequence\tComments\tAvg. Phred Quality\tBit Score\tQuery Length (bp)\tQuery Coverage\tE-value\tPercent Identity\tAccession Number\n")

    file = io.StringIO()
    #Write output fasta header with comments
    file.write("\n;>****************************************************************************")
    file.write("\n;>gp60 SEQUENCE ANALYSIS INPUT PARAMETERS:")
    if customdatabasename:
        file.write("\n  ;>Reference File: " + customdatabasename)
    else:
        file.write("\n  ;>Reference File: " + "gp60_ref.fa (default)") #debug this might not be always true, actually it is blast_gp60.fa as default
    file.write("\n  ;>Program mode: " + typeSeq)
    if typeSeq == 'forward':
        file.write("\n  ;>Forward Primer: " + str(fPrimer))
    elif typeSeq == 'reverse':
        file.write("\n  ;>Reverse Primer: " + str(rPrimer))
    elif typeSeq == 'contig':
        file.write("\n  ;>Forward Primer: " + str(fPrimer))
        file.write("\n  ;>Reverse Primer: " + str(rPrimer))
    file.write("\n;>****************************************************************************")
    file.write("\n;>Program Results:\n")
    #**************************************************************


    if contig:
        LOG.info(f"{typeSeq.upper()} only read input only mode started ...")
        if len(pathlist)%2 == 0:
            for idx, path in enumerate(pathlist):
                filetype = utilities.getFileType(path)
                LOG.info(f"\n\n*** {idx+1}: Running sample {os.path.basename(pathlist[idx])} and {os.path.basename(pathlist[idx+1])} as {filetype} in contig mode ***")
                forward = analyzingGp60() #init forward read object
                reverse = analyzingGp60() #init reverse read object


                filetypeF = utilities.getFileType(path)
                utilities.setFileType(forward,filetypeF)
                filetypeR = utilities.getFileType(pathlist[idx+1])
                utilities.setFileType(forward,filetypeR)
             
                
                forwSeqbool = forward.readFiles(path, True, file, tabfile, filetypeF,customdatabasename)
                revSeqbool = reverse.readFiles(pathlist[idx+1], False, file, tabfile,filetypeR, customdatabasename)
                pathlist.remove(pathlist[idx+1])   


                forwardPhred = forward.averagePhredQuality
                reversePhred = reverse.averagePhredQuality



                forward.averagePhredQuality = str(forwardPhred) + "(F); " + str(reversePhred) + "(R)"
                
                
                if forwSeqbool and revSeqbool:
                    goodTrimF = forward.trimSeq(filetypeF)
                    goodTrimR = reverse.trimSeq(filetypeR)
                    
                    #print(dir(forward),dir(reverse),customdatabasename)
                    #print(f"forward seq top hit blast {forward.species} reverse {reverse.species}")
                    
                    if not goodTrimF and not goodTrimR:
                        forward.repeats="Could not classify repeat region. Check manually."
                        reverse.repeats="Could not classify repeat region. Check manually."

                    else:
                        #if customdatabasename:
                        forward.determineFamily(customdatabasename)
                        reverse.determineFamily(customdatabasename)

                        #if the crypto subfamily was found, find repeat region    
                        if forward.species.split('|')[0] == reverse.species.split('|')[0] and forward.species!= "":
                            LOG.info(f"Determining repeats for {forward.name} and {reverse.name} ...")
                            forward.determineRepeats()
                            reverse.determineRepeats()

                    
                        
                    # Form a contig if species of both reads is the same 
                    # e.g. C.mortiferum == C.mortiferum from forward.species=C.mortiferum|XIVa|a(KP099082)            
                    if forward.species != "" and forward.species.split('|')[0] == reverse.species.split('|')[0]:# and forward.repeats == reverse.repeats and forward.repeats != "":
                        #forward.determineFamily(customdatabasename)
                        #reverse.determineFamily(customdatabasename)
                      
                        Fbitscore,Fevalue,Fquery_coverage,Fquery_length,Fpercent_identity, Faccession, Fnewseq, Fspecies, Fsubject_coverage = forward.blast(str(forward.seq), False, customdatabasename)
                        Rbitscore,Revalue,Rquery_coverage,Rquery_length,Rpercent_identity, Raccession, Rnewseq, Rspecies, Rsubject_coverage = reverse.blast(str(reverse.seq), False, customdatabasename)
                        LOG.info(f"Build contig from forward and reverse extracted sequences of {len(Fnewseq)}bp and {len(Rnewseq)}bp")
                        contig = buildContig(Fnewseq, Rnewseq)
                        
                        sampleName = forward.name.split(".ab1")[0] + ", " + reverse.name.split(".ab1")[0]
                        forward.printFasta(contig, "contig", sampleName, customdatabasename)

                    else:
                        #forward.determineFamily(customdatabasename)
                        #reverse.determineFamily(customdatabasename)
                        forward.printFasta("", "forward", forward.name.split(".ab1")[0])
                        reverse.printFasta("", "reverse", reverse.name.split(".ab1")[0])
                    
                elif forwSeqbool and not revSeqbool:
                    goodTrim = forward.trimSeq(filetypeF)

                    if not goodTrim:
                        forward.repeats="Could not classify repeat region. Check manually."

                    else:
                        forward.determineFamily(customdatabasename)

                    forward.printFasta("", "forward", forward.name.split(".ab1")[0], customdatabasename)


                    reverse.seq = "Poor Sequence Quality"
                    reverse.printFasta("", "reverse", reverse.name.split(".ab1")[0], customdatabasename)

                elif revSeqbool and not forwSeqbool:
                    forward.seq = "Poor Sequence Quality"
                    forward.printFasta("", "forward", forward.name.split(".ab1")[0], customdatabasename)

                    goodTrim = reverse.trimSeq(filetypeR)

                    if not goodTrim:
                        reverse.repeats="Could not classify repeat region. Check manually."

                    else:
                        reverse.determineFamily(customdatabasename)

                    reverse.printFasta("", "reverse", reverse.name.split(".ab1")[0])

                else:
                    forward.seq = "Poor Sequence Quality"
                    forward.printFasta("", "contig", forward.name.split(".ab1")[0] + ", " + reverse.name.split(".ab1")[0])

                writeResults(expName, file, tabfile)        

        else:
            print("ERROR: Uneven number of input files ({}). "
                  "Cannot find all paired forward and reverse files. Aborting ...".format(len(pathlist)))
            tabfile.write("\nCannot find all paired forward and reverse files.  Make sure all files are included to produce the contig.")
            file.write("\nCannot find all paired forward and reverse files.  Make sure all files are included to produce the contig.")

    else:
        LOG.info(f"{typeSeq.upper()} only read input only mode started ...")
        for idx, path in enumerate(pathlist):
            filetype = utilities.getFileType(path)
            LOG.info(f"\n\n*** {idx+1}: Running sample {os.path.basename(path)} as {filetype} file type ***")
            forward = analyzingGp60()
            utilities.setFileType(forward,filetype)  

            #check if a given path is indeed in forward or reverse orientation
            #foward only mode 
            if onlyForwards:
                read_ok = forward.readFiles(path, True, file, tabfile, filetype, customdatabasename)
            #reverse only mode
            elif onlyReverse:
                read_ok = forward.readFiles(path, False, file, tabfile, filetype, customdatabasename)
            
            if read_ok == False:
                forward.seq = "Poor Sequence Quality"
                forward.printFasta("", typeSeq, forward.name.split(f".{filetype}")[0], customdatabasename)
            else:
                
                goodTrim = forward.trimSeq(filetype)

                if goodTrim == False:
                    forward.repeats="Could not classify repeat region. Check manually."
                else:
                   
                    forward.determineFamily(customdatabasename)  
                    forward.determineRepeats()
                    
                  
                forward.printFasta("", typeSeq, forward.name.split(f".{filetype}")[0], filetype, customdatabasename)
        writeResults(expName, file, tabfile)    
    
    if verbose == False:
        LOG.info("Cleaning the temporary FASTA and BLAST database files (if any)")
        utilities.cleanTempFastaFilesDir("tmp_fasta_files_"+expName)
    LOG.info("The gp60 run completed successfully")
    
   # os.system("rm gp60result.xml gp60result2.xml")
    #os.system("rm query.txt")
    #os.system("rm align.fa align.dnd align.aln")

if __name__ == "__main__":
	gp60_main()
