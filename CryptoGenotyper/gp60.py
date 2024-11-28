#******************************************************************************
#******************************************************************************
#                 Cryptosporidium Gp60 Analyzer (Version 1.1)
#             Written by: Christine Yanta (christine.yanta@canada.ca)
#                           Date: May 2, 2017
#   Public Health Agency of Canada - National Microbiology Laboratory Guelph
#******************************************************************************
#******************************************************************************


import io
import os
import re
from CryptoGenotyper import logging, definitions

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


global TESTING
TESTING = False


# setup the application logging
LOG = logging.create_logger(__name__)

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





    def readFiles(self, dataFile, forw, output, tabFile, filetype = "abi"):
        #setting whether the sequence is the forward or reverse strand
        #(used later on know whether the reverse complement is needed)
        if forw:
            self.forwardSeq = True

        else:
            self.forwardSeq = False

        self.file = output
        self.tabfile = tabFile

        #**********Beginning to work with the ab1 file given:**********#

        #retrieves the sample file name (removes directory pathway)
        self.name = dataFile.split("/")[len(dataFile.split("/"))-1]

        LOG.debug("\nSequence: ", self.name) #Lets user know which sequence the program is on

        #opens the ab1 file
        if  filetype == "abi":
            handle=open(dataFile,"rb")
            record=SeqIO.read(handle, filetype)
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
        elif filetype == "fasta" or filetype == "fa":
            handle=open(dataFile,"r")   
            record=SeqIO.read(handle, filetype)
            raw_seq = list(record.seq)
            self.seq = raw_seq
            self.phred_qual = [60] * len(raw_seq)

        
        self.seqLength = len(self.seq)


        self.repeatStarts = 0
        self.repeatEnds = 0

        self.b = 0
        self.e = 0

        


        #if any([isinstance(element,str) == False for element in self.seq]):
        #    self.seq = list(record.annotations['abif_raw']['PBAS2'].decode('UTF-8'))
            #self.oldseq = list(record.annotations['abif_raw']['PBAS2'].decode('UTF-8'))
       #     print(self.seq);exit()

        b = 0
        e = 0

        filename = f"query_vs_blast_gp60.txt"
        myfile = open(filename, 'w')

        LOG.debug(f"Sequence FASTA {self.name} "+''.join(self.seq)) #debug
        myfile.write(''.join(self.seq))

        # Close the file
        myfile.close()
        
        dbpath = os.path.join(os.path.dirname(__file__),"reference_database/blast_gp60.fasta")
        blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn',query=filename, dust='yes',
                                             db= dbpath,
                                             reward=1, penalty=-2,gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="gp60result2.xml")
        LOG.info(f"Querying {self.name} seq of {len(self.seq)}bp against the {os.path.basename(blastn_cline.db)} gp60 database ...")
        stdout, stderr = blastn_cline()
        LOG.debug(f"BLASTN stdout={stdout} and stderr={stderr}")
        

        if (os.stat("gp60result2.xml").st_size == 0):
            LOG.error(f"Generated an empty gp60 BLAST result for {self.name}")
            return False
        else:
            result_handle = open("gp60result2.xml", 'r')
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)

            if len(blast_record.alignments) == 0:
                return False

            br_alignment = blast_record.alignments[0]
            hsp = br_alignment.hsps[0]

            self.species = br_alignment.hit_id
            b = hsp.query_start-1
            e = hsp.query_end
            
            if b > (len(self.seq)*0.25):
                dbpath = os.path.join(os.path.dirname(__file__),"reference_database/gp60_ref.fa")
                LOG.info(f"Querying {self.name} against the {dbpath} gp60 database ...")

                filename = f"query_vs_gp60_ref.txt"

                # Open the file with writing permission
                myfile = open(filename, 'w')

                t = round(len(self.seq)*0.5)
                myfile.write(''.join(self.seq[0:t]))

                # Close the file
                myfile.close()

                blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn',query=filename, dust='yes',
                                                     db=os.path.dirname(__file__)+"/reference_database/gp60_ref.fa", reward=1, penalty=-2,gapopen=5, gapextend=2,evalue=1e-5, outfmt=5, out="gp60result.xml")
                LOG.info(f"Querying {self.name} against the {blastn_cline.db} gp60 database ...")
                stdout, stderr = blastn_cline()

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
        LOG.info(f"Average PHRED {self.name} sequence quality is {self.averagePhredQuality}")
      
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
            LOG.info(f"Final sequence length {self.seqLength}bp selected from raw {len(raw_seq)} bp sequence ({round(self.seqLength/len(raw_seq)*100,1)}%) with selection range from {self.b} to {self.e}")
            #calling trimSeq to trim the ends of the Sanger sequence
            #goodTrim = self.trimSeq()
            return True







    #trimSeq() trims both ends of the Sanger sequence since the sequencing is
    #   never clear (usually messy) at the beginning and the end.  Trims based
    #   on phred_quality (10 bases in a row with 99% base calling accuracy)
    def trimSeq(self, filetype="ab1"):
        if filetype == "fasta" or filetype == "fa":
            self.seq = "".join(self.seq)
            return True
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
        seq = list(self.seq)


        repeatStart = 0
        maxCount = 0

        tempCounter = 0
        inRepeat = False

        repeatForw = ['T','C']
        repeatRev = ['G', 'A']

        index=0


        #loops through entire sequence once
        while True:
            #once entire sequence has been looked at, exit loop
            if index >= self.seqLength:
                #if the last part of the sequence was a repeat region, update stats
                #  if appropriate
                if inRepeat and tempCounter > maxCount:
                    maxCount = tempCounter
                    repeatStart = index - (maxCount*3)

                break


            #if its a TC_ or AG_ part of the sequence (in a potential repeat region)
            if (self.forwardSeq and seq[index:index+2]==repeatForw) or \
                (not self.forwardSeq and seq[index:index+2]==repeatRev):

                #if found previous consecutive TC_ and AG_, keep adding onto how
                #  many and go to the next three bases (next maybe repeat)
                if inRepeat:
                    tempCounter += 1
                    index += 3

                #meaning start of a new potential repeat region
                else:
                    tempCounter = 1
                    index+=3
                    inRepeat=True

            #if it's the end of the repeat region
            elif inRepeat:
                inRepeat = False

                #if that was the longest region found, keep track of where it was
                if tempCounter > maxCount:
                    maxCount = tempCounter
                    repeatStart  = index - (maxCount*3)


                index+=1

            #if not in a repeat region
            else:
                index+=1


        self.repeatStarts = repeatStart
        self.repeatEnds = repeatStart + (maxCount*3)


        #if it's the reverse sequence, want to subtract one since the reverse sequence
        #  reads _GA and looking for the 'GA' portions (meaning you want to start at the
        #  _ part)
        if not self.forwardSeq:
            self.repeatStarts -= 1
            self.repeatEnds -= 1







    #determineRepeats() finds the number of repeats found in the gp60 sequence
    #   Looks for TCA, TCG, TCT, and R (defined below).  Will need to modify
    #   this function if different repeat regions are required to be found.
    def determineRepeats(self):
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
            numAcatca = 0   #C. parvum R
            numHomR = 0     #C. hominis R= A(A/G)(A/G)ACGGTGGTAAGG (15bp minisatellite)
            numHomR2 = 0    #C. hominis R=A(A/G)(A/G)ACGGTGAAGG (13bp minisatellite)
            numHomfR = 0    #C. hominis If R=AAGAAGGCAAAGAAG

            #after repeat region sequence, looking for R

            afterRepeat = copy.copy(self.seq[self.repeatEnds:])
            if (afterRepeat.find("ACATCA")==0) : #returns -1 if substring is not found
                numAcatca = afterRepeat.count("ACATCA")
                self.repeatLen = self.repeatLen+numAcatca*len("ACATCA")

            elif len(self.species.split("|")) >= 2:
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
                    numHomfR = afterRepeat.count("AAGAAGGCAAGAGAAG") + afterRepeat.count("CAGAAGGCAAAGAAG")+ afterRepeat.count("AAGAGGGCAAAGAAG") + afterRepeat.count("AAGAGGGCAGTGAAG")
                self.repeatLen = self.repeatLen+numHomR*len("AAAACGGTGGTAAGG")+numHomR2*len("AAAACGGTGAAGG")    

            #repeat region sequence
            seq = copy.copy(self.seq[self.repeatStarts:self.repeatEnds])
            LOG.info(f"Found repeat {self.repeatLen} bp sequence in {len(self.seq)} bp {self.name} with coordinates {self.repeatStarts+1}-{self.repeatStarts+self.repeatLen}: {self.seq[self.repeatStarts:self.repeatStarts+self.repeatLen]}")

            numTCA = seq.count("TCA")
            numTCG = seq.count("TCG")
            numTCT = seq.count("TCT")

            repeat=""
            family = ""
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
                if family == "Ib" or family == "Ie" or family == "IV" or family == "VIIa" or family == "XIa" or family == "XIVa":
                    repeat += "T" + str(numTCT)
                else:
                    goodRepeat = False
            
            if (numAcatca!=0) and family == "IIa":
                repeat += "R" + str(numAcatca)

            if (numHomR != 0) and family =="Ia":
                repeat += "R" + str(numHomR)

            if (numHomfR != 0) and family == "If":
                repeat += "R" + str(numHomfR)

            if not goodRepeat:
                repeat = ""

            self.repeats = repeat
            LOG.info(f"Final repeat region standard nomenclature encoded value is {self.repeats}")
            return True

        else:
            self.repeats = ""
            return False


    def determineFamily(self,customdatabsename):
        LOG.info(f"Determine family and species of {self.name} of {len(self.seq)}bp sequence ...")
        # Filename to write
        filename = "query.txt"

        # Open the file with writing permission
        myfile = open(filename, 'w')

        # Write a line to the file
        if self.repeatEnds != 0:
            myfile.write(''.join(self.seq[self.repeatEnds:]))
        else:
            myfile.write(''.join(self.seq))


        # Close the file
        myfile.close()

        if customdatabsename:
            blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn', query="query.txt", dust='yes',
                                                 db="custom_db",
                                                 reward=1, penalty=-2, gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="gp60result.xml")
        else:
            blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn',query="query.txt", dust='yes',
                                             db=os.path.dirname(__file__)+"/reference_database/gp60_ref.fa", 
                                             reward=1, penalty=-2,gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="gp60result.xml")
        LOG.info(f"Querying {self.name} against the {os.path.basename(blastn_cline.db)} gp60 database ...")
        stdout, stderr = blastn_cline()

        if (os.stat("gp60result.xml").st_size == 0):
            #self.tabfile.write("Could not determine sequence. Check manually")
            self.species="No blast hits."

        else:
            result_handle = open("gp60result.xml", 'r')
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)

            if len(blast_record.alignments) > 0:
                if len(blast_record.alignments) >= 2:

                    br_alignment1 = blast_record.alignments[0]
                    hsp1 = br_alignment1.hsps[0]

                    br_alignment2 = blast_record.alignments[1]
                    hsp2 = br_alignment2.hsps[0]

                    if hsp1.score == hsp2.score or hsp2.align_length-hsp1.align_length >=25:
                        query_from1 = hsp1.query_start - 1
                        query_to1 = hsp1.query_end - 1
                        sbjctseq1 = hsp1.sbjct
                        sbjct_from1 = hsp1.sbjct_start - 1
                        sbjct_to1 = hsp1.sbjct_end - 1
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
                #exit("No BLAST hits found for the input sequence from {}".format(file))
                self.species="No blast hits."
                #return


    def blast(self, sequence, contig, filetype):
        # Filename to write
        filename = "query.txt"

        # Open the file with writing permission
        myfile = open(filename, 'w')

        # Write a line to the file
        myfile.write(sequence)

        # Close the file
        myfile.close()

        
        blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn',query="query.txt", dust='yes',
                                             db=os.path.dirname(__file__)+"/reference_database/blast_gp60.fasta", 
                                             reward=1, penalty=-2,gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="gp60result.xml")
        
        LOG.info(f"Running BLASTN on {len(sequence)}bp sequence from {self.name} on {os.path.basename(blastn_cline.db)} and contig value is {contig}")
        stdout, stderr = blastn_cline()

        if (os.stat("gp60result.xml").st_size == 0):
            return "","","","","","",""

        elif not contig:
            result_handle = open("gp60result.xml", 'r')
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)

            if len(blast_record.alignments) == 0:
                return "","","","","","",""

            br_alignment = blast_record.alignments[0]
            hsp = br_alignment.hsps[0]


            subjectSeq = hsp.sbjct

            file = open("align.fa", 'w')
            file.write("\n>Seq1\n")
            file.write(sequence)
            file.write("\n>Seq2\n")
            file.write(subjectSeq)
            file.close()

            clustalw_cline = ClustalwCommandline("clustalw", infile="align.fa")
            stdout, stderr = clustalw_cline()

            align = AlignIO.read("align.aln", "clustal")

            for record in align:
                if record.id == "Seq1":
                    query = record.seq
                else:
                    subject = record.seq

            

            if filetype == "ab1" or filetype == "abi":
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
            myfile.write(''.join(self.seq[self.repeatEnds:]))


            # Close the file
            myfile.close()

            blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn',query="query.txt", dust='yes',
                                                 db=os.path.dirname(__file__)+"/reference_database/gp60_ref.fa", reward=1, penalty=-2,gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="gp60result.xml")

            stdout, stderr = blastn_cline()

            if (os.stat("gp60result.xml").st_size == 0):
                return "","","","","","",""

            else:
                result_handle = open("gp60result.xml", 'r')
                blast_records = NCBIXML.parse(result_handle)
                blast_record = next(blast_records)

                if len(blast_record.alignments) == 0:
                    return "","","","","","",""

                br_alignment = blast_record.alignments[0]
                hsp = br_alignment.hsps[0]

                percent_identity = round(hsp.identities/hsp.align_length,3)*100
                evalue = hsp.expect
                bitscore = hsp.score

                if hsp.align_length <= blast_record.query_length :
                    query_coverage = round(hsp.align_length/blast_record.query_length*100)
                else:
                    query_coverage = 100
                query_length = blast_record.query_length
                accession = blast_record.alignments[0].hit_id

                if "|" in accession:
                    accession = accession.split("|")[-1].split("(")[1].split(")")[0]

                return bitscore,evalue,query_coverage,query_length,percent_identity, accession, self.seq

                '''if evalue > 1e-60 or query_coverage < 80:
                    s = newseq[hsp.query_start:hsp.query_end]

                    # Filename to write
                    filename = "query.txt"

                    # Open the file with writing permission
                    myfile = open(filename, 'w')

                    # Write a line to the file
                    myfile.write(s)
                    print("S")
                    print(s)

                    # Close the file
                    myfile.close()

                    blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn',query="query.txt", dust='yes',
                                                         db=os.path.dirname(__file__)+"/reference_database/gp60_ref.fa", reward=1, penalty=-2,gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="gp60result.xml")

                    stdout, stderr = blastn_cline()

                    if (os.stat("gp60result.xml").st_size == 0):
                        return "","","","","","",""

                    else:
                        result_handle = open("gp60result.xml", 'r')
                        blast_records = NCBIXML.parse(result_handle)
                        blast_record = next(blast_records)

                        if len(blast_record.alignments) == 0:
                            return "","","","","","",""

                        br_alignment = blast_record.alignments[0]
                        hsp = br_alignment.hsps[0]

                        percent_identity = round(hsp.identities/hsp.align_length,3)*100
                        evalue = hsp.expect
                        bitscore = hsp.score

                        query_coverage = round(hsp.align_length/blast_record.query_length*100)
                        query_length = blast_record.query_length
                        accession = blast_record.alignments[0].hit_id

                        if "|" in accession:
                            accession = accession.split("|")[1].split("(")[1].split(")")[0]

                        return bitscore,evalue,query_coverage,query_length,percent_identity, accession, newseq
                else:
                    return bitscore,evalue,query_coverage,query_length,percent_identity, accession, newseq'''



        else:
            # Filename to write
            filename = "query.txt"

            # Open the file with writing permission
            myfile = open(filename, 'w')

            self.seq = sequence
            self.findRepeatRegion()
            # Write a line to the file
            myfile.write(''.join(sequence[self.repeatEnds:]))


            # Close the file
            myfile.close()

            blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn',query="query.txt", dust='yes',
                                                 db=os.path.dirname(__file__)+"/reference_database/gp60_ref.fa", reward=1, penalty=-2,gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="gp60result.xml")

            LOG.info(f"Running BLASTN on filtered sequence {len(self.seq)}bp from {self.name} on {blastn_cline.db}")
            stdout, stderr = blastn_cline()

            if (os.stat("gp60result.xml").st_size == 0):
                return "","","","","","",""

            result_handle = open("gp60result.xml", 'r')
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)

            if len(blast_record.alignments) == 0:
                return "","","","","","",""

            br_alignment = blast_record.alignments[0]
            hsp = br_alignment.hsps[0]

            percent_identity = round(hsp.identities/hsp.align_length,3)*100
            evalue = hsp.expect
            bitscore = hsp.score
            if hsp.align_length <= blast_record.query_length :
                query_coverage = round(hsp.align_length/blast_record.query_length*100)
            else:
                query_coverage = 100
            query_length = blast_record.query_length
            accession = blast_record.alignments[0].hit_id

            if "|" in accession:
                accession = accession.split("|")[-1].split("(")[1].split(")")[0]
            LOG.info(f"{self.name} BLAST results on {len(sequence[self.repeatEnds:])} bp sequence after repeat region {query_coverage}% coverage, {percent_identity}% identity, top hit accession {accession}")    

            return bitscore,evalue,query_coverage,query_length,percent_identity, accession, sequence



    #printFasta() prints the results of the alignment and the repeat
    #   region for each sample. It then prints the fasta sequence for
    #   each sample that can successfully be analyzed
    def printFasta(self, contig, mode, sampleName, filetype="ab1"):
        LOG.info(f"Running printFasta() on {self.name} and current {self.species} and repeat encoding {self.repeats} ...")
        #IF WANTING TO PRINT TO SCREEN INSTEAD
        if TESTING and contig!="":
            print(">Contig: ", end="")
            print(self.name, " | ", end="")
            print(self.species.split("(")[0], end="")
            print(" ", self.repeats, end="")
            print()
            print(contig)

        self.file.write("\n>" + sampleName)

        #Output sample name
        #self.tabfile.write(self.name.split(".ab1")[0] + "\t")
        self.tabfile.write(sampleName + "\t")

        #Output type of sequences
        self.tabfile.write(mode + "\t")
        
        #for forward or reverse inputs
        if contig == "":
            sequence=str(self.seq)
            if sequence != "Poor Sequence Quality":
                bitscore,evalue,query_coverage,query_length,percent_identity, accession, seq = self.blast(sequence,False, filetype)

        else:
            sequence = contig

            if sequence != "Poor Sequence Quality":
                LOG.info("Running BLAST on contig of acceptable quality and getting top hit accession ...")
                bitscore,evalue,query_coverage,query_length,percent_identity, accession, seq = self.blast(sequence, True)

        if self.seq == "Poor Sequence Quality":
            self.tabfile.write("\t\t\tPoor Sequence Quality. Check manually.\t" + str(self.averagePhredQuality) + "\t\t\t\t\t\t\n")
            self.file.write(" | Poor Sequence Quality (Average Phred Quality = " + str(self.averagePhredQuality) + "). Check manually.")

        elif seq == "" or self.species == "No blast hits.":
            self.tabfile.write("\t\t\tNo blast hits. Check manually.\t" + str(self.averagePhredQuality)+ "\t\t\t\t\t\t\n")
            self.file.write("\n;No blast hits. Check manually. Sequence Phred Quality is " + str(self.averagePhredQuality))

        elif evalue > 1e-75 or query_coverage < 50:
            self.tabfile.write("\t\t\tCould not determine species from chromatogram. Check manually.\t" + str(self.averagePhredQuality)+ "\t\t\t\t\t\t\n")
            self.file.write("\n;Could not determine species from chromatogram. Check manually. Sequence Phred Quality is "+ str(self.averagePhredQuality))

        else:
            self.seq = seq
            self.determineFamily("")


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

                self.findRepeatRegion()
                if self.checkRepeatManually and self.doublePeaksinRepeat: #self.repeats == "Could not classify repeat region. Check manually." or (
                    foundRepeat = False
                elif self.repeatStarts < 3:
                    foundRepeat = False
                else:
                    foundRepeat = self.determineRepeats()

                if foundRepeat:
                    self.tabfile.write(self.repeats)
                    self.file.write(self.repeats)



                #Still outputting repeats if there's a subfamily
                if len(speciesName.split("|")) == 3 and foundRepeat:
                    self.file.write(speciesName.split("|")[2])
                    self.tabfile.write(speciesName.split("|")[2])

            #Output sequence
            self.tabfile.write("\t" + seq)
            self.file.write("\n" + seq)


            #Quality check in comments
            if not foundRepeat:
                self.tabfile.write("\tCould not classify repeat region. Check manually.\t")

            elif len(self.species.split("|")) > 1 and "If" in self.species.split("|")[1]:
                    begin = seq.find("CACCCCAACTC")
                    end = seq.find("AAGCTACTCCAAAGG")
                    if begin != -1 and end != -1:
                        begin += 10

                    insert = end - begin


                    if insert > 1 and self.checkRepeatManually:
                        self.tabfile.write("\t" + "Note: (1) Not all bases in repeat region had phred quality >= 20. (2) A " + str(insert) + "bp insert was found.\t")
                    elif insert > 1:
                        self.tabfile.write("\t" + "Note: A " + str(insert) + "bp insert was found.\t")
                    elif percent_identity < 99.1:
                        self.tabfile.write("\t" + "BLAST percent identity less than 99.1%. Check manually in case of new gp60 family.\t")
                    else:
                        self.tabfile.write("\tN/A\t")



            elif self.checkRepeatManually:
                self.tabfile.write("\t" + "Note: Not all bases in repeat region had phred quality >= 20.\t")
            elif percent_identity < 99.1:
                self.tabfile.write("\t" + "BLAST percent identity less than 99.1%. Check manually in case of new gp60 family.\t")
            else:
                self.tabfile.write("\tN/A\t")



            #write the average phred quality of the sequence
            #if its a contig, will output as #(F), #(R) to show phred of both
            self.tabfile.write(str(self.averagePhredQuality) + "\t")
            self.tabfile.write(str(bitscore)+"\t")
            self.tabfile.write(str(query_length)+"\t")
            self.tabfile.write(str(query_coverage) + "%\t")
            self.tabfile.write(str(evalue) +"\t")
            self.tabfile.write(str(percent_identity) + "%\t")
            self.tabfile.write(str(accession)+"\n")



#Function to build the contig of forward and reverse sequences
def buildContig(s1, s2):
    LOG.info("Building contig of forward and reverse sequences ...")
    if s1 == "":
        return s2
    elif s2 == "":
        return s1
    else:

        file = open("align.fa", 'w')
        file.write("\n>Seq1\n")
        file.write(s1)
        file.write("\n>Seq2\n")
        file.write(s2)
        file.close()

        clustalw_cline = ClustalwCommandline("clustalw", infile="align.fa")
        stdout, stderr = clustalw_cline()

        align = AlignIO.read("align.aln", "clustal")

        for record in align:
            if record.id == "Seq1":
                seq1 = record.seq
            else:
                seq2 = record.seq

        contig=""
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

        return contig

def getFileType(path):
    #pathtypes = [filetype  for path in pathlist for filetype in FILETYPES if path.endswith(filetype)]
    #print(pathtypes)
    filetype = [filetype for filetype in definitions.FILETYPES if path.endswith(filetype) ]
    if filetype and len(filetype) == 1:
        filetype = filetype[0]
        if 'ab1' == filetype:
            filetype = "abi"
        return filetype 
    else:
        return None


def gp60_main(pathlist, fPrimer, rPrimer, typeSeq, expName, customdatabsename, noheader):
    fPrimer = fPrimer.replace(' ', '')
    rPrimer = rPrimer.replace(' ', '')

    pathlist = [path for path in pathlist if re.search("|".join(definitions.FILETYPES),path)]

    if fPrimer and rPrimer:
        pathlist = [path for path in pathlist if re.search(fPrimer, path) or re.search(rPrimer, path)]  # select only files matching the primers
    elif fPrimer:
        pathlist = [path for path in pathlist if re.search(fPrimer, path)]
    elif rPrimer:
        pathlist = [path for path in pathlist if re.search(rPrimer, path)]

    
   
    if pathlist == []:
        LOG.error(f"No supported input files found in {pathlist}. Supported input filetypes are ab1, fasta, fastq")
        exit()
    
    pathlist.sort()
    pathlistEnumerated = [f"{idx+1}: {i}" for idx, i in enumerate(pathlist)]
    list_of_files_str = "\n".join(pathlistEnumerated)
    LOG.info(f"Processing {len(pathlist)} file(s):\n{list_of_files_str}")

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
    #Write output fasta with comments
    file.write("\n;>****************************************************************************")
    file.write("\n;>gp60 SEQUENCE ANALYSIS INPUT PARAMETERS:")
    if customdatabsename:
        file.write("\n  ;>Reference File: " + customdatabsename)
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
        if len(pathlist)%2 == 0:
            for idx, path in enumerate(pathlist):
                LOG.info(f"Analyzing {os.path.basename(pathlist[idx])} and {os.path.basename(pathlist[idx+1])} in contig mode")
                forward = analyzingGp60() #forward read object
                reverse = analyzingGp60() #reverse read object


                
                #for i in range(0, len(fPrimers)):
                #if fPrimer in path:
                forwSeqbool = forward.readFiles(path, True, file, tabfile)
                revSeqbool = reverse.readFiles(pathlist[idx+1], False, file, tabfile)
                pathlist.remove(pathlist[idx+1])    

                forwardPhred = forward.averagePhredQuality
                reversePhred = reverse.averagePhredQuality



                forward.averagePhredQuality = str(forwardPhred) + "(F); " + str(reversePhred) + "(R)"

                if forwSeqbool and revSeqbool:
                    goodTrimF = forward.trimSeq()
                    goodTrimR = reverse.trimSeq()
                    #print(dir(forward),dir(reverse),customdatabsename)
                    #print(f"forward seq top hit blast {forward.species} reverse {reverse.species}")
                    #'a', 'averagePhredQuality', 'b', 'beginSeq', 'blast', 'c', 'checkRepeatManually', 'determineFamily', 'determineRepeats', 'doublePeaksinRepeat', 'e', 'endSeq', 'file', 'findRepeatRegion', 'fixN', 'forwardSeq', 'g', 'name', 'oldseq', 'peakLoc', 'phred_qual', 'printFasta', 'readFiles', 'repeatEnds', 'repeatStarts', 'repeats', 'seq', 'seqLength', 'species', 't', 'tabfile', 'trimSeq'

                    if not goodTrimF and not goodTrimR:
                        forward.repeats="Could not classify repeat region. Check manually."
                        reverse.repeats="Could not classify repeat region. Check manually."

                    else:
                        
                        if customdatabsename:
                            forward.determineFamily(customdatabsename)
                            reverse.determineFamily(customdatabsename)


                        if forward.species == reverse.species and forward.species!= "":
                            LOG.info(f"Determining repeats for {forward.name} and {reverse.name} ...")
                            forward.determineRepeats()
                            reverse.determineRepeats()

                    
                    #if the crypto subfamily was found, find repeat region
                    if forward.species == reverse.species and forward.species != "" and forward.repeats == reverse.repeats and forward.repeats != "":
                        #forward.determineFamily(customdatabsename)
                        #reverse.determineFamily(customdatabsename)
                        
                        Fbitscore,Fevalue,Fquery_coverage,Fquery_length,Fpercent_identity, Faccession, Fnewseq = forward.blast(str(forward.seq), False)
                        Rbitscore,Revalue,Rquery_coverage,Rquery_length,Rpercent_identity, Raccession, Rnewseq = reverse.blast(str(reverse.seq), False)
                        LOG.info(f"Build contig from forward and reverse extracted sequences of {len(Fnewseq)}bp and {len(Rnewseq)}bp")
                        contig = buildContig(Fnewseq, Rnewseq)

                        sampleName = forward.name.split(".ab1")[0] + ", " + reverse.name.split(".ab1")[0]

                        forward.printFasta(contig, "contig", sampleName)

                    else:
                        #forward.determineFamily(customdatabsename)
                        #reverse.determineFamily(customdatabsename)
                        forward.printFasta("", "forward", forward.name.split(".ab1")[0])
                        reverse.printFasta("", "reverse", reverse.name.split(".ab1")[0])
                    
                elif forwSeqbool and not revSeqbool:
                    goodTrim = forward.trimSeq()

                    if not goodTrim:
                        forward.repeats="Could not classify repeat region. Check manually."

                    else:
                        forward.determineFamily(customdatabsename)

                    forward.printFasta("", "forward", forward.name.split(".ab1")[0])


                    reverse.seq = "Poor Sequence Quality"
                    reverse.printFasta("", "reverse", reverse.name.split(".ab1")[0])

                elif revSeqbool and not forwSeqbool:
                    forward.seq = "Poor Sequence Quality"
                    forward.printFasta("", "forward", forward.name.split(".ab1")[0])

                    goodTrim = reverse.trimSeq()

                    if not goodTrim:
                        reverse.repeats="Could not classify repeat region. Check manually."

                    else:
                        reverse.determineFamily(customdatabsename)

                    reverse.printFasta("", "reverse", reverse.name.split(".ab1")[0])

                else:
                    forward.seq = "Poor Sequence Quality"
                    forward.printFasta("", "contig", forward.name.split(".ab1")[0] + ", " + reverse.name.split(".ab1")[0])    

        else:
            print("ERROR: Uneven number of input files ({}). "
                  "Cannot find all paired forward and reverse files. Aborting ...".format(len(pathlist)))
            tabfile.write("\nCannot find all paired forward and reverse files.  Make sure all files are included to produce the contig.")
            file.write("\nCannot find all paired forward and reverse files.  Make sure all files are included to produce the contig.")

    else:
        for path in pathlist:
            filetype = getFileType(path)
            LOG.info(f"Running {path} as {filetype} file type")
            forward = analyzingGp60()
                
            if onlyForwards:
                read_ok = forward.readFiles(path, True, file, tabfile, filetype)
            elif onlyReverse:
                read_ok = forward.readFiles(path, False, file, tabfile, filetype)
         
            if read_ok == False:
                forward.seq = "Poor Sequence Quality"
                forward.printFasta("", typeSeq, forward.name.split(f".{filetype}")[0])
            else:
                goodTrim = forward.trimSeq(filetype)
                if goodTrim == False:
                    forward.repeats="Could not classify repeat region. Check manually."
                else:
                    forward.determineFamily(customdatabsename)


                forward.printFasta("", typeSeq, forward.name.split(f".{filetype}")[0], filetype)

        LOG.info(f"Finished analyzing sequence {path} ...")

    experimentName = expName + "_"

    output_report_file_name = experimentName+'cryptogenotyper_report.fa'
    filename = os.path.join('.', output_report_file_name)
    
    with open(filename, 'w') as resultFile:
        resultFile.write(file.getvalue())

    output_tabreport_file_name = experimentName+'cryptogenotyper_report.txt'
    tab_filename = os.path.join('.', output_tabreport_file_name)

    with open(tab_filename, 'w') as tabFile:
        tabFile.write(tabfile.getvalue())

    print("Fasta report written to " + os.getcwd()+"/"+output_report_file_name + ".")
    print("Tab-delimited report written to " + os.getcwd() + "/" + output_tabreport_file_name + ".\nThe gp60 run completed successfully")

   # os.system("rm gp60result.xml gp60result2.xml")
    #os.system("rm query.txt")
    #os.system("rm align.fa align.dnd align.aln")

if __name__ == "__main__":
	gp60_main()
