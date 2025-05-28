import io
import os
import re
import logging
from CryptoGenotyper import definitions, utilities
from CryptoGenotyper.logging import create_logger

from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO


from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

import numpy as np
import math
import copy


global minAmp
minAmp = 0


global TESTING
TESTING = False

# setup the application logging
#LOG = create_logger(__name__)
LOG = logging.getLogger(__name__)

##############################################################################
# IUPAC Code:
#	R = A or G
#	Y = C or T
#	S = G or C
#	W = A or T
#	K = G or T
#	M = A or C
#	B = C or G or T (not A)
# 	D = A or G or T (not C)
#	H = A or C or T (not G)
# 	V = A or C or G (not T)
#	N = any
##############################################################################


# revcomp() returns the reverse complement of the DNA sequence given
def revcomp(sequence):
    l = len(sequence)
    rev = ""

    for i in range(l-1, -1, -1):
        if sequence[i] == 'G':
            rev += 'C'
        elif sequence[i] == 'C':
            rev += 'G'
        elif sequence[i] == 'T':
            rev += 'A'
        elif sequence[i] == 'A':
            rev += 'T'

    return rev


#iupac() determines the IUPAC base when there are multiple possibilities.  By
#   following the IUPAC code, it returns the new possible base depending on the
#   current base and the base that is wanting to be added
def addIupac(currBase, base):
    if base == 'A':
        if currBase == 'C':
            return 'M'
        elif currBase == 'G':
            return 'R'
        elif currBase == 'T':
            return 'W'
        elif currBase == 'Y':
            return 'H'
        elif currBase == 'S':
            return 'V'
        elif currBase == 'K':
            return 'D'
        elif currBase == 'B':
            return 'N'

    elif base == 'C':
        if currBase == 'G':
            return 'S'
        elif currBase =='T':
            return 'Y'
        elif currBase == 'A':
            return 'M'
        elif currBase == 'R':
            return 'V'
        elif currBase == 'W':
            return 'H'
        elif currBase == 'K':
            return 'B'
        elif currBase == 'D':
            return 'N'

    elif base == 'G':
        if currBase == 'A':
            return 'R'
        elif currBase == 'C':
            return 'S'
        elif currBase == 'T':
            return 'K'
        elif currBase == 'Y':
            return 'B'
        elif currBase == 'W':
            return 'D'
        elif currBase == 'M':
            return 'V'
        elif currBase == 'H':
            return 'N'

    elif base == 'T':
        if currBase == 'A':
            return 'W'
        elif currBase == 'C':
            return 'Y'
        elif currBase == 'G':
            return 'K'
        elif currBase == 'R':
            return 'D'
        elif currBase == 'S':
            return 'B'
        elif currBase == 'M':
            return 'H'
        elif currBase == 'V':
            return 'N'



#subIupac() subtracts the base from the current IUPAC base and returns
#  the new IUPAC base
def subIupac(currBase, base):
    if currBase == "A" or currBase == "T" or currBase == "G" or currBase == "C":
        return currBase

    elif currBase == 'R':
        if base == 'A':
            return 'G'
        elif base == 'G':
            return 'A'
        else:
            return currBase

    elif currBase == 'Y':
        if base == 'C':
            return 'T'
        elif base == 'T':
            return 'C'
        else:
            return currBase

    elif currBase == 'S':
        if base == 'G':
            return 'C'
        elif base == 'C':
            return 'G'
        else:
            return currBase

    elif currBase == 'W':
        if base == 'A':
            return 'T'
        elif base == 'T':
            return 'A'
        else:
            return currBase

    elif currBase == 'K':
        if base == 'G':
            return 'T'
        elif base == 'T':
            return 'G'
        else:
            return currBase

    elif currBase == 'M':
        if base == 'A':
            return 'C'
        elif base == 'C':
            return 'A'
        else:
            return currBase

    elif currBase == 'N':
        if base == 'A':
            return 'B'
        elif base == 'C':
            return 'D'
        elif base == 'G':
            return 'H'
        elif base == 'T':
            return 'V'

    elif currBase == 'B':
        if base == 'T':
            return 'R'
        elif base == 'G':
            return 'Y'
        elif base == 'C':
            return 'K'

    elif currBase == 'D':
        if base == 'A':
            return 'K'
        elif base == 'G':
            return 'W'
        elif base == 'T':
            return 'R'

    elif currBase == 'H':
        if base == 'A':
            return 'Y'
        elif base == 'C':
            return 'W'
        elif base == 'T':
            return 'M'

    elif currBase == 'V':
        if base == 'A':
            return 'S'
        elif base == 'C':
            return 'R'
        elif base == 'G':
            return 'M'

    else:
        return currBase


# indeligent() algorithm converted into Python Code
#   Developed by: Dmitry Dmitriev & Roman Rakitov (2008)
#   Citation: Dmitriev DA, Rakitov RA. Decoding of superimposed traces produced
#             by direct sequencing of heterozygous indels. Plos Computational
#             Biology. 2008 ;4(7):e1000113. DOI: 10.1371/journal.pcbi.1000113.
def indelligent(gen1):

    sec1 = 0
    sec2 = 2
    sec3 = 1
    l = 10
    ind = 1
    al = "1"
    indM = 15
    subs = 0
    rep0 = 0
    fixind = ""
    ambig = "1"
    ambig1 = 0
    ambig2=0
    align2 = "1"
    indpen=2
    longind = "0"
    bd=0

    w, h = 10000, 3
    MX1 = [["" for x in range(h)] for y in range(w)]

    w, h, d = 10000, 200,3
    MX2 = [[[0 for x in range(d)] for y in range(h)]for z in range(w)]

    w, h, d = 10000, 200,3
    MX3 = [[[0 for x in range(d)] for y in range(h)]for z in range(w)]

    w, h = 10000, 5
    MX4 = [[0 for x in range(h)] for y in range(w)]

    w, h = 5, 3
    mx5 = [["" for x in range(h)] for y in range(w)]


    iMax = len(gen1)+1
    SC1 = 1
    SC2 = 0
    SC3 = indpen


    #mx1 matrix stores pairs of superimposed base calls at each position of analyzed fragment
    for i in range(1, iMax):
        j=gen1[i-1]
        if j == "R":
            MX1[i][1] = "A"
            MX1[i][2] = "G"
        elif j == "Y":
            MX1[i][1] = "C"
            MX1[i][2] = "T"
        elif j == "S":
            MX1[i][1] = "G"
            MX1[i][2] = "C"
        elif j == "W":
            MX1[i][1] = "A"
            MX1[i][2] = "T"
        elif j == "M":
            MX1[i][1] = "A"
            MX1[i][2] = "C"
        elif j == "K":
            MX1[i][1] = "G"
            MX1[i][2] = "T"
        else:
            MX1[i][1] = j
            MX1[i][2] = j

        if j == "B" or j == "D" or j == "H" or j == "V" or j == "N":
            bd = 1
        if j!="A" and j!="C" and j!="G" and j!="T":
            ambig1 += 1

    estimate = 0
    m3 = 0
    m4 = 0

    #limit maximum phase shift magnitude to 1/2 length of sequence
    if len(gen1) <= indM * 2:
        n1 = int(len(gen1) / 2)
    else:
        n1 = indM

    n1 += 1
    #matrices of scores][initial conditions
    for j in range (0,n1):
        for i in range(0,j):
            MX2[i][j][1] = i
            MX2[i][j][2] = i

    for j in range (0,n1):
        for i in range (0,j):
            MX3[iMax+1-i][j][1] = i
            MX3[iMax+1-i][j][2] = i

    for i in range(1,iMax):
        jMax = n1
        if i - 1 < n1 :
            jMax = i - 1
        for j in range(0,jMax):
            for z in range(1,3):
                SCMax = -111111
                for x in range(0,jMax): #or n1
                    if j == 0 and x == 0 and MX1[i][1] != MX1[i][2] :
                        SC = max(MX2[i - 1][x][1],MX2[i - 1][x][2]) - SC2
                    elif j == 0 and MX1[i][1] == MX1[i][2] :
                        SC = max(MX2[i - 1][x][1],MX2[i - 1][x][2]) + SC1
                    elif j != x :
                        if j < x and ((MX1[i][3 - z] == MX1[i - j][1] and MX2[i - j][x][1] >= MX2[i - j][x][2]) or (MX1[i][3 - z] == MX1[i - j][2] and MX2[i - j][x][1] <= MX2[i - j][x][2])) :
                            SC = max(MX2[i - 1][x][1],MX2[i - 1][x][2]) + SC1
                        elif j > x :
                            SC = max(MX2[i - 1][x][1],MX2[i - 1][x][2]) + SC1
                        else:
                            SC = max(MX2[i - 1][x][1],MX2[i - 1][x][2]) - SC2
                    elif (MX1[i][3 - z] == MX1[i - j][1] and MX2[i - j][j][1] >= MX2[i - j][j][2]) or (MX1[i][3 - z] == MX1[i - j][2] and MX2[i - j][j][1] <= MX2[i - j][j][2]) :
                        SC = max(MX2[i - 1][j][1],MX2[i - 1][j][2]) + SC1
                    else:
                        SC = max(MX2[i - 1][j][1],MX2[i - 1][j][2]) - SC2

                    if x != j :
                        SC = SC - SC3 - abs(x - j)
                    if SC > SCMax :
                        SCMax = SC

                MX2[i][j][z] = SCMax

    #Calculate scores in the reverse direction (MX3)

    for i in range(iMax, 1, -1):
        jMax = n1
        if iMax - i < n1 :
            jMax = iMax - i
        for j in range(0, jMax):
            for z in range(1,3):
                SCMax = -111111
                for x in range(0,jMax):
                    if j == 0 and x == 0 and MX1[i][1] != MX1[i][2]:
                        SC = max(MX3[i + 1][x][1],MX3[i + 1][x][2]) - SC2
                    elif j == 0 and MX1[i][1] == MX1[i][2] :
                        SC = max(MX3[i + 1][x][1],MX3[i + 1][x][2]) + SC1
                    elif j != x :
                        if j < x and (MX1[i][z] == MX1[i + j][2] and MX3[i + j][x][1] >= MX3[i + j][x][2]) or (MX1[i][z] == MX1[i + j][1] and MX3[i + j][x][1] <= MX3[i + j][x][2]) :
                            SC = max(MX3[i + 1][x][1],MX3[i + 1][x][2]) + SC1
                        elif j > x :
                            SC = max(MX3[i + 1][x][1],MX3[i + 1][x][2]) + SC1
                        else:
                            SC = max(MX3[i + 1][x][1],MX3[i + 1][x][2]) - SC2

                    elif (MX1[i][z] == MX1[i + j][2] and MX3[i + j][j][1] >= MX3[i + j][j][2]) or (MX1[i][z] == MX1[i + j][1] and MX3[i + j][j][1] <= MX3[i + j][j][2]) :
                        SC = max(MX3[i + 1][j][1],MX3[i + 1][j][2]) + SC1

                    else:
                        SC = max(MX3[i + 1][j][1],MX3[i + 1][j][2]) - SC2

                    if x != j :
                        SC = SC - SC3 - abs(x - j)
                    if SC > SCMax :
                        SCMax = SC

                MX3[i][j][z] = SCMax

    #Sum MX2 and MX3][store results in MX3

    for i in range (1, iMax):
    	for j in range(0, n1):
    		for z in range (1, 3):
    			MX3[i][j][z] = MX3[i][j][z] + MX2[i][j][z]

    temp1 = ""
    temp2 = ""
    indel = "-1"

    indel = "," + indel + ","
    if indel == ",,-1,," :
        indel = ",-1,"

    #Make the list of possible indels

    x = ""
    for i in range (1, iMax):
        SC = -111111
        SCc = 0
        for j in range (0, n1):
            b="," + str(j) + "j"
            if b in indel or indel ==",-1," : # fixind = "" :
                if MX3[i][j][1] > SC :
                    SC = MX3[i][j][1]
                    SCc = j
                if MX3[i][j][2] > SC :
                    SC = MX3[i][j][2]
                    SCc = j
            b="a" + str(SCc) + "a"
            if b not in x:
                x = x + "a" + str(SCc) + "a"

    # Apply scores from MX3 and resolve positions in MX1. Results are stored in MX4.
    # MX4[i][1] the value indicating which base from MX1 goes into the first string. If MX4[i][1] = 0 : the position is ambiguous.
    # MX4[i][2] the value indicating the phase shift which resolves the position.

    SCd = 0

    for i in range (1, iMax):
        SC = -111111
        SCa = 0 # possibility to solve in position 1
        SCb = 0 # possibility to solve in position 2
        SCc = 0 # shift size
        SCMax = ""
        for j in range(0,n1):
            if ("," + str(j) + "," in indel and SCd == 0) or ("a" + str(j) + "a" in x and SCd == 1) or (indel ==",-1," and SCd == 0) :
                if MX3[i][j][1] > SC :
                    SC = MX3[i][j][1]
                    SCa = 0
                    SCb = 0
                    SCc = j
                if MX3[i][j][2] > SC :
                    SC = MX3[i][j][2]
                    SCa = 0
                    SCb = 0
                    SCc = j

                if (MX3[i][j][1] == SC or MX3[i][j][2] == SC) and 0 + j == MX4[i - 1][2] :
                    SCc = j
                if MX3[i][j][1] == SC and "a" + str(j) + "a" in x:
                    SCa = 1
                if MX3[i][j][2] == SC and "a" + str(j) + "a" in x:
                    SCb = 1

        if MX1[i][1] == MX1[i][2] :
            MX4[i][1] = 1
        elif SCa == 1 and SCb == 0 :
            MX4[i][1] = 1
        elif SCb == 1 and SCa == 0 :
            MX4[i][1] = 2
        else:
            MX4[i][1] = 0

        MX4[i][2] = SCc

    # remove phase shifts recovered at the number of consecutive positions smaller than the phase shift magnitude

    u = ""
    SCa = 0

    for i in range (1, iMax):
        exitFor = False
        if MX4[i][2] != MX4[i - 1][2] :
            for z in range (1, MX4[i][2]): # max shift size that can be removed
                if i + z > iMax :
                    exitFor = True
                elif MX4[i + z][2] != MX4[i][2] :
                    SCa = 1
                    exitFor = True
                if (exitFor):
                    break

        if MX4[i][2] != MX4[1][2] :
            SCc = 1
        if z > MX4[i][2] and ("a" + str(MX4[i][2]) + "a") not in u :
            u = u + "a" + str(MX4[i][2]) + "a"

    if u != x and SCd == 0 :
        SCa = 1
    x = u
    if SCa == 1 or SCd == 1 :
        SCd = SCd + 1

    # calculate and mark the ambiguities that could potentially be resolved
    for i in range(1, iMax):
        if MX4[i][1] == 0 :
            j = abs(MX4[i][2])
            if i <= j :
                MX4[i][4] = 1 # 1 - could be resolved][0 - cannot be resolved
            elif j == 0 :
                MX4[i][4] = 1
            elif MX4[i - j][4] == 0 and abs(MX4[i - j][2]) == j and MX4[i - j][1] == 0 :
                MX4[i][4] = 0
            elif MX4[i - j][4] == 1 and abs(MX4[i - j][2]) == j and MX4[i - j][1] == 0 :
                MX4[i][4] = 1
            else:
                SCc = 1
                while MX4[j * SCc + i][1] == 0 and abs(MX4[j * SCc + i][2]) == j and j * SCc + i <= iMax:
                    SCc = SCc + 1

                if abs(MX4[i - j][2]) < j and MX4[i][2] > 0 :
                    MX4[i][4] = 1
                elif j * SCc + i > iMax :
                    MX4[i][4] = 1
                elif abs(MX4[j * SCc + i][2]) != j :
                    MX4[i][4] = 1
                elif abs(MX4[j * SCc + i + MX4[i][3]][2]) != j :
                    MX4[i][4] = 1
                elif MX4[i][2] > 0 and SCc / 2 != int(SCc / 2) and (MX1[i - j][MX4[i - j][1]] == MX1[i][2] and MX1[i][1] == MX1[j * SCc + i][3 - MX4[j * SCc + i][1]]) or (MX1[i - j][MX4[i - j][1]] == MX1[i][1] and MX1[i][2] == MX1[j * SCc + i][3 - MX4[j * SCc + i][1]]) :
                    MX4[i][4] = 1
                elif MX4[i][2] > 0 and SCc / 2 == int(SCc / 2) and (MX1[i - j][MX4[i - j][1]] == MX1[i][2] and MX1[i][2] == MX1[j * SCc + i][3 - MX4[j * SCc + i][1]]) or (MX1[i - j][MX4[i - j][1]] == MX1[i][1] and MX1[i][1] == MX1[j * SCc + i][3 - MX4[j * SCc + i][1]]) :
                    MX4[i][4] = 1
                elif abs(MX4[i - j][2]) != j and MX4[i][2] < 0 :
                    MX4[i][4] = 1
                elif MX4[i - j][1] == 0 :
                    MX4[i][4] = 0
                elif MX4[i][2] < 0 and SCc / 2 != int(SCc / 2) and (MX1[i - j][3 - MX4[i - j][1]] == MX1[i][1] and MX1[i][2] == MX1[j * SCc + i][MX4[j * SCc + i][1]]) or (MX1[i - j][3 - MX4[i - j][1]] == MX1[i][2] and MX1[i][1] == MX1[j * SCc + i][MX4[j * SCc + i][1]]) :
                    MX4[i][4] = 1
                elif MX4[i][2] < 0 and SCc / 2 == int(SCc / 2) and (MX1[i - j][3 - MX4[i - j][1]] == MX1[i][1] and MX1[i][1] == MX1[j * SCc + i][MX4[j * SCc + i][1]]) or (MX1[i - j][3 - MX4[i - j][1]] == MX1[i][2] and MX1[i][2] == MX1[j * SCc + i][MX4[j * SCc + i][1]]) :
                    MX4[i][4] = 1
                else:
                    MX4[i][4] = 0

    j = 0
    for i in range(1, iMax):
        if MX4[i][1] == 0 and MX4[i][4] == 1 :
            for z in range( 1 , n1 + MX4[i][3]):
                if z >= i or z > iMax - i :
                    break
                if MX4[i - z][2] != MX4[i + z][2] :
                    break
            if z < i and z - MX4[i][3] <= n1 and z <= iMax - i :
                ind1 = abs(MX4[i - z][2])
                ind2 = abs(MX4[i + z][2])
                if longind + "" != "1" :
                    if i >= ind1 and i >= ind2 and iMax - i > ind1 and iMax - i > ind2 :
                        if ((MX1[i][2] == MX1[i - ind1][1] and MX3[i - ind1][ind1][1] >= MX3[i - ind1][ind1][2]) or (MX1[i][2] == MX1[i - ind1][2] and MX3[i - ind1][ind1][1] <= MX3[i - ind1][ind1][2]))  and  ((MX1[i][1] == MX1[i + ind1][2] and MX3[i + ind1][ind1][1] >= MX3[i + ind1][ind1][2]) or (MX1[i][1] == MX1[i + ind1][1] and MX3[i + ind1][ind1][1] <= MX3[i + ind1][ind1][2])) :
                            mx5[1][1] = 1
                        else:
                            mx5[1][1] = 0
                        if ((MX1[i][1] == MX1[i - ind1][1] and MX3[i - ind1][ind1][1] >= MX3[i - ind1][ind1][2]) or (MX1[i][1] ==MX1[i - ind1][2] and MX3[i - ind1][ind1][1] <= MX3[i - ind1][ind1][2]))  and  ((MX1[i][2] == MX1[i + ind1][2] and MX3[i + ind1][ind1][1] >= MX3[i + ind1][ind1][2]) or (MX1[i][2] == MX1[i + ind1][1] and MX3[i + ind1][ind1][1] <= MX3[i + ind1][ind1][2])) :
                            mx5[1][2] = 1
                        else:
                            mx5[1][2] = 0
                        if ((MX1[i][2] == MX1[i - ind2][1] and MX3[i - ind2][ind2][1] >= MX3[i - ind2][ind2][2]) or (MX1[i][2] == MX1[i - ind2][2] and MX3[i - ind2][ind2][1] <= MX3[i - ind2][ind2][2]))  and  ((MX1[i][1] == MX1[i + ind2][2] and MX3[i + ind2][ind2][1] >= MX3[i + ind2][ind2][2]) or (MX1[i][1] == MX1[i + ind2][1] and MX3[i + ind2][ind2][1] <= MX3[i + ind2][ind2][2])) :
                            mx5[2][1] = 1
                        else:
                            mx5[2][1] = 0
                        if ((MX1[i][1] == MX1[i - ind2][1] and MX3[i - ind2][ind2][1] >= MX3[i - ind2][ind2][2]) or (MX1[i][1] == MX1[i - ind2][2] and MX3[i - ind2][ind2][1] <= MX3[i - ind2][ind2][2]))  and  ((MX1[i][2] == MX1[i + ind2][2] and MX3[i + ind2][ind2][1] >= MX3[i + ind2][ind2][2]) or (MX1[i][2] == MX1[i + ind2][1] and MX3[i + ind2][ind2][1] <= MX3[i + ind2][ind2][2])) :
                            mx5[2][2] = 1
                        else:
                            mx5[2][2] = 0
                        if ind1 < ind2 and ((MX1[i][2]== MX1[i - ind1][1] and MX3[i - ind1][ind1][1] >= MX3[i - ind1][ind1][2]) or (MX1[i][2] == MX1[i - ind1][2] and MX3[i - ind1][ind1][1] <= MX3[i - ind1][ind1][2]))  and  ((MX1[i][1] == MX1[i + ind2][2] and MX3[i + ind2][ind2][1] >= MX3[i + ind2][ind2][2]) or (MX1[i][1] == MX1[i + ind2][1] and MX3[i + ind2][ind2][1] <= MX3[i + ind2][ind2][2])) :
                            mx5[4][1] = 1
                        else:
                            mx5[4][1] = 0
                        if ind1 < ind2 and ((MX1[i][1] == MX1[i - ind1][1] and MX3[i - ind1][ind1][1] >= MX3[i - ind1][ind1][2]) or (MX1[i][1] == MX1[i - ind1][2] and MX3[i - ind1][ind1][1] <= MX3[i - ind1][ind1][2]))  and  ((MX1[i][2] == MX1[i + ind2][2] and MX3[i + ind2][ind2][1] >= MX3[i + ind2][ind2][2]) or (MX1[i][2] == MX1[i + ind2][1] and MX3[i + ind2][ind2][1] <= MX3[i + ind2][ind2][2])) :
                            mx5[4][2] = 1
                        else:
                            mx5[4][2] = 0
                    elif i <= ind1 :
                        if ((MX1[i][1] == MX1[i + ind1][2] and MX3[i + ind1][ind1][1] >= MX3[i + ind1][ind1][2]) or (MX1[i][1] == MX1[i + ind1][1] and MX3[i + ind1][ind1][1] <= MX3[i + ind1][ind1][2])) :
                            mx5[1][1] = 1
                        else: mx5[1][1] = 0
                        if ((MX1[i][2] == MX1[i + ind1][2] and MX3[i + ind1][ind1][1] >= MX3[i + ind1][ind1][2]) or (MX1[i][2] == MX1[i + ind1][1] and MX3[i + ind1][ind1][1] <= MX3[i + ind1][ind1][2])) :
                            mx5[1][2] = 1
                        else: mx5[1][2] = 0
                        mx5[2][1] = 0
                        mx5[2][2] = 0
                    elif iMax - i < ind2 :
                        mx5[1][1] = 0
                        mx5[1][2] = 0
                        if ((MX1[i][2] == MX1[i - ind2][1] and MX3[i - ind2][ind2][1] >= MX3[i - ind2][ind2][2]) or (MX1[i][2] == MX1[i - ind2][2] and MX3[i - ind2][ind2][1] <= MX3[i - ind2][ind2][2])) :
                            mx5[2][1] = 1
                        else: mx5[2][1] = 0
                        if ((MX1[i][1] == MX1[i - ind2][1] and MX3[i - ind2][ind2][1] >= MX3[i - ind2][ind2][2]) or (MX1[i][1] == MX1[i - ind2][2] and MX3[i - ind2][ind2][1] <= MX3[i - ind2][ind2][2])) :
                            mx5[2][2] = 1
                        else: mx5[2][2] = 0

                    if MX4[i - z][2] < MX4[i + z][2] or i <= ind1 :
                        if ((MX1[i][1] == MX1[i + ind2][2] and MX3[i + ind2][ind2][1] >= MX3[i + ind2][ind2][2]) or (MX1[i][1] == MX1[i + ind2][1] and MX3[i + ind2][ind2][1] <= MX3[i + ind2][ind2][2])) :
                            mx5[3][1] = 1
                        else: mx5[3][1] = 0
                        if ((MX1[i][2] ==MX1[i + ind2][2] and MX3[i + ind2][ind2][1] >= MX3[i + ind2][ind2][2]) or (MX1[i][2] == MX1[i + ind2][1] and MX3[i + ind2][ind2][1] <= MX3[i + ind2][ind2][2])) :
                            mx5[3][2] = 1
                        else: mx5[3][2] = 0
                    else:
                        if ((MX1[i][2] == MX1[i - ind1][1] and MX3[i - ind1][ind1][1] >= MX3[i - ind1][ind1][2]) or (MX1[i][2] == MX1[i - ind1][2] and MX3[i - ind1][ind1][1] <= MX3[i - ind1][ind1][2])) :
                            mx5[3][1] = 1
                        else: mx5[3][1] = 0
                        if ((MX1[i][1] == MX1[i - ind1][1] and MX3[i - ind1][ind1][1] >= MX3[i - ind1][ind1][2]) or (MX1[i][1] == MX1[i - ind1][2] and MX3[i - ind1][ind1][1] <= MX3[i - ind1][ind1][2])) :
                            mx5[3][2] = 1
                        else: mx5[3][2] = 0
                    if ind1 == 0 or ind2 == 0 :
                        mx5[4][1] = 0
                        mx5[4][2] = 0

                else: # long ind
                    if MX4[i - z][2] > 0 :
                        if (MX1[i][2] == MX1[i - ind1][1] and MX4[i - ind1][1] != 2) or (MX1[i][2] == MX1[i - ind1][2] and MX4[i - ind1][1] != 1):
                            mx5[3][1] = 1
                        else:
                            mx5[3][1] = 0
                        if ((MX1[i][1] == MX1[i - ind1][1] and MX4[i - ind1][1] != 2) or (MX1[i][1] == MX1[i - ind1][2] and MX4[i - ind1][1] != 1)) :
                            mx5[3][2] = 1
                        else:
                            mx5[3][2] = 0
                    else:
                        if ((MX1[i][1] == MX1[i - ind1][2] and MX4[i - ind1][1] != 2) or (MX1[i][1] == MX1[i - ind1][1] and MX4[i - ind1][1] != 1)) :
                            mx5[3][1] = 1
                        else: mx5[3][1] = 0
                        if ((MX1[i][2] == MX1[i - ind1][2] and MX4[i - ind1][1] != 2) or (MX1[i][2] == MX1[i - ind1][1] and MX4[i - ind1][1] != 1)) :
                            mx5[3][2] = 1
                        else: mx5[3][2] = 0

                    if MX4[i + z][2] > 0 :
                        if ((MX1[i][1] == MX1[i + ind2][2] and MX4[i + ind2][1] != 2) or (MX1[i][1] == MX1[i + ind2][1] and MX4[i + ind2][1] != 1)) :
                            mx5[4][1] = 1
                        else: mx5[4][1] = 0
                        if ((MX1[i][2] == MX1[i + ind2][2] and MX4[i + ind2][1] != 2) or (MX1[i][2] == MX1[i + ind2][1] and MX4[i + ind2][1] != 1)) :
                            mx5[4][2] = 1
                        else: mx5[4][2] = 0
                    else:
                        if (MX1[i][2] == MX1[i + ind2][1] and MX4[i + ind2][1] != 2) or (MX1[i][2] == MX1[i + ind2][2] and MX4[i + ind2][1] != 1) :
                            mx5[4][1] = 1
                        else: mx5[4][1] = 0
                        if (MX1[i][1] == MX1[i + ind2][1] and MX4[i + ind2][1] != 2) or (MX1[i][1] == MX1[i + ind2][2] and MX4[i + ind2][1] != 1) :
                            mx5[4][2] = 1
                        else: mx5[4][2] = 0

                    if i >= ind1 and i >= ind2 and iMax - i > ind1 and iMax - i > ind2 :
                        if MX4[i - z][2] > 0 :
                            if (MX1[i][2] == MX1[i - ind1][1] and MX4[i - ind1][1] != 2) or (MX1[i][2] == MX1[i - ind1][2] and MX4[i - ind1][1] != 1)  and  (MX1[i][1] == MX1[i + ind1][2] and MX4[i + ind1][1] != 2) or (MX1[i][1] == MX1[i + ind1][1] and MX4[i + ind1][1] != 1) :
                                mx5[1][1] = 1
                            else: mx5[1][1] = 0
                            if (MX1[i][1] == MX1[i - ind1][1] and MX4[i - ind1][1] != 2) or (MX1[i][1] == MX1[i - ind1][2] and MX4[i - ind1][1] != 1)  and  (MX1[i][2] == MX1[i + ind1][2] and MX4[i + ind1][1] != 2) or (MX1[i][2] == MX1[i + ind1][1] and MX4[i + ind1][1] != 1) :
                                mx5[1][2] = 1
                            else: mx5[1][2] = 0
                        else:
                            if (MX1[i][1] == MX1[i - ind1][2] and MX4[i - ind1][1] != 2) or (MX1[i][1] == MX1[i - ind1][1] and MX4[i - ind1][1] != 1)  and  (MX1[i][2] == MX1[i + ind1][1] and MX4[i + ind1][1] != 2) or (MX1[i][2] == MX1[i + ind1][2] and MX4[i + ind1][1] != 1) :
                                mx5[1][1] = 1
                            else: mx5[1][1] = 0
                            if (MX1[i][2] == MX1[i - ind1][2] and MX4[i - ind1][1] != 2) or (MX1[i][2] == MX1[i - ind1][1] and MX4[i - ind1][1] != 1)  and  (MX1[i][1] == MX1[i + ind1][1] and MX4[i + ind1][1] != 2) or (MX1[i][1] == MX1[i + ind1][2] and MX4[i + ind1][1] != 1) :
                                mx5[1][2] = 1
                            else: mx5[1][2] = 0
                        if MX4[i + z][2] > 0 :
                            if (MX1[i][2] == MX1[i - ind2][1] and MX4[i - ind2][1] != 2) or (MX1[i][2] == MX1[i - ind2][2] and MX4[i - ind2][1] != 1)  and  (MX1[i][1] == MX1[i + ind2][2] and MX4[i + ind2][1] != 2) or (MX1[i][1] == MX1[i + ind2][1] and MX4[i + ind2][1] != 1) :
                                mx5[2][1] = 1
                            else: mx5[2][1] = 0
                            if (MX1[i][1] == MX1[i - ind2][1] and MX4[i - ind2][1] != 2) or (MX1[i][1] == MX1[i - ind2][2] and MX4[i - ind2][1] != 1)  and  (MX1[i][2] == MX1[i + ind2][2] and MX4[i + ind2][1] != 2) or (MX1[i][2] == MX1[i + ind2][1] and MX4[i + ind2][1] != 1) :
                                mx5[2][2] = 1
                            else: mx5[2][2] = 0
                        else:
                            if (MX1[i][1] == MX1[i - ind2][2] and MX4[i - ind2][1] != 2) or (MX1[i][1] == MX1[i - ind2][1] and MX4[i - ind2][1] != 1)  and  (MX1[i][2] == MX1[i + ind2][1] and MX4[i + ind2][1] != 2) or (MX1[i][2] == MX1[i + ind2][2] and MX4[i + ind2][1] != 1) :
                                mx5[2][1] = 1
                            else: mx5[2][1] = 0
                            if (MX1[i][2] == MX1[i - ind2][2] and MX4[i - ind2][1] != 2) or (MX1[i][2] == MX1[i - ind2][1] and MX4[i - ind2][1] != 1)  and  (MX1[i][1] == MX1[i + ind2][1] and MX4[i + ind2][1] != 2) or (MX1[i][1] == MX1[i + ind2][2] and MX4[i + ind2][1] != 1) :
                                mx5[2][2] = 1
                            else: mx5[2][2] = 0
                    elif i <= ind1 :
                        if (MX1[i][1] == MX1[i + ind1][2] and MX4[i + ind1][1] != 2) or (MX1[i][1] == MX1[i + ind1][1] and MX4[i + ind1][1] != 1) :
                            mx5[1][1] = 1
                        else: mx5[1][1] = 0
                        if (MX1[i][2] == MX1[i + ind1][2] and MX4[i + ind1][1] != 2) or (MX1[i][2] == MX1[i + ind1][1] and MX4[i + ind1][1] != 1) :
                            mx5[1][2] = 1
                        else: mx5[1][2] = 0
                        mx5[2][1] = 0
                        mx5[2][2] = 0
                    elif iMax - i < ind2 :
                        mx5[1][1] = 0
                        mx5[1][2] = 0
                        if (MX1[i][2] == MX1[i - ind2][1] and MX4[i - ind1][1] != 2) or (MX1[i][2] == MX1[i - ind2][2] and MX4[i - ind1][1] != 1) :
                            mx5[2][1] = 1
                        else: mx5[2][1] = 0
                        if (MX1[i][1] == MX1[i - ind2][1] and MX4[i - ind1][1] != 2) or (MX1[i][1] == MX1[i - ind2][2] and MX4[i - ind1][1] != 1) :
                            mx5[2][2] = 1
                        else: mx5[2][2] = 0

                    if ind1 == 0 :
                        mx5[1][1] = 0
                        mx5[1][2] = 0
                    if ind2 == 0 :
                        mx5[2][1] = 0
                        mx5[2][2] = 0


    # resolve 3-fold degenerate bases

    if bd == 1 :
        for i in range(1, iMax):
            if MX1[i][1] == "B" or MX1[i][1] == "D" or MX1[i][1] == "H" or MX1[i][1] == "V" or MX1[i][1] == "N" :
                if i > MX4[i][2] and i <= iMax - MX4[i][2] :
                    if MX4[i - MX4[i][2]][1] > 0 and MX4[i + MX4[i][2]][1] > 0 :
                        if MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]] == "A" or MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]] == "C" or MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]] == "G" or MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]] == "T" :
                            if MX1[i+ MX4[i][2]][3 - MX4[i + MX4[i][2]][1]] == "A" or MX1[i+ MX4[i][2]][3 - MX4[i + MX4[i][2]][1]] == "C" or MX1[i+ MX4[i][2]][3 - MX4[i + MX4[i][2]][1]] == "G" or MX1[i+ MX4[i][2]][3 - MX4[i + MX4[i][2]][1]] == "T":
                                if MX1[i][1] == "B" and MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]] != "A" and MX1[i+ MX4[i][2]][3 - MX4[i + MX4[i][2]][1]] != "A" and MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]] != MX1[i + MX4[i][2]][3 - MX4[i + MX4[i][2]][1]] :
                                    MX1[i][1] = MX1[i + MX4[i][2]][3 - MX4[i + MX4[i][2]][1]]
                                    MX1[i][2] = MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]]
                                elif  MX1[i][1] == "D" and MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]] != "C" and MX1[i + MX4[i][2]][3 - MX4[i + MX4[i][2]][1]] != "C" and MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]] != MX1[i + MX4[i][2]][3 - MX4[i + MX4[i][2]][1]] :
                                    MX1[i][1] = MX1[i + MX4[i][2]][3 - MX4[i + MX4[i][2]][1]]
                                    MX1[i][2] = MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]]
                                elif  MX1[i][1] == "H" and MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]] != "G" and MX1[i+ MX4[i][2]][3 - MX4[i + MX4[i][2]][1]] != "G" and MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]] != MX1[i + MX4[i][2]][3 - MX4[i + MX4[i][2]][1]] :
                                    MX1[i][1] = MX1[i +  MX4[i][2]][3 - MX4[i + MX4[i][2]][1]]
                                    MX1[i][2] = MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]]
                                elif  MX1[i][1] == "V" and MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]] != "T" and MX1[i+ MX4[i][2]][3 - MX4[i + MX4[i][2]][1]] != "T" and MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]] != MX1[i+ MX4[i][2]][3 - MX4[i + MX4[i][2]][1]] :
                                    MX1[i][1] = MX1[i +  MX4[i][2]][3 - MX4[i + MX4[i][2]][1]]
                                    MX1[i][2] = MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]]
                                elif  MX1[i][1] == "N" and MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]] != MX1[i+ MX4[i][2]][3 - MX4[i + MX4[i][2]][1]] :
                                    MX1[i][1] = MX1[i + MX4[i][2]][3 - MX4[i + MX4[i][2]][1]]
                                    MX1[i][2] = MX1[i - MX4[i][2]][MX4[i - MX4[i][2]][1]]


    # based on MX4 and MX1 reconstruct two allelic sequences
    for i in range(1, iMax):
        if MX4[i][1] == 1 :
            temp1 = temp1 + MX1[i][1]
            temp2 = temp2 + MX1[i][2]
        elif  MX4[i][1] == 2 :
            temp1 = temp1 + MX1[i][2]
            temp2 = temp2 + MX1[i][1]
        else:
            ambig2 = ambig2 + 1
            if MX1[i][1] == "A" and MX1[i][2] == "G" :
                temp1 = temp1 + "R"
                temp2 = temp2 + "R"
            elif  MX1[i][1] == "C" and MX1[i][2] == "T" :
                temp1 = temp1 + "Y"
                temp2 = temp2 + "Y"
            elif  MX1[i][1] == "G" and MX1[i][2] == "C" :
                temp1 = temp1 + "S"
                temp2 = temp2 + "S"
            elif  MX1[i][1] == "A" and MX1[i][2] == "T" :
                temp1 = temp1 + "W"
                temp2 = temp2 + "W"
            elif  MX1[i][1] == "A" and MX1[i][2] == "C" :
                temp1 = temp1 + "M"
                temp2 = temp2 + "M"
            elif  MX1[i][1] == "G" and MX1[i][2] == "T" :
                temp1 = temp1 + "K"
                temp2 = temp2 + "K"

    return temp1, temp2

#Function to build the contig of forward and reverse sequences
#forewardseq = forward sequence in 5'->3'; reverseseq = reverse sequence in 5'->3'
def buildContig(forwardseq, reverseseq, forwardObj=None, reverseObj=None): 
    LOG.info(f"Building contig from forward {len(forwardseq)}bp and reverse {len(reverseseq)}bp sequences")
    if forwardObj and reverseObj:
        with open("forward_original_tmp.fasta", "w") as fp:
            fp.write("".join(forwardObj.seq))
        with open("forward_tmp.fasta", "w") as fp:
            fp.write(forwardseq)

        blastn_cline = NcbiblastnCommandline(query="forward_tmp.fasta", subject="forward_original_tmp.fasta",
                                            reward=1, penalty=-2,gapopen=5, gapextend=2,evalue=0.00001, outfmt=5,
                                            out="blastn_contig_results.xml")  
        blastn_cline()   
        blast_records = list(NCBIXML.parse(open("blastn_contig_results.xml")))
        start_pos = min([(a.hsps[0].sbjct_end, a.hsps[0].sbjct_start) for a in blast_records[0].alignments][0])-1 
        #start_pos = re.search(forwardseq[0:10], "".join(forwardObj.seq)).start() #use 10bp of to find original seq position
        phred_scores_forewardseq = forwardObj.phred_qual[start_pos:len(forwardseq)]
        phred_scores_forewardseq_avg = sum(phred_scores_forewardseq)/len(phred_scores_forewardseq)
        forwardObj.avgPhredQuality = phred_scores_forewardseq_avg
        
        #start_pos = re.search(revcomp(reverseseq[-6:]), "".join(reverseObj.seq)).start()
        #find start position in the reverse sequence via  BLAST as it is more reliable when IUPAC bases  are used (YCCTAM)
        with open("reverseseq_original_tmp.fasta", "w") as fp:
            fp.write("".join(reverseObj.seq))
        with open("reverseseq_tmp.fasta", "w") as fp:
            fp.write(reverseseq)

        blastn_cline = NcbiblastnCommandline(query="reverseseq_tmp.fasta", subject="reverseseq_original_tmp.fasta",
                                            reward=1, penalty=-2,gapopen=5, gapextend=2,evalue=0.00001, outfmt=5,
                                            out="blastn_contig_results.xml")  
        blastn_cline()   
        blast_records = list(NCBIXML.parse(open("blastn_contig_results.xml")))

        start_pos = min([(a.hsps[0].sbjct_end, a.hsps[0].sbjct_start) for a in blast_records[0].alignments][0])-1 
        phred_scores_reverseseq = reverseObj.phred_qual[start_pos:len(reverseseq)][::-1] #reverse order
        phred_scores_reverseseq_avg = sum(phred_scores_reverseseq)/len(phred_scores_reverseseq)
        reverseObj.avgPhredQuality = phred_scores_reverseseq_avg

           
        with open("forwardseq_tmp.fasta", "w") as fp:
            fp.write(forwardseq)
        with open("reverseseq_tmp.fasta", "w") as fp:
            fp.write(reverseseq)

        blastn_cline = NcbiblastnCommandline(query="forwardseq_tmp.fasta", subject="reverseseq_tmp.fasta",
                                            reward=1, penalty=-2,gapopen=5, gapextend=2,evalue=0.00001, outfmt=5,
                                            out="blastn_contig_results.xml")  
        blastn_cline()   
        blast_records = list(NCBIXML.parse(open("blastn_contig_results.xml"))) 
        for record in blast_records:
            if record.alignments:
                hsp = record.alignments[0].hsps[0]
                if phred_scores_forewardseq_avg > phred_scores_reverseseq_avg:
                    #use forward aligned overlap sequence instead of reverse (forward+forward_overlap_alignment(blast)+reverse)
                    contig = forwardseq[0:hsp.query_end]+reverseseq[hsp.sbjct_end:]
                else:
                    #use reverse aligned overlap sequence instead of reverse  (forward+reverse_overlap_alignment(blast)+reverse)
                    contig = forwardseq[0:hsp.query_start]+reverseseq[hsp.sbjct_start-1:]
                LOG.info(f"A contig of {len(contig)}bp was successfully formed ...")    
    else:
        substring = reverseseq[0:10]
    
        p = forwardseq.find(substring)
        contig = ""
    
        LOG.info(f"Trying to build contig from forward {len(forwardseq)}bp and reverse {len(reverseseq)}bp sequence based on first reverse {len(substring)}bp overlap ({substring}) ...")
        if p != -1: #if 10bp overlap is found in forward sequence
            contig += forwardseq[:p]
            contig += reverseseq
        else:
            LOG.warning(f"Could not form contig as was not able to find overlap based on {substring} substring")
    



                                

    return contig


#MixedSeq Class Description:
#   This class keeps all the information of each Cryptosporidium sample's
#   18S Sanger sequence and analyzes it to classify.
#   The data is obtained from reading the given ab1 file and is then aligned
#   with the reference sequences to determine family and subtype.  It then
#   counts the number of repeats within the sequence, depending on the
#   family determined.
class MixedSeq(object):

    def __init__(self, File, tabFile, analysisType):
        self.g = []
        self.a = []
        self.t = []
        self.c = []

        self.file = File
        self.tabfile = tabFile
        self.analysistype = analysisType

        self.name=""  #name of the file
        self.species1=""
        self.species2=""

        self.species1Seq=""
        self.species2Seq=""

        self.species1RefSeq = ""
        self.species2RefSeq = ""

        self.species1queryTo = 0
        self.species1queryFrom = 0
        self.species1sbjtTo = 0
        self.species1sbjtFrom = 0

        self.species2queryTo = 0
        self.species2queryFrom = 0
        self.species2sbjtTo = 0
        self.species2sbjtFrom = 0

        self.mixed=False

        self.oldseq = []  #before trimming the ends of the sequence
        self.seq = []    #after trimming the ends of the sequence
        self.seqLength = 0

        self.peakLoc = [] #peak locations
        self.phred_qual = [] #phred quality: 10 = 90% base call accuracy, 20=99% (log scale)

        self.forwardSeq = True

        self.avgPhredQuality = 0
        self.avgPhred = 0



    def readFiles(self, dataFile, forw, filetype="abi"):
        #setting whether the sequence is the forward or reverse strand
        #(used later on know whether the reverse complement is needed)
        #print("Datafile:",dataFile)
        if forw:
            self.forwardSeq = True

        else:
            self.forwardSeq = False



        #**********Beginning to work with the ab1 file given:**********#

        #retrieves the sample file name (removes directory pathway)
        self.name = dataFile.split("/")[len(dataFile.split("/"))-1]

        if TESTING:
            print("\nSequence: ", self.name) #Lets user know which sequence the program is on

        #opens the ab1 file
        if filetype == "abi" or filetype == "ab1":
            handle=open(dataFile,"rb")
            record=SeqIO.read(handle, "abi")

        elif filetype == "fasta" or filetype == "fa":
            handle=open(dataFile,"r")   
            record=SeqIO.read(handle, "fasta")
            raw_seq = list(record.seq.upper()) #making sure all bases converted to upper case so matching works
            self.seq = raw_seq
            self.phred_qual = [60] * len(raw_seq)
            return True


        #retrieving base amplitude data
        self.g=np.array(record.annotations['abif_raw']['DATA9'])
        self.a=np.array(record.annotations['abif_raw']['DATA10'])
        self.t=np.array(record.annotations['abif_raw']['DATA11'])
        self.c=np.array(record.annotations['abif_raw']['DATA12'])

        self.aOrig = copy.copy(self.a)
        self.gOrig = copy.copy(self.g)
        self.tOrig = copy.copy(self.t)
        self.cOrig = copy.copy(self.c)


        #obtaining primary bases, includes N's
        self.seq=list(record.annotations['abif_raw']['PBAS2'])
        self.oldseq=list(record.annotations['abif_raw']['PBAS2'])
        self.seqLength = len(self.seq)
        self.origLength = len(self.seq)

        #in case the class is of a binary nature
        if any([isinstance(element,str) == False for element in self.seq]):
            self.seq = list(record.annotations['abif_raw']['PBAS2'].decode('UTF-8'))
            self.oldseq = list(record.annotations['abif_raw']['PBAS2'].decode('UTF-8'))
            #exit("Sequence extracted is of wrong data type. Expected lisf of stings, got {}".format(self.seq))

        #peak locations
        self.peakLoc=np.array(record.annotations['abif_raw']['PLOC2'])

        #phred quality of the bases
        self.phred_qual=record.letter_annotations['phred_quality']

        self.mixed=False


        gTemp = []
        aTemp = []
        tTemp = []
        cTemp = []



        for i in range(0, len(self.peakLoc)):
            gTemp.append(self.g[self.peakLoc[i]])
            aTemp.append(self.a[self.peakLoc[i]])
            tTemp.append(self.t[self.peakLoc[i]])
            cTemp.append(self.c[self.peakLoc[i]])

        self.g=gTemp
        self.a=aTemp
        self.c=cTemp
        self.t=tTemp



        length = len(self.phred_qual)

        if length > 800:
            self.a=self.a[:800]
            self.g=self.g[:800]
            self.c=self.c[:800]
            self.t=self.t[:800]
            self.phred_qual=self.phred_qual[:800]
            self.peakLoc=self.peakLoc[:800]
            self.seq=self.seq[:800]
            self.oldseq=self.oldseq[:800]
            length=len(self.phred_qual)
            self.seqLength=800
            self.origLength=800

        sum=0
        for i in range(0, length):
            sum+=self.phred_qual[i]
            #print(self.phred_qual[i])

        self.avgPhredQuality = int(sum/length)
        self.avgPhred = int(sum/length)



        self.majorSeq = ""
        self.minorSeq = ""

        self.counter = 0
        self.avgLRI = 0


        return True


    #fixN() finds the best base when an 'N' is encountered in the sequence by finding
    #   the maximum amplitude of all four bases at that position
    #   **Not yet implemented double peaks
    def fixN(self):
        if self.avgPhredQuality > 20:
            qual_cutOff = self.avgPhredQuality
        else:
            qual_cutOff = 20 #99% base calling certainty
        foundBegin = False

        for x in range(0, self.seqLength):
            if self.phred_qual[x] >= qual_cutOff and self.oldseq[x]!='N':
                self.beginSeq = x
                break
        #print(self.beginSeq)

        if x < 30:
            self.beginSeq = 30
            x=30

        self.oldseq = copy.copy(self.oldseq[x:self.seqLength-30])
        self.phred_qual = copy.copy(self.phred_qual[x:self.seqLength-30])
        self.a = copy.copy(self.a[x:self.seqLength-30])
        self.t = copy.copy(self.t[x:self.seqLength-30])
        self.g = copy.copy(self.g[x:self.seqLength-30])
        self.c = copy.copy(self.c[x:self.seqLength-30])
        self.peakLoc = copy.copy(self.peakLoc[x:self.seqLength-30])
        self.seqLength = len(self.oldseq)
        self.fixedSeq = self.oldseq

        self.backgroundG = 0
        self.backgroundC = 0
        self.backgroundT = 0
        self.backgroundA = 0

        avgA = 0
        avgC = 0
        avgT = 0
        avgG = 0

        cA = 0
        cC = 0
        cG = 0
        cT = 0



        for i in range(0, self.seqLength):
            a = self.a[i]
            g = self.g[i]
            c = self.c[i]
            t = self.t[i]

            maxAmp = max(a,g,c,t)

            if a == maxAmp:
                if self.oldseq[i]!='A' and self.oldseq[i]!='G' and self.oldseq[i]!='C' and self.oldseq[i]!='T':
                    self.oldseq[i] = 'A'
                secondMaxAmp = max(g,c,t)

                avgG += self.g[i]
                cG += 1
                avgT += self.t[i]
                cT += 1
                avgC += self.c[i]
                cC += 1

            elif g == maxAmp:
                if self.oldseq[i]!='A' and self.oldseq[i]!='G' and self.oldseq[i]!='C' and self.oldseq[i]!='T':
                    self.oldseq[i] = 'G'
                secondMaxAmp = max(a,c,t)

                avgA += self.a[i]
                cA += 1
                avgT += self.t[i]
                cT += 1
                avgC += self.c[i]
                cC += 1

            elif c == maxAmp:
                if self.oldseq[i]!='A' and self.oldseq[i]!='G' and self.oldseq[i]!='C' and self.oldseq[i]!='T':
                    self.oldseq[i] = 'C'
                secondMaxAmp = max(g,a,t)

                avgG += self.g[i]
                cG += 1
                avgT += self.t[i]
                cT += 1
                avgA += self.a[i]
                cA += 1

            elif t == maxAmp:
                if self.oldseq[i]!='A' and self.oldseq[i]!='G' and self.oldseq[i]!='C' and self.oldseq[i]!='T':
                    self.oldseq[i] = 'T'
                secondMaxAmp = max(g,c,a)

                avgG += self.g[i]
                cG += 1
                avgA += self.a[i]
                cA += 1
                avgC += self.c[i]
                cC += 1

            if secondMaxAmp != 0:
                self.avgLRI += math.log((maxAmp/secondMaxAmp),2)
                #print(i, self.oldseq[i], math.log((maxAmp/secondMaxAmp),2), a, g, c,t)

        if cA != 0:
            self.backgroundA = int((avgA/cA) *2)
        if cC != 0:
            self.backgroundC = int((avgC/cC) *2)
        if cG != 0:
            self.backgroundG = int((avgG/cG) *2)
        if cT != 0:
            self.backgroundT = int((avgT/cT) *2)

        if self.seqLength != 0:
            self.avgLRI = self.avgLRI/self.seqLength
        else:
            self.avgLRI = 0

        #print(self.backgroundA,self.backgroundC, self.backgroundG,self.backgroundT)




    # findHeteroBases() finds the IUPAC bases within the chromatogram
    def findHeteroBases(self, lri):
        maxLRI1 = lri
        maxLRI = lri
        if self.avgLRI != 0:
            maxLRI1 = math.log((self.avgLRI/2),2)
            maxLRI = ((self.avgLRI - maxLRI1) + 2.0) /2

        if maxLRI1 > 1.0:
            maxLRI = self.avgLRI/2


        if maxLRI < lri:
            maxLRI = lri

        self.avgLRI = maxLRI


        counter = 0
        newSeq = []

        lriBefore1 = maxLRI
        lriBefore2 = maxLRI
        lriBefore3 = maxLRI
        lriBefore4 = maxLRI
        lriBefore5 = maxLRI


        for i in range(0, self.seqLength):
            a = self.a[i]
            g = self.g[i]
            c = self.c[i]
            t = self.t[i]

            maxAmp = max(a,g,c,t)

            if a == maxAmp:
                secondMaxAmp = max(g,c,t)
            elif g == maxAmp:
                secondMaxAmp = max(a,c,t)
            elif c == maxAmp:
                secondMaxAmp = max(g,a,t)
            elif t == maxAmp:
                secondMaxAmp = max(g,c,a)

            if secondMaxAmp != 0:
                LRI = math.log((maxAmp/secondMaxAmp),2)

            else:
                LRI = maxLRI + 0.1



            if LRI <= maxLRI:
                counter += 1
                if a == maxAmp:
                    firstBase = "A"
                    self.majorSeq += "A"
                    self.fixedSeq[i] = "A"
                elif g == maxAmp:
                    firstBase = "G"
                    self.majorSeq += "G"
                    self.fixedSeq[i] = "G"
                elif c == maxAmp:
                    firstBase = "C"
                    self.majorSeq += "C"
                    self.fixedSeq[i] = "C"
                elif t == maxAmp:
                    firstBase = "T"
                    self.majorSeq += "T"
                    self.fixedSeq[i] = "T"

                if g == secondMaxAmp and "G"!= firstBase:
                    secondBase = "G"
                    self.minorSeq += "G"
                elif a == secondMaxAmp and "A"!=firstBase:
                    secondBase = "A"
                    self.minorSeq += "A"
                elif c == secondMaxAmp and "C" != firstBase:
                    secondBase = "C"
                    self.minorSeq += "C"
                elif t ==  secondMaxAmp and "T" != firstBase:
                    secondBase = "T"
                    self.minorSeq += "T"


                iupacBase = addIupac(firstBase, secondBase)
                newSeq.append(iupacBase)

                if i > 50 and i < self.seqLength-50 and LRI < 1.0:
                    self.mixed = True

            else:
                aBefore = 0
                cBefore = 0
                gBefore = 0
                tBefore = 0

                aAfter = 0
                cAfter = 0
                gAfter = 0
                tAfter = 0

                if i-1 >= 0:
                    lowA = self.a[i-1]
                    lowC = self.c[i-1]
                    lowG = self.g[i-1]
                    lowT = self.t[i-1]

                    for j in range(self.peakLoc[i-1], self.peakLoc[i]):
                        if self.aOrig[j] < lowA:
                            lowA = self.aOrig[j]
                        if self.gOrig[j] < lowG:
                            lowG = self.gOrig[j]
                        if self.tOrig[j] < lowT:
                            lowT = self.tOrig[j]
                        if self.cOrig[j] < lowC:
                            lowC = self.cOrig[j]

                    aBefore = lowA
                    cBefore = lowC
                    gBefore = lowG
                    tBefore = lowT

                if i + 1 < len(self.a):
                    lowA = self.a[i]+1
                    lowC = self.c[i]+1
                    lowG = self.g[i]+1
                    lowT = self.t[i]+1

                    for j in range(self.peakLoc[i]+1, self.peakLoc[i+1]):
                        if self.aOrig[j] < lowA:
                            lowA = self.aOrig[j]
                        if self.gOrig[j] < lowG:
                            lowG = self.gOrig[j]
                        if self.tOrig[j] < lowT:
                            lowT = self.tOrig[j]
                        if self.cOrig[j] < lowC:
                            lowC = self.cOrig[j]

                    aAfter = lowA
                    cAfter = lowC
                    gAfter = lowG
                    tAfter = lowT


                notBaseline = False
                secondBase = 'N'
                if t == secondMaxAmp:
                    if t > self.backgroundT and t>=tBefore and t>=tAfter:
                        notBaseline = True
                        secondBase = 'T'
                elif a == secondMaxAmp:
                    if a > self.backgroundA and a>=aBefore and a>=aAfter:
                        notBaseline = True
                        secondBase = 'A'
                elif c == secondMaxAmp:
                    if c > self.backgroundC and c>=cBefore and c>=cAfter:
                        notBaseline = True
                        secondBase = 'C'
                elif g == secondMaxAmp:
                    if g > self.backgroundG and g>=gBefore and g>=gAfter:
                        notBaseline = True
                        secondBase = 'G'

                avgLRI = (lriBefore1 + lriBefore2 + lriBefore3 + lriBefore4 + lriBefore5)/5



                if (abs(LRI-avgLRI) >= 1.5 or abs(LRI-lriBefore1) >= 1.5)  and notBaseline:
                    counter += 1
                    if a == maxAmp:
                        firstBase = "A"
                        self.majorSeq += "A"
                        self.fixedSeq[i] = "A"
                    elif g == maxAmp:
                        firstBase = "G"
                        self.majorSeq += "G"
                        self.fixedSeq[i] = "G"
                    elif c == maxAmp:
                        firstBase = "C"
                        self.majorSeq += "C"
                        self.fixedSeq[i] = "C"
                    elif t == maxAmp:
                        firstBase = "T"
                        self.majorSeq += "T"
                        self.fixedSeq[i] = "T"

                    if g == secondMaxAmp and "G"!= firstBase:
                        secondBase = "G"
                        self.minorSeq += "G"
                    elif a == secondMaxAmp and "A"!=firstBase:
                        secondBase = "A"
                        self.minorSeq += "A"
                    elif c == secondMaxAmp and "C" != firstBase:
                        secondBase = "C"
                        self.minorSeq += "C"
                    elif t ==  secondMaxAmp and "T" != firstBase:
                        secondBase = "T"
                        self.minorSeq += "T"
                    iupacBase = addIupac(firstBase, secondBase)
                    newSeq.append(iupacBase)

                else:
                    if a == maxAmp and self.oldseq[i] == 'A':
                        newSeq.append('A')
                        self.majorSeq += "A"
                        self.minorSeq += "A"
                        self.fixedSeq[i] = "A"
                    elif g == maxAmp and self.oldseq[i] == 'G':
                        newSeq.append('G')
                        self.majorSeq += "G"
                        self.minorSeq += "G"
                        self.fixedSeq[i] = "G"
                    elif c == maxAmp and self.oldseq[i] == 'C':
                        newSeq.append('C')
                        self.majorSeq += "C"
                        self.minorSeq += "C"
                        self.fixedSeq[i] = "C"
                    elif t == maxAmp and self.oldseq[i] == 'T':
                        newSeq.append('T')
                        self.majorSeq += "T"
                        self.minorSeq += "T"
                        self.fixedSeq[i] = "T"
                    else:
                        newSeq.append(self.oldseq[i])
                        self.majorSeq += self.oldseq[i]
                        self.minorSeq += self.oldseq[i]

            lriBefore5 = lriBefore4
            lriBefore4 = lriBefore3
            lriBefore3 = lriBefore2
            lriBefore2 = lriBefore1
            lriBefore1 = LRI


        self.seq = newSeq
        self.counter = counter

        if counter > len(newSeq)*0.33:
            self.avgPhredQuality = 0


        return newSeq


    # determineAllType() finds the multiple mixed species, if present
    def determineAllTypes(self, customdatabsename):
        reverse1 = False
        reverse2 = False

        #get file type information
        filetype = utilities.getFileType(self.name)
        # Filename to write
        filename_seq1 = "query1.txt"
        filename_seq2 = "query2.txt"

        # Open the file with writing permission
        file1 = open(filename_seq1, 'w')

        seq = ''.join(self.seq)

       
        seq1, seq2 = indelligent(seq)


        l1 = len(seq1)
        l2 = len(seq2)
        
       
        if self.counter < int(len(self.seq)*0.2) and filetype == "abi":
            s1=""
            s2=""

            for i in range(0, l1):
                if seq1[i] != 'A' and seq1[i] != 'C' and seq1[i] != 'G' and seq1[i]!='T':
                    s1 += self.majorSeq[i]
                else:
                    s1 += seq1[i]

            for i in range(0, l2):
                if seq2[i]!='A' and seq2[i] != 'C' and seq2[i] != 'G' and seq2[i]!='T':
                    s2 += self.minorSeq[i]
                else:
                    s2 += seq2[i]


            seq1 = s1
            seq2 = s2

        self.species1Seq=seq1
        self.species2Seq=seq2
        #print(self.species1Seq, self.species2Seq)


        # Write a line to the file
        file1.write(seq1)
        #file1.write(self.majorSeq)


        # Close the file
        file1.close()


        if customdatabsename:
            blastn_cline = NcbiblastnCommandline(cmd='blastn',  query="query1.txt",
                                                 dust='yes',
                                                 db="custom_db",
                                                 reward=1, penalty=-2, gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="result1.xml")
        else:
            blastn_cline = NcbiblastnCommandline(query="query1.txt", db= os.path.dirname(__file__)+"/reference_database/msr_ref.fa",
                                              outfmt=5, out="result1.xml")

        stdout, stderr = blastn_cline()


        if (os.stat("result1.xml").st_size == 0):
            self.species1=";>No blast hits."
            self.species1RefSeq=""
            seq1=""

        else:
            result_handle = open("result1.xml", 'r')
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)
            
            blast_record.alignments = utilities.sort_blast_hits_by_id_and_bitscore(blast_record)


            if len(blast_record.alignments) > 0:
                hit = blast_record.alignments[0].hit_id
                data = hit.split(" ")[0]
                self.species1 = blast_record.alignments[0].hit_id
                self.species1RefSeq = blast_record.alignments[0].hsps[0].sbjct

                if blast_record.alignments[0].hsps[0].sbjct_start > blast_record.alignments[0].hsps[0].sbjct_end:
                    reverse1 = True
            else:
                self.species1=";>No blast hits."
                self.species1RefSeq=""
                seq1=""


        file2 = open(filename_seq2, 'w')
        file2.write(seq2)
        file2.close()


        if customdatabsename:
            blastn_cline = NcbiblastnCommandline(cmd='blastn',
                                                 query="query2.txt",
                                                 dust='yes',
                                                 db="custom_db",
                                                 reward=1, penalty=-2, gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="result2.xml")
        else:
            blastn_cline = NcbiblastnCommandline(
                                             query="query2.txt",

                                             db=os.path.dirname(__file__)+"/reference_database/msr_ref.fa",
                                              outfmt=5, out="result2.xml")
        blastn_cline
        stdout, stderr = blastn_cline()

        if (os.stat("result2.xml").st_size == 0):
            self.species2=";>No blast hits."
            self.species2RefSeq=""
            seq2 = ""

        else:
            result_handle = open("result2.xml")
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)
            blast_record.alignments = utilities.sort_blast_hits_by_id_and_bitscore(blast_record)

            if len(blast_record.alignments) > 0:
                hit = blast_record.alignments[0].hit_id
                data = hit.split(" ")[0]
                self.species2 = blast_record.alignments[0].hit_id
                self.species2RefSeq = blast_record.alignments[0].hsps[0].sbjct

                if blast_record.alignments[0].hsps[0].sbjct_start > blast_record.alignments[0].hsps[0].sbjct_end:
                    reverse2 = True
            else:
                self.species2=";>No blast hits."
                self.species2RefSeq=""
                seq2 = ""


        if self.species1 == ";>No blast hits." and self.species2 == ";>No blast hits.":
            file2 = open(filename_seq2, 'w')
            file2.write(''.join(self.oldseq))
            file2.close()


            if customdatabsename:
                blastn_cline = NcbiblastnCommandline(cmd='blastn',
                                                     query="query2.txt",
                                                     dust='yes',
                                                     db="custom_db",
                                                     reward=1, penalty=-2, gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="result2.xml")
            else:
                blastn_cline = NcbiblastnCommandline(
                                                 query="query2.txt",
                                                 db=os.path.dirname(__file__)+"/reference_database/msr_ref.fa",
                                                  outfmt=5, out="result2.xml")
            blastn_cline
            stdout, stderr = blastn_cline()

            if (os.stat("result2.xml").st_size == 0):
                self.species2=";>No blast hits."
                self.species2RefSeq=""
                seq2 = ""

                self.species1=";>No blast hits."
                self.species1RefSeq=""
                seq1 = ""

            else:
                result_handle = open("result2.xml")
                blast_records = NCBIXML.parse(result_handle)
                blast_record = next(blast_records)
                blast_record.alignments = utilities.sort_blast_hits_by_id_and_bitscore(blast_record)

                if len(blast_record.alignments) > 0:
                    hit = blast_record.alignments[0].hit_id
                    data = hit.split(" ")[0]
                    self.species2 = blast_record.alignments[0].hit_id
                    self.species2RefSeq = blast_record.alignments[0].hsps[0].sbjct
                    self.species1 = blast_record.alignments[0].hit_id
                    self.species1RefSeq = blast_record.alignments[0].hsps[0].sbjct
                    seq1 = ''.join(self.oldseq)
                    seq2 = ''.join(self.oldseq)


        elif self.species1 != ";>No blast hits.":

            os.system("blastdbcmd -db " + os.path.dirname(__file__)+"/reference_database/msr_ref.fa -entry " + "'{}'".format(self.species1) + " -out refseq.fa")

            file = open("refseq.fa", 'r')
            lines = file.readlines()

            refseq = ""

            for i in range (1, len(lines)):
                refseq += lines[i]

            if reverse1:
                reverse = ""
                l = len(refseq) -1

                for i in range(l, -1 , -1):
                    if refseq[i] == 'A':
                        reverse += 'T'
                    elif refseq[i] == 'C':
                        reverse += 'G'
                    elif refseq[i] == 'G':
                        reverse += 'C'
                    elif refseq[i] == 'T':
                        reverse += 'A'
                self.species1RefSeq = reverse
            else:
                self.species1RefSeq = refseq

        elif self.species2 != ";>No blast hits.":

            os.system("blastdbcmd -db " + os.path.dirname(__file__)+"/reference_database/msr_ref.fa -entry " + "'{}'".format(self.species2) + " -out refseq.fa")

            file = open("refseq.fa", 'r')
            lines = file.readlines()

            refseq2 = ""

            for i in range(1, len(lines)):
                refseq2 += lines[i]

            if reverse2:
                reverse = ""
                l = len(refseq2) -1

                for i in range(l, -1 , -1):
                    if refseq2[i] == 'A':
                        reverse += 'T'
                    elif refseq2[i] == 'C':
                        reverse += 'G'
                    elif refseq2[i] == 'G':
                        reverse += 'C'
                    elif refseq2[i] == 'T':
                        reverse += 'A'
                self.species2RefSeq = reverse
            else:
                self.species2RefSeq = refseq2


        if self.species1RefSeq == "" and self.species2RefSeq != "":
            self.species1RefSeq = self.species2RefSeq

        elif self.species1RefSeq != "" and self.species2RefSeq == "":
            self.species2RefSeq = self.species1RefSeq

        elif self.species1RefSeq == "" and self.species2RefSeq == "":
            return

        if seq1 == "":
            seq1 = ''.join(self.oldseq)
        if seq2 == "":
            seq2 = ''.join(self.oldseq)


        file = open("align.fa", 'w')
        file.write(">Ref1\n")
        file.write(self.species1RefSeq)
        file.write("\n>Seq1\n")
        file.write(seq1)
        file.write("\n>Ref2\n")
        file.write(self.species2RefSeq)
        file.write("\n>Seq2\n")
        file.write(seq2)
        file.write("\n>IUPAC\n")
        file.write(seq)
        file.close()


        clustalw_cline = ClustalwCommandline("clustalw", infile="align.fa")
        stdout, stderr = clustalw_cline()

        align = AlignIO.read("align.aln", "clustal")


        for record in align:
            if record.id == "Ref1":
                ref_align= record.seq
            elif record.id == "Ref2":
                ref2_align=record.seq
            elif record.id == "Seq1":
                seq1_align = record.seq
            elif record.id == "Seq2":
                seq2_align = record.seq
            elif record.id == "IUPAC":
                iupac_align = record.seq
            elif record.id == "major":
                majorseq = record.seq
            elif record.id == "minor":
                minorseq = record.seq

        length1 = len(ref_align)
        length2 = len(ref2_align)

        begin = 0

        for i in range(0, length1):
            if ref_align[i]!='-' and ref2_align[i]!='-' and seq1_align[i]!='-' and seq2_align[i]!='-' and iupac_align[i]!='-': #and majorseq[i] !='-' and minorseq[i]!='-':
                begin = i
                break

        amp_cutoff = 0

        for i in range(0, begin+1):
            if seq1_align[i] != '-' and seq2_align[i] != '-':
                amp_cutoff = i
                break

        a = self.a[begin-amp_cutoff:]
        c = self.c[begin-amp_cutoff:]
        g = self.g[begin-amp_cutoff:]
        t = self.t[begin-amp_cutoff:]


        newseq1 = ""
        newseq2 = ""
        seq1temp = ""
        seq2temp = ""
        ref1temp = ""
        ref2temp = ""
        iupactemp = ""
        majortemp = ""
        minortemp=""

        l1 = len(seq1_align)
        l2 = len(seq2_align)
        r1 = len(ref_align)
        r2 = len(ref2_align)

        for i in range(begin, length1):
            s1 = "N"
            s2 = "N"
            if i < r1:
                s1 = ref_align[i]
            if i < r2:
                s2 = ref2_align[i]

            if seq1_align[i] != '-' and (s1 != '-' or s2 != '-'):
                seq1temp += seq1_align[i]

        for i in range(begin, length2):
            s1 = "N"
            s2 = "N"
            if i < r1:
                s1 = ref_align[i]
            if i < r2:
                s2 = ref2_align[i]

            if seq2_align[i] != '-' and (s1 != '-' or s2 != '-'):
                seq2temp += seq2_align[i]

        for i in range(begin, r1):
            if ref_align[i] != '-':
                s1 = "N"
                s2 = "N"
                s3 = "N"
                if i < l1:
                    s1 = seq1_align[i]
                if i < l2:
                    s2 = seq2_align[i]
                if i < r2:
                    s3 = ref2_align[i]

                if s3 == '-':

                    if (s1 != '-' or s2 != '-') or iupac_align[i] != '-':
                        ref1temp += ref_align[i]

                elif ref_align[i] == s3 and s1 == '-' and self.species1 == self.species2:
                    pass

                else:
                    if (s1 != '-' and s2 != '-') or (iupac_align[i]!='-' and s1 != '-') or (iupac_align[i]!='-' and s2 != '-'):
                        ref1temp += ref_align[i]

        for i in range(begin, r2):
            if ref2_align[i] != '-':
                s1 = "N"
                s2 = "N"
                s3 = "N"
                if i < l1:
                    s1 = seq1_align[i]
                if i < l2:
                    s2 = seq2_align[i]
                if i < r1:
                    s3 = ref_align[i]

                if s3 == '-':

                    if (s1 != '-' or s2 != '-') or iupac_align[i] != '-':
                        ref2temp += ref2_align[i]

                elif ref2_align[i] == s3 and s2 == '-' and self.species1 == self.species2:
                    pass


                else:
                    if (s1 != '-' and s2 != '-') or (iupac_align[i]!='-' and s1 != '-') or (iupac_align[i]!='-' and s2 != '-'):
                        ref2temp += ref2_align[i]


        for i in range(begin, len(iupac_align)):
            s1 = "N"
            s2 = "N"
            if i < r1:
                s1 = ref_align[i]
            if i < r2:
                s2 = ref2_align[i]
            if iupac_align[i] != '-' and (s1 != '-' or s2 != '-'):
                iupactemp += iupac_align[i]


        ref_align = ref1temp
        ref2_align = ref2temp
        seq1_align = seq1temp
        seq2_align = seq2temp
        iupac_align = iupactemp


        length = min(len(ref_align), len(ref2_align), len(seq1_align), len(seq2_align), len(iupac_align))


        for i in range(0, length):

            if ref_align[i] == seq1_align[i] and ref2_align[i] == seq2_align[i]:
                newseq1 += seq1_align[i]
                newseq2 += seq2_align[i]

            elif ref_align[i] == seq2_align[i] and ref2_align[i] == seq1_align[i]:
                newseq2 += seq1_align[i]
                newseq1 += seq2_align[i]

            else:
                Gamp = g[i]
                Tamp = t[i]
                Aamp = a[i]
                Camp = c[i]
                maxAmp = max(Gamp, Tamp, Camp, Aamp)
                maxBase = ''
                minorBase = ''

                if Gamp == maxAmp:
                    maxBase = 'G'
                    minorAmp = max(Tamp,Camp,Aamp)

                elif Tamp == maxAmp:
                    maxBase = 'T'
                    minorAmp = max(Gamp, Aamp, Camp)
                elif Aamp == maxAmp:
                    maxBase = 'A'
                    minorAmp = max(Gamp, Camp, Tamp)
                else:
                    maxBase = 'C'
                    minorAmp = max(Gamp, Aamp, Tamp)

                if Gamp == minorAmp:
                    if Gamp > self.backgroundG:
                        minorBase = 'G'
                    else:
                        minorBase = maxBase
                elif Tamp == minorAmp:
                    if Tamp > self.backgroundT:
                        minorBase = 'T'
                    else:
                        minorBase = maxBase
                elif Aamp == minorAmp:
                    if Aamp > self.backgroundA:
                        minorBase='A'
                    else:
                        minorBase = maxBase
                else:
                    if Camp > self.backgroundC:
                        minorBase='C'
                    else:
                        minorBase = maxBase

                if maxBase == ref_align[i] and minorBase == ref2_align[i]:
                    newseq1 += maxBase
                    newseq2 += minorBase
                elif minorBase==ref_align[i] and maxBase == ref2_align[i]:
                    newseq1 += minorBase
                    newseq2 += maxBase
                elif maxBase == ref_align[i] and maxBase == ref2_align[i]:
                    newseq1 += maxBase
                    newseq2 += maxBase
                else:
                    p = self.avgLRI/2

                    if p < 2.0:
                        p = 2.0

                    if ref2_align[i] == 'G':
                        if Gamp != 0:
                            p = math.log((maxAmp/Gamp),2)

                        if p < 2.0:
                            newseq2 += 'G'
                        else:
                            newseq2 += maxBase

                    elif ref2_align[i] == 'A':
                        if Aamp!=0:
                            p = math.log((maxAmp/Aamp),2)


                        if p < 2.0:
                            newseq2 += 'A'
                        else:
                            newseq2 += maxBase

                    elif ref2_align[i] == 'C':
                        if Camp != 0:
                            p = math.log((maxAmp/Camp),2)

                        if p < 2.0:
                            newseq2 += 'C'
                        else:
                            newseq2 += maxBase

                    elif ref2_align[i] == 'T':
                        if Tamp != 0:
                            p = math.log((maxAmp/Tamp),2)

                        if p < 2.0:
                            newseq2 += 'T'
                        else:
                            newseq2 += maxBase

                    p = 2.1

                    if ref_align[i] == 'G':
                        if Gamp != 0:
                            p = math.log((maxAmp/Gamp),2)

                        if p < 2.0:
                            newseq1 += 'G'
                        else:
                            newseq1 += maxBase

                    elif ref_align[i] == 'A':
                        if Aamp!=0:
                            p = math.log((maxAmp/Aamp),2)


                        if p < 2.0:
                            newseq1 += 'A'
                        else:
                            newseq1 += maxBase

                    elif ref_align[i] == 'C':
                        if Camp != 0:
                            p = math.log((maxAmp/Camp),2)

                        if p < 2.0:
                            newseq1 += 'C'
                        else:
                            newseq1 += maxBase

                    elif ref_align[i] == 'T':
                        if Tamp != 0:
                            p = math.log((maxAmp/Tamp),2)

                        if p < 2.0:
                            newseq1 += 'T'
                        else:
                            newseq1 += maxBase

        bitscore,evalue,query_coverage,query_length,percent_identity1, accession, species1,sequence = self.blast(newseq1, False)
        bitscore,evalue,query_coverage,query_length,percent_identity2, accession, species2,sequence = self.blast(newseq2, False)


        if percent_identity1 < 99 or percent_identity2 < 99:
            stop1 = len(newseq1)
            stop2 = len(newseq2)

            counter = 0
            counter2 = 0

            breakloop1 = False
            breakloop2 = False

            for i in range(length-1, 0, -1):
                if i < length/2:
                    stop1 = i
                    stop2 = i
                    break
                if breakloop1 and breakloop2:
                    break

                if i < len(newseq1) and i < len(ref_align):
                    if newseq1[i] == ref_align[i] and not breakloop1:
                        if counter == 30:
                            breakloop1 = True
                        else:
                            counter += 1
                    elif newseq1[i] != ref_align[i] and not breakloop1:
                        counter = 0
                        stop1 = i-1

                if i < len(newseq2) and i < len(ref2_align):

                    if newseq2[i] == ref2_align[i] and not breakloop2:
                        if counter2 == 30:
                            breakloop2 = True
                        else:
                            counter2 += 1

                    elif newseq2[i] != ref2_align[i] and not breakloop2:
                        stop2 = i-1
                        counter2 = 0


            self.species1Seq = newseq1[:stop1]
            self.species2Seq = newseq2[:stop2]

        else:
            self.species1Seq = newseq1
            self.species2Seq = newseq2


    #blast() performs a blast search and returns the results against the msr_ref.fa reference database
    def blast(self, sequence, customdatabsename):
        # Filename to write
        filename = "query.txt"

        # Open the file with writing permission
        myfile = open(filename, 'w')

        # Write a line to the file
        myfile.write(sequence)
        #print(sequence)


        # Close the file
        myfile.close()
        
        if customdatabsename:
            LOG.info(f"Running BLAST on {len(sequence)}bp query from {self.name} on {customdatabsename}")
            blastn_cline = NcbiblastnCommandline(cmd='blastn',query="query.txt", dust='yes',
                                             db="custom_db", reward=1, penalty=-2,gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="SSUresult.xml")
        else:
            LOG.info(f"Running BLAST on {len(sequence)}bp query from {self.name} on {os.path.dirname(__file__)+'/reference_database/msr_ref.fa'}")    
            blastn_cline = NcbiblastnCommandline(cmd='blastn',query="query.txt", dust='yes',
                                             db=os.path.dirname(__file__)+"/reference_database/msr_ref.fa", reward=1, penalty=-2,gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="SSUresult.xml")

        stdout, stderr = blastn_cline()

        if (os.stat("SSUresult.xml").st_size == 0):
            return "","",0,"",0,"","",sequence

        else:

            result_handle = open("SSUresult.xml", 'r')
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)

            if len(blast_record.alignments) == 0:
                LOG.error("No BLAST hits were found! No species will be identified!")
                return "","",0,"",0,"","",""

            blast_record.alignments = utilities.sort_blast_hits_by_id_and_bitscore(blast_record)

            maxBitScore = 0; identicalAlignHits = [] #BLAST hits that have identical top score (if any). Usually only single hit with unique top score
            for idx, alignment in enumerate(blast_record.alignments):
                hsp = alignment.hsps[0]
                if idx == 0:
                   maxBitScore =  alignment.hsps[0].score
                if maxBitScore == alignment.hsps[0].score: #
                    identicalAlignHits.append(alignment)
                

            if len(identicalAlignHits) >= 2:   
                identical_score_hits_ids_str = '\n'.join([f"{a.hit_id}\t{a.hsps[0].score}" for a in identicalAlignHits]) 
                LOG.warning(f"!!! Found {len(identicalAlignHits)} identically scored candidate BLAST hits in reference database with bitscores:\n{identical_score_hits_ids_str}.\nBe careful with species ID!") 
                #min_gaps = min([align.hsps[0].gaps for align in identicalAlignHits])
                #min_gaps_alignments = [align for align in identicalAlignHits if align.hsps[0].gaps == min_gaps]
                #if len(min_gaps_alignments) > 1:
                #    LOG.warning(f"Could not resolve candidate BLAST top hits based on min gaps. Will pick the first top hit ({min_gaps_alignments[0].hit_id})")
                #else:
                #    LOG.info(f"Successfully resolved the tie and picked {min_gaps_alignments[0].hit_id} hit")
                #br_alignment = min_gaps_alignments[0]
                #hsp = br_alignment.hsps[0]
            #else:
            br_alignment = blast_record.alignments[0]
            hsp = br_alignment.hsps[0]    
           
        
            percent_identity = round(hsp.identities/hsp.align_length,3)*100
            evalue = hsp.expect
            bitscore = hsp.score
            query_coverage = min(round(hsp.align_length/blast_record.query_length,3)*100,100) #make sure coverage not higher than 100 due to alignment gaps, etc.
            query_length = blast_record.query_length
            species = br_alignment.hit_id
            accession = blast_record.alignments[0].hit_id
            

            LOG.debug(f"TOP hit species={species} query_length={query_length} originalLength={self.origLength} percent_identity={percent_identity} evalue {evalue}")
            if query_length < int(0.6*self.origLength) or percent_identity < 85 or evalue > 1e-200:
                LOG.warning(f"Query length {query_length}bp was either reduced to less than 60% of its original length of {self.origLength}bp or top hit {percent_identity}% identity < 85% or e-value {evalue} > 1e-200")
            #    return "","",0,"",0,"","",""
            if percent_identity < 95:
                LOG.warning(f"The %identity of the top hit ({accession}) is less than 95% ({percent_identity}%) which may lead to incorrect species identification. Check reference database and input.")
                return "","",0,"",0,"","",""

            
            if "|" in accession:
                accession = accession.split("|")[1]

            filename = "query.txt"

            # Open the file with writing permission
            myfile = open(filename, 'w')

            # Write a line to the file
            myfile.write(sequence)


            # Close the file
            myfile.close()

            #blastn_cline = NcbiblastnCommandline(cmd='blastn',query="query.txt", dust='yes',
            #                                     db=os.path.dirname(__file__)+"/reference_database/msr_ref.fa", reward=1, penalty=-2,gapopen=5, gapextend=2,evalue=0.00001, outfmt=5, out="SSUresult.xml")

            #stdout, stderr = blastn_cline()

            #if (os.stat("SSUresult.xml").st_size == 0):
            #    return "","",0,"",0,"","",""

            #result_handle = open("SSUresult.xml", 'r')
            #blast_records = NCBIXML.parse(result_handle)
            #blast_record = next(blast_records)

            #if len(blast_record.alignments) == 0:
            #    return "","",0,"",0,"","",""

            #species = blast_record.alignments[0].hit_id
            

            return bitscore,evalue,query_coverage,query_length,percent_identity, accession, species,sequence


    #outputResults() outputs the results in .txt and .fa file formats
    def outputResults(self, contig, customdatabsename, mode, filetype="abi"):
        #print(self.avgLRI)
       
        self.file.write("\n>Sequence: " + self.name.split(f".{filetype}")[0] + " | ")
        self.tabfile.write(self.name.split(f".{filetype}")[0] + "\t" + mode + "\t")
        
        if (self.species1 == ";>No blast hits." and self.species2 == ";>No blast hits."):
            if self.avgPhredQuality < 10:
                self.tabfile.write("\t\t\t" + "Could not analyze input file. Please check manually." + "\t\t\t\t\t\t\t\n")
            else:
                self.tabfile.write("\t\t\t" + "No blast hits." + "\t\t\t\t\t\t\t\n")
            return
       
        if filetype == "abi" or filetype == "ab1":
            bitscore,evalue,query_coverage,query_length,percent_identity, accession, species, seq = self.blast(self.species1Seq,False)
            bitscore2,evalue2,query_coverage2,query_length2,percent_identity2, accession2, species2, seq2 = self.blast(self.species2Seq,False)
        elif filetype == "fasta" or filetype == "fna"  or filetype == "fa":
            self.avgPhredQuality = 60
            bitscore,evalue,query_coverage,query_length,percent_identity, accession, species, seq = self.blast("".join(self.seq), customdatabsename)
            query_coverage2=query_coverage; species2 = species; percent_identity2 = percent_identity
        else:
            raise TypeError(f"Unsupported filetype '{filetype}' for {self.name}. Aborting")

        if self.avgPhred >= 20 and query_coverage < 50 and query_coverage2 < 50:
            bitscore,evalue,query_coverage,query_length,percent_identity, accession, species, seq = self.blast(''.join(self.fixedSeq),False)


        if mode == 'reverse':
            seq = revcomp(seq)
            seq2 = revcomp(seq2)
        
        if (species == "" and species2 == ""): #or (query_coverage < 50 and query_coverage2 < 50):
            self.tabfile.write("\t\t\t" + "Could not analyze input file. Please check manually." + "\t\t\t\t\t\t\t\n")



        elif species == "":# or query_coverage < 50:
            self.file.write(species2.split("|")[0])

            self.tabfile.write("No\t")
            self.tabfile.write(species2.split("|")[0] + "\t")
            self.tabfile.write(seq2 + "\t")

            if self.avgPhredQuality < 10 and "C.hominis" not in species2:
                self.tabfile.write("Average Phred Quality < 10, could be other potential mixed seqs. Check manually.\t")
                self.file.write("  (Note: Average Phred Quality < 10, could be other potential mixed seqs. Check manually.)")

            elif self.mixed:
                self.tabfile.write("Could be other potential mixed seqs. Check manually.\t")
                self.file.write("  (Note: Could be other potential mixed seqs. Check manually.)")
            else:
                self.tabfile.write(" \t")

            self.file.write("\n" + seq2)

            self.tabfile.write(str(bitscore2)+"\t")
            self.tabfile.write(str(query_length2)+"\t")
            self.tabfile.write(str(query_coverage2) + "%\t")
            self.tabfile.write(str(evalue2) +"\t")
            self.tabfile.write(str(percent_identity2) + "%\t")
            self.tabfile.write(str(accession2)+"\n")

        elif species2 == "":# or query_coverage2 < 50:
            self.file.write(species.split("|")[0])

            self.tabfile.write("No\t")
            self.tabfile.write(species.split("|")[0] + "\t")
            self.tabfile.write(seq + "\t")

            if self.avgPhredQuality < 10 and "C.hominis" not in species:
                self.tabfile.write("Average Phred Quality < 10, could be other potential mixed seqs. Check manually.\t")
                self.file.write("  (Note: Average Phred Quality < 10, could be other potential mixed seqs. Check manually.)")
            elif self.mixed:
                self.tabfile.write("Could be other potential mixed seqs. Check manually.\t")
                self.file.write("  (Note: Could be other potential mixed seqs. Check manually.)")
            else:
                self.tabfile.write(" \t")

            self.file.write("\n" + seq)

            self.tabfile.write(str(bitscore)+"\t")
            self.tabfile.write(str(query_length)+"\t")
            self.tabfile.write(str(query_coverage) + "%\t")
            self.tabfile.write(str(evalue) +"\t")
            self.tabfile.write(str(percent_identity) + "%\t")
            self.tabfile.write(str(accession)+"\n")

        elif species == species2 or ("C.hominis" in species and "C.hominis" in species2):
            #if query_coverage < 50 and query_coverage2 < 50:
            #    self.tabfile.write("\t\t\t" + "Could not analyze input file. Please check manually." + "\t\t\t\t\t\t\t\n")
            #    return
            self.file.write(species.split("|")[0])
            self.tabfile.write("No\t")
            self.tabfile.write(species.split("|")[0] + "\t")

            if percent_identity >= percent_identity2:
                self.tabfile.write(seq + "\t")
             
                if self.avgPhredQuality < 10 and "C.hominis" not in species:
                    self.tabfile.write("Average Phred Quality < 10, could be other potential mixed seqs. Check manually.\t")
                    self.file.write("  (Note: Average Phred Quality < 10, could be other potential mixed seqs. Check manually.)")
                elif self.mixed:
                    self.tabfile.write("Could be other potential mixed seqs. Check manually.\t")
                    self.file.write("  (Note: Could be other potential mixed seqs. Check manually.)")
                else:
                    self.tabfile.write(" \t")

                self.file.write("\n" + seq)
                self.tabfile.write(str(bitscore)+"\t")
                self.tabfile.write(str(query_length)+"\t")
                self.tabfile.write(str(query_coverage) + "%\t")
                self.tabfile.write(str(evalue) +"\t")
                self.tabfile.write(str(percent_identity) + "%\t")
                self.tabfile.write(str(accession)+"\n")
            else:
                self.tabfile.write(seq2 + "\t")

                if self.avgPhredQuality < 10 and "C.hominis" not in species:
                    self.tabfile.write("Average Phred Quality < 10, could be other potential mixed seqs. Check manually.\t")
                    self.file.write("  (Note: Average Phred Quality < 10, could be other potential mixed seqs. Check manually.)")
                elif self.mixed:
                    self.tabfile.write("Average Phred Quality < 10, could be other potential mixed seqs. Check manually.\t")
                    self.file.write("  (Note: Average Phred Quality < 10, could be other potential mixed seqs. Check manually.)")
                else:
                    self.tabfile.write(" \t")

                self.file.write("\n" + seq2)
                self.tabfile.write(str(bitscore2)+"\t")
                self.tabfile.write(str(query_length2)+"\t")
                self.tabfile.write(str(query_coverage2) + "%\t")
                self.tabfile.write(str(evalue2) +"\t")
                self.tabfile.write(str(percent_identity2) + "%\t")
                self.tabfile.write(str(accession2)+"\n")

        else:
            if ("C.parvum" in species2 and "C.hominis" in species and seq.find("TCACAATTAATG") == -1):
                os.system("blastdbcmd -db " + os.path.dirname(__file__)+"/reference_database/msr_ref.fa -entry " + "'{}'".format("C.parvum|KT948751.1") + " -out refseq.fa")

                file = open("refseq.fa", 'r')
                lines = file.readlines()

                refseq = ""

                for i in range (1, len(lines)):
                    refseq += lines[i].strip('\n')

                bitscore,evalue,query_coverage,query_length,percent_identity, accession, species, seq = self.blast(refseq,False)



            elif ("C.parvum" in species and "C.hominis" in species2 and seq2.find("TCACAATTAATG") == -1):
                os.system("blastdbcmd -db " + os.path.dirname(__file__)+"/reference_database/msr_ref.fa -entry " + "'{}'".format("C.parvum|KT948751.1") + " -out refseq.fa")

                file = open("refseq.fa", 'r')
                lines = file.readlines()

                refseq = ""

                for i in range (1, len(lines)):
                    refseq += lines[i].strip("\n")

                bitscore2,evalue2,query_coverage2,query_length2,percent_identity2, accession2, species2, seq2 = self.blast(refseq,False)



            self.tabfile.write("Yes\t")

            self.tabfile.write(species.split("|")[0] + "\t")
            self.tabfile.write(seq + "\t \t")
            self.tabfile.write(str(bitscore)+"\t")
            self.tabfile.write(str(query_length)+"\t")
            self.tabfile.write(str(query_coverage) + "%\t")
            self.tabfile.write(str(evalue) +"\t")
            self.tabfile.write(str(percent_identity) + "%\t")
            self.tabfile.write(str(accession)+"\n")

            self.tabfile.write("\t\tYes\t")
            self.tabfile.write(species2.split("|")[0] + "\t")
            self.tabfile.write(seq2 + "\t \t")
            self.tabfile.write(str(bitscore2)+"\t")
            self.tabfile.write(str(query_length2)+"\t")
            self.tabfile.write(str(query_coverage2) + "%\t")
            self.tabfile.write(str(evalue2) +"\t")
            self.tabfile.write(str(percent_identity2) + "%\t")
            self.tabfile.write(str(accession2)+"\n")


def msr_main(pathlist_unfiltered, forwardP, reverseP, typeSeq, expName, customdatabsename, noheader, verbose):       
     
    tabfile = io.StringIO()

    forwardP = forwardP.replace(' ', '')
    reverseP= reverseP.replace(' ', '')

    pathlist = [path for path in pathlist_unfiltered if re.search("$|".join(definitions.FILETYPES),path)]

    if forwardP and reverseP:
        pathlist = [path for path in pathlist if re.search(forwardP, path) or re.search(reverseP, path)]  # select only files matching the primers
    elif reverseP:
        pathlist = [path for path in pathlist if re.search(forwardP, path)]
    elif reverseP:
        pathlist = [path for path in pathlist if re.search(reverseP, path)]

    #if multi-FASTA file is present in the list slice it up into individual files https://www.metagenomics.wiki/tools/fastq/multi-fasta-format 
    fasta_paths = [path for path in pathlist for fasta_extension in definitions.FASTA_FILETYPES if path.endswith(fasta_extension)]
    for fasta_path in fasta_paths:
        with open(fasta_path,"r") as handle:
            records=list(SeqIO.parse(handle, "fasta"))
            LOG.info(f"File {handle.name} has {len(records)} sequences {[r.name for r in records]}")
            if len(records) > 1:
                for fasta_record in records:
                    tempFastaFilePath = utilities.createTempFastaFiles(expName,fasta_record)
                    pathlist.append(tempFastaFilePath)
                pathlist.remove(fasta_path)     
    
    
    
    if pathlist == []:
        msg = f"Not supported input file(s) found in pathlist {pathlist_unfiltered}. Supported input filetypes are {definitions.FILETYPES}"
        LOG.error(msg)
        raise Exception(msg)
   
    pathlist.sort()

    if not noheader:
        tabfile.write("Sample Name\tType of Sequences\tMixed?\tSpecies\tSequence\tComments\tBit Score\tQuery Length (bp)\tQuery Coverage\tE-value\tPercent Identity\tAccession Number\n")


    file = io.StringIO()
    file.write("\n;>****************************************************************************")
    file.write("\n;>MIXED SEQUENCE ANALYSIS INPUT PARAMETERS:")
    if customdatabsename:
        file.write("\n  ;>Reference File: " + customdatabsename)
    else:
        file.write("\n  ;>Reference File: " + "msr_ref.fa")
    file.write("\n  ;>Forward Primer: " + str(forwardP))
    file.write("\n  ;>Reverse Primer: " + str(reverseP))
    file.write("\n;>****************************************************************************")
    file.write("\n;>Program Results:\n")
    #**************************************************************




    if len(pathlist) == 0:
        exit("No files in the input are matching the forward or reverse primer. Aborting.")


    if len(pathlist)%2 != 0 and typeSeq == 'contig':
        print("ERROR: Uneven number of input files ({}). "
              "Cannot find all paired forward and reverse files. Aborting ...".format(len(pathlist)))
        file.write("\nError: Need to include both forward and reverse sequences of ALL samples to produce contig.")
        file.write("\nSequences files given:\n")

        for i in range(len(pathlist)):
            file.write(pathlist[i])
            file.write("\n")

        tabfile.write("Error: Need to include both forward and reverse sequences of ALL samples to produce contig.\t\t\t\t\t\t\t\t\t\t\t\n")

    elif typeSeq == "contig":
        LOG.info(f"Processing {len(pathlist)} file(s):\n{pathlist}")
        for idx in range(0,len(pathlist),2):
            forward = MixedSeq(file, tabfile, 'contig')
            reverse = MixedSeq(file, tabfile, 'contig')

            f_goodSeq = forward.readFiles(pathlist[idx], True)
            r_goodSeq = reverse.readFiles(pathlist[idx+1], False)

            forward.fixN()
            reverse.fixN()

            f_bitscore,f_evalue,f_query_coverage,f_query_length,f_percent_identity, f_accession, f_species, f_seq = forward.blast(str(''.join(forward.oldseq)),False)
            r_bitscore,r_evalue,r_query_coverage,r_query_length,r_percent_identity, r_accession, r_species, r_seq = reverse.blast(str(''.join(reverse.oldseq)),False)

            if f_goodSeq and r_goodSeq:
                
                if f_species == "" and r_species == "":
                    forward.species1=";>No blast hits."
                    forward.species2=";>No blast hits."
                    forward.outputResults("", customdatabsename, typeSeq)

                else:
                    forward.findHeteroBases(2.0)
                    reverse.findHeteroBases(2.0)

                    reverse.determineAllTypes(customdatabsename)
                    forward.determineAllTypes(customdatabsename)

                    f_bitscore,f_evalue,f_query_coverage,f_query_length,f_percent_identity, f_accession, f_species, f_seq = forward.blast(forward.species1Seq,False)
                    f_bitscore2,f_evalue2,f_query_coverage2,f_query_length2,f_percent_identity2, f_accession2, f_species2, f_seq2 = forward.blast(forward.species2Seq,False)

                    if f_species == "" and f_species2 == "":
                        f_bitscore,f_evalue,f_query_coverage,f_query_length,f_percent_identity, f_accession, f_species, f_seq = forward.blast(str(''.join(forward.fixedSeq)),False)
                        forward.species1Seq = str(''.join(forward.oldseq))
                        forward.species2Seq = str(''.join(forward.oldseq))

                    r_bitscore,r_evalue,r_query_coverage,r_query_length,r_percent_identity, r_accession, r_species, r_seq = reverse.blast(reverse.species1Seq,False)
                    r_bitscore2,r_evalue2,r_query_coverage2,r_query_length2,r_percent_identity2, r_accession2, r_species2, r_seq2 = reverse.blast(reverse.species2Seq,False)

                    if r_species == "" and r_species2 == "":
                        r_bitscore,r_evalue,r_query_coverage,r_query_length,r_percent_identity, r_accession, r_species, r_seq = reverse.blast(str(''.join(reverse.fixedSeq)),False)
                        reverse.species1Seq = str(''.join(reverse.oldseq))
                        reverse.species2Seq = str(''.join(reverse.oldseq))

                    if f_species == "":
                        f_species = f_species2
                        forward.species1Seq = forward.species2Seq
                    if f_species2 == "":
                        f_species2 = f_species
                        forward.species2Seq = forward.species1Seq
                    if r_species == "":
                        r_species = r_species2
                        reverse.species1Seq = reverse.species2Seq
                    if r_species2 == "":
                        r_species2 = r_species
                        reverse.species2Seq = reverse.species1Seq

                    if f_species == "" and r_species == "" and f_species2 == "" and r_species2 == "":
                        f_species = "|."
                        r_species = "|."
                        f_species2 = "|."
                        r_species2 = "|."

                    elif f_species == "" and f_species2 == "":
                        f_species = "|."
                        f_species2 = "|."

                    elif r_species == "" and r_species2 == "":
                        r_species = "|."
                        r_species2 = "|."

                    if ((r_species.split('|')[0] == r_species2.split('|')[0] and f_species.split('|')[0] ==f_species2.split('|')[0]) and f_species==r_species) or ("C.hominis" in r_species and "C.hominis" in r_species2 and "C.hominis" in f_species and "C.hominis" in f_species2):
                        if "C.hominis|GQ183513.11" in f_species and "C.hominis|GQ183513.8" in f_species2:
                            forwardseq = forward.species1Seq
                        elif "C.hominis|GQ183513.11" in f_species2 and "C.hominis|GQ183513.8" in f_species:
                            forwardseq = forward.species2Seq
                        elif f_percent_identity >= f_percent_identity2:
                            forwardseq = forward.species1Seq
                        else:
                            forwardseq = forward.species2Seq


                        if "C.hominis|GQ183513.11" in r_species and "C.hominis|GQ183513.8" in r_species2:
                            reverseseq = reverse.species1Seq
                        elif "C.hominis|GQ183513.11" in r_species2 and "C.hominis|GQ183513.8" in r_species:
                            reverseseq = reverse.species2Seq
                        if r_percent_identity >= r_percent_identity2:
                            reverseseq = reverse.species1Seq
                        else:
                            reverseseq = reverse.species2Seq

                        
                        
                        reverseseq = revcomp(reverseseq)
                        
                        contig = buildContig(forwardseq, reverseseq, forward, reverse)
                        
                        if contig != "":
                            forward.species1Seq = contig
                            forward.species2Seq = contig
                        else:
                            forward.species1Seq = forwardseq
                            forward.species2Seq = reverseseq

                        forward.avgPhredQuality = round((forward.avgPhredQuality + reverse.avgPhredQuality), 2)
                        forward.outputResults(contig, customdatabsename, typeSeq)

                    elif (r_species.split('|')[1].split('.')[0] == f_species.split('|')[1].split('.')[0] and r_species2.split('|')[1].split('.')[0] == f_species2.split('|')[1].split('.')[0]) :
                        reverseseq = revcomp(r_seq)
                        reverseseq2 = revcomp(r_seq2)

                        contig1 = buildContig(f_seq, reverseseq)
                        contig2 = buildContig(f_seq2, reverseseq2)

                        if contig1 != "":
                            forward.species1Seq = contig1
                        if contig2 != "":
                            forward.species2Seq = contig2

                        forward.avgPhredQuality = round((forward.avgPhredQuality + reverse.avgPhredQuality), 2)
                        forward.outputResults(contig1, customdatabsename, typeSeq)

                    elif (r_species.split('|')[1].split('.')[0] == f_species2.split('|')[1].split('.')[0] and r_species2.split('|')[1].split('.')[0] == f_species.split('|')[1].split('.')[0]) :
                        reverseseq = revcomp(r_seq)
                        reverseseq2 = revcomp(r_seq2)

                        if ("C.parvum" in r_species and "C.hominis" in r_species2 and reverseseq2.find("TCACAATTAATG") == -1 and f_seq.find("TCACAATTAATG") ==-1):
                            os.system("blastdbcmd -db " + os.path.dirname(__file__)+"/reference_database/msr_ref.fa -entry " + "'{}'".format("C.parvum|KT948751.1") + " -out refseq.fa")

                            file = open("refseq.fa", 'r')
                            lines = file.readlines()

                            refseq = ""

                            for i in range (1, len(lines)):
                                refseq += lines[i].strip("\n")

                            reverseseq2 = refseq
                            f_seq = refseq


                        elif ("C.parvum" in r_species2 and "C.hominis" in f_species2 and reverseseq.find("TCACAATTAATG") == -1 and f_seq2.find("TCACAATTAATG") ==-1):
                            os.system("blastdbcmd -db " + os.path.dirname(__file__)+"/reference_database/msr_ref.fa -entry " + "'{}'".format("C.parvum|KT948751.1") + " -out refseq.fa")

                            file = open("refseq.fa", 'r')
                            lines = file.readlines()

                            refseq = ""

                            for i in range (1, len(lines)):
                                refseq += lines[i].strip("\n")

                            reverseseq = refseq
                            f_seq2 = refseq

                        contig1 = buildContig(f_seq, reverseseq2)
                        contig2 = buildContig(f_seq2, reverseseq)

                        if contig1 != "":
                            forward.species1Seq = contig1
                        if contig2 != "":
                            forward.species2Seq = contig2

                        forward.outputResults(contig1, customdatabsename, typeSeq)

                    else:
                        LOG.info(f"Did not build a contig but analyzed forward and reverse as separate sequences as species are not identical (forward species:{f_species} {f_species2}, reverse species: {r_species} {r_species2})")
                        forward.outputResults("", customdatabsename, "forward")
                        reverse.outputResults("", customdatabsename, "reverse")
            else:
                LOG.info("Contig would not be built because of quality issues")            


    elif typeSeq=="forward":
        LOG.info(f"Processing {len(pathlist)} file(s):\n{pathlist}")
        for idx, path in enumerate(pathlist):
            filetype = utilities.getFileType(path)
            forward = MixedSeq(file, tabfile, 'forward')
            goodSeq = forward.readFiles(pathlist[idx], True, filetype)
            
            if goodSeq and filetype == "abi":
                forward.fixN()
                forward.findHeteroBases(2.0)
                forward.determineAllTypes(customdatabsename)
                forward.outputResults("", customdatabsename, typeSeq)
            elif filetype == "fasta" or filetype == "fna"  or filetype == "fa":
                forward.origLength = len(forward.seq)
                forward.outputResults("", customdatabsename, typeSeq, filetype)
            
    elif typeSeq=="reverse":
        for idx in range(0, len(pathlist)):
            reverse=MixedSeq(file,tabfile, 'reverse')
            goodSeq = reverse.readFiles(pathlist[idx], False)
            reverse.fixN()

            if goodSeq:
                reverse.findHeteroBases(2.0)
                reverse.determineAllTypes(customdatabsename)
                reverse.outputResults("", customdatabsename, typeSeq)

    experimentName = expName + "_"
    output_report_file_name = experimentName + 'cryptogenotyper_report.fa'
    filename = os.path.join('.', output_report_file_name )

    with open(filename, 'w') as resultFile:
        resultFile.write(file.getvalue())

    output_tabreport_file_name = experimentName + 'cryptogenotyper_report.txt'
    tabfilename = os.path.join('.', output_tabreport_file_name )

    with open(tabfilename, 'w') as resultFile:
        resultFile.write (tabfile.getvalue())
    
    print("\n>>>18S RESULTS  REPORT (only first 10 lines are printed)")
    tabfile.seek(0)
    for idx, line in enumerate(tabfile.read().split("\n")):
        print(line)
        if idx == 10:
            LOG.warning(f"Please check the {output_tabreport_file_name} for the completed output ...")
            break    

    print(">>> FASTA report written to " + os.getcwd()+"/"+output_report_file_name)
    print(">>> Tab-delimited report written to " + os.getcwd() + "/" + output_tabreport_file_name + "\nThe 18S run completed successfully")

    
    #remove files that were made during the analysis
    
    if verbose == False:
        LOG.info("Cleaning the temporary FASTA and BLAST database files (if any)")
        utilities.cleanTempFastaFilesDir()
    LOG.info("The 18S run completed successfully")

if __name__ == "__main__":
    msr_main()
