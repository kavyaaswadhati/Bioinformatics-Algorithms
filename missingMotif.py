#!/usr/bin/env python3

########################################################################
# File: missingMotif.py
#  executable: python3 missingMotif.py -minMotif 3 -maxMotif 8 -cutoff 0.0 < Zm4-genomic.fna > out
# Purpose: Generate a list of motifs in a genome and quantify how under-represented they are. 
#   stderr: None.
#   stdout: List of sequence:reverse    count   expect  z-score
#          
# Author: Kavya Aswadhati
# History: kaswadha 08/20/2011 Created
# group: Danilo Dubocanin, Alex Zee, David Parks
# SCAFFOLD FROM DAVID BERNICK
#               
########################################################################

from fastaReader import FastaReader
import argparse
import math
import sys

########################################################################
# Main
########################################################################

class FindMissingMotifs(object):
    """
    Generate a list of potential motifs in a given sequence, ordered by underrepresentation.
    INPUT: list of fasta sequences, min motif length, max motif length, z-score cutoff
    """
    def __init__(self,fastas,min,max,zCutOff):
        # re-include ,min,max,zCutOff
        self.min = min
        self.max = max
        self.zCutOff = zCutOff
        self.genomeLength = 0
        self.potentialMotifs = {}
        for head,fasta in fastas:
            self.findKmers(fasta)
            self.genomeLength += len(fasta)
        self.getZAndExpected()
        self.sortAndPrint()

    def findKmers(self,sequence):
        '''
        Slide a 'window' along the sequence, where the window size is each length in range(min,max) to find and store counts of all kmers
        ARGS:
            sequence: A sequence in the inputted fasta file
        UPDATES:
            self.potentialMotifs = {(kmer, reverse kmer): [count]}
            - if empty, populated is with kmers and counts found
            - otherwise updates self.potentialMotifs to include kmers from new sequence
        '''

        ## At each base position, this loop finds kmers of all acceptable lengths. 
        for position,base in enumerate(sequence):
            # Window considers 'kmers' down to length 1 and 2, since these counts are required
            # for z-score and expected value calculations. They are not actually valid motifs
            for window in range(self.min-2,self.max+1):
                kmer  = sequence[position:position+window]

                # Kmers that are 'hanging off the end' of the sequence are not considered. 
                ## Note previously used to define the acceptable region to check in based on kmer length
                ## i.e. (len(sequence)-window +1), changed based on comments by Prof. Bernick on 10/7/19
                ## to define region by the largest kmer size
                if position > (len(sequence)-window +1):
                    continue
                # Kmers that include N's are not considered.
                if 'N' in kmer:
                    continue
                # Generate the unique key associated with each kmer and its complement
                this = self.getTuple(kmer)
                if this in self.potentialMotifs:
                        self.potentialMotifs[this][0] += 1
                else:
                    self.potentialMotifs[this] = [1]
           
    def getTuple(self,kmerIn):
        ''' 
        Helper funnction to generate tuple of kmer and its complement.
        ARGS:
            kmerIn: a kmer (string) to generate into a tuple. 
        RETURNS: 
            Tuple of (kmer, reverse kmer). Note kmer vs reverse kmer is determined by alpha heirarchy 
            of each string, where kmer alphabetically preceeds reverse kmer. This way each potential kmer 
            and its complement has only one possible key. 
        '''

        # Generate the complement.
        complement = self.generateComplement(kmerIn)

        # Order and return tuple based on alpha-ordering. 
        if complement < kmerIn:
            return(complement,kmerIn)
        else:
            return(kmerIn,complement)


    def generateComplement(self,kmer):
        ''' 
        Helper function to generate the reverse complement of a found kmer.
        ARGS: 
            kmer: a kmer (string) to generate into a reverse complement
        RETURNS: 
            out: reverse complement of the kmer.
        '''
        temp = []
        for n in kmer:
            if n  in 'A':
                temp.append('T')
            if n in 'T':
                temp.append('A')
            if n in 'G':
                temp.append('C')
            if n in 'C':
                temp.append('G')
        tempii = ''
        tempi = temp[::-1]
        out = tempii.join(tempi)
        return out


    def getZAndExpected(self):
        ''' 
        Calculate z-score and expected value for each found kmer
        UPDATES:
            potentialMotifs[kmer] : [count, expected, zscore]
         '''
        smallKmers = []

        # Iterate through all found kmers and calculate expected counts and zscores
        for k_mer in self.potentialMotifs:
            # Only consider actual motifs (counts for 'kmers' of size 1 and 2 are used in calculations
            # but are not in fact actual motifs)
            if len(k_mer[0]) > 2:

                # Since kmers are stored in tuples, get the actual kmer string
                kmer = k_mer[0]

                # Get keys for the sub-kmers, their counts are used in calculating p
                ksuff = self.getTuple(kmer[1:])
                kpref = self.getTuple(kmer[:-1])
                kmid = self.getTuple(kmer[1:-1])

                # Some sub-kmers don't exist, in these cases expected and zscore cannot be calculated
                if (ksuff not in self.potentialMotifs):
                    continue
                if (kmid not in self.potentialMotifs):
                    continue
                if (kpref not in self.potentialMotifs):
                    continue
                else:
                    # Get count of the actual kmer.
                    s = self.potentialMotifs[self.getTuple(kmer)][0]

                    # Get counts of sub kmers.
                    ksuff_count = self.potentialMotifs[ksuff][0]
                    kpref_count = self.potentialMotifs[kpref][0]
                    kmid_count = self.potentialMotifs[kmid][0]

                    # Calculate probability of the kmer being in the genome.
                    p = float((ksuff_count * kpref_count)/(kmid_count*self.genomeLength))
                    # Calculate mean and standard deviation for occurance of the kmer. (mean is expected)
                    mean = self.genomeLength * p
                    sd = math.sqrt((self.genomeLength*p)*(1-p))
                    zScore = ((s - mean)/sd)

                    # Add calculated values of expected and zscore to dictionary of potential motifs
                    self.potentialMotifs[k_mer].append(mean)
                    self.potentialMotifs[k_mer].append(zScore)

            # Find keys of size 1,2 that were used in calculations, but are invalid motifs.
            else:
                smallKmers.append(k_mer)
        # Remove invalid motifs.
        for k in smallKmers:
                del self.potentialMotifs[k]


    def sortAndPrint(self):
        ''' Print motifs of greatest length and lowest z-score first, kmer pairs printed in alpha order'''

        # Sort dictionary of motifs by zscore (greatest negative first). Note: sorted outputs a list
        zSorted = sorted(self.potentialMotifs.items(), key = lambda k: k[1][2])
        # Sort motifs by length (longest first). 
        lenSorted = sorted(zSorted, key = lambda k: len(k[0][0]), reverse = True)
        
        # Print header followed by final self.potentialMotifs (Note: cases don't match program description,
        # changed for consistency in formatting)
        print('{0}:{1}\t{2}\t{3}\t{4}'.format('sequence', 'reverse', 'count', 'expect','z-score'))
        for mer in lenSorted:
            if mer[1][2] < self.zCutOff:
                # Prints kmer:reversekmer   count   expected    zscore
                print('{0:8}:{1:8}\t{2:0d}\t{3:0.2f}\t{4:0.2f}'.format(mer[0][0], mer[0][1], mer[1][0],mer[1][1],mer[1][2]))
        
def getInput():
    ''' 
    Get input from command line
    ARGS:
        command line input: --minMotif int in range(3,8) --maxMotif int in range(3,8) --cutoff float
    RETURNS:
        minMotif (int), maxMotif(int), cutoff (float)
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-min", "--minMotif", type = int, choices = range(3,9), help = "input min motif length")
    parser.add_argument("-max", "--maxMotif", type = int, choices = range(3,9), help = "input max motif length")
    parser.add_argument("-z","--cutoff", type = float, help = "input z-score")
    args = parser.parse_args()
    return args.minMotif, args.maxMotif,args.cutoff


def main():

    myFastas = []
    # remove the in file parameter to just have input from stdin
    myReader = FastaReader()
    for head, seq in myReader.readFasta():
        myFastas.append([head,seq])
   
    if getInput() is None:
        min = 3
        max = 8
        zScoreCutoff = 0.0
    else:
        min,max, z = getInput()
        # If input is not valid, exit the program
        if max < min:
            print('USAGE ERROR: maxMotif must be greater than minMotif', file=sys.stderr)
            sys.exit(0)
    FindMissingMotifs(myFastas, min,max, z)


if __name__ == "__main__":
    main()

