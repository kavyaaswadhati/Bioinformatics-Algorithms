#!/usr/bin/env python3

########################################################################
# File: randomizedMotifSearch.py
#   executable: python3 randomizedMotifSearch.py -i int -p int -k int -r (opt) -g (opt) -m (opt) < fasta > out
#   purpose: Find the most likely promoter sequence motif (consensus sequence) 
#          from a collection of 50bp sequences upstream from CRISPR arrays. 
#   stderr: Usage warnings
#   stdout: score, consensus motif.
#          
# Author: Kavya Aswadhati
# History: kaswadha 10/8/19 Created
# group: Danilo Dubocanin, Alex Zee
# SCAFFOLD FROM DAVID BERNICK
#               
########################################################################

from fastaReader import FastaReader # by David Bernick
import argparse
import math
import sys
import random
import operator

########################################################################
# Main
########################################################################

class FindConsensusMotif(object):
    """
    Find the consensus motif amongst a list of sequences.
    """
    def __init__(self,i,p,k,randomize,gibbs,printMotif):
        self.numseqs = 0
        self.iterationNum = i
        self.pseudocounts = p
        self.kmerLen = k
        self.randomized = randomize # bool
        self.gibbs = gibbs # bool
        self.printMotif = printMotif # bool
        self.dnaSequences = {} # {name:seq,...}
        self.shuffledSeqs = []
    ######################################################################################################################

    def shuffleSequence(self,sequence):
        '''Helper function to randomize input sequences'''
        seqList = list(sequence)
        random.shuffle(seqList)
        shuffledSeq = ''.join(seqList)
        return shuffledSeq

    def generateProfile(self,motifs):
        '''
        Generate a profile based on input of list of motifs
        INPUT:
            list of motifs from each sequence in self.dnaSequences
        RETURNS:
            dict: key = base, value =  list of len(self.kmerLen) of count of each base at each position.
        '''
        # Make a new profile each time this method is called
        profile = []
        # Iterate through each position in the motif
        for pos in range(0,self.kmerLen):
            # In each motif, get the counts of each nucleotide at each position
            # Initialized with the pseudo counts and set uo so that base counts can be added as probabilities
            posProfile = {
                        'A' : self.pseudocounts/((4*self.pseudocounts)+len(self.dnaSequences)),
                        'C' : self.pseudocounts/((4*self.pseudocounts)+len(self.dnaSequences)),
                        'T' : self.pseudocounts/((4*self.pseudocounts)+len(self.dnaSequences)),
                        'G' : self.pseudocounts/((4*self.pseudocounts)+len(self.dnaSequences))}
            # Get the counts of each base at the position in question in each motif
            for motif in motifs:
                posProfile[motif[pos]] += 1/((4*self.pseudocounts)+len(self.dnaSequences))
            # Profile is now a list of dictionaries.
            profile.append(posProfile)
        return(profile)


    def getMostProbableMotifs(self,profile):
        ''' 
        Find the most probable motif in each sequence, based on a given probability profile (derived from a set of motifs.
        INPUT:
            A probability profile.
        RETURNS: 
            List of the most probable motif for each sequence
        '''
        motifs = []
        # Go through each sequence
        for seq in self.dnaSequences:
            # Each sequence will have its own "best" motif
            bestKmerProb = [0,'']
            # Move a "window" of size kmer across the length of the sequence to check each potential kmer.
            for i, base in enumerate(self.dnaSequences[seq]):
                # Ensure that the "window" is valid (not passing the end of the sequence)
                if i < (len(self.dnaSequences[seq]) - (self.kmerLen)):
                    kmer =  self.dnaSequences[seq][i:i+self.kmerLen]
                    prob = 1
                    # Find the product of the probabilities that each base in the kmer will appear at the position it is in.
                    for position, nuc in enumerate(kmer):
                        #  Get the nucleotide probability for each position in each generated kmer from the profile. 
                        prob = (profile[position][nuc]) * prob
                    # Update the "best" kmer based on this probability (highest probability is best)
                    if prob > bestKmerProb[0]:
                        bestKmerProb = [prob, kmer]
            # Add the found best kmer in each seq to te list of motifs, and return.
            motifs.append(bestKmerProb[1])
        return motifs

    def consensus(self,motifs):
        '''
        Find the consensus motif from a given list of motifs.
        INPUT:
            List of motifs.
        RETURNS:
            Consensus motif (string)
        '''
        consensus = ''
        # Generate profile to base consensus off of.
        profile = self.generateProfile(motifs)
        for pos in profile:
            # Get the base with the highest probability at each position.
            bestBase = max(pos.items(), key=operator.itemgetter(1))[0]
            #  Add the best base at the position to the consensus
            consensus += bestBase
        return consensus

    def score(self,motifs):
        ''' 
        Get score (entropy/encoded cost) from a motif set.
        INPUT:
            List of motifs.
        RETURNS:
            Entropy score of the motif set (float).
        '''
        profile = self.generateProfile(motifs)
        entropy = 0
        for pos in profile:
            probs = []
            # sum of entropies for each nuc at this position
            for nuc in pos:
                probs.append(pos[nuc])
            entropy += self.getShannonsEntropy(probs)
        return entropy


    def getShannonsEntropy(self,probs):
        ''' Helper function to calculate entropy at a single position'''
        # entropy is calculated based on every base at that position, not just the best one. 
        sum = 0
        for prob in probs:
            # entropy describes encoding costs
            entropy = (prob) * math.log2(prob)
            sum+= entropy
        return (-1 * sum)
  

    def getKmersRandomly(self,seq):
        '''
        Choose a random kmer from inputted sequence.
        INPUT:
            A single sequence from the set of sequences.
        RETURNS:
            A randomly chosen kmer (String).
        '''
        pos = random.randint(1,len(seq)-self.kmerLen)
        # Return the kmer the begins at that position
        return seq[pos:pos+self.kmerLen]

    

    ######################################################################################################################
    ######################################################################################################################

    def runRandomized(self,motifs):
        '''
        Execute the randomized motif search algorithm
        INPUT:
            A set of motifs to start with.
        RETURNS:
            The set of motifs with the lowest entropy score
        '''
        # Algorithm from Pevzner and Compeau
        # Set the best motifs to be the initial ones
        bestMotifs = motifs
        # Keep searching for motifs until you find a set better than the best.
        while True:
            # Generate a proability profile from the motifs at hand
            profile = self.generateProfile(motifs)
            # Based on this profile, find the most likely kmer from each sequence.
            motifs = self.getMostProbableMotifs(profile)
            # Compare the encoding cost of this set of motifs from the best set (their dissimilarity).
            # As long as you can find better sets of motifs, keep searching.
            if self.score(motifs) < self.score(bestMotifs):
                bestMotifs = motifs
            else:
                return bestMotifs

    def runGibbs(self):
        ''' Place holder for future implementation of gibbs sampling'''
        pass


    ##################################################################################################################
    ##################################################################################################################

    def execute(self):
        '''
        Execute the main task of the program
        INPUT:
            None
        OUTPUT:
            stdout: the consensus motif found and its associated score
        '''
        if self.gibbs:
            # This option is not yet available
            self.runGibbs()
        else: 
            bestMotifs = []
            for i in range(0,self.iterationNum):
                # start with a random set of motifs
                motifs = []
                for seq in self.dnaSequences:
                    # kmer = self.getKmersRandomly(self.dnaSequences[seq])
                    motifs.append(self.getKmersRandomly(self.dnaSequences[seq]))
                bestMotifs = motifs
                theseMotifs = self.runRandomized(motifs)
                if self.score(theseMotifs) < self.score(bestMotifs):
                    bestMotifs = theseMotifs
            print(self.consensus(bestMotifs),self.score(bestMotifs))


######################################################################################################################    
        
def getInput():
    ''' 
    Get input from command line
    ARGS:
        command line input: -i int -p int -k int -r (opt) -g (opt) -m (opt)
    RETURNS:  terationNum, pseudocountNum, motifLength, randomize, gibbs, printMotif

    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--iterationNum", type = int, help = "enter number of time to iterate")
    parser.add_argument("-p", "--pseudocountNum", type = float, help = "enter number of pseudocounts")
    parser.add_argument("-k","--motifLength", type = int, help = "enter motif length to search for")
    parser.add_argument("-r","--randomize", action = 'store', help = "include to enable sequence randomization")
    parser.add_argument("-g","--gibbs", action = 'store', help = "include to perform search with Gibbs Sampling")
    parser.add_argument("-m","--printMotif", action = 'store', help = "print motif and names of contributing sequences")
    args = parser.parse_args()
    return args.iterationNum, args.pseudocountNum, args.motifLength, args.randomize, args.gibbs, args.printMotif

###############################################################################################################################
def main():

    if getInput() is None:
        # Hard code default params to run if no user input is given
        i = 10000
        p = 1
        k = 13
        randomize = False
        gibbs = False
        printMotif = False
    else:
        i, p, k, randomize, gibbs, printMotif = getInput()
        if randomize or gibbs or printMotif:
            print('Your selected options are not available yet', file=sys.stderr)
            sys.exit(0)
    # Instantiate class in which motif finding will occur.
    findMotif = FindConsensusMotif(i,p,k,randomize,gibbs,printMotif)

    # Provide DNA sequences to findMotif.
    myReader = FastaReader()
    for head, seq in myReader.readFasta():
        findMotif.dnaSequences[head] = seq
    
    # Find consensus motif.    
    findMotif.execute()

if __name__ == "__main__":
    main()

