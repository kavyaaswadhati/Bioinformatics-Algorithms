#!/usr/bin/env python3

########################################################################
# File: problem26.py
#		executable: python3 problem26.py < input.txt > output.txt
# Purpose: Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum.
#		stderr: None.
#   	stdout: List of amino acid strings "peptide", space delimited
#          
# Author: Kavya Aswadhati
# History: kaswadha 12/6/19 Created
# Group: None
#        
########################################################################


import sys

aaMasses = [57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186]

aaIntegerMass = {
	'G':57,'A':71,'S':87,'P':97,'V':99,
	'T':101,'C':103,'I':113, 'L':113,'N':114,
	'D':115, 'K':128,'Q':128, 'E':129,'M':131, 
	'H':137,'F':147, 'R':156,'Y':163,'W':186}


class FindCyclicPeptides(object):
	"""
	Given an ideal experimental spectrum, find a cyclic peptide whose theoretical spectrum
	matches the ideal spectrum.
	INPUT:
		(list) idealSpectrum: List of ints.
	OUTPUTS:
		-STDOUT: List of amino acid strings "peptide", space delimited
	ATTRIBUTES:
		- (str) peptide: A sequence of amino acids.
	"""

	def __init__(self, idealSpectrum):
		'''
		Initialize class-wide objects. Call functions that will execute program function.
		INPUT:
			(list) idealSpectrum: List of ints.
		OUTPUT:
			NONE
		'''

		# List of ints
		self.idealSpectrum = idealSpectrum
		self.massFreq = self.listToFreqDict(idealSpectrum)
		self.validProteins = []
		self.parentMass = self.idealSpectrum[-1]
		self.findCyclicPeptides()
		# Sort and output results.
		self.output()


	def findCyclicPeptides(self):
		'''
		Find all potential cyclic peptides for a given spectrum.
		INPUT: 
			NONE
		RETURNS
			NONE
		'''

		# Stored as a list of strings.
		peptides = []
		
		# Initially populate peptides with peptides of size one.
		for mass in self.idealSpectrum:
			if (mass in aaMasses) and (mass not in peptides):
				peptides.append(str(mass))

		# Keep searching as long as there are valid peptides.
		while len(peptides) > 0:
			# Branch step, generate possible next aas in the peptide.
			peptides = self.expandPeptides(peptides)
			# list.remove() was giving me problems so I pruned via inclusion rather than
			# removal. newPeps saves these potential peptides, needs to be wiped clean each iteration.
			newPeps = []

			# Check each peptide
			for pep in peptides:

				#Gather mass
				massPep = self.getMass(pep)
				# If it has the right mass it could be a valid protein.
				if massPep == self.parentMass:
					# Generate the spectrum it would emit if it was circularized.
					pepCycloSpec = self.generateCyclicSpectrum(pep)
					# If its spectrum matches, it needs to be saved.
					if pepCycloSpec == self.idealSpectrum:
						self.validProteins.append(pep)
					peptides.remove(pep)
				# If the protein was not the right mass, but was still consistent with the Spectrum,
				# save it. One of its expanded forms may be a valid protein.
				elif self.isConsistent(pep) == True:
					newPeps.append(pep)
			# Update peptides with the collected peptides.
			peptides = newPeps


	def isConsistent(self,peptide):
		'''
		Return True if the peptide is consistent with the spectrum, False otherwise
		INPUT:
			peptide
		RETURNS
			boolean: consistent?
		'''
		# Peptide is stored as a string.
		masses = peptide.split('-')

		# Find the mass of the peptide.
		pepMass = self.getMass(peptide)
		
		# Generate the theoretical spectrum.
		linearSpectrum = self.generateLinearSpectrum(masses)
		pepFreq = self.listToFreqDict(linearSpectrum)

		# All masses in the peptide need to be consistent with the spectrum,so store how each mass fared.
		consistent = []

		for mass in pepFreq:

			# All sub-masses in the peptide spectrum must be in the experimental spectrum.
			if mass in self.massFreq:
				# Each mass may only occur in the protein as many times or less than it appeards in the experimental spectrum.
				if pepFreq[mass] <= self.massFreq[mass]:
					if pepMass < self.parentMass:
						consistent.append(True)
					else:
						return False
				else:
					consistent.append(False)
			else:
				return False
		# If all of the masses in the peptide are valid, return True.
		return all(consistent)



	def getMass(self,peptide):
		'''
		Return the mass of a peptide. (Sum its parts)
		INPUT: 
			peptide
		RETURNS
			mass
		'''
		masses = peptide.split('-')
		stm = [int(m) for m in masses]
		mass = sum(stm)
		return mass

	def expandPeptides(self,peptides):
		'''
		Generate new peptides by adding all 18 amino acid masses to all previous peptides.
		INPUT: 
			peptide
		RETURNS
			mass
		'''
		expandedPeps = []
		for p in peptides:
			# newPeps = []
			for aa in aaMasses:
				newP = p + '-' + str(aa)
				# newp = p.append(aa)
				expandedPeps.append(newP)
		return expandedPeps

	def listToFreqDict(self,spectrum):
		'''
		Collapse a spectrum in list form to a dictionary of the frequencies of its subunits.
		INPUT: 
			spectrum
		RETURNS
			dict of spectrum
		'''
		massFreq = {}
		for mass in spectrum:
			if mass in massFreq:
				massFreq[mass] +=1
			else:
				massFreq[mass] = 1
		return massFreq

	def generateCyclicSpectrum(self,pep):
		'''
		Generate the spectrum for the given cyclic protein. Follows algorithm (Cyclic Spectrum) on pg. 214 of Compeau & Pevzner.
		INPUT:
			pep: string of masses. 
		RETURNS:
			NONE
		'''
		# Array that stores the mass of all prefix subsequences so that masses of subpeptides can be easily calculated.
		# Offset by an entry of zero so that calculations can reach "backwards" from the start.
		prefixMass = [0]
		spectrum = [0]
		peptide = pep.split('-')

		# First traverse the protein to populate prefixMass.
		for i in range(0,len(peptide)):
			# The mass of the prefix at a given pos is the preceeding mass and the mass of the amino acid at the position.
			prevMass = prefixMass[i]
			aaMass = int(peptide[i])
			nextMass = prevMass + aaMass
			# Add the calculated prefix.
			prefixMass.append(nextMass)

		# Store the mass of the total protein.
		peptideMass = prefixMass[len(peptide)]

		# Now all subprotein masses can be calculated.
		for i in range(0,len(peptide)):
			# Same as linear spectrum, at each position, calculate the mass of all possible subproteins 
			# begining at this position. 
			# sub-protein mass = (mass of prefix at the end position) - (mass of prefix at the startposition).
			for j in range(i+1,len(peptide)+1):
				subSpectrum = prefixMass[j] - prefixMass[i]
				spectrum.append(subSpectrum)
				# To circularize this, the mass of the subproteins that wrap across the start/stop of this linearized 
				# representation are calculated via:
				# mass wrapped subprotein = mass entire protein - mass subprotein excluded from wrapped section 
				if i>0 and (j<len(peptide)):
					wrapSpectrum = peptideMass - (prefixMass[j] - prefixMass[i])
					spectrum.append(wrapSpectrum)
		return sorted(spectrum)

	def generateLinearSpectrum(self,peptide):
		'''
		Generate the spectrum a linearized protein.
		INPUT:
			pep: string of masses. 
		RETURNS:
			NONE
		'''
		# Array that stores the mass of all prefix subsequences so that masses of subpeptides can be easily calculated.
		# Offset by an entry of zero so that calculations can reach "backwards" from the start.
		prefixMass = [0]
		spectrum = [0]

		# First traverse the protein to populate prefixMass.
		for i in range(0,len(peptide)):
			# The mass of the prefix at a given pos is the preceeding mass and the mass of the amino acid at the position.
			prevMass = prefixMass[i]
			aaMass = int(peptide[i])
			nextMass = prevMass + aaMass
			# Add the calculated prefix.
			prefixMass.append(nextMass)

		# Store the mass of the total protein.
		peptideMass = prefixMass[len(peptide)]

		# Now all subprotein masses can be calculated.
		for i in range(0,len(peptide)):
			# Same as linear spectrum, at each position, calculate the mass of all possible subproteins 
			# begining at this position. 
			# sub-protein mass = (mass of prefix at the end position) - (mass of prefix at the startposition).
			for j in range(i+1,len(peptide)+1):
				subSpectrum = prefixMass[j] - prefixMass[i]
				spectrum.append(subSpectrum)
		return sorted(spectrum)


	def output(self):
		'''
		Sort and print spectrum.
		INPUT:
			NONE
		OUTPUT:
			STDOUT: Spectrum, space delimited
		'''
		print(' '.join([peptide for peptide in self.validProteins]))
			
		

def getInput():
	'''
	Get input from standard in (rosalind file) and format it into a dictionary.
	INPUT:
		Stdin(file): (str) peptide
	RETURNS:
		- (str) peptide: aa string
	'''

	# Get input from standard in.
	inp = sys.stdin
	inList = inp.readline().strip().split(' ')
	idealSpectrum = list(map(int,inList))
	return idealSpectrum


#################################


def main():
	'''
	Call functions to gather input and execute main program task.
	INPUT: 
		NONE
	RETURNS:
		NONE
	'''
	# print(len(aaMasses))

	# Get input from stdin.
	idealSpectrum = getInput()
	
	# Find possible cyclic proteins given an experimental spectrum.
	FindCyclicPeptides(idealSpectrum)



if __name__ == "__main__":
    main()
