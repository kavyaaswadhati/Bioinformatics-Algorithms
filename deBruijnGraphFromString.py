# problem11.py
# Kavya Aswadhati
# kaswadha
# group: none
# Construct the DeBruijn Graph of a string.

import sys

def getInput():
	'''
	Get user input.
	INPUT:
		stdin
	RETURNS:
		desired kmer length (int), text (string)
	'''
	inp = sys.stdin
	for i,line in enumerate(inp):
		if i == 0:
			kmerLen = int(line.strip())
		if i == 1:
			text = line.strip()
	return kmerLen, text

def generateKmerNodes(kmerLen,text):
	'''
	Generate the nodes in the De Bruijn graph (kmer prefixes), as well as a list of kmers (length of kmerLen)
	INPUT:
		kmerLen,text
	RETURNS:
		kNodes (nodes from kmers in dict form), kmers (list)
	'''
	kmers = []
	kNodes = {}
	for position, base in enumerate(text):
		# Slide a window across the sequence of size kmer.
		if position < (len(text) -kmerLen+1):
			kmer = text[position:position+kmerLen]
			prefix = kmer[:-1]
			# Add the prefixes as keys in the dict kNodes
			if prefix not in kNodes:
				kNodes[prefix] = []
			kmers.append(kmer)
	return kNodes, kmers

def makeDeBruijn(kNodes,kmers):
	'''
	Connect kmer nodes with their overlaps
	INPUT:
		kNodes (dict), kmers (list)
	RETURNS:
		kNodes (nodes from kmers in dict form, values are suffix overlaps)
	'''
	for node in kNodes:
		for kmer in kmers:
			# If the prefix is the same as the node, it is an overlap, add the suffix.
			if kmer[:-1] in node:
				kNodes[node].append(kmer[1:])
	return kNodes

def main():

	# Get input from stdin.
	kmerLen, text = getInput()

	# Generate Kmers.
	kNodes,kmers = generateKmerNodes(kmerLen,text)
	# Make the De Bruijn graph.
	outNodes = makeDeBruijn(kNodes,kmers)
	
	# Output results.
	for kmer in sorted(outNodes):
		if len(outNodes[kmer]) > 1:
			sortedOut = sorted(outNodes[kmer])
			print(kmer,'->',','.join([x for x in sortedOut]))
		elif len(outNodes[kmer]) == 1:
			print(kmer,'->',outNodes[kmer][0])
		else:
			continue

if __name__ == "__main__":
    main()