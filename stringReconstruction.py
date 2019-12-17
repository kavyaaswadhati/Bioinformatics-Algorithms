# problem14.py
# Kavya Aswadhati
# kaswadha
# group: none

'''
Reconstruct a String from its k-mer Composition by generating a Eulerian path from a DeBruijn graph
executable : python3 rosalind_14.txt < inFile.txt
'''

import sys
from numpy import random

################################################################################################################
################################################################################################################


class EulerianPath(object):
	"""
	Find an Eulerian Path from a directed graph.
	INPUT:
		Directed graph in dictionary form.
	OUTPUT:
		Stdout: Eulerian path connected via '->'
	"""
	def __init__(self,directedGraph):
		'''Initialize class-wide objects'''
		# Will be used to store and update the resultant path
		self.finalPath = []
		# Will be updated for a continuous reference of nodes with untraversed edges
		self.availableNodes = []

		# Created class object with user-inputted directed graph
		self.directedGraph = directedGraph
		# Number of edges in directed graph, necessary for determining the head and tail of the graph
		self.inpointers = self.getNumInpointers()
		# Call to method that does the actual work of making the Eulerian Path
		# self.makeEulerian()
	
	def makeEulerian(self):
		'''
		Generate sub-Eulerian paths and piece them together to include the whole dataset.
		OUTPUT:
			Stdout: Eulerian path connected via '->'
		'''
		# Find the values of the first and last path values
		head = self.findStart()
		tail = self.findEnd()

		# Add a key that points the tail to the head to circularize the graph 
		if tail not in self.directedGraph:                                              
			self.directedGraph[tail] = [head]
		else:
			self.directedGraph[tail].append(head)

		# Get an initial path that starts at the start
		self.finalPath = self.findPath(head)
		# Keep generating and inserting all available sublists as long as there are nodes with 
		# untraversed edges.
		while (len(self.availableNodes) > 0):
			# Look for subpaths branching off of nodes that have already been traversed
			# (often there are more untraversed nodes than traversed, checking the nodes already
			# in the list prevents you from creating unattached Eulerian cycles).
			for node in self.finalPath:
				if node in self.availableNodes:
					# Path based off of an available node.
					subpath = self.findPath(node)
					# Insert the subpath at the first instance of the node in the larger path
					pos = self.finalPath.index(node)
					self.finalPath[pos:pos] = subpath[:-1]
					# Update which nodes are available, as some new edges have been traversed in the course of 
					# generating a subpath
					self.updateAvailableNodes()
		# The last index is a repeat of the start, so chop it off.
		toOut = self.finalPath[:-1]
		# Output result to standard out
		# print('->'.join([x for x in toOut]))
		return(toOut)

	def updateAvailableNodes(self):
		'''
		Helper function to check for additional available nodes. 
		'''
		# Reset to empty every time
		self.availableNodes = []
		# Check each node for available edges.
		for node in self.directedGraph:
			if len(self.directedGraph[node]) > 0:
				# Add found nodes to list
				self.availableNodes.append(node)

	def findPath(self,start):
		'''
		Traverse nodes in directed graph to create and output a Eulerian Path.
		INPUT:
			The value of a node to start at.
		RETURNS:
			A Eulerian cycle in the form of a list.
		'''
		node = start
		# Initialize a path, and add the first node.
		path = []
		path.append(node)
		# While the current node has an available edge...
		while(len(self.directedGraph[node])> 0):
			if len(self.directedGraph[node]) == 1:
				# If it only has one edge, at this point the edge will be used up, and the node is no longer viable.
				if node in self.availableNodes:
					self.availableNodes.remove(node)
				# There is only one possible edge, so traverse it.
				nextVal = self.directedGraph[node][0]
				# Once the edge has been traversed, it is no longer viable, so remove it.
				self.directedGraph[node].remove(nextVal)
				# Add the value to the path (this is "traversal").
				path.append(nextVal)
				# "Point" to this current node.
				node = nextVal
			else:
				# If it has more than one available edge, it will have available edges and need to be traversed again.
				self.availableNodes.append(node)
				# choose a random next edge to traverse
				choice = random.randint(0, len(self.directedGraph[node]))
				nextVal = self.directedGraph[node][choice]
				# Once the edge has been traversed, it is no longer viable, so remove it.
				self.directedGraph[node].remove(nextVal)
				# Add the value to the path (this is "traversal")
				path.append(nextVal)
				# "Point" to this current node
				node = nextVal
		return(path)

	def getNumInpointers(self):
		'''
		Find the number of inpointers a given path has, populate global dict so this
		this information can be easily queried.
		'''
		inpointers = {}
		# Check the edges of each node in the directed graph
		for node in self.directedGraph:
			if node not in inpointers:
				# A node may never appear in the edges (the start), even so it must have a key
				inpointers[node] = 0
			for edge in self.directedGraph[node]:
				# Update the dict with the number of inpointers a given value has.
				if edge not in inpointers:
					inpointers[edge] = 1
				else:
					inpointers[edge] += 1
		return inpointers

	def findEnd(self):
		'''
		Find the last node in the Eulerian path. The tail is defined as having more in-edges than out-edges.
		In some cases the tail may not even have any out-edges. 
		RETURNS:
			value of end node (int)
		'''
		for node in self.inpointers:
			# A value could be the tail because it has no outpointers
			if node not in self.directedGraph:
				tail = node
			# Or because it has more inpointers than outpointers.
			elif self.inpointers[node] > len(self.directedGraph[node]):
				tail = node
		return tail


	def findStart(self):
		'''
		Find the first node in the Eulerian path. The head is defined as having more out-edges than in-edges. 
		RETURNS:
			value of start node (int)
		'''
		for path in self.directedGraph:
			if (self.inpointers[path]) < len(self.directedGraph[path]):
				start = path
		return start

################################################################################################################
################################################################################################################

class DeBruijnGraph(object):
	"""Construct the De Bruijn Graph of a collection of kmers"""
	def __init__(self):
		# Generate Kmers
		self.kmers = self.getInput()
		self.kNodes = {}
		self.generateKmerNodes()
		# Construct De Bruijn Graph
		# self.makeDeBruijn(kNodes,kmers)


	def getInput(self):
		'''
		Get user input.
		INPUT:
			stdin
		RETURNS:
			list of kmers
		'''
		inp = sys.stdin
		kmers = []
		rawInp = []
		for line in inp:
			kmer = line.strip()
			rawInp.append(kmer)
		kmerLen = rawInp[0]
		kmers = rawInp[1:]
		return kmers

	def generateKmerNodes(self):
		'''
		Generate the nodes in the De Bruijn graph (kmer prefixes), as well as a list of kmers (length of kmerLen)
		INPUT:
			kmers (list)
		RETURNS:
			kNodes (nodes from kmers in dict form), kmers (list)
		'''
		for kmer in self.kmers:
			prefix = kmer[:-1]
			# Add the prefixes as keys in the dict kNodes
			if prefix not in self.kNodes:
				self.kNodes[prefix] = []

	def makeDeBruijn(self):
		'''
		Connect kmer nodes with their overlaps
		INPUT:
			kNodes (dict), kmers (list)
		RETURNS:
			kNodes (nodes from kmers in dict form, values are suffix overlaps)
		'''
		for node in self.kNodes:
			for kmer in self.kmers:
				# If the prefix is the same as the node, it is an overlap, add the suffix.
				if kmer[:-1] in node:
					self.kNodes[node].append(kmer[1:])
		return self.kNodes

################################################################################################################
################################################################################################################ 
def getStringFromPath(kmers):
	'''
	Generate a consensus string from an ordered list of kmers.
	INPUT:
		list of kmers
	RETURNS:
		consensus string
	'''
	# Start by adding the first kmer to the consensus
	consensus = kmers[0]
	for i, kmer in enumerate(kmers):
		# Skip the first kmer.
		if i == 0:
			continue
		# Add the last character in the kmer to the consensus.
		consensus += kmer[(len(kmer)-1)]
	return consensus
################################################################################################################
################################################################################################################ 

def main():

	# Get input from stdin
	myGraph = DeBruijnGraph()
	directedGraph = myGraph.makeDeBruijn()
	# Generate a Eulerian Path
	myPath = EulerianPath(directedGraph)
	path = myPath.makeEulerian()
	# Convert it to a consensus string
	consensus = getStringFromPath(path)
	print(consensus)



if __name__ == "__main__":
    main()
