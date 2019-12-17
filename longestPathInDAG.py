# Kavya Aswadhati
# problem15.py
# Find the Length, Longest Path in a Directed Acyclic Graph via topological sorting, scoring, and backtracking
# executable:
#	python3 problem15.py < input.txt > output.txt
# Bioinformatics Models and Algorithms, Professor David Bernick
# Group: None.

import sys
import random
import math


class LongestPathInDAG(object):
	'''
	Find the longest path through a DAG and its length.
	ATTRIBUTES:
		- dag: dictionary of Nodes in DAG stored by value
		- source: key for the source Node
		- sink: key for the sink Node
		- setIncomingEdges(): Helper function to set all the incoming edges for each node in the dag.
		- topologicalOrder: Generate a valid topological ordering of the DAG in list form.
		- generateLengthLongestPath(): Initialize Node scores and generate the length/longest path through the DAG.
		- backtrack(): Traverse the scored, ordered DAG from sink to source to determine the Longest Path.
		- output(): Output results to stdout.


	'''
	def __init__(self,dag,source,sink):
		'''
		Initialize class attributes and call functions to execute main class function of determining length and identity of the 
		longest path through a given DAG.
		INPUT:
			- (dict) dag: dictionary of Nodes stored by value
			- (str) source: key for the source Node
			- (str) sink: key for the sink Node
		OUTPUT:
			STDOUT: (line 1) Length of longest path. 
					(line [n:]) Longest path from source to sink. 
		'''
		# Initialize global variables (dag, source, sink)
		self.dag = dag
		self.source = source
		self.sink = sink
		self.invalidNodes = []

		# Update Nodes in DAG to include their inedges.
		self.setIncomingEdges()

		# Get a valid topological order for the DAG. (There may be multiple, but any will work).
		topo = self.topologicalOrder(dag)

		# Reset the Incoming edges as the topologicalOrder() algorithm is dependent of removing inedges.
		self.setIncomingEdges()

		# Initialize Node scores and determine the length of the longest path through the DAG. 
		length = self.generateLengthLongestPath(topo)

		# Follow the longest path from sink to source to save and outout the longest path.
		longestPath = self.backtrack()

		# Output results to stdout.
		self.output(length,longestPath)


	def output(self,length,longestPath):
		'''
		Print program results to stdout.
		INPUT:
			- (int) Length of longest path.
			- (list) Longest path.
		OUTPUTS:
			STDOUT: (line 1) Length of path
					(line [n:]) Nodes in forward order longest path -- '->' delimited
		'''
		print(length)
		# List comprehension to format output as '->' delimited.
		print('->'.join([x for x in longestPath]))


	def setIncomingEdges(self):
		'''

		Helper function to set all the incoming edges for each node in the dag.
		INPUT: 
			NONE
		RETURNS:
			NONE

		'''

		# Map the outedges for all nodes in the dag to their corresponding inedges.
		toPrune = []
		for d in self.dag:
			# Each node may have multiple outedges.
			for outedge in self.dag[d].outedges:
				# Here fwd should be formatted 'nextNode:edgeWeight'. Deconstruct the outedge to 
				# generate the inedge that it would correspond to.
				node, weight = outedge.split(':')

				# The node to with the outedge points my not be in the DAG dict (if it were a dead-end)
				if node in self.dag:
					# Reconstruct the corresponding inedge, 'currentNode:edgeWeight'
					inedge = d + ':' + weight
					# Add this edge to Node attribute inedges for the node the outedge was pointing to. 
					self.dag[node].addPrevEdge(inedge)

	def topologicalOrder(self,dag):
		'''

		Generate a list that is a valid topological order (each node comes after all of the nodes that point to it
		and before all of the nodes it points to) for the DAG. Multiple such orderings exist.
		INPUT: 
			A dictionary that contains all of the node objects in the Directed Acyclic Graph. 
				Each node in the DAG has a value (these are the keys in the dict), a list of incoming edges, 
				a list of outgoing edges, a pointer to its preceding node (uninitialized), and a score (uninitialized).
		RETURNS:
			A valid topological order of the nodes in the DAG (list). The first value of the list is always the source,
			while the last the value in the list is not always the sink.

		'''

		# The topological list that will be populated and returned.
		topoList = []
		# List of nodes that are available for consideration as the next node in the topological order.
		# These are node with no inedges (no nodes precede them so they can be the next node).
		candidates = []
		# Boolean to ensure that the first node considered and added is the source.
		first = True

		# Populate the list of candidates (node objects).
		for node in dag:
			# All nodes with no inedges are candidates. Though on the first one we will only consider the 
			# source node, there are other nodes that start with no inedges and they need to be saved.
			if (dag[node].inedges == []): 
				candidates.append(node)

		# While there are still candidates (nodes that have not been ordered) keep adding to the topological list
		while len(candidates) > 0:

			# Make sure the source is the first node in the list.
			if first == True:
				# currentVal is the next value to be added to the topological list, in this case it is the source node.
				currentVal = dag[self.source].val
				# The first node has been considered so first is false.
				first = False
			else:
				# If this is not the first node, any node in candidates is a valid choice as the next node in the topoList.
				# Choose an arbitrary node in candidates (choice is the index in the list of candidates).

				# The source node is inviable after the first pass.
				if self.source in candidates:
					candidates.remove(self.source)

				choice = random.randint(0,len(candidates)-1)
				# currentVal is the integer value of the Node object at position choice in candidates.
				currentVal = candidates[choice]
				# The chosen node is no longer a candidate.
				candidates.remove(candidates[choice])
			
			# Add the node value to the topological list.
			topoList.append(currentVal)
			
			# In order to move on to the next set of candidate nodes, the outedges of the current node (connections
			# to the nodes following it) need to be removed, and the candidates list needs to be updated appropriately.

			# fwds is the list out outgoing edges of the current node in consideration (currentVal)
			fwds = dag[currentVal].outedges

			# For each outgoing edge from currentVal, remove the edge. If the node to which it pointed as no other inedges,
			# add it to the list of candidates.
			for i, outedge in enumerate(fwds):
				# Some nodes in the DAG may have inedges but no outedges. These are dead-ends and do not have key values in self.dag.
				# Assume that nodes are not dead-ends until they are proven to be.
				deadEnd = False

				# Outedges are stored in nodes as a string of 'nextNode:edgeWeight', parse this so it's more easy to access.
				node, weight = outedge.split(':')

				# The inpointer for the node pointed to by the current outedge would be the value of the current node 
				# plus the weight of the edge connecting them. Contstruct the string this object is stored as in
				# the Node object attribute "inedges"
				inpointer = currentVal+':'+weight

				# Check if the node is a dead-end
				if node not in dag:
					# Node has inpointers but no outpointers (is a dead end), give it a key. None of the Node attributes are 
					# initialized, but it is added to both the global and local DAG for downstream use. 
					dag[node] = Node(node)
					self.dag[node] = Node(node)
					# It is a dead end.
					deadEnd = True
				
				# If the node is not a dead-end and it has inedges, we need to remove the edge between the current node and its 
				# next one.
				if(deadEnd == False) and (len(dag[node].inedges) > 0):
					# Remove the inpointer from the local dag.
					if inpointer in dag[node].inedges:
						dag[node].inedges.remove(inpointer)

				# If after removing the edge connecting the current node and the corresponding outedge node the outedge node
				# does not have any inedges, it is a candidate.
				if (dag[node].inedges == []): 
					candidates.append(node)
		# print(topoList)

		# Return the topological ordering of the DAG.
		return topoList


	def generateLengthLongestPath(self, topoList):
		'''
		Populate the scores of all of the nodes in the DAG by following the topological ordering of the graph. The value of 
		score at the sink is then the length of the longest path through the DAG from source to sink.
		INPUT:
			topoList: A list containing a topological ordering of the DAG
		RETURNS: 
			(int) The length of the longest path through the DAG (score of sink + 1)
		'''

		# Initialize the scores of all nodes in the DAG
		for node in self.dag:
			# All nodes other than the source are initialized as negative infinity to preserve
			# the topological order. Since there are nodes other than the source with an indegree of 0, this ensures that 
			# only paths begining at the source have valid paths. Other paths will have scores of -inf as the addition of 
			# an int to -inf is still -inf.
			self.dag[node].score = -math.inf
		# Set the source node score to 0.
		self.dag[self.source].score = 0

		#Traverse the DAG following the topological order and assign scores.
		for i,b in enumerate(topoList):
			# Temporary variables to store the max inpointer. These are the values required for 
			# determining the score of the current node.
			maxScore = -math.inf

			# Save the value of the best node so that we can populate the prev pointer in the Node object. 
			bestNode = None
			if i == 0:
				# If i == 0 we are at the source node, it has no predecessors, its score is 0
				maxScore = 0
				self.dag[b].score = maxScore 
				# The first node has no previous node.
			else:
				# For every other node (up until the sink) its score is equal to the score
				# of its predecessor with the largest score plus the weight of the edge from that 
				# node to the current one.
				for inpointer in self.dag[b].inedges:
					# a is the id of node a, weight is the weight of the edge from a to b.
					a, weight  = inpointer.split(':')
					currentScore = self.dag[a].score + int(weight)
					# Check if the score of this inpointing node a is the largest of all the inpointers.
					if currentScore > maxScore:
						# If it is the largest score found so far, save its score and the weight of 
						# its edge leading to b
						maxScore = currentScore
						bestNode = a
				# The score of node b is = to the score of its predecessor with the largest score + the weight
				# of the edge from that predecessor (a) to the current node (b).
				self.dag[b].score = maxScore  
				self.dag[b].previous = bestNode

			if b == self.sink:
				# When we reach the sink we dont care about nodes that may fall after it in the topo graph,
				# so we return the score of the sink to traverse backwards and find the longest path via backtrack()
				# The length of the longest path through the DAG is the score of the sink plus one since
				# the score of sink is 0 but it must be counted as the first node.
				return(self.dag[self.sink].score)
	
	def backtrack(self):
		'''
		Backtrack from the sink to the source to find the ID of a longest path through the dag.
		INPUT:
			NONE
		RETURNS:
			(list) A ordered list of the node values comprising the longest path through the dag from source to sink.
		'''

		# Saves the path as it is first encountered in reverse order.
		revLongestPath = []

		# Start at the sink.
		currentNode = self.sink
		revLongestPath.append(currentNode)

		# Keep traversing backwards through the dag until you read the source node.
		while True:
			# The previous node is the node amongst the inedges whose score is equal to 
			# score of the current node minus the weight of the edge connecting them.
			# (There could potentially be more than one if there is a symmetrical bubble in the
			# DAG, this implementation simply considers the last correct previous node it encounters.)
			# print('current',currentNode)
			for inedge in self.dag[currentNode].inedges:
				# Deconstruct the inedge object for ease of use.
				inNode, weight = inedge.split(':')
				if self.dag[inNode].score == (self.dag[currentNode].score - int(weight)):
					# If so add it to the path...
					revLongestPath.append(self.dag[inNode].val)
					# ... and move the pointer currentNode to this previous node.
					currentNode = inNode

				# Once you reach the source Node, return to exit this loop.
				if currentNode == self.source:
					# Flip revLongestPath to get the forward order longest path through the DAG.
					longestPath = revLongestPath[::-1]
					# Return this longest path
					return longestPath
		

############################################################################################
############################################################################################

def getInput():
	'''
	Get input from standard in (rosalind file) and format it into a dictionary.
	INPUT:
		Stdin(file): (line 1) source (int), (line 2) sink(int), (lines [n:]) adjacency list of [nodeA --> nodeB:edgeWeight, ...]
	RETURNS:
		Dictionary of nodes: {node value (int): Node (object), ...}
	'''

	# Get input from standard in.
	inp = sys.stdin

	# Directed Acyclic Graph in dictionary form to populate.
	dag = {}
	toPrune = []
	# Parse through each line in the input file. Needs to be enumerated because the first two lines are treated differently.
	for i,line in enumerate(inp):

		# The first two lines in the file are always the source then the sink (rosalind file format).
		if i == 0:
			# Save the value of the source, a Node object for it will be created once it is encountered in the adjacency list.
			source = line.strip()
		elif i ==1:
			# Save the value of the sink.
			sink = line.strip()
			# The sink has no outedges so it would not be initialized but it needs to be entered in the DAG.
			dag[sink] = Node(sink)
			# Some rosalind data sets have the same source and sink node, these are invalid data sets. (A DAG cannot be a single node).
			if source == sink:
				# Print an error and exit.
				print('The source and sink node are the same, there is no path', file=sys.stderr)
				print('0')
				print(source)
				sys.exit(0)

		# All other lines in the file are an adjacency list of [nodeA --> nodeB:edgeWeight, ...]
		# Note that Nodes may have many outedges and inedges, so they may be encountered multiple times.
		# Inedges are not initialized here.
		else:
			# Splitting along the arrow gives us a list (a) of [node,'nextNode:edgeWeight'].
			a = line.strip().split('->')
			nextNode, edgeWeight = a[1].split(':')

			# Save all nodes that preceed the source to prune later.
			if nextNode == source:
				toPrune.append(a[0])

			# Conserve the colon designation between the node pointer and the weight.
			# If the node has multiple outpointers (has been found before) append the second 
			# pointer to the same Node object's list of outpointers (Node.outedges)
			if a[0] in dag:
				# Call Node object helper function that will add the edge to the list outedges.
				dag[a[0]].addNextEdge(a[1])

			# If this is the first time the node has been seen, initialize the dictionary 
			# with its value and create a Node object with its value.
			else:
				dag[a[0]] = Node(a[0])
				# Initialize the new Node object with its outpointer.
				dag[a[0]].addNextEdge(a[1])
		if source not in dag:
			dag[source] = Node(source)

	# Prune out nodes preceeding the source.
	[dag.pop(n) for n in toPrune]

	## Note: this input parsing method does not deal with DAGs with dead-end branches (nodes with inedges but no outedges).

	# Return the dag dictionary and the values of the source and sink.
	return dag,source,sink
############################################################################################

############################################################################################

class Node(object):
	'''
	Create a data structure that will be used as the node units of a Directed Acyclic Graph.
	INPUT:
		(str) value
	RETURNS:
		NONE
	ATTRIBUTES:
		- (str) value
		- (int) score
		- (Node Object) previous
		- (list) outedges
		- (list) inedges
		- addNextEdge(edge): Helper function to add an edge to self.outedges
		- addPrevEdge(edge): Helper function to add an edge to self.inedges
		- getStatus(): Helper function for debugging that prints the status of all Object 
					   attributes to standard out.
	'''
	def __init__(self,val):
		'''
		Initialize Node Object attributes.
		INPUT:
			(int) val
		RETURNS:
			NONE
		'''
		# Value is the only attribute required to be initialized upon object creation.
		self.val = val

		# The rest of the attributes are initialized with place-holder or empty values of the correct data-type.
		self.score = 0
		self.previous = None # Will be a Node Object of the Node that preceeds it.(Maybe this should just be an int value of the prev).
		self.outedges = []
		self.inedges = []

	def addNextEdge(self,edge):
		''' 
		Add a given edge to the outedges list of this Object instance.
		INPUT:
			(str) 'nextNode:edgeWeight'
		RETURNS:
			NONE
		'''
		if edge in self.outedges:
			pass
		else:
			self.outedges.append(edge)

	def addPrevEdge(self,edge):
		''' 
		Add a given edge to the inedges list of this Object instance.
		INPUT:
			(str) 'prevNode:edgeWeight'
		RETURNS:
			NONE
		'''
		if edge in self.inedges:
			pass
		else:
			self.inedges.append(edge)

	def getStatus(self):
		'''
		Print the status of all attributes of this Node object.
		INPUT:
			NONE
		RETURNS:
			NONE
		'''
		print('val: ', self.val)
		print('score: ', self.score)
		print('previous: ', self.previous)
		# List comprehension to print all outedge/inedge values on one line -- comma-delimited.
		print('outedges: ',','.join([x for x in self.outedges]))
		print('inedges: ',','.join([x for x in self.inedges]))


############################################################################################

def main():
	'''
	Call functions to gather input and execute main program tast of finding the longest path in a DAG
	INPUT: 
		NONE
	RETURNS:
		NONE
	'''

	# Get input from stdin.
	directedGraph, source, sink = getInput()

	# Find the length and path of the longest path in the given directed acyclic graph (DAG)
	LongestPathInDAG(directedGraph,source,sink)



if __name__ == "__main__":
    main()
