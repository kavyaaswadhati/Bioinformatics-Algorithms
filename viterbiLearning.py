#!/usr/bin/env python3
#!/usr/bin/env python3

########################################################################
# File: problem18.py
#		executable: python3 problem18.py < input.txt > output.txt
# Purpose: Implement the Viterbi Algorithm to find the most likely hidden path given a string x,
#			its alphabet, the possible states, a transition matrix, and an emission matrix.
#		stderr: None.
#   	stdout: A path that maximizes the unconditional probability Pr(x,pi) over all possible paths pi.
#          
# Author: Kavya Aswadhati
# Group: None
#        
########################################################################


import sys
import numpy as np
import math



class PathFromViterbi(object):
	"""
	Given a string x, the alphabet from which it was constructed, a hidden path pi, the states,
	and the transition matrix Transition of an HMM, find the most likely hidden path.
	ATTRIBUTES:
		- self.viterbi: List of dicts, saves state information at each position
		- generateViterbi(): Creates the viterbi data structure
		- setScores(): Populate score values for each state position.
	"""

	def __init__(self, x, transitionMatrix, emissionMatrix):
		'''
		Initialize class-wide objects. Call functions that will execute program function.
		INPUT:
			- x: String of emissions.
			- transitionMatrix: Dict of dicts, transition probabilities to a current state given a prev.
			- emissionMatrix: Dict of dicts, emission probabilities of a letter given a state.
		'''
		# Save class-wide variables.
		self.x = x
		self.transitionMatrix = transitionMatrix
		self.emissionMatrix = emissionMatrix

		# Create the viterbi graph data structure.
		self.viterbi = self.generateViterbi()

		# Set scores for each position in the emission string
		bestLast, highestProb = self.setScores()

		# Backtrack through the scores to save the best hidden path.
		statePath = self.backtrackViterbi(bestLast)
		# Output state path in rosalind format.
		print(''.join(x for x in statePath))
		

	def generateViterbi(self):
		'''
		Create the data structure that will store the viterbi graph
		RETURNS:
			viterbi: A list of dictionaries. Each list index is a position along string x
					 and contains a dictionary of Node objects, one for each state, keyed by state.
		'''
		numStates = len(self.transitionMatrix)
		viterbi = []

		# Viterbi will need to be of length string x.
		for i in range(0,len(self.x)):
			# Initialize empty dict to populate and add to viterbi. One per pos.
			nodesAtPos = {}
			# Each dict will contain all states.
			for state in self.transitionMatrix:
				nodesAtPos[state] = Node(state)
			viterbi.append(nodesAtPos)

		# Return result.
		return viterbi

	def setScores(self):
		'''
		Set scores for each state at each position.
		INPUT:
			NONE
		RETURNS:
			State, probability of the state with highest probability at the last
			index of the viterbi graph.
		'''

		# All states at each position need to be populates with scores.
		for i, position in enumerate(self.viterbi):
			# Go through each state at this position.
			for node in position:
				# Emission is only state and position dependent.
				pEmission = self.emissionMatrix[node][self.x[i]]

				# All positions but the start have the same considerations.
				if i > 0:
					# Temporary objects to save the highest scoring state.
					highScore = -math.inf
					save = ''
					# Calculate the potential score of the current position for each possible predecessor.
					for prev in self.viterbi[i-1]:
						# Gather the score of the current predecessor in consideration.
						prevScore = self.viterbi[i-1][prev].score
						# Gather the probability of the transitioning from the previous state to the current one.
						pTransition = self.transitionMatrix[node][prev]
						# Calculate the weight in log form so that scores can be sums rather than products.
						weight = math.log10(pTransition * pEmission)

						# Check if this current score is the highest so far.
						tempScore = prevScore + weight
						if tempScore > highScore:
							# If it is the highest, save the score and the state.
							save = prev
							highScore = tempScore
					# The score of the current state is the highest of the previous ones.
					position[node].score = highScore

				# The first set of states have unique score calculations.
				elif i == 0:
					# Each state has an equal probability of starting.
					pTransition = np.float128(1/(len(self.transitionMatrix)))
					weight = math.log10(pTransition * pEmission)
					position[node].score = weight

		# Return the highest scoring state at the last position and its score.
		return self.getMax(self.viterbi[-1])



	def backtrackViterbi(self, lastState):
		'''
		Traverse the viterbi graph backwards to gather the state path based on scores.
		INPUT:
			lastState: The highest scoring node at the last position in the viterbi graph.
		RETURNS:
			statePath: The resultant string of states in forward order.
		'''

		# List in which to save the states.
		statePath = []
		# The last state is the first state as we are traveling backwards.
		statePath.append(lastState)

		# Pointer to current state under consideration.
		currentState = lastState

		# Traverse the viterbi graph in reverse order.
		for i, nodes in reversed(list(enumerate(self.viterbi))):

			# The state is already known for the last state, start here.
			if i == len(self.viterbi)-1:
				currentState = lastState

			# Gather the score of the current state at the current position.
			currentScore = self.viterbi[i][currentState].score

			# Emission is state and position dependent only.
			pEmission = self.emissionMatrix[currentState][self.x[i]]

			# Collect and check all the states at the preceeding position.
			previousNodes = self.viterbi[i-1]
			for prevNode in previousNodes:
				# Probability of transitioning from this previous state to the current one.
				pTransition = self.transitionMatrix[currentState][prevNode]
				# Edge weight connecting this previous state to the current one.
				weight = math.log10(pEmission * pTransition)
				prevScore = previousNodes[prevNode].score

				# If the previous score is the current score - their edge weight, it is the
				# state in the state path.
				if currentScore == (prevScore + weight):
					statePath.append(prevNode)
					currentState = prevNode
		# Since the state path was generated in reverse, flip it and return it.
		return statePath[::-1]
			


	def getMax(self,nodes):
		'''
		Find the state with the highest score in a dict of states.
		INPUT:
			(dict) nodes: Dictionary of node objects for each state.
		RETURNS:
			state: The highest scoring state
			maxScore: The score of that state.
		'''
		# Temporary variables.
		maxScore = -math.inf
		state = ''
		# Check each state in the dict.
		for node in nodes:
			# Check if it has the highest score so far.
			if nodes[node].score > maxScore:
				maxScore = nodes[node].score
				state = node
		# Return the results.
		return state, maxScore

##################################################################
##################################################################

class Node(object):
	'''
	Create a data structure that will be used as the node units of a Viterbi Graph.
	INPUT:
		(str) value
	RETURNS:
		NONE
	ATTRIBUTES:
		- (str) value
		- (int) score

	'''
	def __init__(self,state):
		'''
		Initialize Node Object attributes.
		INPUT:
			(int) val
		RETURNS:
			NONE
		'''
		# State is the only attribute required to be initialized upon object creation.
		self.state = state

		# The rest of the attributes are initialized with place-holder or empty values of the correct data-type.
		self.score = 0


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

##################################################################
##################################################################

def getInput():
	'''
	Get input from standard in (rosalind file) and format it.
	INPUT:
		Stdin(file): 
	RETURNS:
		- x: String of emissions.
		- transitionMatrix: Dict of dicts, transition probabilities to a current state given a prev.
		- emissionMatrix: Dict of dicts, emission probabilities of a letter given a state.
	'''

	# Get input from standard in.
	
	# Temporary objects.
	alphabet = []
	states = []
	tempTransitions = []
	tempEmissions = []

	# Objects to populate and return.
	x = ''
	transitionMatrix = {}
	emissionMatrix = {}

	# Messy parsing for rosalind's file format.
	# Keep in consideration that number of states and number of 
	# letters are flexible inputs.
	inp = sys.stdin
	for i,line in enumerate(inp):
		if i == 0:
			x = line.strip()
		elif i == 2:
			alphabet = line.strip().split()
		elif i == 4:
			states = line.strip().split()
		elif 6 < i < (7+len(states)):
			transition = line.strip().split()
			tempTransitions.append(transition[1:])
		elif (8+len(states)) < i :
			emission = line.strip().split()
			tempEmissions.append(emission)

	# Probability of an emission given a state
	# ie emissionMatrix[currentState][letter]
	for emission in tempEmissions:
		state = emission[0]
		emissionMatrix[state] = {}
		emission.remove(state)
		for i, e in enumerate(emission):
			emissionMatrix[state][alphabet[i]] = np.float128(e)

	# Transition to current state from previous state
	# ie transitionMatrix[currentState][previousState]
	for i,state in enumerate(states):
		transitionMatrix[state] = {}
		for j,temp in enumerate(tempTransitions):
			transitionMatrix[state][states[j]] = np.float128(temp[i])

	# Return pertinant inputs.
	return x, transitionMatrix, emissionMatrix


##################################################################
##################################################################


def main():
	'''
	Call functions to gather input and execute main program task.
	INPUT: 
		NONE
	RETURNS:
		NONE
	'''
	# Get input from stdin.
	x, transitionMatrix, emissionMatrix = getInput()

	# Implement the Viterbi algorithm to find the most probable path of an HMM, string.
	PathFromViterbi(x, transitionMatrix, emissionMatrix)	



if __name__ == "__main__":
    main()
