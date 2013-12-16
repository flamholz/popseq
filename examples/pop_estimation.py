#!/usr/bin/python

import numpy as np
import pylab

from scikits.statsmodels.discrete.discrete_model import Poisson


def GenCaptureHistories(true_pop_size,
 					    capture_rounds,
 					    capture_prob):
	"""Generates a record of captures

	TODO: the returned matrix may be huge because it contains
	individuals never captured. Can save memory for very large populations.

	Args:
		true_pop_size: the size of the population.
		capture_rounds: the number of rounds of capture.
		capture_prob: the probability that any individual is captured 
		  in a given round.

	Returns:
		A binary matrix true_pop_size x capture_rounds. 1 at
		index i,j implies individual i was captured in round j.
	"""
 	capture_prs = np.random.random((true_pop_size, capture_rounds))
	capture_histories = np.zeros((true_pop_size, capture_rounds))
	capture_histories[np.where(capture_prs < capture_prob)] = 1
	return capture_histories
 

def CalcHistoryFreqs(capture_histories, capture_rounds):
	"""Calculates the observed frequencies of all capture histories.
	
	Args:
		capture_histories: the binary matrix of capture histories
			for all individuals.
		capture_rounds: the number of rounds of capture.
	
	Returns:
		A 1-d array of length 2^capture_rounds - 1 with entry i being
		the number of times history i was observed. History i is the 
		history whose binary representation = i + 1.
	"""
	n_possible_histories = 2**capture_rounds - 1
	obs_history_freq = np.zeros(n_possible_histories)
	for row in capture_histories:
	  int_form = int('0b%s' % ''.join(['%d' % i for i in row]), 2)
	  if int_form == 0:
		continue

	  # We have observed this non-zero history once more
	  obs_history_freq[int_form - 1] += 1
	return obs_history_freq


def InitIndependentVariables(capture_rounds):
	"""Initializes ind. variables for poisson regression.
	
	Args:
		capture_rounds: the number of rounds of capture.
	
	Returns:
		A (2^capture_rounds - 1) x 2 matrix.
	"""
	n_possible_histories = 2**capture_rounds - 1
	X = np.ones((n_possible_histories, 2))
	for i_plus in xrange(1, n_possible_histories + 1):
	  index = i_plus - 1
	  X[index, 1] = bin(i_plus).count("1")
	return X
	

def EstimatePopnSize(obs_history_freq, n_caught, capture_rounds):
	"""Estimate the population size from observed frequencies.
	
	Args:
		obs_history_freq: the observed frequencies of the capture histories.
		n_caught: the total number of individuals caught at least once.
		capture_rounds: the number of captures.
	
	Returns:
		Estimated population size.
	"""
	X = InitIndependentVariables(capture_rounds)
	fitted = Poisson(obs_history_freq, X).fit(disp=False)
	intercept = fitted.params[0]
	slope = fitted.params[1]
	expected_num_uncaught = np.exp(intercept)
	return expected_num_uncaught + n_caught
	

def Main():
	# Parameters and sample generation
	popn_size = 1e4
	capture_rounds = 10
	capture_prob = 0.05
	
	capture_probs = np.power(2.0, np.arange(-20, 1, 1))
	estimated_popn_sizes = np.zeros(capture_probs.size)
	
	for i, capture_prob in enumerate(capture_probs):
		capture_histories = GenCaptureHistories(
			popn_size, capture_rounds, capture_prob)
		obs_history_freq = CalcHistoryFreqs(
			capture_histories, capture_rounds)

		n_caught = np.sum(np.sum(capture_histories, 1) > 0)
		estimated_popn_size = EstimatePopnSize(
			obs_history_freq, n_caught, capture_rounds)
		estimated_popn_sizes[i] = estimated_popn_size
	
	error = np.abs(estimated_popn_sizes - popn_size)
	print capture_probs
	print error
	
	pylab.figure()
	pylab.xlabel('Capture Probability')
	pylab.ylabel('Error in Size Estimate')
	pylab.loglog(capture_probs, error, 'r.')
	pylab.ylim((1e-3, np.max(error) * 10))
	pylab.show()
	
	
if __name__ == '__main__':
	Main()

