# python version = 3.11.4
from scipy import optimize # 1.11.2
import patsy # 0.5.3
import numpy as np # 1.25.2

def stabilize(expression_matrix):
	''' Use Anscombes approximation to variance stabilize Negative Binomial data

	See https://f1000research.com/posters/4-1041 for motivation.

	Assumes columns are samples, and rows are genes
	'''
	phi_hat, _ = optimize.curve_fit(lambda mu, phi: mu + phi * mu ** 2, expression_matrix.mean(1), expression_matrix.var(1))

	return np.log(expression_matrix + 1. / (2 * phi_hat[0]))

def regress_out(sample_info, expression_matrix, covariate_formula, design_formula='1', rcond=-1):
	''' Implementation of limma's removeBatchEffect function
	'''
	# Ensure intercept is not part of covariates
	covariate_formula += ' - 1'

	covariate_matrix = patsy.dmatrix(covariate_formula, sample_info)
	design_matrix = patsy.dmatrix(design_formula, sample_info)

	design_batch = np.hstack((design_matrix, covariate_matrix))

	coefficients, res, rank, s = np.linalg.lstsq(design_batch, expression_matrix.T, rcond=rcond)
	beta = coefficients[design_matrix.shape[1]:]
	regressed = expression_matrix - covariate_matrix.dot(beta).T

	return regressed