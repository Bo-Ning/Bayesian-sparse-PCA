# VBsparsePCA: The Variational Bayesian Method for Sparse PCA

This package contains functions for a variational Bayesian approach
for sparse PCA. The algorithm is the PX-CAVI algorithm proposed by Ning (2020) if assuming
the jointly row-sparsity assumption for the loadings matrix and the batch PX-CAVI algorithm
if otherwise. The outputs of the main function, VBsparsePCA, include the mean and covariance
of the loadings matrix, the score functions, the variable selection results, and the estimated
variance of the random noise. 

See *VBsparsePCA.pdf* file for details.
