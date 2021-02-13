# VBsparsePCA: The Variational Bayesian Method for Sparse PCA

This package contains functions for a variational Bayesian approach
for sparse PCA. The algorithm is the PX-CAVI algorithm proposed by Ning (2021) [(arXiv:2102.00305)](https://arxiv.org/abs/2102.00305) if assuming
the jointly row-sparsity assumption for the loadings matrix and the batch PX-CAVI algorithm
if otherwise. The outputs of the main function, VBsparsePCA, include the mean and covariance
of the loadings matrix, the score functions, the variable selection results, and the estimated
variance of the random noise. 

See *VBsparsePCA.pdf* for details.
