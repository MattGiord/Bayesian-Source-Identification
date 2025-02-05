# Bayesian-Source-Identification

MATLAB code for nonparametric Bayesian inference on the source function in an elliptic PDE.

Author: Matteo Giordano, https://matteogiordano.weebly.com.

This repository is associated with the article "A Bayesian approach with Gaussian priors to the inverse problem 
of source identification in elliptic PDEs" by Matteo Giordano. The abstract of the paper is:

"We consider the statistical linear inverse problem of making inference on an  unknown source function in an elliptic 
partial differential equation from noisy observations of its solution. We employ nonparametric Bayesian procedures based on Gaussian priors, 
leading to convenient conjugate formulae for posterior inference. We review recent results providing theoretical guarantees on the 
quality of the resulting posterior-based estimation and uncertainty quantification, and we discuss the application of the theory to the 
important classes of Gaussian series priors defined on the Dirichlet-Laplacian eigenbasis and Matérn process priors. 
We provide an implementation of posterior inference for both classes of priors, and investigate its performance in a numerical simulation study."

This repository contains the MATLAB code to replicate the numerical simulation study presented in Section 3 of the article. 
It contains the following:

1. GenerateObservations.m, code to generate the observations (discrete point evaluations of the elliptic PDE solution corrupted by
   additive Gaussian measurement errors).
2. BayesSeries.m, code to implement posterior inference with Gaussian series priors defined on the Dirichlet-Laplacian eigenbasis.
   It requires the output of GenerateObservations.m.
3. BayesSeriesOneDim.m, code to implement and investigate the performance of semiparametric inference for one dimensional functionals
   of the source function. It requires the output of GenerateObservations.m.
4. BayesMatern.m, code to implement posterior inference with Matérn process priors.
   It requires the output of GenerateObservations.m and the auxiliary function in K_mat.m.
5. K_mat.m, auxiliary code for the Matérn covariance kernel, required by BayesMatern.m.

For questions or for reporting bugs, please e-mail Matteo Giordano (matteo.giordano@unito.it).

Please cite the following article if you use this repository in your research: Giordano, M (2024). A Bayesian approach with Gaussian priors 
to the inverse problem of source identification in elliptic PDEs.
