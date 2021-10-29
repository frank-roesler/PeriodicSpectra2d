# PeriodicSpectra2d
Matlab routine to compute the spectrum of a Schrödinger operators with periodic potential defined on **R**<sup>2</sup>. The implementation is based on a method recently developed in https://arxiv.org/abs/2104.09575. A Floquet-Bloch transform is used to replace **R**<sup>2</sup> by a bounded domain, then eigenvalues are identified as the zeros of a certain operator-valued function of the form `I-K(z)`. This operator is discretized by choosing a Fourier basis. As proved in the article, the method is guaranteed to converge.
The code takes the potential function as an input and returns a set of complex numbers, which approximates the spectrum of the Schrödinger operator.

The main building blocks of the routine are
* `Main.m` - Main script that returns approximation of the spectrum;
* `potential(x,y)` - returns value of potential function for input x,y in [0,1]<sup>2</sup>;
* `compute_potential_matrix(a,N)` - Computes operator matrix of the potential in Fourier representation. Input: vector of frequencies `k`, Fourier coefficients `a` and size of potential matrix `N`.
* `GD(z_start, stepsize, maxiter, tol)` - performs gradient descent minimization until `|det(I-K(z))|<tol` or maxiter steps have been taken.

Any comments or queries are welcome at https://frank-roesler.github.io/contact/
