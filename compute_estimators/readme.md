# Compute estimator

## Content

We provide functions for computing estimators of the hyperuniformity exponent, visualisation of the scaling curve, and asymptotic confidence intervals.

Theoritical garentee and implementation details are provided in the companion paper [''Estimating the hyperunformity exponent''](https://arxiv.org).

### Content of ``compute_alpha_hat.py``

- `compute_alpha(Phi, J,  i_min, i_max)`: compute the estimator of the hyperuniformity exponent with the point pattern Phi, using the set of scale J and centered Hermites wavelets of indexes with components between i_min and i_max -1 (refer to Section 4.1 of the companion paper). 

- `compute_wavelet_transforms(Phi, J,  i_min, i_max)`: compute the curve $`j \in J \mapsto \mathcal{C}(j)`$ used for estimating the hyperuniformity exponent (refer to Section 4.1 of the companion paper) with the point pattern Phi, using the set of scale J and centered Hermites wavelets of indexes with components between i_min and i_max -1.

- `psi(x, n)` : compute the value at x of the n-th one-dimensional Hermite wavelet.
