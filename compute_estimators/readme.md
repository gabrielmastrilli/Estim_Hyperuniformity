# compute_estimator

We provide functions for computing estimators of the hyperuniformity exponent, for plotting the scaling curve, and computing the asymptotic confidence interval.

## Content of ``alpha_hat.py``

- `alpha_hat(Phi, J,  i_min, i_max)`: computes the estimator of the hyperuniformity exponent with the point pattern Phi, using the set of scales J and centered Hermite wavelets of indexes with components between i_min and i_max -1 (refer to Section 4.1 of the companion paper). 

- `wavelet_transform(Phi, J,  i_min, i_max)`: compute the curve $`j \in J \mapsto \mathcal{C}(j)`$ used for estimating the hyperuniformity exponent (refer to Section 4.1 of the companion paper) with the point pattern Phi, using the set of scales J and centered Hermite wavelets of indexes with components between i_min and i_max -1.

- `psi(x, n)` : compute the value at x of the n-th one-dimensional Hermite wavelet.

## Content of ``confident_interval.py``

- `cov_matrix(r, alpha, J, i_min, i_max)`: computes the covariance matrix of Proposition 3.14 of the companion paper, for an observation window $`[-r, r]^2`$, using the set of scales J and centered Hermite wavelets of indexes with components between i_min and i_max -1 (refer to Section 4.1 of the companion paper). For details about computation of this matrix, we refer to Section 5.7.

- `confidence_interval_cov_R(q_1, q_2, r, alpha_hat, J, i_min, i_max)`: computes the confidence interval of coverage level q_2 - q_1 (by default q_1 =0.025 and q_2 = 0.975) of Proposition 3.14 of the companion paper, for an observation window  $`[-r, r]^2`$, using the set of scale J and centered Hermite wavelets of indexes with components between i_min and i_max -1 (refer to Section 4.1 of the companion paper).

- `confidence_interval_cov_R(q_1, q_2, r, alpha, J, i_min, i_max)` : estimate the quantile of the hyperuniformity exponent with the asymptotic covariance matrix of Theorem 3.9 of the companion paper. Note that this function has been implemented only for comparison with `compute_confident_interval` and should not be used as it could lead to non accurate confidence intervals for non really large point patterns (refer to the discussion of the end of Section 4.2 of the companion paper)
