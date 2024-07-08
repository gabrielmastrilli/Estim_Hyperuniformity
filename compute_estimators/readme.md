# Compute estimator

## Content

We provide functions for computing estimators of the hyperuniformity exponent, visualisation of the scaling curve, and asymptotic confidence intervals.

Theoritical garentee and implementation details are provided in the companion paper 

### Content of ``compute_alpha_hat.py``

- `compute_alpha(Phi, J,  i_min, i_max)`: compute the estimator of the hyperuniformity exponent with the point pattern Phi, using the set of scale J and centered Hermites wavelets of indexes with components between i_min and i_max -1 (refer to Section 4.1 of the companion paper). 

- 
- one for estimating the hyperuniformity exponent of a cloacked and perturbed lattice by stable distribution in dimension two
- one for computing the asymptotic confidence interval associated to one estimation of the hyperuniformity exponent for a Poisson point process.

They use the multi-tapers, multi-scale estimators (with Hermite wavelets tapers) introduced in [''Estimating the hyperunformity exponent''](https://arxiv.org). 

## Companion paper 

For theoretical properties of the estimators, we refer to Sections 3 of the aforementioned paper, and concerning the choices of its parameters, we refer to Section 4.1.
