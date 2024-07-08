
# Companion paper

## Content

We provide scripts that have been used for creating the figure of Section 4 of the companion paper [''Estimating the hyperunformity exponent''](https://arxiv.org).

### Content of the scripts

- ``figure_explanations.py``correspond to Figure 1, that provides examples of estimation of the hyperuniformity exponent of Ginibre point process and RSA point process.

-   ``figure_box_perturbed_lattice.py`` correspond to Figure 2, that study the distribution of the estimator of the hyperuniformity exponent when considering perturbed lattice. 

- `compute_wavelet_transforms(Phi, J,  i_min, i_max)`: compute the curve $`j \in J \mapsto \mathcal{C}(j)`$ used for estimating the hyperuniformity exponent (refer to Section 4.1 of the companion paper) with the point pattern Phi, using the set of scale J and centered Hermites wavelets of indexes with components between i_min and i_max -1.

- `psi(x, n)` : compute the value at x of the n-th one-dimensional Hermite wavelet.

