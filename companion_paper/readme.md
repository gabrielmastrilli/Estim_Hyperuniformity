
# Companion paper

## Content

We provide scripts that have been used for creating the figure of Section 4 of the companion paper [''Estimating the hyperunformity exponent''](https://arxiv.org).

- ``figure_explanations.py``correspond to Figure 1, that provides examples of estimation of the hyperuniformity exponent of Ginibre point process and RSA point process. It use pre-computed realisation of these two point process, that have been stored in the [./data/point_patterns](./data/point_patterns) folder.

-   ``figure_box_perturbed_lattice.py`` correspond to Figure 2, that studies the distribution of the estimator of the hyperuniformity exponent when considering perturbed lattice. First part of the script generate the data that have been stored in the [./data](./data) folder. Second part create the Figure 2.

-   ``table_coverage.py`` correspond to Table 1, that computes the coverage rate for perturbed lattices. Data concerning the computed covariance matrices are not on the [./data/cov](./data/cov) due to their size, but  are available uppon request. 

  - ``figure_match.py`` correspond to Figure 3, that study the distribution of the estimator of the hyperuniformity exponent when considering matching point process. First part of the script generate the data that have been stored in the [./data](./data) folder. Second part create the Figure 3. 

- ``figure_algae.py`` correspond to Figure 4, that estimate the hyperuniformity exponent for a system of marine algae studied in [''Circular swimming motility and disordered hyperuniform state in an algae system''](https://www.pnas.org/doi/full/10.1073/pnas.2100493118). Position of the points are in [./data/point_patterns](./data/point_patterns) folder.

