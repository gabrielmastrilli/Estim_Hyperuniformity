
# Companion paper

We provide scripts and data that have been used for creating the figures of Section 4 of the companion paper [''Estimating the hyperunformity exponent''](https://arxiv.org).

## Content

- ``figure_explanations.py``corresponds to Figure 1, that provides examples of estimation of the hyperuniformity exponent for a Ginibre point process and for a RSA point process. It uses pre-computed realisation of these two point processes, that have been stored in the [./data/point_patterns](./data/point_patterns) folder.

-   ``figure_box_perturbed_lattice.py`` corresponds to Figure 2, that plots the distribution of the estimator of the hyperuniformity exponent when considering perturbed lattices. First part of the script generates the data that have been stored in the [./data](./data) folder. The second part creates Figure 2.

-   ``table_coverage.py`` corresponds to Table 1, that computes the coverage rates for perturbed lattices. Data concerning the computed covariance matrices are not on the [./data/cov](./data/cov) due to their size, but  are available uppon request. 

  - ``figure_match.py`` corresponds to Figure 3, that studies the distribution of the estimator of the hyperuniformity exponent when considering matching point processes. First part of the script generates the data that have been stored in the [./data](./data) folder. Second part creates the Figure 3. 

- ``figure_algae.py`` corresponds to Figure 4, that estimates the hyperuniformity exponent for a system of marine algae studied in [''Circular swimming motility and disordered hyperuniform state in an algae system''](https://www.pnas.org/doi/full/10.1073/pnas.2100493118). Positions of the points are in [./data/point_patterns](./data/point_patterns) folder.

