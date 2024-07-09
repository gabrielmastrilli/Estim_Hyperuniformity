
# Companion paper

We provide scripts and data that have been used for creating the figures of Section 4 of the companion paper [''Estimating the hyperunformity exponent''](https://arxiv.org).

## Content

- ``figure_explanations.py``corresponds to Figure 1, that provides examples of estimation of the hyperuniformity exponent for a Ginibre point process and for a RSA point process. It uses pre-computed realisation of these two point processes, that have been stored in the [./data/point_patterns](./data/point_patterns) folder.

-   ``figure_box_perturbed_lattice.py`` corresponds to Figure 2, that plots the distribution of the estimator of the hyperuniformity exponent when considering perturbed lattices. First part of the script generates the data. These data have been pre-computed and stored in the [./data](./data) folder so that the first part of the script is not necessary to run. The second part creates Figure 2 based on these data.

-   ``table_coverage.py`` provides the script that has been used to create Table 1, that computes the coverage rates for perturbed lattices. Due to the computation of the covariance matrices it takes several hours to tun. 

  - ``figure_match.py`` corresponds to Figure 3, that shows the distribution of the estimator of the hyperuniformity exponent when considering matched point processes. First part of the script generates the data. These data have been pre-computed and stored in the [./data](./data) folder so that the first part of the script is not necessary to run. Second part creates Figure 3. 

- ``figure_algae.py`` corresponds to Figure 4, that estimates the hyperuniformity exponent for a system of marine algae studied in [''Circular swimming motility and disordered hyperuniform state in an algae system''](https://www.pnas.org/doi/full/10.1073/pnas.2100493118). Positions of the points are in [./data/point_patterns](./data/point_patterns) folder.

