
# Companion paper

We provide scripts and data that have been used for creating the figures of Section 4 of the companion paper [''Estimating the hyperunformity exponent''](https://arxiv.org/abs/2407.16797).

## Content

- ``figure_explanations.py``corresponds to Figure 1, that provides examples of estimation of the hyperuniformity exponent for a Ginibre point process and for a RSA (also called Matérn III hard-core) point process. It uses pre-computed realisation of these two point processes, that have been stored in the [./data/point_patterns](./data/point_patterns) folder.

-   ``figure_box_perturbed_lattice.py`` corresponds to Figure 2, that plots the distribution of the estimator of the hyperuniformity exponent when considering cloaked perturbed lattices [[1]](#1). First part of the script generates the data. These data have been pre-computed and stored in the [./data/estim_perturbed_lattices](./data/estim_perturbed_lattices) folder so that the first part of the script is not necessary to run. The second part creates Figure 2 based on these data.

-   ``table_coverage.py`` provides the script that has been used to create Table 1, that computes the coverage rates for perturbed lattices. Due to the computation of the covariance matrices it takes several hours to tun. 

  - ``figure_match.py`` corresponds to Figure 3, that shows the distribution of the estimator of the hyperuniformity exponent when considering matched point processes [[2]](#2). First part of the script generates the data. These data have been pre-computed and stored in the [./data/estim_matched](./data/estim_matched) folder so that the first part of the script is not necessary to run. Second part creates Figure 3. 

- ``figure_algae.py`` corresponds to Figure 4, that estimates the hyperuniformity exponent for a system of marine algae studied in [[3]](#3). Positions of the points are in [./data/point_patterns/algae](./data/point_patterns/algae) folder.

- ## References

<a id="1">[1]</a> 
Klatt, M. A., Kim, J., & Torquato, S. (2020). 
Cloaking the underlying long-range order of randomly perturbed lattices. 
Physical Review E, 101(3), 032118

<a id="2">[2]</a> 
Andreas Klatt, M., Last, G., & Yogeshwaran, D. (2020).
Hyperuniform and rigid stable matchings. 
Random Structures & Algorithms, 57(2), 439-473.
  
<a id="3">[3]</a> 
Huang, Mingji and Hu, Wensi and Yang, Siyuan and Liu, Quan-Xing and Zhang, HP (2021). 
Circular swimming motility and disordered hyperuniform state in an algae system. 
Proceedings of the National Academy of Sciences, 118(18), e2100493118.

