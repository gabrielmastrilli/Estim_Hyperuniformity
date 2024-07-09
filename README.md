# Estimation of the Hyperuniformity Exponent (Estim_Hyperuniformity)

[![Python >=3.8,<3.10](https://img.shields.io/badge/python->=3.8,<3.10-blue.svg)](https://www.python.org/downloads/release/python-371/)

## Introduction

`Estim_Hyperuniformity` currently contains Python functions that allows for computing the estimator of the hyperuniformity exponent of point processes, and its associated asymptotic confidence interval. It also contains Python functions for sampling hyperuniform perturbed lattices with prescribed hyperuniformity exponent. 

Finally, it provides the scripts that have been used for creating the figures of the associated research paper titled [''Estimating the hyperunformity exponent''](https://arxiv.org).

## Dependencies

- `EstimHyperuniformityAlpha` works with [![Python >=3.8,<3.10](https://img.shields.io/badge/python->=3.8,<3.10-blue.svg)](https://www.python.org/downloads/release/python-371/).

- ``figure_match.py`` requires the Python package `structure-factor` (https://github.com/For-a-few-DPPs-more/structure-factor) for sampling matching point processes. 

## Getting started

### Companion paper

This project serves as a companion code for the research paper titled [''Estimating the hyperunformity exponent''](https://arxiv.org).

In this paper, we develop a multi-scale, multi-taper estimator of the hyperuniformity exponent and we analyze its asymptotic behavior such as consistency and asymptotic confidence intervals. We also provide insights into the influence of tapers on the bias-variance trade-off of the estimator. Finally, we investigate its performance through simulations, and we apply our method to the analysis of hyperuniformity in a real dataset of marine algae.

### Content

- In [./compute_estimators](./compute_estimators), we provide functions to compute the estimator of the hyperuniformity exponent, along with its associated asymptotic confidence intervals, given a point pattern.

- In [./tutorial](./tutorial), we provide tutorial for using the functions of the [./compute_estimators](./compute_estimators) folder. We also provide functions for sampling the point processes of the tutorial.

- In [./companion_paper](./companion_paper), we provide the scripts and the data that have been used for creating the figures of the companion paper.

## How to cite this work

  ```bash
  @preprint{?,
    arxivid = {},
    journal = {arXiv},
    title = {Estimating the hyperuniformity exponent of a point process},
    author = {Mastrilli, Gabriel and B\l{}aszczyszyn, Bartlomiej and Lavancier, Fr\'ed\'eric}, 
    year    = {2024},
  }
  ```






