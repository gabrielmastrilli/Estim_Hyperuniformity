# Estimation of the Hyperuniformity Exponent (EstimHyperuniformityAlpha)

[![Python >=3.8,<3.10](https://img.shields.io/badge/python->=3.8,<3.10-blue.svg)](https://www.python.org/downloads/release/python-371/)

> Estimating hyperuniformity exponent via multi-tapers (Hermites) spectral estimator, computing asymptotic confident intervals, sampling cloacked and alpha-stable perturbed lattices 

## Introduction

`EstimHyperuniformityAlpha` currently contains Python functions that allows for computing estimators of the hyperuniformity exponent of point processes, and associated asymptotic confidence intervals. It also contains Python functions for sampling hyperuniform perturbed lattices with prescribed hyperuniformity exponent. 

Finally, it provides the scripts that has been used for creating the figures of the associated research paper titled [''Estimating the hyperunformity exponent''](https://arxiv.org).

## Dependencies

- `EstimHyperuniformityAlpha` works with [![Python >=3.8,<3.10](https://img.shields.io/badge/python->=3.8,<3.10-blue.svg)](https://www.python.org/downloads/release/python-371/).

- ``figure_match.py`` requires the Python package `structure-factor` (https://github.com/For-a-few-DPPs-more/structure-factor) for sampling matching point processes. 

## Getting started

### Companion paper

This project serves as a companion code for the research paper titled [''Estimating the hyperunformity exponent''](https://arxiv.org).

In this paper, we develop multi-scale, multi-taper estimators the hyperuniformity exponent and analyze its asymptotic behavior : consistency, asymptotic confidence intervals  We provide insights into the influence of tapers on the bias-variance trade-off of the estimator. Finally, we investigate its performance through simulations, and we apply our method to the analysis of hyperuniformity in a real dataset of marine algae.

### Content

1. In [./compute_estimators](./compute_estimators), we provide functions to compute estimators of the hyperuniformity exponent introduced in the paper [''Estimating the hyperunformity exponent''](https://arxiv.org) and its associated asymptotic confident intervals.

- ``compute_alpha_hat.py``: compute the estimators of the hyperuniformity exponent for point processes in dimension two, using Hermites wavelets, and compute the regression curve leading to its estimator (refer to Section 4 of the companion paper). 
- ``compute_confident_interval``: compute the covariance matrices used for the confident intervals and estimate the quantile of the asymptotic distribution (refer to Section 3.4  of the companion paper).

2. In [./tutorial](./tutorial), we provide tutorial for using ``compute_alpha_hat.py`` and ``compute_confident_interval``. We also provide function for sampling the point processes of the tutorial.

- ``tutorial.py``: estimate the hyperuniformity exponent for a Poisson point process, a perturbed lattice and compute an asymptotic confident interval for associated to one estimation of the hyperuniformity exponent of a Poisson point process.

- ``generate_pp.py`` : sample (two dimensional) Poisson point processes and cloaked perturbed lattice by stable distributions.

3. In [./companion_paper](./companion_paper), we provide the scripts and the data that has been used for creating the figures of the companion paper.

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






