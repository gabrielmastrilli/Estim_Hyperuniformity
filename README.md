# Estimation of the Hyperuniformity Exponent (EstimHyperuniformityAlpha)

[![Python >=3.8,<3.10](https://img.shields.io/badge/python->=3.8,<3.10-blue.svg)](https://www.python.org/downloads/release/python-371/)

> Estimating hyperuniformity exponent via multi-tapers (Hermites) spectral estimator, computing asymptotic confident intervals, sampling cloacked and alpha-stable perturbed lattices 

## Introduction

`EstimHyperuniformityAlpha` currently contains Python functions that allows for computing an estimator of the hyperuniformity exponent of a point process, and an associated asymptotic confidence interval. It also contains code for sampling hyperuniform perturbed lattices with prescribed hyperuniformity exponent. 

Finally, it provides scripts that has been used in the associated research paper titled [''Estimating the hyperunformity exponent''](https://arxiv.org).

## Dependencies

- `EstimHyperuniformityAlpha` works with [![Python >=3.8,<3.10](https://img.shields.io/badge/python->=3.8,<3.10-blue.svg)](https://www.python.org/downloads/release/python-371/).

- ``figure_macth.py`` require the Python package `structure-factor` (https://github.com/For-a-few-DPPs-more/structure-factor) for sampling matching point processes. 

## Getting started

### Companion paper

This project serves as a companion code for the research paper titled [''Estimating the hyperunformity exponent''](https://arxiv.org).

In this paper, we develop a multi-scale, multi-taper estimator the hyperuniformity exponent and analyze its asymptotic behavior : consistency, asymptotic confidence intervals  We provide insights into the influence of tapers on the bias-variance trade-off of the estimator. Finally, we investigate its performance through simulations (whose scripts are in the [./companion_paper](./companion_paper) folder) , and we apply our method to the analysis of hyperuniformity in a real dataset of marine algae.





