# Tutorial

## Content
The python file ``tutorial.py`` contains three examples:

- the first one estimates the hyperuniformity exponent of a Poisson point process in dimension two and plot the associated regression curve,
- the second one focuses on a cloacked [[1]](#1) and perturbed lattice by stable distribution in dimension two for which a prescribed hyperuniformity exponent (denoted alpha in the code) can be set,
- the third one computes the asymptotic confidence interval associated to the estimation of the hyperuniformity exponent for a Poisson point process.
  
The choice of the parameters (scales, number of tapers) is described in Section 4.1. of the companion paper.

The python file ``generate_pp.py`` contains auxiliary functions to generate the point patterns used in ``tutorial.py``.

## References
<a id="1">[1]</a> 
Klatt, M. A., Kim, J., & Torquato, S. (2020). 
Cloaking the underlying long-range order of randomly perturbed lattices. 
Physical Review E, 101(3), 032118

