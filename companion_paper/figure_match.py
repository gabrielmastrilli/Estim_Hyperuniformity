import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams.update({'font.size': 13})
from structure_factor.point_processes import (
    HomogeneousPoissonPointProcess,
    mutual_nearest_neighbor_matching,
)
from structure_factor.spatial_windows import BoxWindow
from structure_factor.utils import meshgrid_to_column_matrix
from structure_factor.point_processes import HomogeneousPoissonPointProcess

from compute_estimator import alpha_hat

### Creating points pattern of matched point process and estimating alpha with 75 tapers and 16 tapers


#Simulation of a matched point process in [-R, R]^2, with Poisson intensity paremeter intens
#The function requires the python library structure-factor. Documentation of this library: https://for-a-few-dpps-more.github.io/structure-factor/point_processes.html

def match(R, intens):
    #Observation windows
    bounds = np.array([[0, 2*R], [0, 2*R]])
    window_m = BoxWindow(bounds)

    # Perturbed grid
    stepsize = 1.0
    ranges = (np.arange(a, b, step=stepsize) for a, b in window_m.bounds)
    X = meshgrid_to_column_matrix(np.meshgrid(*ranges))
    shift = np.random.uniform(0.0, stepsize, size=window_m.dimension)
    X += shift

    #Poisson process
    rho = (1 + intens) / (stepsize ** window_m.dimension)
    ppp = HomogeneousPoissonPointProcess(rho)
    Y = ppp.generate_sample(window_m)

    #Matching
    boxsize = window_m.bounds[:, 1]
    matching = mutual_nearest_neighbor_matching(X, Y, boxsize=boxsize)

    Phi = np.zeros((len(Y[matching, 0]),2))
    for i in range(len(Phi)):
        Phi[i, 0] = Y[matching, 0][i] -R
        Phi[i, 1] = Y[matching, 1][i] -R
    return Phi

#Estimation with 75 tapers

#Parameters related to the intensity of the Poisson point process
intenss= [1, 0.5, 0.2]
#Lower bound for the scales used for estimating alpha
j_min = [0.7, 0.75, 0.8]

#Side lenght of the observation window [-r, r]^2
R = 40

#Number of  estimation of alpha
n_sim = 500

#Max/min degree of the hermites polynomials
i_min = 0
i_max = 10

for i_intens in range(len(intenss)):
    intens = intenss[i_intens]

    #Set of scales used used for estimating alpha
    J = np.linspace(j_min[i_intens], 1)

    #If we haven't computed the estimated alphas for this set of parameter
    if os.path.exists("./companion_paper/data/estim_matched/match_intens_"+str(intens)+".txt") == False:
        alpha_est = []
        T_j_mean = np.zeros(len(J))
        for i in range(n_sim):
            print(i, intens)

            #Generate a realization of a matched point process
            Phi = match(R, intens)

            #Estimate alpha for this realization
            alpha_est.append(alpha_hat.alpha_hat(Phi, J, i_min, i_max))


            print(alpha_est[-1], np.mean(alpha_est), np.std(alpha_est))
            print("\n")
        np.savetxt("./companion_paper/data/estim_matched/match_intens_"+str(intens)+".txt", alpha_est)

#Estimation with 16 tapers

#Max/min degree of the hermites polynomials
i_min = 0
i_max = 5

intens = 0.2
j_min = 0.75
R = 40

J = np.linspace(j_min, 1)
if os.path.exists("./companion_paper/data/estim_matched/match_intens_"+str(intens)+"_16.txt") == False:
    alpha_est = []
    T_j_mean = np.zeros(len(J))
    for i in range(n_sim):
        print(i, intens, i_max)
        Phi = match(R, intens)
        alpha_est.append(alpha_hat.alpha_hat(Phi, J, i_min, i_max))
        print(alpha_est[-1], np.mean(alpha_est), np.std(alpha_est))
        print("\n")
    np.savetxt("./companion_paper/data/estim_matched/match_intens_"+str(intens)+"_16.txt", alpha_est)

### Box plot of the estimated alphas

## Estimation with 75 tapers
intenss= [0.5, 1]
colors = ["tab:blue","tab:orange", "tab:green"]

x = 0.5
y = 3.5
for i_intens in range(len(intenss)):
    intens = intenss[i_intens]
    alpha_hats = np.loadtxt("./companion_paper/data/estim_matched/match_intens_"+str(intens)+".txt")
    bplot = plt.boxplot(alpha_hats, positions=[x],  patch_artist=True, widths = 0.4)
    plt.text(x-0.03*(intens==0.5), 0.7, str(1+intens))
    x+=0.5
    for patch in bplot['boxes']:
        patch.set_facecolor(colors[i_intens])
    for median in bplot['medians']:
        median.set_color('black')


plt.hlines(2, -1, x, color = 'k', linewidth = 0.2)
ax= plt.gca()
ax.set(xticklabels=[])
ax.tick_params(bottom=False)
plt.ylim(0.9, 3.5)
plt.xlim(0.25, x-0.25)
plt.show()

##Estimation with 16 tapers
intens= 0.2

x = 0.5
y = 3.5
alpha_hats = np.loadtxt("./companion_paper/data/estim_matched/match_intens_"+str(intens)+"_16.txt")
bplot = plt.boxplot(alpha_hats, positions=[x],  patch_artist=True, widths = 0.4)
plt.text(x-0.1, -1.5, "16 tapers")
x+=0.5
for patch in bplot['boxes']:
    patch.set_facecolor("tab:green")
for median in bplot['medians']:
    median.set_color('black')

alpha_hats = np.loadtxt("./companion_paper/data/estim_matched/match_intens_"+str(intens)+".txt")
bplot = plt.boxplot(alpha_hats, positions=[x],  patch_artist=True, widths = 0.4)
plt.text(x-0.1, -1.5, "75 tapers")
x+=0.5
for patch in bplot['boxes']:
    patch.set_facecolor("tab:red")
for median in bplot['medians']:
    median.set_color('black')


plt.hlines(2, -1, x, color = 'k', linewidth = 0.2)
ax= plt.gca()
ax.set(xticklabels=[])
ax.tick_params(bottom=False)
plt.ylim(-1, 5)
plt.xlim(0.25, x-0.25)
plt.show()


