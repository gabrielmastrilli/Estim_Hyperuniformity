import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})
import os

from alpha_hat  import alpha_hat
from tutorial import generate_pp

### Creating points pattern of perturbed lattices and estimating alpha

# Parameters of  the perturbed lattices
# alpha is the hyperuniformity exponent
# sig is a parameter related to the strengh of the perturbation
alphas = [0.5, 1, 1.5]
#sigs = [0.15, 0.25, 0.35]
sigs = [0.15, 0.3, 0.35]

# Side length of the observation windows [-R, R]^2
Rs = [40, 35, 30, 25, 20, 15]

#Lower bound for the scales used for estimating alpha
j_min = [0.45, 0.5, 0.5, 0.55, 0.6, 0.65]

# Number of estimations of alpha
n_sim = 500

#Max/min degree of the hermites polynomials
i_min = 0
i_max = 10


for i_R in range(len(Rs)):
    R = Rs[i_R]
    for i_alpha in range(len(alphas)):
        alpha = alphas[i_alpha]

        #Set of scales used used for estimating alpha
        J = np.linspace(j_min[i_R], 1)

        #If we haven't computed the estimated alphas for this set of parameter
        if os.path.exists("data/alpha_"+str(alpha)+"_R_"+str(R)+".txt") == False:
            alpha_hat = []
            for i in range(n_sim):
                print(i, alpha, R)

                #Generate a realization of a perturbed lattice
                Phi = generate_pp.perturbed_latt(R, alpha, sigs[i_alpha])

                #Estimate alpha for this realization
                alpha_hat.append(alpha_hat.alpha_hat(Phi, J, i_min, i_max))

                print(alpha_hat[-1], np.mean(alpha_hat), np.std(alpha_hat))
                print("\n")

            np.savetxt("data/alpha_"+str(alpha)+"_R_"+str(R)+".txt", np.array(alpha_hat))


### Box plot of the estimated alphas

alphas = [0.5, 1, 1.5]
Rs = [15, 20, 25, 30, 35, 40]

colors = ["tab:blue","tab:orange", "tab:green"]
x = 1
is_legend = False
ax= plt.gca()

for R in Rs:
    alpha_hats = []
    for alpha in alphas:
        alpha_hats.append(np.loadtxt("data/alpha_"+str(alpha)+"_R_"+str(R)+".txt"))
        print(alpha, R, np.quantile(alpha_hats, 0.05), np.quantile(alpha_hats, 0.95))

    bplot = plt.boxplot(alpha_hats, positions=[x, x+0.5, x+1],  patch_artist=True, widths = 0.4)
    plt.text(x+0.1, -0.4, str(4*R**2))

    x+=2.5

    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
    for median in bplot['medians']:
        median.set_color('black')
    if is_legend ==False:
        ax.legend([bplot["boxes"][0], bplot["boxes"][1], bplot["boxes"][2]], ['0.5', '1', '1.5'], loc='upper left', title = "$\\alpha$")
        is_legend = True

plt.hlines(alphas, -1, x, color = 'k', linewidth = 0.2)
ax.set(xticklabels=[])
ax.tick_params(bottom=False)
plt.ylim(-0.2, 2.3)
plt.show()
