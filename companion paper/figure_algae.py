import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})
import csv
import os
import compute_alpha_hat

import compute_alpha_hat
import generate_pp


#Mean of the curves C associated to each frame
L_j_mean = np.zeros(50)

#Sets of scales used for visualization
J = np.linspace(0.2, 1.2)
#Min and max degree of the 1d-hermites polynomials involved in the hermites wavelets

i_min = 0
i_max = 10

nb_frame = 0
for id in range(98):
    print(id)
    if os.path.exists("data/point_patterns/algae_"+str(id)+".txt"):
        nb_frame+=1
        #Algae point pattern associated to the frame id
        Phi = np.loadtxt("data/point_patterns/algae_"+str(id)+".txt")
        #Compute and plot the curve C for the frame id
        L_j, r = compute_alpha_hat.compute_wavelet_transforms(Phi, J, i_min, i_max)
        plt.plot(J, np.log(L_j)/np.log(r), linewidth = 0.1, color = 'k')

        L_j_mean+=np.array(np.log(L_j)/np.log(r))
#Plot the mean of the curves C
plt.plot(J, L_j_mean/nb_frame, color = 'r')

#Small scale portion of the curve
j_0 = 0.2
j_1 = 0.45

#Get the indexes in J corresponding to j_0 and j_1
k_0 = -1
k_1 = -1
for i_j in range(len(J)):
    if J[i_j] >= j_0 and k_0==-1:
        k_0= i_j
    if J[i_j] >= j_1:
        k_1 = i_j
        break

#Estimate the local slope on [0.2, 0.45]
fit = np.polyfit(J[k_0:k_1], L_j_mean[k_0:k_1]/nb_frame, 1)
plt.plot(J[k_0:k_1], J[k_0:k_1]*fit[0]+fit[1]-0.1, color = 'orange', linestyle = '--')
plt.plot(J[k_0:k_1], J[k_0:k_1]*fit[0]+fit[1]+0.2, color = 'orange', linestyle = '--')
x_txt = 0.4
plt.text(0.35, 0.8, "$\widehat{\\alpha} = $"+str(2 - fit[0])[:3], color = 'orange')

#Intermediate scale portion of the curve
j_0 = 0.45
j_1 = 0.7
#Get the indexes in J corresponding to j_0 and j_1
k_0 = -1
k_1 = -1
for i_j in range(len(J)):
    if J[i_j] >= j_0 and k_0==-1:
        k_0= i_j
    if J[i_j] >= j_1:
        k_1 = i_j
        break

#Estimate the local slope on [0.45, 0.7]
fit = np.polyfit(J[k_0:k_1], L_j_mean[k_0:k_1]/nb_frame, 1)
plt.plot(J[k_0:k_1], J[k_0:k_1]*fit[0]+fit[1]-0.1, color = 'violet', linestyle = '--')
plt.plot(J[k_0:k_1], J[k_0:k_1]*fit[0]+fit[1]+0.2, color = 'violet', linestyle = '--')
x_txt = 0.4
plt.text(0.55, 0.8, "$\widehat{\\alpha} = $"+str(2 - fit[0])[:3], color = 'violet')

#Large scale portion of the curve
j_0 = 0.7
j_1 = 0.95
k_0 = -1
k_1 = -1
for i_j in range(len(J)):
    if J[i_j] >= j_0 and k_0==-1:
        k_0= i_j
    if J[i_j] >= j_1:
        k_1 = i_j
        break

#Estimate the local slope on [0.7, 0.95]
fit = np.polyfit(J[k_0:k_1], L_j_mean[k_0:k_1]/nb_frame, 1)
plt.plot(J[k_0:k_1], J[k_0:k_1]*fit[0]+fit[1]-0.1, color = 'blue', linestyle = '--')
plt.plot(J[k_0:k_1], J[k_0:k_1]*fit[0]+fit[1]+0.2, color = 'blue', linestyle = '--')
x_txt = 0.7
plt.text(0.8, 0.8, "$\widehat{\\alpha} = $"+str(2 - fit[0])[:3], color = 'blue')

plt.xlabel("j", loc='right', fontsize = 22)
plt.ylabel("$\mathcal{C}$", loc = "top", rotation=0, fontsize = 22)

ax = plt.gca()
ax.xaxis.set_label_coords(1.05, 0)
ax.yaxis.set_label_coords(0,1)


#Compute the curve C for a Poisson

#Number of realization
n_sim_poisson = 20
L_j_poisson_mean = np.zeros(len(J))

#We average the curves corresponding to independents realization.
for i in range(n_sim_poisson):
    print(i)
    L_j_poisson, r_poisson = compute_alpha_hat.compute_wavelet_transforms(generate_pp.PPP(r), J, i_min, i_max)
    L_j_poisson_mean += np.log(L_j_poisson)/np.log(r_poisson)

plt.plot(J, L_j_poisson_mean/n_sim_poisson, color = 'g', linestyle = "dotted")
J  = np.array(J)
plt.text(0.55, 1.8, "Poisson.", color = 'g')
plt.vlines(0.95, np.min(L_j_poisson_mean/n_sim_poisson), np.max(L_j_poisson_mean/n_sim_poisson), color = 'k')

plt.show()