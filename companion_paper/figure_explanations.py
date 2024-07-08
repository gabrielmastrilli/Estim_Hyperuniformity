import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})

import compute_alpha_hat
import generate_pp

### Explication RSA

Phi = np.loadtxt("data/point_patterns/matern.txt")

#Plots the point pattern
#Rescale the point pattern to intensity 1, in order to compare the points configurations of RSA and Ginibre
lambda_hat = np.sqrt(len(Phi)/(4*100**2))
plt.scatter(Phi[:,0]*lambda_hat, Phi[:, 1]*lambda_hat, s = 0.25)
ax = plt.gca()

#Print only the point inside [-20, 20]^2
plt.xlim(-20, 20)
plt.ylim(-20, 20)
#plt.axis('off')
ax = plt.gca()
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])

plt.show()

#Plot the curve C

#The set of scales for which the curve C will be computed
J = np.linspace(0.2, 1.4)

#Min and max degree of the 1d-hermites polynomials involved in the hermites wavelets
i_min = 0
i_max = 10


L_j, r = compute_alpha_hat.compute_wavelet_transforms(Phi, J, i_min, i_max)

#Small scale portion of the curve
j_0 = 0
j_1 = 0.33

#Get the indexes in J corresponding to j_0 and j_1
k_0 = -1
k_1 = -1
for i in range(len(J)):
    if J[i] >= j_0 and k_0 == -1:
        k_0 = i
    if J[i] >= j_1:
        k_1 = i
        break
plt.plot(J[k_0:k_1], np.log(L_j[k_0:k_1])/np.log(r), color = 'b', linewidth = 2, label = "Small scale.")

#Intermediate scale portion of the curve
j_0 = 0.3
j_1 = 0.5

#Get the indexes in J corresponding to j_0 and j_1
k_0 = -1
k_1 = -1
for i in range(len(J)):
    if J[i] >= j_0 and k_0 == -1:
        k_0 = i
    if J[i] >= j_1:
        k_1 = i
        break
plt.plot(J[k_0:k_1], np.log(L_j[k_0:k_1])/np.log(r), color = 'y', linewidth = 2, label = "Intermediate scale.")

#Large scale portion of the curve
j_0 = 0.46
j_1 = 1.02

#Get the indexes in J corresponding to j_0 and j_1
k_0 = -1
k_1 = -1
for i in range(len(J)):
    if J[i] >= j_0 and k_0 == -1:
        k_0 = i
    if J[i] >= j_1:
        k_1 = i
        break
plt.plot(J[k_0:k_1], np.log(L_j[k_0:k_1])/np.log(r), color = 'r', linewidth = 2, label = "Large scale.")

#Estimate alpha with the log-linear regression
fit = np.polyfit(J[k_0:k_1], np.log(L_j[k_0:k_1])/np.log(r), 1)
plt.plot(np.array(J[k_0:k_1]), np.array(J[k_0:k_1])*fit[0]+ fit[1], color = 'r', linestyle = "--", linewidth = 1)
plt.text(0.68, 0.95, "$\widehat{\\alpha} = "+str(2 -fit[0])[:6]+"$",color = 'red')
print(2 - np.polyfit(J[k_0:k_1], np.log(L_j[k_0:k_1])/np.log(r), 1)[0])

#Border effects portion of the curve
j_0 = 1
j_1 = 1.4

#Get the indexes in J corresponding to j_0 and j_1
k_0 = -1
k_1 = -1
for i in range(len(J)):
    if J[i] >= j_0 and k_0 == -1:
        k_0 = i
    if J[i] >= j_1:
        k_1 = i
        break
plt.plot(J[k_0:k_1], np.log(L_j[k_0:k_1])/np.log(r), color = 'c', linewidth = 2, label = "Border effects.")

ax = plt.gca()
plt.xlabel("j", fontsize = 22)
plt.ylabel("$\mathcal{C}$", rotation=0, fontsize = 22)
plt.xlim(0.15, 1.25)
plt.ylim(0.6, 2.65)
ax.xaxis.set_label_coords(1.05, -0.025)
ax.yaxis.set_label_coords(-0.025,1)

#Compute the curve C for a Poisson

#Number of realization
n_sim_poisson = 5
L_j_poisson_mean = np.zeros(len(J))

#We average the curves corresponding to independents realization.
for i in range(n_sim_poisson):
    print(i)
    L_j_poisson, r_poisson = compute_alpha_hat.compute_wavelet_transforms(generate_pp.PPP(r), J, i_min, i_max)
    L_j_poisson_mean += np.log(L_j_poisson)/np.log(r)

plt.plot(J, L_j_poisson_mean/n_sim_poisson, color = 'g', linestyle = "dotted")
J  = np.array(J)
plt.text(0.4, 1.65, "Poisson.", color = 'g')
plt.vlines(1, np.min(L_j_poisson_mean/n_sim_poisson), np.max(L_j_poisson_mean/n_sim_poisson), color = 'k')

plt.legend()
plt.show()

### Explication Ginibre

Phi = np.loadtxt("data/point_patterns/ginibre.txt")

#Plots the point pattern

#Rescale the point pattern to intensity 1, in order to compare the points configurations of RSA and Ginibre
lambda_hat = np.sqrt(len(Phi)/(4*52**2))

plt.scatter(Phi[:,0]*lambda_hat, Phi[:, 1]*lambda_hat, s = 0.25)
ax = plt.gca()
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
#Print only the point inside [-20, 20]^2
plt.xlim(-20, 20)
plt.ylim(-20, 20)

plt.show()


#Plot the curve C
#The set of scales for which the curve C will be computed
J = np.linspace(0.2, 1.4)

#Min and max degree of the 1d-hermites polynomials involved in the hermites wavelets
i_min = 0
i_max = 10


L_j, r = compute_alpha_hat.compute_wavelet_transforms(Phi, J, i_min, i_max)

#Small scale portion of the curve
j_0 = 0
j_1 = 0.3

#Get the indexes in J corresponding to j_0 and j_1
k_0 = -1
k_1 = -1
for i in range(len(J)):
    if J[i] >= j_0 and k_0 == -1:
        k_0 = i
    if J[i] >= j_1:
        k_1 = i
        break
plt.plot(J[k_0:k_1], np.log(L_j[k_0:k_1])/np.log(r), color = 'b', linewidth = 2, label = "Small scale.")

#Intermediate scale portion of the curve
j_0 = 0.26
j_1 = 0.62

#Get the indexes in J corresponding to j_0 and j_1
k_0 = -1
k_1 = -1
for i in range(len(J)):
    if J[i] >= j_0 and k_0 == -1:
        k_0 = i
    if J[i] >= j_1:
        k_1 = i
        break
plt.plot(J[k_0:k_1], np.log(L_j[k_0:k_1])/np.log(r), color = 'y', linewidth = 2, label = "Intermediate scale.")


#Large scale portion of the curve
j_0 = 0.6
j_1 = 1.02

#Get the indexes in J corresponding to j_0 and j_1
k_0 = -1
k_1 = -1
for i in range(len(J)):
    if J[i] >= j_0 and k_0 == -1:
        k_0 = i
    if J[i] >= j_1:
        k_1 = i
        break
plt.plot(J[k_0:k_1], np.log(L_j[k_0:k_1])/np.log(r), color = 'r', linewidth = 2, label = "Large scale.")

fit = np.polyfit(J[k_0:k_1], np.log(L_j[k_0:k_1])/np.log(r), 1)
plt.plot(np.array(J[k_0:k_1]), np.array(J[k_0:k_1])*fit[0]+ fit[1], color = 'r', linestyle = "--", linewidth = 1)
plt.text(0.66, 1.05, "$\widehat{\\alpha} = "+str(2 -fit[0])[:6]+"$",color = 'red')
print(2 - np.polyfit(J[k_0:k_1], np.log(L_j[k_0:k_1])/np.log(r), 1)[0])

#Border effect portion of the curve
j_0 = 1
j_1 = 1.4

#Get the indexes in J corresponding to j_0 and j_1
k_0 = -1
k_1 = -1
for i in range(len(J)):
    if J[i] >= j_0 and k_0 == -1:
        k_0 = i
    if J[i] >= j_1:
        k_1 = i
        break
plt.plot(J[k_0:k_1], np.log(L_j[k_0:k_1])/np.log(r), color = 'c', linewidth = 2, label = "Border effects.")

ax = plt.gca()
plt.xlabel("j", fontsize = 22)
plt.ylabel("$\mathcal{C}$", rotation=0, fontsize = 22)
plt.xlim(0.15, 1.25)
plt.ylim(0.6, 2.65)
ax.xaxis.set_label_coords(1.05, -0.025)
ax.yaxis.set_label_coords(-0.025,1)

#Compute the curve C for a Poisson

#Number of realization
n_sim_poisson = 10
L_j_poisson_mean = np.zeros(len(J))

#We average the curves corresponding to independents realization.
for i in range(n_sim_poisson):
    print(i)
    L_j_poisson, r_poisson = compute_alpha_hat.compute_wavelet_transforms(generate_pp.PPP(r), J, i_min, i_max)
    L_j_poisson_mean += np.log(L_j_poisson)/np.log(r)

plt.plot(J, L_j_poisson_mean/n_sim_poisson, color = 'g', linestyle = "dotted")
J  = np.array(J)
plt.text(0.4, 1.65, "Poisson.", color = 'g')
plt.vlines(1, np.min(L_j_poisson_mean/n_sim_poisson), np.max(L_j_poisson_mean/n_sim_poisson), color = 'k')


plt.legend()
plt.show()