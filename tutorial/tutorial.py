import numpy as np
import matplotlib.pyplot as plt

from compute_estimator import confidence_interval
from compute_estimator  import alpha_hat
from tutorial import generate_pp

#Make sure to have the console on the folder Estim_Hyperuniformity

### Example: estimation for a poisson point process.

#Length side of the observation window [-R, R]^2
R = 20
#Generate a Poisson pattern of intensity 1
Phi = generate_pp.PPP(R)

#Plot the curve C
plt.subplot(1, 2, 2)
#The set of scales for which the curve C will be computed
J = np.linspace(0.2, 1.4, num= 100)

#Min and max degree of the 1d-hermites polynomials involved in the 2d-hermites wavelets
i_min = 0
i_max = 10

#Compute and plot of the curve C
L_j, r = alpha_hat.wavelet_transform(Phi, J, i_min, i_max)
C_j =  np.log(L_j)/np.log(r)
plt.plot(J, C_j)

#The set of scales used for estimating alpha with a linear regression
J_alpha= np.linspace(0.5, 1)

#Compute and print of alpha hat.
alpha_est= alpha_hat.alpha_hat(Phi, J_alpha, i_min, i_max)
print(alpha_est)

plt.vlines(np.min(J_alpha), np.min(C_j), np.max(C_j),color = 'k')
plt.vlines(np.max(J_alpha), np.min(C_j), np.max(C_j), color = 'k')

plt.xlabel("j")
plt.ylabel("$\mathcal{C}$")
plt.title("$\widehat{\\alpha} = $" +str(alpha_est)[:8])

#Plots the point pattern
plt.subplot(1, 2, 1)
plt.scatter(Phi[:,0], Phi[:, 1], s = 0.2)
ax = plt.gca()
ax.set_aspect('equal', adjustable="datalim")
plt.axis('off')
plt.title("Poisson point pattern in $[$-"+str(R)+","+str(R)+"$]^2$. \n $\\alpha$ =0")
print("Number of points :"+str(len(Phi)))

plt.show()


### Example: estimation for a perturbed lattice.

#Length side of the observation window [-R, R]^2
R = 35
alpha = 1

#Generate a perturbed lattice of intensity 1 with hyperuniformity exponent alpha
Phi = generate_pp.perturbed_latt(R, alpha, 0.3)

#Plots the point pattern
plt.subplot(1, 2, 1)
plt.scatter(Phi[:,0], Phi[:, 1], s = 0.25)
ax = plt.gca()
ax.set_aspect('equal', adjustable="datalim")
plt.axis('off')
plt.title("Perturbed lattice in $[$-"+str(R)+","+str(R)+"$]^2$ \n $\\alpha$ = "+str(alpha))
print("Number of points :"+str(len(Phi)))

#Plot the curve C
plt.subplot(1, 2, 2)

#The set of scales for which the curve C will be computed
J = np.linspace(0.2, 1.4, num = 100)

#Min and max degree of the 1d-hermites polynomials involved in the hermites wavelets
i_min = 0
i_max = 10
L_j, r = alpha_hat.wavelet_transform(Phi, J, i_min, i_max)
C_j = np.log(L_j)/np.log(r)
plt.plot(J, C_j)

#The set of scales used for estimating alpha with a linear regression
J_alpha = np.linspace(0.45, 1)
alpha_est = alpha_hat.alpha_hat(Phi, J_alpha, i_min, i_max)
print(alpha_est)

plt.vlines(np.min(J_alpha), np.min(np.log(L_j)/np.log(r)), np.max(np.log(L_j)/np.log(r)),color = 'k')
plt.vlines(np.max(J_alpha), np.min(np.log(L_j)/np.log(r)), np.max(np.log(L_j)/np.log(r)), color = 'k')

plt.xlabel("j")
plt.ylabel("$\mathcal{C}$")

plt.title("$\widehat{\\alpha} = $" +str(alpha_est)[:8])
plt.show()

### Example: computation of the asymptotic confidence interval for a poisson point process.

#Fix the seed in order to possibly re-use the precomputed convariance matrix in the /data/cov folder
np.random.seed(1)

#Length side of the observation window [-R, R]^2
R = 25
#Generate a Poisson pattern of intensity 1
Phi = generate_pp.PPP(R)

#Min and max degree of the 1d-hermites polynomials involved in the 2d-hermites wavelets
#To reduce computation time, the use i_max = 7 instead of i_max = 10.
i_min = 0
i_max = 7

#The set of scales used for estimating alpha with a linear regression
#To reduce computation time, only 20 scales have been used.
J = np.linspace(0.5, 1, num = 20)

#Compute and print of alpha hat.
alpha_est = alpha_hat.alpha_hat(Phi, J, i_min, i_max)
print(alpha_est)

q_1 = 0.025
q_2 = 0.975
a, b = confidence_interval.confidence_interval_cov_R(q_1, q_2, R, alpha_est, J, i_min, i_max)
print("Asymptotic confidence intervals \n")
print("["+str(a)[:6]+", "+str(b)[:6]+"]")
