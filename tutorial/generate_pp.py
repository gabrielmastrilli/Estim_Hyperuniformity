import numpy as np


#Poisson point process
def PPP(r):
    """
    Generate one realization of a poisson point process of intensity 1 in the observation window [-r, r]^2
    """
    n = np.random.poisson(4*r**2)
    Phi = np.zeros((n,2))
    for i in range(n):
        Phi[i, 0] = (2*np.random.uniform()-1)*r
        Phi[i, 1] = (2*np.random.uniform()-1)*r
    return Phi

#Perturbed lattice
def one_side_stable(a):
    """
    Condition:
    0 < a<1

    Simulate a one sided stable distribution with parameter a, i.e. a random variable X >= 0 with E[e^{-sX}] = e^{-s^{a}}.

    Justification of the algorithm is provied by Corollary 4.1. of Stable Densities Under Change of Scale and Total Variation Inequalities, Marek Kanter, 1975, Ann. Probab. 3(4): 697-707.
    """
    u_1 = np.random.uniform()
    u_2 = np.random.uniform()
    t = np.pi*u_1
    v = np.power(np.sin(a*t)/np.sin(t), a)*np.power(np.sin((1-a)*t)/np.sin(t), 1-a)
    return np.sqrt(np.power(v/np.abs(np.log(u_2)), (1 - a)/a))

def perturbed_latt(r, alpha, std):
    """

    Generate a realization of a perturbed lattice (p + U + U_p + sig V_p| p in Z^2)  in [-r, r]^2 with hyperuniformity exponent alpha.

    U is a uniform random variable on [-1/2, 1/2]^2
    (U_p) are i.i.d. uniform random variables on [-1/2, 1/2]^2
    (V_p) are i.i.d. with characteristic function phi(k) = exp(-|k|^(alpha))

    Conditions:
    0 < alpha < 2

    Recommendations:
    |std| < 1
    """

    x = []
    y = []
    U_x = np.random.uniform()-1/2
    U_y = np.random.uniform()-1/2

    #To handle border effects we simulate in [-2r, 2r]^2 and then keep the points in [-r, r]^2
    for i in range(-int(2*r), int(2*r)):
        for j in range(-int(2*r), int(2*r)):
            S = std*one_side_stable(alpha/2)
            n_1 = np.random.normal()
            n_2 = np.random.normal()
            u_x = np.random.uniform()-1/2
            u_y = np.random.uniform()-1/2
            if np.abs(i+ U_x+ u_x+ S*n_1) < r and np.abs(j+ U_y + u_y+ S*n_2) < r:
                x.append(i+ U_x + u_x+ S*n_1)
                y.append(j+ U_y + u_y+ S*n_2)
    Phi = np.zeros((len(x),2))
    for i in range(len(x)):
        Phi[i, 0] = x[i]
        Phi[i, 1] = y[i]

    return Phi

