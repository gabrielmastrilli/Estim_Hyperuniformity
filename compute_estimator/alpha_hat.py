import numpy as np

#Coefficients of the hermites polynomials.

"""

The 1d hermites polynomials are  stored as:
H_{2i}(x) = cs_h[i]*sum_{j = 0}^i co_h[j] x^{2j}
H_{2i+1}(x) = cs_h[i]*sum_{j = 0}^i co_h[j] x^{2j+1}

Warning: coefficients for i <= 20. For i > 20, errors starts to appears with this implementation.

"""

co_h= [[1], [2], [-2, 4], [-12, 8], [12, -48, 16], [120, -160, 32], [-120, 720, -480, 64], [-1680, 3360, -1344, 128], [1680, -13440, 13440, -3584, 256], [30240, -80640, 48384, -9216, 512], [-30240, 302400, -403200, 161280, -23040, 1024], [-665280, 2217600, -1774080, 506880, -56320, 2048], [665280, -7983360, 13305600, -7096320, 1520640, -135168, 4096], [17297280, -69189120, 69189120, -26357760, 4392960, -319488, 8192], [-17297280, 242161920, -484323840, 322882560, -92252160, 12300288, -745472, 16384], [-518918400, 2421619200, -2905943040, 1383782400, -307507200, 33546240, -1720320, 32768], [518918400, -8302694400, 19372953600, -15498362880, 5535129600, -984023040, 89456640, -3932160, 65536], [17643225600, -94097203200, 131736084480, -75277762560, 20910489600, -3041525760, 233963520, -8912896, 131072], [-17643225600, 317578060800, -846874828800, 790416506880, -338749931520, 75277762560, -9124577280, 601620480, -20054016, 262144], [-670442572800, 4022655436800, -6436248698880, 4290832465920, -1430277488640, 260050452480, -26671841280, 1524105216, -44826624, 524288]]

cs_h = np.zeros(len(co_h))
u = 1/np.power(np.pi, 1/4)
cs_h[0]= u
for i in range(1, len(co_h)):
    u /=(np.power(i, 1/2)*np.power(2, 1/2))
    cs_h[i] = u


#Computation of alpha_hat

def alpha_hat(Phi, J,  i_min, i_max):
    """
    Entries:

    Phi is a numpy array of size (number of points, 2); Phi[:, 0] is the x coordinate of the points
    J is a numpy array of size (number of scale); J is the set of scale

    The hermites functions that will be used are of index i = (i_1, i_2) in I with i_min <= min(i_1,i_2) <= max(i_1, i_2) <= i_max and i_1 or i_2 odd.

    Return

    The value of the estimator alpha_hat:
    alpha_hat =  2 - slope of {j -> sum_{i in I} |sum_{x \in Phi} 1_{[-r, r]}(x) psi_(5*x/r^j)|^2} on J

    Recommandations:

    Visualize before the curve C with compute_wavelet_transforms(Phi, J, i_min, i_max) in order to select the set of scale J.

    J < 1

    i_min = 0
    For large number of points:
    i_max = 10

    For small number of point, the estimator may be biased with i_max = 10. Biais can be reduced with smaller value of i_max (for example i_max = 7 or i_max = 5), but this could increase the variance of the estimation.

    """

    L_j, r = wavelet_transform(Phi, J, i_min, i_max)
    denom_W = (len(J)*np.sum(np.power(np.array(J),2)) - np.power(np.sum(np.array(J)), 2))
    sum_J = np.sum(np.array(J))
    slope_hat = 0
    for j in range(len(L_j)):
        slope_hat += (len(J)*J[j] - sum_J)/denom_W*np.log(L_j[j])/np.log(r)
    return 2 - slope_hat


def wavelet_transform(Phi, J, i_min, i_max):
    """
    Entries:

    Phi is a numpy array of size (number of points, 2); Phi[:, 0] is the x-coordinate of the points
    J is a numpy array of size (number of scale); J is the set of scale

    The Hermite functions that will be used are of index i = (i_1, i_2) in I with i_min <= min(i_1,i_2) <= max(i_1, i_2) <= i_max and i_1 or i_2 odd.


    Return:
    L_j = sum_{i in I} |sum_{x \in Phi} 1_{[-r, r]}(x) psi_(5*x/r^j)|^2
    r, the side  length of the windows [-r, r]^2 of the rescaled points

    Recommandations:

    J contains scales between 0.1 and 1.3
    i_min = 0
    i_max = 10

    """

    #Ensure that phi contains points in that are centered with respect to 0
    #Rescale the point pattern to get a realization of intensity of approximately 1

    p_1 = np.array(Phi[:, 0])- np.mean(Phi[:, 0])
    p_2 = np.array(Phi[:, 1]) - np.mean(Phi[:, 1])
    R = np.max(p_1)
    lambda_hat = np.sqrt(len(p_1)/(4*R**2))
    p_1 = p_1*lambda_hat
    p_2 = p_2*lambda_hat
    r = np.max(p_1)

    l= 5
    #Compute L_j = sum_{i in I} |sum_{x \in Phi} 1_{[-r, r]}(x) psi_(l*x/r^j)|^2
    #We compute Q[k_1, k_2] = sum_{x \in Phi} 1_{[-r, r]}(x) exp(-|l x/r^j|^2/2) (l x_1/r^j)^{k_1} (l x_2/r^j)^{k_2} if the computation has not been done before. Then, we use the expression of the hermites wavelets as polynomials to compute T_j using S[k_1, k_2]
    L_j = []
    for j in J:
        Q = np.ones((i_max, i_max))
        Q_is_computed = np.zeros((i_max, i_max))
        m = 0
        for i_1 in range(i_min, i_max):
            for i_2 in range(i_min, i_max):
                #To ensure that the hermites function is of zero integral
                if i_1%2 ==1 or i_2%2==1:
                    s = 0
                    u_1 = int(i_1%2==1)
                    u_2 = int(i_2%2==1)
                    for k_1 in range(len(co_h[i_1])):
                        for k_2 in range(len(co_h[i_2])):
                            #If S[2*k_1+u_1, 2*k_2+u_2] has not been computed
                            if Q_is_computed[2*k_1+u_1, 2*k_2+u_2]==0:
                                Q_is_computed[2*k_1+u_1, 2*k_2+u_2] = 1
                                a = np.power(l*p_1/r**j, 2*k_1+u_1)*np.exp(-np.power(l*p_1/r**j, 2)/2)
                                b = np.power(l*p_2/r**j, 2*k_2+u_2)*np.exp(-np.power(l*p_2/r**j, 2)/2)
                                Q[2*k_1+u_1, 2*k_2+u_2] = np.dot(a, b)
                            s+=cs_h[i_1]*cs_h[i_2]*co_h[i_1][k_1]*co_h[i_2][k_2]*Q[2*k_1+u_1, 2*k_2+u_2]
                    m+=s*s
        L_j.append(m)
    return L_j, r
