import numpy as np
import scipy.linalg as lg
import os
import scipy as sp
import scipy.special as spp


#Coefficients of the hermites polynomials.

"""
H_i(x) = cs_h[i]*sum_{j = 0}^i co_h[j] x^j
Warning: coefficients for i <= 20. For i > 20, errors starts to appears with this implementation.
"""

co_h= [[1], [2], [-2, 4], [-12, 8], [12, -48, 16], [120, -160, 32], [-120, 720, -480, 64], [-1680, 3360, -1344, 128], [1680, -13440, 13440, -3584, 256], [30240, -80640, 48384, -9216, 512], [-30240, 302400, -403200, 161280, -23040, 1024], [-665280, 2217600, -1774080, 506880, -56320, 2048], [665280, -7983360, 13305600, -7096320, 1520640, -135168, 4096], [17297280, -69189120, 69189120, -26357760, 4392960, -319488, 8192], [-17297280, 242161920, -484323840, 322882560, -92252160, 12300288, -745472, 16384], [-518918400, 2421619200, -2905943040, 1383782400, -307507200, 33546240, -1720320, 32768], [518918400, -8302694400, 19372953600, -15498362880, 5535129600, -984023040, 89456640, -3932160, 65536], [17643225600, -94097203200, 131736084480, -75277762560, 20910489600, -3041525760, 233963520, -8912896, 131072], [-17643225600, 317578060800, -846874828800, 790416506880, -338749931520, 75277762560, -9124577280, 601620480, -20054016, 262144], [-670442572800, 4022655436800, -6436248698880, 4290832465920, -1430277488640, 260050452480, -26671841280, 1524105216, -44826624, 524288]]

cs_h = np.zeros(len(co_h))
u = 1/np.power(np.pi, 1/4)
cs_h[0]= u
for i in range(1, len(co_h)):
    u /=(np.power(i, 1/2)*np.power(2, 1/2))
    cs_h[i] = u



#I_1 contains store value of P(p, q) for the range of (p, q) that will be used for the computation of Sigma_R.
I_1_is_computed = np.zeros((len(cs_h),len(cs_h)))
I_1 = np.ones((len(cs_h),len(cs_h)))

def P(p, q):
    """
    Compute int_0^{2pi} cos(x)^p sin(x)^q dx using the recursion formula of Proposition 5.13 of the companion papers.
    """
    if p%2==1 or q%2==1:
        return 0
    else:
        if I_1_is_computed[p, q] ==1:
            return I_1[p, q]
        if I_1_is_computed[q, p] ==1:
            return I_1[q, p]
        if p == 0 and q==0:
            return 2*np.pi
        elif q==0 or p== 0:
            o = 1
            k = int(p/2)+int(q/2)
            for i in range(1, k+1):
                o*=(2*i-1)/(2*i)
            I_1[p, q] = 2*np.pi*o
            I_1[q, p] = 2*np.pi*o
            I_1_is_computed[p, q] = 1
            I_1_is_computed[q, p] = 1
            return 2*np.pi*o
        elif I_1_is_computed[p-2, q-2] ==1:
            return (p-1)*(q-1)/((p+q)*(p+q-2))*I_1[p-2, q-2]
        elif I_1_is_computed[q-2, p-2] ==1:
            return (p-1)*(q-1)/((p+q)*(p+q-2))*I_1[q-2, p-2]
        else:
            a = (p-1)*(q-1)/((p+q)*(p+q-2))*P(p-2, q-2)
            I_1[p, q] = a
            I_1[q, p] = a
            I_1_is_computed[q, p] = 1
            I_1_is_computed[p, q] = 1
            return a


def cov_matrix(r, alpha, J, i_min, i_max):
    """
    Using the Proposition 5.13 of the companion paper, compute the coefficient of the matrix:
    Sig[i_1,i_2, j_1, j_2] = r^{(d+alpha)(j_1 +j_2)/2} int psi_{i_1}(r^{j_1} x) conjugate{psi_{i_2}(r^{j_2} x)} |x|^{alpha}, (j_1, j_2) in J^2 and (i_1, i_2) in I, where I = \{i = (i_x, i_y) in N^2| i_min =< i_x, i_y < i_max and i_x or i_y is odd.

    Warning:

    Computation time can be long when i_max - i_min and/or |J| is large. For i_min = 0 and i_max = 10, and |J|= 50, the computation took 1.5 hours on a standard computer. For the previous setting the matrix contains approximately 14 000 000 millions coefficients.

    To reduce the computation time, we recommand to use fewer scales (for example |J| = 20).
    """

    #Computed the indexes of the involved hermites functions.
    I_ = []
    for i_1 in range(i_min, i_max):
        for i_2 in range(i_min, i_max):
            if (i_1%2==1) or (i_2%2==1):
                I_.append([i_1, i_2])

    #For the computation the covariance matrix is implented as a 4-dimensionnal numpy array
    S = np.zeros((len(I_), len(I_), len(J), len(J)))

    #In order to reduce the computation time, we keep track of the already computed coefficient, to exploit the sympetries of the matrix Sigma_R.
    S_is_computed = np.zeros((len(I_), len(I_), len(J), len(J)), dtype = complex)

    #Save the computed values of the gamma function in order to reduce computation time. If gamma_save[p] = -1, then gamma_save[p] has not been yet computed.
    gamma_save = -np.ones(5*i_max)

    n_progress = 0
    for i_j_1 in range(len(J)):
        for i_j_2 in range(len(J)):
            j_1 = J[i_j_1]
            j_2 = J[i_j_2]
            for i_1 in range(len(I_)):
                for i_2 in range(len(I_)):
                    if n_progress%2500 == 0:
                        print(str(n_progress/(len(J)**2*len(I_)**2)*100)[:4]+"%")
                    n_progress +=1
                    if S_is_computed[i_2, i_1, i_j_2, i_j_1] == 0:
                        i_1_1 = I_[i_1][0]
                        i_1_2 = I_[i_1][1]
                        i_2_1 = I_[i_2][0]
                        i_2_2 = I_[i_2][1]
                        s = 0
                        u_1 = int(i_1_1%2==1)
                        u_2 = int(i_1_2%2==1)
                        u_3 = int(i_2_1%2==1)
                        u_4 = int(i_2_2%2==1)
                        #If (i_1_1 - i_2_1)%2==1 or (i_1_2 - i_2_2)%2==1, by parity argument S[i_2, i_1, i_j_2, i_j_1] = 0.
                        if (i_1_1 - i_2_1)%2==0 and (i_1_2 - i_2_2)%2==0:
                            #Range over the coefficient of the hermites polynomials.
                            for l_1 in range(len(co_h[i_1_1])):
                                for l_2 in range(len(co_h[i_1_2])):
                                    for l_3 in range(len(co_h[i_2_1])):
                                        for l_4 in range(len(co_h[i_2_2])):
                                            a = 0
                                            if I_1_is_computed[2*l_1+u_1+2*l_3+u_3, 2*l_2+u_2+2*l_4+u_4] == 0:
                                                a = P(2*l_1+u_1+2*l_3+u_3, 2*l_2+u_2+2*l_4+u_4)
                                                I_1[2*l_1+u_1+2*l_3+u_3, 2*l_2+u_2+2*l_4+u_4] = a
                                            else:
                                                a = I_1[2*l_1+u_1+2*l_3+u_3, 2*l_2+u_2+2*l_4+u_4]
                                            if np.abs(a) > 0:
                                                b = 0
                                                if gamma_save[2*(l_1+l_2+l_3+l_4)+u_1+u_2+u_3+u_4] == -1:
                                                    b  = 0.5*spp.gamma((2*(l_1+l_2+l_3+l_4)+u_1+u_2+u_3+u_4 + alpha +2)/2)
                                                    gamma_save[2*(l_1+l_2+l_3+l_4)+u_1+u_2+u_3+u_4] = b
                                                else:
                                                    b = gamma_save[2*(l_1+l_2+l_3+l_4)+u_1+u_2+u_3+u_4]
                                                c= r**((1+alpha/2)*(j_2 + j_1)+j_1*(2*l_1+u_1+2*l_2+u_2)+j_2*(2*l_3+u_3+2*l_4+u_4))/(np.power(np.sqrt((r**(2*j_1)+r**(2*j_2))/2), 2+alpha+2*l_1+u_1+2*l_2+u_2+2*l_3+u_3+2*l_4+u_4))
                                                s+= (-1)**((i_1_1 - i_2_1)/2+(i_1_2 - i_2_2)/2)*co_h[i_1_1][l_1]*co_h[i_1_2][l_2]*co_h[i_2_1][l_3]*co_h[i_2_2][l_4]*cs_h[i_1_1]*cs_h[i_1_2]*cs_h[i_2_1]*cs_h[i_2_2]*a*b*c
                            S[i_1, i_2, i_j_1, i_j_2] = s
                            S[i_2, i_1, i_j_2, i_j_1] = s
                            S_is_computed[i_1, i_2, i_j_1, i_j_2] = 1
                            S_is_computed[i_2, i_1, i_j_2, i_j_1] = 1
    #Store the 4-dimensionnal array S intro a matrix S
    S_ = np.zeros((len(I_)*len(J),len(I_)*len(J)))
    for i_j_1 in range(len(J)):
        for i_j_2 in range(len(J)):
            for i_1 in range(len(I_)):
                for i_2 in range(len(I_)):
                    S_[i_j_1*len(I_)+i_1, i_j_2*len(I_)+i_2] = S[i_1, i_2, i_j_1, i_j_2]
    return S_


def quantile_non_asymptotic(q_1, q_2, r, alpha, J, i_min, i_max):
    """
    Compute the theorical quantiles of order q_1 and q_2 of alpha_hat, using the non asymptotic covariance matrix Sigma_R defined in Propostion 3.14.
    """
    print("Computation cov matrix: start")
    if os.path.exists("./companion_paper/data/cov/sig_alpha_"+str(alpha)+"_r_"+str(r)+"_jmin_"+str(np.min(J))+"_jmax"+str(np.max(J))+"_imin_"+str(i_min)+"_imax_"+str(i_max)+".txt") == False:
        Sig = cov_matrix(r, alpha, J, i_min, i_max)
        np.savetxt("./companion_paper/data/cov/sig_alpha_"+str(alpha)+"_r_"+str(r)+"_jmin_"+str(np.min(J))+"_jmax"+str(np.max(J))+"_imin_"+str(i_min)+"_imax_"+str(i_max)+".txt", Sig)
    else:
        print("Already computed. Loading...")
        Sig = np.loadtxt("./companion_paper/data/cov/sig_alpha_"+str(alpha)+"_r_"+str(r)+"_jmin_"+str(np.min(J))+"_jmax"+str(np.max(J))+"_imin_"+str(i_min)+"_imax_"+str(i_max)+".txt")
    print("Computation cov matrix: ok")
    print("Cholesky decomposition: start")
    if os.path.exists("./companion_paper/data/cov/sig12_alpha_"+str(alpha)+"_r_"+str(r)+"_jmin_"+str(np.min(J))+"_jmax"+str(np.max(J))+"_imin_"+str(i_min)+"_imax_"+str(i_max)+".txt") == False:
        #Due to possibly large size of the matrice numerical errors could appear leading to small negative eigenvalues of the covariance matrix.
        #Theses numerical errors leads to small imaginary part of the square root of the covariance matrix
        #We thus consider the real part of the square root.
        Sig12 = np.real(lg.sqrtm(Sig))
        np.savetxt("./companion_paper/data/cov/sig12_alpha_"+str(alpha)+"_r_"+str(r)+"_jmin_"+str(np.min(J))+"_jmax"+str(np.max(J))+"_imin_"+str(i_min)+"_imax_"+str(i_max)+".txt", Sig12)
    else:
        print("Already computed. Loading...")
        Sig12 = np.loadtxt("./companion_paper/data/cov/sig12_alpha_"+str(alpha)+"_r_"+str(r)+"_jmin_"+str(np.min(J))+"_jmax"+str(np.max(J))+"_imin_"+str(i_min)+"_imax_"+str(i_max)+".txt")
    print("Cholesky decomposition: ok")

    #Number of independent realization of the theoritical limit distribution used to  estimate the quantiles of alpha_hat
    N_sim = 10000
    T_ = []

    denom_W = (len(J)*np.sum(np.power(np.array(J),2)) - np.power(np.sum(np.array(J)), 2))
    sum_J = np.sum(np.array(J))
    W = (len(J)*J - sum_J)/denom_W

    I_ = []
    for i_1 in range(i_min, i_max):
        for i_2 in range(i_min, i_max):
            if (i_1%2==1) or (i_2%2==1):
                I_.append([i_1, i_2])

    print("Estimation quantile: start")
    for i in range(N_sim):
        T = 0
        if i%1000==0:
            print(str(i/N_sim*100)+"%")
        Z = np.random.normal(0, 1, size = len(I_)*len(J))
        N = np.dot(Sig12, Z)
        for j in range(len(J)):
            T+= W[j]*np.log(np.sum(np.power(N[j*len(I_):(j+1)*len(I_)], 2)))
        T_.append(T)

    return np.quantile(T_, q_1)/np.log(r), np.quantile(T_, q_2)/np.log(r)

def quantile_asymptotic(q_1, q_2, r, alpha, J, i_min, i_max):
    """
    Compute the theorical quantiles of order q_1 and q_2 of alpha_hat - alpha, using the asymptotic covariance matrix Sigma.
    """
    print("Computation cov matrix: start")
    if os.path.exists("./companion_paper/data/cov/sig_alpha_"+str(alpha)+"_r_"+str(r)+"_jmin_"+str(np.min(J))+"_jmax"+str(np.max(J))+"_imin_"+str(i_min)+"_imax_"+str(i_max)+".txt") == False:
        Sig = cov_matrix(r, alpha, J, i_min, i_max)
        np.savetxt("./companion_paper/data/cov/sig_alpha_"+str(alpha)+"_r_"+str(r)+"_jmin_"+str(np.min(J))+"_jmax"+str(np.max(J))+"_imin_"+str(i_min)+"_imax_"+str(i_max)+".txt", Sig)
    else:
        print("Already computed. Loading...")
        Sig = np.loadtxt("./companion_paper/data/cov/sig_alpha_"+str(alpha)+"_r_"+str(r)+"_jmin_"+str(np.min(J))+"_jmax"+str(np.max(J))+"_imin_"+str(i_min)+"_imax_"+str(i_max)+".txt")
    print("Computation cov matrix: ok")
    print("Cholesky decomposition: start")
    I_ = []
    for i_1 in range(i_min, i_max):
        for i_2 in range(i_min, i_max):
            if (i_1%2==1) or (i_2%2==1):
                I_.append([i_1, i_2])

    #Due to possibly large size of the matrice numerical errors could appear leading to small negative eigenvalues of the covariance matrix.
    #Theses numerical errors leads to small imaginary part of the square root of the covariance matrix
    #We thus consider the real part of the square root.

    #We use the already computed covariance matrix Sigma_R, in which the asymptotic covariance matrix appears in the first top-left block
    Sig12_ind = np.real(lg.sqrtm(Sig[:len(I_), :len(I_)]))
    print("Cholesky decomposition: ok")

    #Number of independent realization of the theoritical limit distribution used to  estimate the quantiles of alpha_hat
    N_sim = 10000
    T_ = []
    denom_W = (len(J)*np.sum(np.power(np.array(J),2)) - np.power(np.sum(np.array(J)), 2))
    sum_J = np.sum(np.array(J))
    W = (len(J)*J - sum_J)/denom_W


    print("Estimation quantile: start")
    for i in range(N_sim):
        T = 0
        if i%1000==0:
            print(str(i/N_sim*100)+"%")
        for j in range(len(J)):
            Z = np.random.normal(0, 1, size = len(I_))
            N = np.dot(Sig12_ind, Z)
            T+= W[j]*np.log(np.sum(np.power(N, 2)))
        T_.append(T)
    return np.quantile(T_, q_1)/np.log(r), np.quantile(T_, q_2)/np.log(r)

def confidence_interval_cov_R(q_1, q_2, r, alpha_hat, J, i_min, i_max):
    """
    Compute the asymptotic confidence intervals of order q_2 - q_1:

    [alpha_hat - F^{-1}(q_2; alpha_hat)/log(r), alpha_hat -F^{-1}(q_1; alpha_hat)/log(r)],

    when alpha_hat has been estimated with scales J and with Hermites functions indexed by I = \{i = (i_x, i_y) in N^2| i_min =< i_x, i_y < i_max and i_x or i_y is odd.

    Warning:

    Computation time can be long when i_max - i_min and/or |J| is large. For i_min = 0 and i_max = 10, and |J|= 50, the computation took 1.5 hours on a standard computer. For the previous setting the matrix contains approximately 14 000 000 millions coefficients.

    To reduce the computation time, we recommand to use fewer scales (for example |J| = 20).
    """

    q_q_1, q_q_2 = quantile_non_asymptotic(q_1, q_2, r, alpha_hat, J, i_min, i_max)
    return alpha_hat - q_q_2, alpha_hat - q_q_1

def confidence_interval_cov_asymptotic(q_1, q_2, r, alpha_hat, J, i_min, i_max):
    """
    Compute the asymptotic confidence intervals of order q_2 - q_1 of the Corollary 3.10 of the companion paper:

    [alpha_hat - F^{-1}(q_2; alpha_hat)/log(r), alpha_hat -F^{-1}(q_1; alpha_hat)/log(r)]
    when alpha_hat has been estimated with scales J and with Hermites functions indexed by I = \{i = (i_x, i_y) in N^2| i_min =< i_x, i_y < i_max and i_x or i_y is odd.
    and where F^{-1}(q_2; alpha_hat) is the quantile of order q_2 of the random variable
    Warning:

    Note that this function has been implemented only for comparison with `confidence_interval_cov_R` and should not be used as it could lead to non accurate confidence intervals for non really large point patterns (refer to the discussion of the end of Section 4.2 of the companion paper)
    Computation time can be long when i_max - i_min and/or |J| is large. For i_min = 0 and i_max = 10, and |J|= 50, the computation took 1.5 hours on a standard computer. For the previous setting the matrix contains approximately 14 000 000 millions coefficients.
    To reduce the computation time, we recommand to use fewer scales (for example |J| = 20).
    """

    q_q_1, q_q_2 = quantile_asymptotic(q_1, q_2, r, alpha_hat, J, i_min, i_max)
    return alpha_hat - q_q_2, alpha_hat - q_q_1

