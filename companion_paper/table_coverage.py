import numpy as np

### Table of the coverage of alpha_hat for perturbed lattices.


compute_coverage_R = False
compute_coverage_ind = False

if compute_coverage_R:
    #Compute the coverage probability for the parameters used for estimating alpha for perturbed lattices using the non asymptotic covariance matrix.

    alphas = [0.5, 1, 1.5]
    # Side length of the observation windows [-R, R]^2
    Rs = [40, 35, 30, 25, 20, 15]
    #Lower bound for the scales used for estimating alpha
    j_min = [0.45, 0.5, 0.5, 0.55, 0.6, 0.65]
    i_min = 0
    i_max = 10
    p_insides = np.zeros((len(alphas), len(Rs)))
    for i_alpha in range(len(alphas)):
        alpha = alphas[i_alpha]
        for i_R in range(len(Rs)):
            R = Rs[i_R]
            J = np.linspace(j_min[i_R], 1)
            print(alpha, R)
            alpha_hats = np.loadtxt("data/alpha_"+str(alpha)+"_R_"+str(R)+".txt")
            #Compute the quantiles of order q_1 = 0.025 and q_2 = 0.975 of log(R)(alpha_hat - alpha).
            q_1, q_2 =  estimate_quantile_cov_R(0.025, 0.975, R, alpha, J, i_min, i_max)
            print(q_1, q_2)
            n_inside = 0
            for alpha_hat in alpha_hats:
                if alpha_hat > alpha + q_1 and alpha_hat < alpha+q_2:
                    n_inside +=1
            p_inside = n_inside/500
            print(p_inside)
            p_insides[i_alpha, i_R]= p_inside
    np.savetxt("data/nb_inside.txt", p_insides)
    print(p_inside)

if compute_coverage_ind:
    #Compute the coverage probability for the parameters used for estimating alpha for perturbed lattices using the asymptotic covariance matrix.
    alphas = [0.5, 1, 1.5]
    # Side length of the observation windows [-R, R]^2
    Rs = [40, 35, 30, 25, 20, 15]
    #Lower bound for the scales used for estimating alpha
    j_min = [0.45, 0.5, 0.5, 0.55, 0.6, 0.65]
    i_min = 0
    i_max = 10
    p_insides = np.zeros((len(alphas), len(Rs)))
    for i_alpha in range(len(alphas)):
        alpha = alphas[i_alpha]
        for i_R in range(len(Rs)):
            R = Rs[i_R]
            J = np.linspace(j_min[i_R], 1)
            print(alpha, R)
            #Compute the quantiles of order q_1 = 0.025 and q_2 = 0.975 of log(R)(alpha_hat - alpha).
            alpha_hats = np.loadtxt("data/alpha_"+str(alpha)+"_R_"+str(R)+".txt")
            q_1, q_2 =  estimate_quantile_ind(0.025, 0.975, R, alpha, J, i_min, i_max)
            print(q_1, q_2)
            n_inside = 0
            for alpha_hat in alpha_hats:
                if alpha_hat > alpha + q_1 and alpha_hat < alpha+q_2:
                    n_inside +=1
            p_inside = n_inside/500
            print(p_inside)
            p_insides[i_alpha, i_R]= p_inside
    np.savetxt("data/nb_inside_ind.txt", p_insides)
    print(p_inside)