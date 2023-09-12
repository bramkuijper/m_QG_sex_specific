#!/usr/bin/env python3


# python script to systematically vary 
# parameters

import numpy as np
import math,sys
import socket

# initial values of hm, hp, m0, m1
init_vals = [ 0.0, 0.0, 0.0, 0.0 ]
mu_hm = [ 0.01 ]
mu_hp = [ 0.01 ]
mu_m0 = [ 0.01 ]
mu_m1 = [ 0.0 ]
sdmu = 0.01

omega_y = math.sqrt(0.1)
omega_hm = math.sqrt(1000)
omega_hp = math.sqrt(1000)
omega_m0 = math.sqrt(1000)
omega_m1 = math.sqrt(1000)

sd_mat_error = 0.0
sigma_e = [0.32]

nsteps = 20
#sinusoidal_frequency = list(np.linspace(0.0, math.pi, num=nsteps))
sinusoidal_frequency = [0]
gauss_autocorr = list(np.linspace(-1,1.0,num=nsteps))
gauss_sd_eps = 1.0

A = 0.0
#B = list(np.linspace(2.0,4.0,num=10))

B = [1.0]

s0 = 0.1

ctr = 0

exe = "./m_matpat.exe"

# number of replicates of each 
nreps = 8

# calculate number of jobs
njobs = len(mu_hm) * len(mu_hp) * len(mu_m0) * len(mu_m1) * len(sinusoidal_frequency) * len(B) * len(gauss_autocorr) * nreps

print(njobs)

for rep_i in range(0, nreps):
    for mu_hm_i in mu_hm:
        for mu_hp_i in mu_hp:
            for mu_m0_i in mu_m0:
                for mu_m1_i in mu_m1:
                    for freq_i in sinusoidal_frequency:
                        for B_i in B:
                            for sigma_e_i in sigma_e:
                                for autocorr_i in gauss_autocorr:
                                    print("echo " + str(ctr))
                                    ctr += 1

                                    sinusval = 0

                                    output_str = (exe + " "
                                                    + " ".join([ str(val) for val in init_vals]) + " "
                                                    + str(mu_hm_i) + " "
                                                    + str(mu_hp_i) + " "
                                                    + str(mu_m0_i) + " "
                                                    + str(mu_m1_i) + " "
                                                    + " ".join([ str(sdmu) for i in range(0,4)]) + " "
                                                    + str(omega_y) + " "
                                                    + str(omega_hm) + " "
                                                    + str(omega_hp) + " "
                                                    + str(omega_m0) + " "
                                                    + str(omega_m1) + " "
                                                    + str(sd_mat_error) + " "
                                                    + str(freq_i) + " "
                                                    + str(autocorr_i) + " "
                                                    + str(gauss_sd_eps) + " "
                                                    + str(A) + " "
                                                    + str(B_i) + " "
                                                    + str(s0) + " "
                                                    + str(sigma_e_i) + " ")

                                    print(output_str)
