import os
import sys
import time
import numpy as np

from liblibra_core import *
from libra_py.workflows.nbra import mapping2
from libra_py import data_conv, data_outs, data_stat, units


dt = 40 # in atomic units
homo_index = 5
fragment = 1
num_occ_states = 6
num_unocc_states = 6
istep = 1
fstep = 2002
# ============================== buidling electron-only excitation basis

# Reading the energies
energy_filename = F'NAD/E_{istep}_{fragment}_1_re'
E = np.loadtxt(energy_filename)


print(E)


average_energy = np.mean(E)
print("Average Energy:", average_energy)


#sys.exit(0)
exc_energies = []
basis = []
GS = []

for occ_state in range(1,1+num_occ_states):
    GS.append(occ_state)
    GS.append(-occ_state)
print(GS)

basis.append(GS)
exc_energies.append(0)

for exc in range(num_occ_states, 0, -1):
    for unocc_state in range(num_occ_states+1,num_occ_states+num_unocc_states+1):
        exc_energy = E[unocc_state-1] - E[exc-1] 
        exc_energies.append(exc_energy)

        excited_state = []
        for occ_state in range(1,1+num_occ_states):
            if occ_state==exc:
                excited_state.append(unocc_state)
                excited_state.append(-occ_state)
            else:
                excited_state.append(occ_state)
                excited_state.append(-occ_state)
        #print(excited_state)
        basis.append(excited_state)

exc_energies = np.array(exc_energies)
new_indices = list(np.argsort(exc_energies))


basis_list = []
for i in range(len(new_indices)):
    basis_list.append(basis[ new_indices[i] ])
print(basis_list)
print(len(basis_list))


print(np.sort(exc_energies)*27.211385)

sys.exit(0)

# ================================ Computing the St of excited states Slater determinants
os.system(F'mkdir SD_basis_{fragment}')
for step in range(istep,fstep):        
    st_filename = F'NAD/S_{step}_{fragment}_1_tre'
    st_numpy = np.loadtxt(st_filename)    
    st_rks = data_conv.nparray2CMATRIX(st_numpy)
    #st_rks.real().show_matrix()
    st_sd = mapping2.ovlp_mat_arb(basis_list, basis_list, st_rks, reduce_det=0)
    st_sd.real().show_matrix()
    nac_sd = (st_sd - st_sd.H() )/(2.0*dt)

    st_sd.real().show_matrix(F'SD_basis_{fragment}/St_sd_{step}_re')
    st_sd.imag().show_matrix(F'SD_basis_{fragment}/St_sd_{step}_im')
    nac_sd.real().show_matrix(F'SD_basis_{fragment}/Hvib_sd_{step}_im')

    energy_filename = F'NAD/E_{step}_{fragment}_1_re'
    
    E_sd = []
    E_sd.append(0.0)
    E = np.loadtxt(energy_filename)
    for j in range(num_occ_states, num_occ_states+num_unocc_states):
        for i in range(0,num_occ_states):
            exc_ener = E[j] - E[i]
            E_sd.append(exc_ener) 

    E_sd = np.array(E_sd)
    E_sd = np.sort(E_sd)
    E_sd = np.diag( np.array(E_sd) )
    E_sd_MATRIX = data_conv.nparray2MATRIX(E_sd) 
    E_sd_MATRIX.show_matrix(F'SD_basis_{fragment}/Hvib_sd_{step}_re')








