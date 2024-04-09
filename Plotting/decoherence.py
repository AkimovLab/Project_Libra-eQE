import os
import time
import multiprocessing as mp
import matplotlib.pyplot as plt   # plots
import numpy as np
import scipy.sparse as sp
from scipy.optimize import curve_fit
import h5py
import warnings

from liblibra_core import *
import util.libutil as comn
from libra_py import units, data_conv, dynamics_plotting
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.dynamics.tsh.plot as tsh_dynamics_plot
import libra_py.data_savers as data_savers
import libra_py.workflows.nbra.decoherence_times as decoherence_times

from recipes import dish_nbra, fssh_nbra, fssh2_nbra, gfsh_nbra, ida_nbra,idf_nbra, mash_nbra, msdm_nbra

#from matplotlib.mlab import griddata
#%matplotlib inline 
warnings.filterwarnings('ignore')

colors = {}
colors.update({"11": "#8b1a0e"})  # red       
colors.update({"12": "#FF4500"})  # orangered 
colors.update({"13": "#B22222"})  # firebrick 
colors.update({"14": "#DC143C"})  # crimson   
colors.update({"21": "#5e9c36"})  # green
colors.update({"22": "#006400"})  # darkgreen  
colors.update({"23": "#228B22"})  # forestgreen
colors.update({"24": "#808000"})  # olive      
colors.update({"31": "#8A2BE2"})  # blueviolet
colors.update({"32": "#00008B"})  # darkblue  
colors.update({"41": "#2F4F4F"})  # darkslategray

clrs_index = ["11", "21", "31", "41", "12", "22", "32", "13","23", "14", "24"]
nthreads = 10
istep = 1    # the first timestep to read
fstep = 1772 # the last timestep to read
dt = 40   # integration time-step [a.u. of time]

nsteps = fstep - istep
NSTEPS = nsteps
print(F"Number of steps = {nsteps}")

x = np.loadtxt(F'SD_basis_1/Hvib_sd_1_im')
nstates = x.shape[0]
NSTATES = nstates                                                            
print(F"Number of states = {nstates}")

#sys.exit(0)
# ================= Reading the data

#================== Read energies =====================
E = []
for step in range(istep,fstep):
    energy_filename = F"SD_basis_1/Hvib_sd_{step}_re"
    energy_mat = np.loadtxt(energy_filename)
    # For data conversion we need to turn np.ndarray to np.array so that 
    # we can use data_conv.nparray2CMATRIX
    E.append( np.array( np.diag( energy_mat ) ) )
E = np.array(E)
NSTATES = E[0].shape[0]
#================== Read time-overlap =====================
St = []
for step in range(istep,fstep):        
    St_filename = F"SD_basis_1/St_sd_{step}_re"
    St_mat = np.loadtxt(St_filename)
    St.append( np.array( St_mat ) )
St = np.array(St)
#================ Compute NACs and vibronic Hamiltonians along the trajectory ============    
NAC = []
Hvib = [] 
for c, step in enumerate(range(istep,fstep)):
    nac_filename = F"SD_basis_1/Hvib_sd_{step}_im"
    nac_mat = np.loadtxt(nac_filename)
    NAC.append( np.array( nac_mat ) )
    Hvib.append( np.diag(E[c, :])*(1.0+1j*0.0)  - (0.0+1j)*nac_mat[:, :] )

NAC = np.array(NAC)
Hvib = np.array(Hvib)

print('Number of steps:', NSTEPS)
print('Number of states:', NSTATES)


# ================= Computing the energy gaps and decoherence times
HAM_RE = []
for step in range(E.shape[0]):
    HAM_RE.append( data_conv.nparray2CMATRIX( np.diag(E[step, : ]) ) )

# Average decoherence times and rates
tau, rates = decoherence_times.decoherence_times_ave([HAM_RE], [0], NSTEPS, 0)

# Computes the energy gaps between all states for all steps
dE = decoherence_times.energy_gaps_ave([HAM_RE], [0], NSTEPS)

# Decoherence times in fs
avg_deco = data_conv.MATRIX2nparray(tau) * units.au2fs

# Zero all the diagonal elements of the decoherence matrix
np.fill_diagonal(avg_deco, 0)

# Saving the average decoherence times
np.savetxt('avg_deco.txt',avg_deco.real)

# Computing the average energy gaps
gaps = MATRIX(NSTATES, NSTATES)
for step in range(NSTEPS):
    gaps += dE[step]
gaps /= NSTEPS

rates.show_matrix()
gaps.show_matrix()
#sys.exit(0)



class tmp:
    pass

def compute_model(q, params, full_id):
    timestep = params["timestep"]
#    icond = params["icond"]
#    maxstep = params["maxstep"]
#    timestep = (timestep + icond) % maxstep
    nst = params["nstates"]
    obj = tmp()

    obj.ham_adi = data_conv.nparray2CMATRIX( np.diag(E[timestep, : ]) )
    obj.nac_adi = data_conv.nparray2CMATRIX( NAC[timestep, :, :] )
    obj.hvib_adi = data_conv.nparray2CMATRIX( Hvib[timestep, :, :] )
    obj.basis_transform = CMATRIX(nst,nst); obj.basis_transform.identity()  #basis_transform
    obj.time_overlap_adi = data_conv.nparray2CMATRIX( St[timestep, :, :] )
    
    return obj




#================== Model parameters ====================
model_params = { "timestep":0, "icond":0,  "model0":0, "nstates":NSTATES }




#=============== Some automatic variables, related to the settings above ===================
NSTEPS = 1771
dyn_general = { "nsteps":NSTEPS, "ntraj":250, "nstates":NSTATES, "dt":dt,                                                 
                "decoherence_rates":rates, "ave_gaps":gaps,                
                "progress_frequency":0.1, "which_adi_states":range(NSTATES), "which_dia_states":range(NSTATES),
                "mem_output_level":2,
                "properties_to_save":[ "timestep", "time","se_pop_adi", "sh_pop_adi" ],
                "prefix":F"NBRA", "prefix2":F"NBRA", "isNBRA":0, "nfiles": nsteps - 1, "maxstep":885
              }

#=========== Select the method ===========
#dish_nbra.load(dyn_general); prf = "DISH"  # DISH
fssh_nbra.load(dyn_general); prf = "FSSH"  # FSSH
#fssh2_nbra.load(dyn_general); prf = "FSSH2"  # FSSH2

#gfsh_nbra.load(dyn_general); prf = "GFSH"  # GFSH
#idf_nbra.load(dyn_general); prf = "IDF"  # IDF

#ida_nbra.load(dyn_general); prf = "IDA"  # IDA
#mash_nbra.load(dyn_general); prf = "MASH"  # MASH
#msdm_nbra.load(dyn_general); prf = "MSDM"  # MSDM

#dyn_general.update( {"prefix":prf, "prefix2":prf } )


#=================== Dynamics =======================
# Nuclear DOF - these parameters don't matter much in the NBRA calculations
nucl_params = {"ndof":1, "init_type":3, "q":[-10.0], "p":[0.0], "mass":[2000.0], "force_constant":[0.01], "verbosity":-1 }

# Amplitudes are sampled
elec_params = {"ndia":NSTATES, "nadi":NSTATES, "verbosity":-1, "init_dm_type":0}
elec_params.update( {"init_type":1,  "rep":1,  "istate":6 } )  # how to initialize: random phase, adiabatic representation

if prf=="MASH":
    istates = list(np.zeros(NSTATES))
    istates[4] = 1.0
    elec_params.update( {"init_type":4,  "rep":1,  "istate":4, "istates":istates } )  # different initialization for MASH

def function1(icond):
    print('Running the calculations for icond:', icond)
    time.sleep( icond * 0.01 )
    rnd=Random()
    mdl = dict(model_params)
   # mdl.update({"icond": icond})  #create separate copy
   # dyn_general = dict(model_params)
    dyn_general.update({"icond": icond})
    dyn_gen = dict(dyn_general)
    dyn_gen.update({"prefix":F"{prf}_new88_icond_{icond}", "prefix2":F"{prf}_new88_icond_{icond}" })
    #res = tsh_dynamics.generic_recipe(dyn_gen, compute_model, elec_params, nucl_params, rnd)
    res = tsh_dynamics.generic_recipe(dyn_gen, compute_model, mdl, elec_params, nucl_params, rnd)


pool = mp.Pool(nthreads)
pool.map(function1, list(range(1,1002,200)))
pool.close()                            
pool.join()


