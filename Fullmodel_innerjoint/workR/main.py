#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import time
"""Script to run a direct dynamic analysis on a multibody system.

Summary
-------
This template loads the data file *.mbs and execute:
 - the coordinate partitioning module
 - the direct dynamic module (time integration of equations of motion).
 - if available, plot the time evolution of the first generalized coordinate.

It may have to be adapted and completed by the user.


Universite catholique de Louvain
CEREM : Centre for research in mechatronics

http://www.robotran.eu
Contact : info@robotran.be

(c) Universite catholique de Louvain
"""

# ============================================================================
# Packages loading
# =============================================================================
start_time = time.time()
try:
    import MBsysPy as Robotran
except:
    raise ImportError("MBsysPy not found/installed."
                      "See: https://www.robotran.eu/download/how-to-install/"
                      )
    
import sys
sys.path.insert(0, "C:\\Users\\ndine\\Documents\\MBProjects\\full_model_plot_test\\Fullmodel_innerjoint\\User_function")
import sys
sys.path.insert(1, "C:\\Users\\ndine\\Documents\\MBProjects\\full_model_plot_test\\Fullmodel_innerjoint\\userfctR")


# ===========================================================================
# Project loading
# =============================================================================
mbs_data = Robotran.MbsData('../dataR/Fullmodel_innerjoint.mbs')

# ===========================================================================
# Partitionning
# =============================================================================
mbs_data.process = 1
mbs_part = Robotran.MbsPart(mbs_data)
mbs_part.set_options(rowperm=1, verbose=1)
mbs_part.run()

# ===========================================================================
# Direct Dynamics
# =============================================================================
mbs_data.process = 3
mbs_dirdyn = Robotran.MbsDirdyn(mbs_data)
mbs_dirdyn.set_options(dt0=250e-7, tf=1.5, save2file=1)
results = mbs_dirdyn.run()


# ===========================================================================
# Plotting results
# =============================================================================
end_time = time.time()
execution_time = end_time - start_time
print("temps d'exécution : ", execution_time, "secondes.")

# import Neural_control_layer as neural
# import user_ExtForces as ext_forces

# left = neural.left_stim
# right = neural.right_stim
# left_phase = neural.left_phase
# right_phase = neural.right_phase




# print(left)

# # x = uf.graphe
# import matplotlib.pyplot as plt



# plt.scatter(left[:,7],left[:,0])
# plt.title("left VAS stim")
# plt.show()

# plt.scatter(left[:,7],left[:,1])
# plt.title("left SOL stim")
# plt.show()

# plt.scatter(left[:,7],left[:,2])
# plt.title("left GAS stim")
# plt.show()

# plt.scatter(left[:,7],left[:,3])
# plt.title("left TA stim")
# plt.show()

# plt.scatter(left[:,7],left[:,4])
# plt.title("left HAM stim")
# plt.show()

# plt.scatter(left[:,7],left[:,5])
# plt.title("left GLU stim")
# plt.show()

# plt.scatter(left[:,7],left[:,6])
# plt.title("left HFL stim")
# plt.show()

# plt.scatter(right[:,7],right[:,0])
# plt.title("right VAS stim")
# plt.show()

# plt.scatter(right[:,7],right[:,1])
# plt.title("right SOL stim")
# plt.show()

# plt.scatter(right[:,7],right[:,2])
# plt.title("right GAS stim")
# plt.show()

# plt.scatter(right[:,7],right[:,3])
# plt.title("right TA stim")
# plt.show()

# plt.scatter(right[:,7],right[:,4])
# plt.title("right HAM stim")
# plt.show()

# plt.scatter(right[:,7],right[:,5])
# plt.title("right GLU stim")
# plt.show()

# plt.scatter(right[:,7],right[:,6])
# plt.title("right HFL stim")
# plt.show()


# import user
# try:
#     import matplotlib.pyplot as plt
# except Exception:
#     raise RuntimeError('Unable to load matplotlib, plotting results unavailable.')



# Figure creation
# fig = plt.figure(num='Example of plot')
# axis = fig.gca()

# # Plotting data's

# axis.plot(results.q[:, 0],results.q[:, 7], label='q ankleL [rad]')
# axis.plot(current_t,ankleL_q,'--',color='black',label='ankleL Simulink [rad]')
# axis.set_title('ankle L position')


# axis.plot(results.q[:, 0], - results.q[:, 8], label='q kneeL [rad]')
# axis.plot(current_t,kneeL_q,'--',color='black',label='kneeL Simulink [rad]')
# axis.set_title('knee L position')


# axis.plot(results.q[:, 0],results.q[:, 10], label='q hipL [rad]')
# axis.plot(current_t,hipL_q,'--',color='black',label='HipL Simulink [rad]')
# axis.set_title('hip L position')


# axis.plot(results.q[:, 0],- results.q[:, 14], label='q ankleR [rad]')
# axis.plot(current_t,ankleR_q,'--',color='black',label='ankleR Simulink [rad]')
# axis.set_title('ankle R position')


# axis.plot(results.q[:, 0], results.q[:, 13], label='q kneeR [rad]')
# axis.plot(current_t,kneeR_q,'--',color='black',label='kneeR Simulink [rad]')
# axis.set_title('knee R position')


#axis.plot(results.q[:, 0],- results.q[:, 11], label='q hipR [rad]')
#axis.plot(current_t,hipR_q,'--',color='black',label='HipR Simulink [rad]')
#axis.set_title('hip R position')




# Figure enhancement
# axis.grid(True)
# axis.set_xlim(left=mbs_dirdyn.get_options('t0'), right=mbs_dirdyn.get_options('tf'))
# axis.set_xlabel('Time (s)')
# axis.set_ylabel('Coordinate value (m or rad)')
# axis.legend()
# plt.show()
