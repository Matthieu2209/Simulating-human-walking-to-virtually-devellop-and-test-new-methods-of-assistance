#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 14:58:27 2023

@author: matthieuxaussems
"""
import pandas as pd
import useful_functions as u_f
import Muscle_actuation_layer as muscle
import matplotlib.pyplot as plt
import numpy as np

joint_limits = pd.read_excel ("/Users/matthieuxaussems/Documents/MBProjects/Fullmodel_innerjoint/User_function/jointlimits.xlsx")
# Tupdate_ankle = pd.read_excel ("/Users/matthieuxaussems/Documents/MBProjects/Fullmodel_innerjoint/User_function/Torque_update_ankle.xlsx")
# Tupdate_knee = pd.read_excel ("/Users/matthieuxaussems/Documents/MBProjects/Fullmodel_innerjoint/User_function/Tupdate_knee.xlsx")
# Tupdate_hip = pd.read_excel ("/Users/matthieuxaussems/Documents/MBProjects/Fullmodel_innerjoint/User_function/Tupdate_hip.xlsx")
#pressure_sheet = pd.read_excel ("/Users/matthieuxaussems/Documents/MBProjects/Fullmodel_innerjoint/User_function/pressure_sheet.xlsx")


VAS =0
SOL =1
GAS = 2
TA = 3
HAM = 4
GLU = 5
HFL = 6

current_t = joint_limits["time"].to_numpy()

ankleL_q = joint_limits["ankleL_q"].to_numpy()
kneeL_q = joint_limits["kneeL_q"].to_numpy()
ankleL_qd = joint_limits["ankleL_qd"].to_numpy()
kneeL_qd = joint_limits["kneeL_qd"].to_numpy()
knee_stop = joint_limits["Knee_stop"].to_numpy()
ankle_stop = joint_limits["Ankle_stop"].to_numpy()

Knee_jl = np.zeros(len(current_t))
Ankle_jl = np.zeros(len(current_t))
Hip_jl=np.zeros(len(current_t))

# current_t = Tupdate_ankle["time"].to_numpy()
# ankleL_q = Tupdate_ankle["ankle_q"].to_numpy()
# Fm_GAS_a = Tupdate_ankle["Fm_GAS"].to_numpy()
# Fm_SOL = Tupdate_ankle["Fm_SOL"].to_numpy()
# Fm_TA = Tupdate_ankle["Fm_TA"].to_numpy()
# T_GAS_a = Tupdate_ankle["TorqueGASankle"].to_numpy()
# T_SOL = Tupdate_ankle["TorqueSOLankle"].to_numpy()
# T_TA = Tupdate_ankle["TorqueTAankle"].to_numpy()
# T_ankle = Tupdate_ankle["torque_ankle"].to_numpy()

# current_t = Tupdate_knee["time"].to_numpy()
# kneeL_q = Tupdate_knee["knee_q"].to_numpy()
# Fm_GAS_k = Tupdate_knee["Fm_GAS"].to_numpy()
# Fm_HAM_k = Tupdate_knee["Fm_HAM"].to_numpy()
# Fm_VAS = Tupdate_knee["Fm_VAS"].to_numpy()
# T_GAS_k = Tupdate_knee["T_GAS"].to_numpy()
# T_VAS = Tupdate_knee["T_VAS"].to_numpy()
# T_HAM_k = Tupdate_knee["T_HAM"].to_numpy()
# T_knee = Tupdate_knee["T_knee"].to_numpy()

# current_t = Tupdate_hip["time"].to_numpy()
# hipL_q = Tupdate_hip["hip_q"].to_numpy()
# Fm_GLU = Tupdate_hip["Fm_GLU"].to_numpy()
# Fm_HAM_h = Tupdate_hip["Fm_HAM"].to_numpy()
# Fm_HFL = Tupdate_hip["Fm_HFL"].to_numpy()
# T_GLU = Tupdate_hip["T_GLU"].to_numpy()
# T_HFL = Tupdate_hip["T_HFL"].to_numpy()
# T_HAM_h = Tupdate_hip["T_HAM"].to_numpy()
# T_hip = Tupdate_hip["T_hip"].to_numpy()

#current_t = pressure_sheet["time"].to_numpy()
#v = pressure_sheet["v"].to_numpy()
#p = pressure_sheet["p"].to_numpy()
#F_final = pressure_sheet["F_final"].to_numpy()


############# test des joints limits 

for i in range(len(current_t)):
    
    Ankle_jl[i]=muscle.joint_limits(0, ankleL_q[i], ankleL_qd[i])
    Knee_jl[i]=muscle.joint_limits(1, kneeL_q[i], kneeL_qd[i])


#Figure creation
fig = plt.figure(num='Example of plot')
axis = fig.gca()
# axis.plot(current_t,Knee_jl, label='Knee joint limits')
# axis.plot(current_t,knee_stop,'--',color='black',label='Knee joint limits simulink')
# axis.set_title('Knee Joint limits')
# axis.scatter(current_t,Knee_jl-knee_stop,s=1,label='Knee joint limits')
# axis.set_title('diff Knee Joint limits')

# axis.plot(current_t,Ankle_jl, label='Ankle joint limits')
axis.plot(current_t,ankle_stop,'--',color='black',label='Ankle joint limits simulink')
axis.set_title('Ankle Joint limits')
#axis.scatter(current_t,Ankle_jl-ankle_stop,s=1,label='Ankle joint limits')
#axis.set_title('diff Ankle Joint limits')

#Figure enhancement
axis.grid(True)
axis.set_xlabel('Time (s)')
axis.set_ylabel('Torque')
axis.legend()
plt.show()

# Torque_TA_ankle = np.zeros(len(current_t))
# Torque_GAS_ankle = np.zeros(len(current_t))
# Torque_SOL_ankle = np.zeros(len(current_t))

# Torque_GAS_knee = np.zeros(len(current_t))
# Torque_VAS_knee = np.zeros(len(current_t))
# Torque_HAM_knee = np.zeros(len(current_t))

# Torque_HAM_hip = np.zeros(len(current_t))
# Torque_GLU_hip = np.zeros(len(current_t))
# Torque_HFL_hip = np.zeros(len(current_t))

# Torque_ankle = np.zeros(len(current_t))
# Torque_knee = np.zeros(len(current_t))
# Torque_hip = np.zeros(len(current_t))

################ test de torque compute / torque update 

# for i in range(len(current_t)):
    
#     Torque_GAS_ankle[i] = muscle.torque_updateANKLE(ankleL_q[i], Fm_GAS_a[i], GAS)
#     Torque_TA_ankle[i] = muscle.torque_updateANKLE(ankleL_q[i], Fm_TA[i], TA)
#     Torque_SOL_ankle[i] = muscle.torque_updateANKLE(ankleL_q[i], Fm_SOL[i], SOL)
    
#     Torque_GAS_knee[i] = muscle.torque_updateKNEE(kneeL_q[i], Fm_GAS_k[i], GAS)
#     Torque_VAS_knee[i] = muscle.torque_updateKNEE(kneeL_q[i], Fm_VAS[i], VAS)
#     Torque_HAM_knee[i] = muscle.torque_updateKNEE(kneeL_q[i], Fm_HAM_k[i], HAM)
    
#     Torque_HAM_hip[i] = muscle.torque_updateHIP(hipL_q[i], Fm_HAM_h[i], HAM)
#     Torque_GLU_hip[i] = muscle.torque_updateHIP(hipL_q[i], Fm_GLU[i], GLU)
#     Torque_HFL_hip[i] = muscle.torque_updateHIP(hipL_q[i], Fm_HFL[i], HFL)
    
    
#     Torque_ankle[i] = Ankle_jl[i] + Torque_GAS_ankle[i] + Torque_SOL_ankle[i] - Torque_TA_ankle[i] 
#     Torque_knee[i] = Knee_jl[i] + Torque_VAS_knee[i] - Torque_GAS_knee[i] - Torque_HAM_knee[i]
#     Torque_hip[i] = Hip_jl[i] + Torque_GLU_hip[i] + Torque_HAM_hip[i] - Torque_HFL_hip[i]

# Figure creation
# fig = plt.figure(num='Example of plot')
# axis = fig.gca()
#axis.plot(current_t,Torque_hip, label='ankle torque')
#axis.plot(current_t,T_ankle,'--',color='black',label='ankle GAS torque simulink')
#axis.set_title('Ankle torque')
# axis.scatter(current_t,Torque_knee-T_knee,s=1,label='knee torque')
# axis.set_title('diff knee HAM torque')

############### test de pressure sheet for inner joint 

# F_pressure = np.zeros(len(current_t))

# for i in range(len(current_t)):
    
#     F_pressure[i] = u_f.pressure_sheet(p[i], v[i])

# # Figure creation
# fig = plt.figure(num='Example of plot')
# axis = fig.gca()
# #axis.plot(current_t,F_pressure, label='F_pressure')
# #axis.plot(current_t,F_final,'--',color='black',label='F_final')
# #axis.set_title('Pressure sheet force')
# axis.scatter(current_t,F_pressure-F_final,s=1,label='knee torque')
# axis.set_title('diff pressure sheet force')

# plt.show()






