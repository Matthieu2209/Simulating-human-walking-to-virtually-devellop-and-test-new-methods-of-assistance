#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 10:28:48 2023

@author: matthieuxaussems
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import useful_functions as u_f
import Muscle_actuation_layer as muscle
from scipy.interpolate import interp1d

#muscle_test = pd.read_excel ("/Users/matthieuxaussems/Documents/MBProjects/Fullmodel_innerjoint/User_function/muscle_test.xlsx")
muscle_test = pd.read_excel ("/Users/matthieuxaussems/Documents/MBProjects/Fullmodel_innerjoint/User_function/muscle_total.xlsx")


current_t = muscle_test["time"].to_numpy()


#Pour test spécifique à VAS
# A_VAS = muscle_test["A_VAS"].to_numpy()
# Stim_VAS = muscle_test["Stim_VAS"].to_numpy()
# L_mtu_VAS = muscle_test["L_mtu_VAS"].to_numpy()
# Lce_VAS = muscle_test["Lce_VAS"].to_numpy()
# Fm_VAS = muscle_test["Fm_VAS"].to_numpy()
# kneeL_q = muscle_test["kneeL_q"].to_numpy()
# Torque_knee_VAS = muscle_test["Torque_knee_VAS"].to_numpy()

# Pour test avec tous les muscles et sur toutes les articulations

kneeL_q = muscle_test["kneeL_q"].to_numpy()
kneeL_qd = muscle_test["kneeL_qd"].to_numpy()
hipL_q = muscle_test["hipL_q"].to_numpy()
hipL_qd = muscle_test["hipL_qd"].to_numpy()
ankleL_q = muscle_test["ankleL_q"].to_numpy()
ankleL_qd = muscle_test["ankleL_qd"].to_numpy()

# Jointlimits_ankle = muscle_test["Jointlimits_ankle"].to_numpy()
# Jointlimits_hip = muscle_test["Jointlimits_hip"].to_numpy()
# jointlimits_knee = muscle_test["jointlimits_knee"].to_numpy()

#Torque_ankle_TA = muscle_test["Torque_ankle_TA"].to_numpy()
Stim_TA = muscle_test["Stim_TA"].to_numpy()
#Torque_ankle_GAS = muscle_test["Torque_ankle_GAS"].to_numpy()
#Torque_knee_GAS = muscle_test["Torque_knee_GAS"].to_numpy()
Stim_GAS = muscle_test["Stim_GAS"].to_numpy()
# Torque_ankle_SOL = muscle_test["Torque_ankle_SOL"].to_numpy()
Stim_SOL = muscle_test["Stim_SOL"].to_numpy()
# Torque_knee_VAS = muscle_test["Torque_knee_VAS"].to_numpy()
Stim_VAS = muscle_test["Stim_VAS"].to_numpy()
#Torque_knee_HAM = muscle_test["Torque_knee_HAM"].to_numpy()
#Torque_hip_HAM = muscle_test["Torque_hip_HAM"].to_numpy()
Stim_HAM = muscle_test["Stim_HAM"].to_numpy()
#Torque_hip_GLU = muscle_test["Torque_hip_GLU"].to_numpy()
Stim_GLU = muscle_test["Stim_GLU"].to_numpy()
#Torque_hip_HFL = muscle_test["Torque_hip_HFL"].to_numpy()
Stim_HFL = muscle_test["Stim_HFL"].to_numpy()

Torque_ankle_total = muscle_test["Torque_ankle_total"].to_numpy()
Torque_knee_total = muscle_test["Torque_knee_total"].to_numpy()
Torque_hip_total = muscle_test["Torque_hip_total"].to_numpy()



# données
VAS =0
SOL =1
GAS = 2
TA = 3
HAM = 4
GLU = 5
HFL = 6

ankle=0
knee=1
hip=2

### Test de vce compute et integration

#n=len(current_t)
n=20000




lce_test = np.zeros(n)
Fm_test = np.zeros(n)
lce_control = np.zeros(n)
Fm_control = np.zeros(n)
t=np.zeros(n)


# for i in range(n):
    
#     print(i)
#     t[i] = current_t[i]
#     lce_control[i] = Lce_VAS[i]
#     Fm_control[i] = Fm_VAS[i]
    
#     if i == 0 :
#         x_0 = L_mtu_VAS[0] - muscle.l_slack_muscle[0]
#         act_memory = np.array(([A_VAS[0],current_t[0]]))
#         lmtc_memory =np.array(([L_mtu_VAS[0],current_t[0]]))
        
#         lce_test[i] = x_0
#         Fm_test[i] = muscle.vce_compute(lce_test[0], L_mtu_VAS[0], A_VAS[0], 0)[1]
#     else:
#         act_memory = np.vstack([act_memory,[A_VAS[i],current_t[i]]])
#         lmtc_memory =np.vstack([lmtc_memory,[L_mtu_VAS[i],current_t[i]]])
#         lce_test[i] = muscle.integrateur2000(lce_test[i-1], 0.00005, 5, lmtc_memory[:,0], act_memory[:,0], current_t[i-1], 0,lmtc_memory[:,1])
#         Fm_test[i] = muscle.vce_compute(lce_test[i], L_mtu_VAS[i], A_VAS[i], 0)[1]


# error_fm  = Fm_test-Fm_control
# error_lce = lce_test - lce_control

# # f_m
# plt.scatter(t, Fm_test,s=1)
# plt.grid()
# plt.title("f_m_test euler N_iteration=5 dt = 0.00005")
# #plt.savefig("f_m_test euler N_iteration=5 dt = 0.00005",dpi=300)
# plt.show()

# #l_ce
# plt.scatter(t, lce_test,s=1)
# plt.grid()
# plt.title("l_ce_test euler N_iteration=5 dt = 0.00005")
# #plt.savefig("l_ce_test euler N_iteration=5 dt = 0.00005",dpi=300)
# plt.show()

# #erreur_fm
# plt.scatter(t, error_fm,s=1)
# plt.grid()
# plt.title("error_fm euler N_iteration=5 dt = 0.00005")
# #plt.savefig("error_fm euler N_iteration=5 dt = 0.00005",dpi=300)
# plt.show()

# #erreur_lce
# plt.scatter(t, error_lce,s=1)
# plt.grid()
# plt.title("error_lce euler N_iteration=5 dt = 0.00005")
# #plt.savefig("error_lce euler N_iteration=5 dt = 0.00005",dpi=300)
# plt.show()

### Test du filtre qui permet de passer de stimulation à activation

A_test = np.zeros(n)
A_control = np.zeros(n)
Stim_control = np.zeros(n)
A_prec = 0
tau =0.01
t=np.zeros(n)

# for i in range(n):
#     print(i)
#     t[i]=current_t[i]
#     A_control[i]=A_VAS[i]
#     Stim_control[i] = Stim_VAS[i]
    
#     if i ==0 : 
#         A_test[i] = u_f.low_filter(Stim_VAS[i],tau,0,0)
#     else :
#         A_test[i] = u_f.low_filter(Stim_VAS[i],tau,0.00005,A_prec)
    
#     A_prec = A_test[i]

# # A
# plt.scatter(t, A_test,s=1)
# #plt.scatter(t, Stim_control,s=1)
# plt.grid()
# plt.title("Activation low filter")
# plt.show()

# #erreur_A
# error_A = A_test - A_control 
# plt.scatter(t, error_A,s=1)
# plt.grid()
# plt.title("error activation low filter")
# plt.show()

### Test lmtu_update

lmtu_test = np.zeros(n)
lmtu_control = np.zeros(n)

# for i in range(n):
#     print(i)
#     t[i] = current_t[i]
#     lmtu_control[i] = L_mtu_VAS[i]
#     lmtu_test[i] = muscle.lmtu_updateVAS(kneeL_q[i], 0)

# # lmtu
# plt.scatter(t, lmtu_test,s=1)
# plt.grid()
# plt.title("lmtu VAS")
# plt.show()

# #erreur_lmtu
# error_lmtu = lmtu_test - lmtu_control 
# plt.scatter(t, error_lmtu,s=1)
# plt.grid()
# plt.title("error lmtu VAS")
# plt.show()

### Test torque update

torqueVAS_test = np.zeros(n)
torqueVAS_control = np.zeros(n)

# for i in range(n):
#     print(i)
#     t[i]=current_t[i]
#     torqueVAS_control[i] = Torque_knee_VAS[i]
#     torqueVAS_test[i] = muscle.torque_updateKNEE(kneeL_q[i], Fm_VAS[i], 0)

# # Torque update 
# plt.scatter(t, torqueVAS_test,s=1)
# plt.grid()
# plt.title("torque update knee VAS")
# plt.show()

# #erreur_torque update
# error_lmtu = torqueVAS_test - torqueVAS_control 
# plt.scatter(t, error_lmtu,s=1)
# plt.grid()
# plt.title("error torque update knee VAS")
# plt.show()

### Test de la mise ensemble

# for i in range(n):
    
#         print(i)
#         t[i] = current_t[i]
#         lce_control[i] = Lce_VAS[i]
#         Fm_control[i] = Fm_VAS[i]
#         torqueVAS_control[i] = Torque_knee_VAS[i]
        
#         lmtu_test[i] = muscle.lmtu_updateVAS(kneeL_q[i], 0)
        
#         if i == 0 :
#             x_0 = lmtu_test[i] - muscle.l_slack_muscle[0]
#             A_test[i] = u_f.low_filter(Stim_VAS[i],tau,0,0)
#             act_memory = np.array(([A_test[i],current_t[i]]))
#             lmtc_memory =np.array(([lmtu_test[i],current_t[i]]))
            
#             lce_test[i] = x_0
#             Fm_test[i] = muscle.vce_compute(lce_test[0], lmtu_test[i], A_test[i], 0)[1]
#         else:
#             A_test[i] = u_f.low_filter(Stim_VAS[i],tau,0.00005,A_prec)
#             act_memory = np.vstack([act_memory,[A_test[i],current_t[i]]])
#             lmtc_memory =np.vstack([lmtc_memory,[lmtu_test[i],current_t[i]]])
#             lce_test[i] = muscle.integrateur2000(lce_test[i-1], 0.00005, 5, lmtc_memory[:,0], act_memory[:,0], current_t[i-1], 0,lmtc_memory[:,1])
#             Fm_test[i] = muscle.vce_compute(lce_test[i], lmtu_test[i], A_test[i], 0)[1]
        
#         A_prec = A_test[i]
#         torqueVAS_test[i] = muscle.torque_updateKNEE(kneeL_q[i], Fm_test[i], 0)

# error_fm  = Fm_test-Fm_control
# error_lce = lce_test - lce_control
# error_torqueVAS = torqueVAS_test - torqueVAS_control

# # f_m
# plt.scatter(t, Fm_test,s=1)
# plt.grid()
# plt.title("f_m_test euler N_iteration=5 dt = 0.00005")
# #plt.savefig("f_m_test euler N_iteration=5 dt = 0.00005",dpi=300)
# plt.show()

# #l_ce
# plt.scatter(t, lce_test,s=1)
# plt.grid()
# plt.title("l_ce_test euler N_iteration=5 dt = 0.00005")
# #plt.savefig("l_ce_test euler N_iteration=5 dt = 0.00005",dpi=300)
# plt.show()

# #Torque update VAS
# plt.scatter(t, torqueVAS_test,s=1)
# plt.grid()
# plt.title("torqueVAS knee_test")
# #plt.savefig("l_ce_test euler N_iteration=5 dt = 0.00005",dpi=300)
# plt.show()

# #erreur_fm
# plt.scatter(t, error_fm,s=1)
# plt.grid()
# plt.title("error_fm euler N_iteration=5 dt = 0.00005")
# #plt.savefig("error_fm euler N_iteration=5 dt = 0.00005",dpi=300)
# plt.show()

# #erreur_lce
# plt.scatter(t, error_lce,s=1)
# plt.grid()
# plt.title("error_lce euler N_iteration=5 dt = 0.00005")
# #plt.savefig("error_lce euler N_iteration=5 dt = 0.00005",dpi=300)
# plt.show()

# #erreur Torque update VAS
# plt.scatter(t, error_torqueVAS,s=1)
# plt.grid()
# plt.title("error torqueVAS knee_test")
# #plt.savefig("l_ce_test euler N_iteration=5 dt = 0.00005",dpi=300)
# plt.show()

####### Test d'obtentions des torques sur articulations à partir des angles des articulations et des stimulations

jl_ankle = np.zeros(n)
jl_knee = np.zeros(n)
jl_hip = np.zeros(n)
# jl_ankle_control = np.zeros(n)
# jl_knee_control = np.zeros(n)
# jl_hip_control = np.zeros(n)
#Torque_ankle_TA_test=np.zeros(n)
#Torque_ankle_TA_control = np.zeros(n)
# Torque_ankle_GAS_test=np.zeros(n)
# Torque_ankle_GAS_control=np.zeros(n)
# Torque_knee_GAS_test=np.zeros(n)
# Torque_knee_GAS_control=np.zeros(n)
# Torque_ankle_SOL_test=np.zeros(n)
# Torque_ankle_SOL_control=np.zeros(n)
# Torque_knee_VAS_test=np.zeros(n)
# Torque_knee_VAS_control=np.zeros(n)
#Torque_knee_HAM_test=np.zeros(n)
#Torque_hip_HAM_test=np.zeros(n)
#Torque_knee_HAM_control=np.zeros(n)
#Torque_hip_HAM_control=np.zeros(n)
#Torque_hip_GLU_test=np.zeros(n)
#Torque_hip_GLU_control=np.zeros(n)
#Torque_hip_HFL_test=np.zeros(n)
#Torque_hip_HFL_control=np.zeros(n)

Torque_ankle=np.zeros(n)
Torque_ankle_control=np.zeros(n)
Torque_knee=np.zeros(n)
Torque_knee_control=np.zeros(n)
Torque_hip=np.zeros(n)
Torque_hip_control=np.zeros(n)

A_prec_TA=0
A_prec_GAS=0
A_prec_SOL=0
A_prec_VAS=0
A_prec_HAM=0
A_prec_GLU=0
A_prec_HFL=0

lce_prec_TA=0
lce_prec_GAS=0
lce_prec_SOL=0
lce_prec_VAS=0
lce_prec_HAM=0
lce_prec_GLU=0
lce_prec_HFL=0


t=np.zeros(n)

for i in range(n):
    
    print(i)
    
    t[i] = current_t[i]
    # jl_ankle_control[i] = Jointlimits_ankle[i]
    # jl_knee_control[i] = jointlimits_knee[i]
    # jl_hip_control[i] = Jointlimits_hip[i]
    # Torque_ankle_TA_control[i] = Torque_ankle_TA[i]
    # Torque_ankle_GAS_control[i] = Torque_ankle_GAS[i]
    #Torque_knee_GAS_control[i]=Torque_knee_GAS[i]
    #Torque_ankle_SOL_control[i]=Torque_ankle_SOL[i]
    #Torque_knee_VAS_control[i]=Torque_knee_VAS[i]
    #Torque_knee_HAM_control[i]=Torque_knee_HAM[i]
    #Torque_hip_HAM_control[i]=Torque_hip_HAM[i]
    #Torque_hip_GLU_control[i]=Torque_hip_GLU[i]
    #Torque_hip_HFL_control[i]=Torque_hip_HFL[i]
    Torque_ankle_control[i] = Torque_ankle_total[i]
    Torque_knee_control[i] = Torque_knee_total[i]
    Torque_hip_control[i] = Torque_hip_total[i]
    
    lmtu_TA = muscle.lmtu_updateTA(ankleL_q[i], TA)
    lmtu_GAS = muscle.lmtu_updateGAS(ankleL_q[i],kneeL_q[i],GAS)
    lmtu_SOL = muscle.lmtu_updateSOL(ankleL_q[i], SOL)
    lmtu_VAS = muscle.lmtu_updateVAS(kneeL_q[i], VAS)
    lmtu_HAM = muscle.lmtu_updateHAM(kneeL_q[i], hipL_q[i], HAM)
    lmtu_GLU = muscle.lmtu_updateGLU(hipL_q[i], GLU)
    lmtu_HFL = muscle.lmtu_updateHFL(hipL_q[i], HFL)
    
    if i == 0 :
        #TA
        x_0_TA = lmtu_TA - muscle.l_slack_muscle[TA]
        A_TA = u_f.low_filter(Stim_TA[i],tau,0,0)
        act_memory_TA = np.array(([A_TA,current_t[i]]))
        lmtc_memory_TA =np.array(([lmtu_TA,current_t[i]]))
        lce_TA = x_0_TA
        Fm_TA = muscle.vce_compute(lce_TA, lmtu_TA, A_TA, TA)[1]
        
        #GAS
        x_0_GAS = lmtu_GAS - muscle.l_slack_muscle[GAS]
        A_GAS = u_f.low_filter(Stim_GAS[i],tau,0,0)
        act_memory_GAS = np.array(([A_GAS,current_t[i]]))
        lmtc_memory_GAS =np.array(([lmtu_GAS,current_t[i]]))
        lce_GAS = x_0_GAS
        Fm_GAS = muscle.vce_compute(lce_GAS, lmtu_GAS, A_GAS, GAS)[1]
        
        #SOL
        x_0_SOL = lmtu_SOL - muscle.l_slack_muscle[SOL]
        A_SOL = u_f.low_filter(Stim_SOL[i],tau,0,0)
        act_memory_SOL = np.array(([A_SOL,current_t[i]]))
        lmtc_memory_SOL =np.array(([lmtu_SOL,current_t[i]]))
        lce_SOL = x_0_SOL
        Fm_SOL = muscle.vce_compute(lce_SOL, lmtu_SOL, A_SOL, SOL)[1]
        
        #VAS
        x_0_VAS = lmtu_VAS - muscle.l_slack_muscle[VAS]
        A_VAS = u_f.low_filter(Stim_VAS[i],tau,0,0)
        act_memory_VAS = np.array(([A_VAS,current_t[i]]))
        lmtc_memory_VAS =np.array(([lmtu_VAS,current_t[i]]))
        lce_VAS = x_0_VAS
        Fm_VAS = muscle.vce_compute(lce_VAS, lmtu_VAS, A_VAS, VAS)[1]
        
        #HAM
        x_0_HAM = lmtu_HAM - muscle.l_slack_muscle[HAM]
        A_HAM = u_f.low_filter(Stim_HAM[i],tau,0,0)
        act_memory_HAM = np.array(([A_HAM,current_t[i]]))
        lmtc_memory_HAM =np.array(([lmtu_HAM,current_t[i]]))
        lce_HAM = x_0_HAM
        Fm_HAM = muscle.vce_compute(lce_HAM, lmtu_HAM, A_HAM, HAM)[1]
        
        #GLU
        x_0_GLU = lmtu_GLU - muscle.l_slack_muscle[GLU]
        A_GLU = u_f.low_filter(Stim_GLU[i],tau,0,0)
        act_memory_GLU = np.array(([A_GLU,current_t[i]]))
        lmtc_memory_GLU =np.array(([lmtu_GLU,current_t[i]]))
        lce_GLU = x_0_GLU
        Fm_GLU = muscle.vce_compute(lce_GLU, lmtu_GLU, A_GLU, GLU)[1]
        
        #HFL
        x_0_HFL = lmtu_HFL - muscle.l_slack_muscle[HFL]
        A_HFL = u_f.low_filter(Stim_HFL[i],tau,0,0)
        act_memory_HFL = np.array(([A_HFL,current_t[i]]))
        lmtc_memory_HFL =np.array(([lmtu_HFL,current_t[i]]))
        lce_HFL = x_0_HFL
        Fm_HFL = muscle.vce_compute(lce_HFL, lmtu_HFL, A_HFL, HFL)[1]
           
        
    else:
        #TA
        A_TA = u_f.low_filter(Stim_TA[i],tau,0.00005,A_prec_TA)
        act_memory_TA = np.vstack([act_memory_TA,[A_TA,current_t[i]]])
        lmtc_memory_TA =np.vstack([lmtc_memory_TA,[lmtu_TA,current_t[i]]])
        lce_TA = muscle.integrateur2000(lce_prec_TA, 0.00005, 5, lmtc_memory_TA[:,0], act_memory_TA[:,0], current_t[i-1], TA,lmtc_memory_TA[:,1])
        Fm_TA = muscle.vce_compute(lce_TA, lmtu_TA, A_TA, TA)[1]
        
        #GAS
        A_GAS = u_f.low_filter(Stim_GAS[i],tau,0.00005,A_prec_GAS)
        act_memory_GAS = np.vstack([act_memory_GAS,[A_GAS,current_t[i]]])
        lmtc_memory_GAS =np.vstack([lmtc_memory_GAS,[lmtu_GAS,current_t[i]]])
        lce_GAS = muscle.integrateur2000(lce_prec_GAS, 0.00005, 5, lmtc_memory_GAS[:,0], act_memory_GAS[:,0], current_t[i-1], GAS,lmtc_memory_GAS[:,1])
        Fm_GAS = muscle.vce_compute(lce_GAS, lmtu_GAS, A_GAS, GAS)[1]
        
        #SOL
        A_SOL = u_f.low_filter(Stim_SOL[i],tau,0.00005,A_prec_SOL)
        act_memory_SOL = np.vstack([act_memory_SOL,[A_SOL,current_t[i]]])
        lmtc_memory_SOL =np.vstack([lmtc_memory_SOL,[lmtu_SOL,current_t[i]]])
        lce_SOL = muscle.integrateur2000(lce_prec_SOL, 0.00005, 5, lmtc_memory_SOL[:,0], act_memory_SOL[:,0], current_t[i-1], SOL,lmtc_memory_SOL[:,1])
        Fm_SOL = muscle.vce_compute(lce_SOL, lmtu_SOL, A_SOL, SOL)[1]
        
        #VAS
        A_VAS = u_f.low_filter(Stim_VAS[i],tau,0.00005,A_prec_VAS)
        act_memory_VAS = np.vstack([act_memory_VAS,[A_VAS,current_t[i]]])
        lmtc_memory_VAS =np.vstack([lmtc_memory_VAS,[lmtu_VAS,current_t[i]]])
        lce_VAS = muscle.integrateur2000(lce_prec_VAS, 0.00005, 5, lmtc_memory_VAS[:,0], act_memory_VAS[:,0], current_t[i-1], VAS,lmtc_memory_VAS[:,1])
        Fm_VAS = muscle.vce_compute(lce_VAS, lmtu_VAS, A_VAS, VAS)[1]
        
        #HAM
        A_HAM = u_f.low_filter(Stim_HAM[i],tau,0.00005,A_prec_HAM)
        act_memory_HAM = np.vstack([act_memory_HAM,[A_HAM,current_t[i]]])
        lmtc_memory_HAM =np.vstack([lmtc_memory_HAM,[lmtu_HAM,current_t[i]]])
        lce_HAM = muscle.integrateur2000(lce_prec_HAM, 0.00005, 5, lmtc_memory_HAM[:,0], act_memory_HAM[:,0], current_t[i-1], HAM,lmtc_memory_HAM[:,1])
        Fm_HAM = muscle.vce_compute(lce_HAM, lmtu_HAM, A_HAM, HAM)[1]
        
        #GLU
        A_GLU = u_f.low_filter(Stim_GLU[i],tau,0.00005,A_prec_GLU)
        act_memory_GLU = np.vstack([act_memory_GLU,[A_GLU,current_t[i]]])
        lmtc_memory_GLU =np.vstack([lmtc_memory_GLU,[lmtu_GLU,current_t[i]]])
        lce_GLU = muscle.integrateur2000(lce_prec_GLU, 0.00005, 5, lmtc_memory_GLU[:,0], act_memory_GLU[:,0], current_t[i-1], GLU,lmtc_memory_GLU[:,1])
        Fm_GLU = muscle.vce_compute(lce_GLU, lmtu_GLU, A_GLU, GLU)[1]
        
        #HFL
        A_HFL = u_f.low_filter(Stim_HFL[i],tau,0.00005,A_prec_HFL)
        act_memory_HFL = np.vstack([act_memory_HFL,[A_HFL,current_t[i]]])
        lmtc_memory_HFL =np.vstack([lmtc_memory_HFL,[lmtu_HFL,current_t[i]]])
        lce_HFL = muscle.integrateur2000(lce_prec_HFL, 0.00005, 5, lmtc_memory_HFL[:,0], act_memory_HFL[:,0], current_t[i-1], HFL,lmtc_memory_HFL[:,1])
        Fm_HFL = muscle.vce_compute(lce_HFL, lmtu_HFL, A_HFL, HFL)[1]
        
    
    #TA
    A_prec_TA = A_TA
    lce_prec_TA = lce_TA
    Torque_ankle_TA = muscle.torque_updateANKLE(ankleL_q[i],Fm_TA,TA)
    
    #GAS
    A_prec_GAS = A_GAS
    lce_prec_GAS = lce_GAS
    Torque_ankle_GAS = muscle.torque_updateANKLE(ankleL_q[i],Fm_GAS,GAS)
    Torque_knee_GAS=muscle.torque_updateKNEE(kneeL_q[i], Fm_GAS, GAS)
    
    #SOL
    A_prec_SOL = A_SOL
    lce_prec_SOL = lce_SOL
    Torque_ankle_SOL= muscle.torque_updateANKLE(ankleL_q[i],Fm_SOL,SOL)
    
    #VAS
    A_prec_VAS = A_VAS
    lce_prec_VAS = lce_VAS
    Torque_knee_VAS= muscle.torque_updateKNEE(kneeL_q[i],Fm_VAS,VAS)
    
    #HAM
    A_prec_HAM = A_HAM
    lce_prec_HAM = lce_HAM
    Torque_knee_HAM= muscle.torque_updateKNEE(kneeL_q[i],Fm_HAM,HAM)
    Torque_hip_HAM=muscle.torque_updateHIP(hipL_q[i], Fm_HAM, HAM)
    
    #GLU
    A_prec_GLU = A_GLU
    lce_prec_GLU = lce_GLU
    Torque_hip_GLU=muscle.torque_updateHIP(hipL_q[i], Fm_GLU, GLU)
    
    #HFL 
    A_prec_HFL = A_HFL
    lce_prec_HFL = lce_HFL
    Torque_hip_HFL=muscle.torque_updateHIP(hipL_q[i], Fm_HFL, HFL)
    
    
    # joint limits 
    jl_ankle[i]=muscle.joint_limits(ankle, ankleL_q[i], ankleL_qd[i])
    jl_knee[i]=muscle.joint_limits(knee, kneeL_q[i], kneeL_qd[i])
    jl_hip[i]=muscle.joint_limits(hip, hipL_q[i], hipL_qd[i])
    
    Torque_ankle[i] = jl_ankle[i] + Torque_ankle_GAS + Torque_ankle_SOL - Torque_ankle_TA
    Torque_knee[i] = jl_knee[i] + Torque_knee_VAS - Torque_knee_GAS - Torque_knee_HAM
    Torque_hip[i] = jl_hip[i] + Torque_hip_GLU + Torque_hip_HAM - Torque_hip_HFL



plt.scatter(t, Torque_hip-Torque_hip_control,s=1)
plt.grid()
plt.title("diff  hip torque")
plt.show()

plt.scatter(t, Torque_ankle-Torque_ankle_control,s=1)
plt.grid()
plt.title("diff ankle torque")
plt.show()

plt.scatter(t, Torque_ankle,s=1)
plt.grid()
plt.title("diff ankle torque")
plt.show()

plt.scatter(t, Torque_knee-Torque_knee_control,s=1)
plt.grid()
plt.title("diff knee torque")
plt.show()





