#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 10:46:13 2023

@author: matthieuxaussems
"""

#TEST muscle and stimulations ensemble
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import useful_functions as u_f
import Muscle_actuation_layer as muscle
import Neural_control_layer as neural


#### Valeurs exportées

muscle_stimulation_test = pd.read_excel ("/Users/matthieuxaussems/Documents/MBProjects/Fullmodel_innerjoint/User_function/muscle_stimulation_test2.xlsx")

# les données que on a vraiment

current_t = muscle_stimulation_test["time"].to_numpy()

Pz_heelL =  - muscle_stimulation_test["Pz_heelL"].to_numpy()
Pz_ballL =  - muscle_stimulation_test["Pz_ballL"].to_numpy()
Pz_heelR =  - muscle_stimulation_test["Pz_heelR"].to_numpy()
Pz_ballR =  - muscle_stimulation_test["Pz_ballR"].to_numpy()

Ldx_nofilter = muscle_stimulation_test["Ldx_nofilter"].to_numpy()
Rdx_nofilter = muscle_stimulation_test["Rdx_nofilter"].to_numpy()

kneeL_q = muscle_stimulation_test["knee_q"].to_numpy()
kneeL_qd = muscle_stimulation_test["knee_qd"].to_numpy()
hipL_q = muscle_stimulation_test["hip_q"].to_numpy()
hipL_qd = muscle_stimulation_test["hip_qd"].to_numpy()
ankleL_q = muscle_stimulation_test["ankle_q"].to_numpy()
ankleL_qd = muscle_stimulation_test["ankle_qd"].to_numpy()

trunk_theta_nodelay = muscle_stimulation_test["trunk_theta_nodelay"].to_numpy()
trunk_dtheta_nodelay = muscle_stimulation_test["trunk_dtheta_nodelay"].to_numpy()

# LVAS_fm = muscle_stimulation_test["LVAS_fm"].to_numpy()
# LSOL_fm = muscle_stimulation_test["LSOL_fm"].to_numpy()
# LTA_fm = muscle_stimulation_test["LTA_fm"].to_numpy()
# LHAM_fm = muscle_stimulation_test["LHAM_fm"].to_numpy()
# LGAS_fm = muscle_stimulation_test["LGAS_fm"].to_numpy()
# LGLU_fm = muscle_stimulation_test["LGLU_fm"].to_numpy()
# LHFL_fm = muscle_stimulation_test["LHFL_fm"].to_numpy()

# LVAS_lce = muscle_stimulation_test["LVAS_lce"].to_numpy()
# LSOL_lce = muscle_stimulation_test["LSOL_lce"].to_numpy()
# LTA_lce = muscle_stimulation_test["LTA_lce"].to_numpy()
# LHAM_lce = muscle_stimulation_test["LHAM_lce"].to_numpy()
# LGAS_lce = muscle_stimulation_test["LGAS_lce"].to_numpy()
# LGLU_lce = muscle_stimulation_test["LGLU_lce"].to_numpy()
# LHFL_lce = muscle_stimulation_test["LHFL_lce"].to_numpy()


# Ldx_test = muscle_stimulation_test["Ldx_filtered"].to_numpy()
# Rdx_test = muscle_stimulation_test["Rdx_filtered"].to_numpy()

# StanceL_test = muscle_stimulation_test["L_stance_Nodelay"].to_numpy()
# StanceR_test = muscle_stimulation_test["R_Stance_nodelay"].to_numpy()

# Stim_TA = muscle_stimulation_test["Stim_TA"].to_numpy()
# Stim_GAS = muscle_stimulation_test["Stim_GAS"].to_numpy()
# Stim_SOL = muscle_stimulation_test["Stim_SOL"].to_numpy()
# Stim_VAS = muscle_stimulation_test["Stim_VAS"].to_numpy()
# Stim_HAM = muscle_stimulation_test["Stim_HAM"].to_numpy()
# Stim_GLU = muscle_stimulation_test["Stim_GLU"].to_numpy()
# Stim_HFL = muscle_stimulation_test["Stim_HFL"].to_numpy()

Torque_ankle_total = muscle_stimulation_test["Torque_ankle_total"].to_numpy()
Torque_knee_total = muscle_stimulation_test["Torque_knee_total"].to_numpy()
Torque_hip_total = muscle_stimulation_test["Torque_hip_total"].to_numpy()

#n=len(current_t)
n=20000
#n=50000

#### données
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

ground_limit = 0

n_muscle = 7

tau =0.01

dt=0.00005 #à calculer avec tsim dans joint forces 


#VARIABLE test

t= np.zeros(n)
Ldx_control = np.zeros(n)
Rdx_control = np.zeros(n)
Ldx_nf_control = np.zeros(n)
StanceL_control = np.zeros(n)
StanceR_control = np.zeros(n)
Stim_VAS_control = np.zeros(n)
Stim_SOL_control = np.zeros(n)
Stim_GAS_control = np.zeros(n)
Stim_TA_control = np.zeros(n)
Stim_HAM_control = np.zeros(n)
Stim_GLU_control = np.zeros(n)
Stim_HFL_control = np.zeros(n)

Torque_ankle=np.zeros(n)
Torque_ankle_control=np.zeros(n)
Torque_knee=np.zeros(n)
Torque_knee_control=np.zeros(n)
Torque_hip=np.zeros(n)
Torque_hip_control=np.zeros(n)


# test1 = np.zeros(n)
# test2 = np.zeros(n)
# test3 = np.zeros(n)
# test4 = np.zeros(n)


##### Initialisation des mémoires et des paramètres globaux

#Ldx/Rdx  / theta-dtheta trunk / stance

LDx_memory = np.array([0,0])
RDX_memory = np.array([0,0])

theta_trunk_memory = np.array([0, 0])
dtheta_trunk_memory = np.array([0, 0])

Stance_memory = np.array([0,0,0])
Stance_memory_delayed = np.zeros((round((0.01/dt)),3))

# Fm, lce

f_m_memory = np.array([0, 0, 0, 0, 0, 0, 0, 0])
lce_memory = np.array([0, 0, 0, 0, 0, 0, 0, 0])

#Variable global

Ldx_prec = 0
Rdx_prec = 0

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


# Stimulations

StimL = np.zeros((n,7))

for i in range(n):
    
    print(i)
    t[i] = current_t[i]
    # Ldx_control[i] = Ldx_test[i]
    # Rdx_control[i] = Rdx_test[i]
    # Ldx_nf_control[i] = Ldx_nofilter[i]
    # StanceL_control[i] = StanceL_test[i]
    # StanceR_control[i] = StanceR_test[i]
    # Stim_VAS_control[i] = Stim_VAS[i]
    # Stim_SOL_control[i] = Stim_SOL[i]
    # Stim_GAS_control[i] = Stim_GAS[i]
    # Stim_TA_control[i] = Stim_TA[i]
    # Stim_HAM_control[i] = Stim_HAM[i]
    # Stim_GLU_control[i] = Stim_GLU[i]
    # Stim_HFL_control[i] = Stim_HFL[i]
    Torque_ankle_control[i] = Torque_ankle_total[i]
    Torque_knee_control[i] = Torque_knee_total[i]
    Torque_hip_control[i] = Torque_hip_total[i]
    
    #global Ldx_prec
    #global Rdx_prec
    #global LDx_memory
    #global RDX_memory
    #global theta_trunk_memory
    #global dtheta_trunk_memory
    #global Stance_memory
    #global Stance_memory_delayed
    #global f_m_memory
    #global lce_memory
    #global A_prec_TA
    # global A_prec_GAS
    # global A_prec_SOL
    # global A_prec_VAS
    # global A_prec_HAM
    # global A_prec_GLU
    # global A_prec_HFL

    # global lce_prec_TA
    # global lce_prec_GAS
    # global lce_prec_SOL
    # global lce_prec_VAS
    # global lce_prec_HAM
    # global lce_prec_GLU
    # global lce_prec_HFL
    
    ### obtention StanceL
    
    # control of contact for BallL,HeelL,BallR,HeelR
    
    BallL_cnt,HeelL_cnt,BallR_cnt,HeelR_cnt = u_f.contact_cnt(Pz_ballL[i], Pz_heelL[i], Pz_ballR[i], Pz_heelR[i],ground_limit)
    
    # Stance or not ?
    
    StanceL,StanceR = u_f.Stance_cnt(BallL_cnt,HeelL_cnt,BallR_cnt,HeelR_cnt)
    
    # test3[i] = StanceL
    # test4[i] = StanceR
    
    ### obtention Ldx/Rdx Filtered
    
    if i == 0 :
        Ldx = u_f.low_filter(Ldx_nofilter[i], 0.02, 0, Ldx_prec)
        Rdx = u_f.low_filter(Rdx_nofilter[i], 0.02, 0, Rdx_prec)
    else :
        Ldx = u_f.low_filter(Ldx_nofilter[i], 0.02, dt, Ldx_prec)
        Rdx = u_f.low_filter(Rdx_nofilter[i], 0.02, dt, Rdx_prec)
    
    Ldx_prec = Ldx
    Rdx_prec = Rdx
    
    # test1[i] = Ldx
    # test2[i] = Rdx
    
    ### update memory Stance, ldx/Rdx filtered, theta-dtheta trunk
    
    if i == 0 :
        
        LDx_memory = np.array([Ldx,current_t[0]])
        RDX_memory = np.array([Rdx,current_t[0]])
        theta_trunk_memory = np.array([trunk_theta_nodelay[0], current_t[0]])
        dtheta_trunk_memory = np.array([trunk_dtheta_nodelay[0], current_t[0]])
        Stance_memory = np.array([StanceL,StanceR, current_t[0]])
        f_m_memory = np.array([0, 0, 0, 0, 0, 0, 0, current_t[0]])
        lce_memory = np.array([0, 0, 0, 0, 0, 0, 0, current_t[0]])
        
        
    else :
        
        LDx_memory =  np.vstack([LDx_memory,[Ldx,current_t[i]]])
        RDX_memory = np.vstack([RDX_memory, [Rdx, current_t[i]]])
        theta_trunk_memory = np.vstack( [theta_trunk_memory, [trunk_theta_nodelay[i], current_t[i]]])
        dtheta_trunk_memory = np.vstack([dtheta_trunk_memory, [trunk_dtheta_nodelay[i], current_t[i]]])
        Stance_memory = np.vstack([Stance_memory, [StanceL, StanceR, current_t[i]]])
        f_m_memory = np.vstack([f_m_memory, [0, 0, 0, 0, 0, 0, 0, current_t[i]]])
        lce_memory = np.vstack([lce_memory, [0, 0, 0, 0, 0, 0, 0, current_t[i]]])
    
    if i >= round((0.01/dt)):
            
        Stance_memory_delayed = np.vstack([Stance_memory_delayed, [Stance_memory[i-round((0.01/dt)),0], Stance_memory[i-round((0.01/dt)),1], current_t[i]]])
    else: 
        Stance_memory_delayed[i,2] = current_t[i]
    
    
    ### obtention des stimulations
    
    StimL[i] = neural.Feedback(current_t[i], dt, Stance_memory_delayed,f_m_memory, LDx_memory, RDX_memory, lce_memory,kneeL_q[i],kneeL_qd[i],theta_trunk_memory, dtheta_trunk_memory,0)
    
    
    ### obtention des Fm,lce 
    
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
        A_TA = u_f.low_filter(StimL[i,TA],tau,0,0)
        act_memory_TA = np.array(([A_TA,current_t[i]]))
        lmtc_memory_TA =np.array(([lmtu_TA,current_t[i]]))
        lce_TA = x_0_TA
        Fm_TA = muscle.vce_compute(lce_TA, lmtu_TA, A_TA, TA)[1]
        
        #GAS
        x_0_GAS = lmtu_GAS - muscle.l_slack_muscle[GAS]
        A_GAS = u_f.low_filter(StimL[i,GAS],tau,0,0)
        act_memory_GAS = np.array(([A_GAS,current_t[i]]))
        lmtc_memory_GAS =np.array(([lmtu_GAS,current_t[i]]))
        lce_GAS = x_0_GAS
        Fm_GAS = muscle.vce_compute(lce_GAS, lmtu_GAS, A_GAS, GAS)[1]
        
        #SOL
        x_0_SOL = lmtu_SOL - muscle.l_slack_muscle[SOL]
        A_SOL = u_f.low_filter(StimL[i,SOL],tau,0,0)
        act_memory_SOL = np.array(([A_SOL,current_t[i]]))
        lmtc_memory_SOL =np.array(([lmtu_SOL,current_t[i]]))
        lce_SOL = x_0_SOL
        Fm_SOL = muscle.vce_compute(lce_SOL, lmtu_SOL, A_SOL, SOL)[1]
        
        #VAS
        x_0_VAS = lmtu_VAS - muscle.l_slack_muscle[VAS]
        A_VAS = u_f.low_filter(StimL[i,VAS],tau,0,0)
        act_memory_VAS = np.array(([A_VAS,current_t[i]]))
        lmtc_memory_VAS =np.array(([lmtu_VAS,current_t[i]]))
        lce_VAS = x_0_VAS
        Fm_VAS = muscle.vce_compute(lce_VAS, lmtu_VAS, A_VAS, VAS)[1]
        
        #HAM
        x_0_HAM = lmtu_HAM - muscle.l_slack_muscle[HAM]
        A_HAM = u_f.low_filter(StimL[i,HAM],tau,0,0)
        act_memory_HAM = np.array(([A_HAM,current_t[i]]))
        lmtc_memory_HAM =np.array(([lmtu_HAM,current_t[i]]))
        lce_HAM = x_0_HAM
        Fm_HAM = muscle.vce_compute(lce_HAM, lmtu_HAM, A_HAM, HAM)[1]
        
        #GLU
        x_0_GLU = lmtu_GLU - muscle.l_slack_muscle[GLU]
        A_GLU = u_f.low_filter(StimL[i,GLU],tau,0,0)
        act_memory_GLU = np.array(([A_GLU,current_t[i]]))
        lmtc_memory_GLU =np.array(([lmtu_GLU,current_t[i]]))
        lce_GLU = x_0_GLU
        Fm_GLU = muscle.vce_compute(lce_GLU, lmtu_GLU, A_GLU, GLU)[1]
        
        #HFL
        x_0_HFL = lmtu_HFL - muscle.l_slack_muscle[HFL]
        A_HFL = u_f.low_filter(StimL[i,HFL],tau,0,0)
        act_memory_HFL = np.array(([A_HFL,current_t[i]]))
        lmtc_memory_HFL =np.array(([lmtu_HFL,current_t[i]]))
        lce_HFL = x_0_HFL
        Fm_HFL = muscle.vce_compute(lce_HFL, lmtu_HFL, A_HFL, HFL)[1]
           
        
    else:
        #TA
        A_TA = u_f.low_filter(StimL[i,TA],tau,dt,A_prec_TA)
        act_memory_TA = np.vstack([act_memory_TA,[A_TA,current_t[i]]])
        lmtc_memory_TA =np.vstack([lmtc_memory_TA,[lmtu_TA,current_t[i]]])
        lce_TA = muscle.integrateur2000(lce_prec_TA, dt, 5, lmtc_memory_TA[:,0], act_memory_TA[:,0], current_t[i-1], TA,lmtc_memory_TA[:,1])
        Fm_TA = muscle.vce_compute(lce_TA, lmtu_TA, A_TA, TA)[1]
        
        #GAS
        A_GAS = u_f.low_filter(StimL[i,GAS],tau,dt,A_prec_GAS)
        act_memory_GAS = np.vstack([act_memory_GAS,[A_GAS,current_t[i]]])
        lmtc_memory_GAS =np.vstack([lmtc_memory_GAS,[lmtu_GAS,current_t[i]]])
        lce_GAS = muscle.integrateur2000(lce_prec_GAS, dt, 5, lmtc_memory_GAS[:,0], act_memory_GAS[:,0], current_t[i-1], GAS,lmtc_memory_GAS[:,1])
        Fm_GAS = muscle.vce_compute(lce_GAS, lmtu_GAS, A_GAS, GAS)[1]
        
        #SOL
        A_SOL = u_f.low_filter(StimL[i,SOL],tau,dt,A_prec_SOL)
        act_memory_SOL = np.vstack([act_memory_SOL,[A_SOL,current_t[i]]])
        lmtc_memory_SOL =np.vstack([lmtc_memory_SOL,[lmtu_SOL,current_t[i]]])
        lce_SOL = muscle.integrateur2000(lce_prec_SOL, dt, 5, lmtc_memory_SOL[:,0], act_memory_SOL[:,0], current_t[i-1], SOL,lmtc_memory_SOL[:,1])
        Fm_SOL = muscle.vce_compute(lce_SOL, lmtu_SOL, A_SOL, SOL)[1]
        
        #VAS
        A_VAS = u_f.low_filter(StimL[i,VAS],tau,dt,A_prec_VAS)
        act_memory_VAS = np.vstack([act_memory_VAS,[A_VAS,current_t[i]]])
        lmtc_memory_VAS =np.vstack([lmtc_memory_VAS,[lmtu_VAS,current_t[i]]])
        lce_VAS = muscle.integrateur2000(lce_prec_VAS, dt, 5, lmtc_memory_VAS[:,0], act_memory_VAS[:,0], current_t[i-1], VAS,lmtc_memory_VAS[:,1])
        Fm_VAS = muscle.vce_compute(lce_VAS, lmtu_VAS, A_VAS, VAS)[1]
        
        #HAM
        A_HAM = u_f.low_filter(StimL[i,HAM],tau,dt,A_prec_HAM)
        act_memory_HAM = np.vstack([act_memory_HAM,[A_HAM,current_t[i]]])
        lmtc_memory_HAM =np.vstack([lmtc_memory_HAM,[lmtu_HAM,current_t[i]]])
        lce_HAM = muscle.integrateur2000(lce_prec_HAM, dt, 5, lmtc_memory_HAM[:,0], act_memory_HAM[:,0], current_t[i-1], HAM,lmtc_memory_HAM[:,1])
        Fm_HAM = muscle.vce_compute(lce_HAM, lmtu_HAM, A_HAM, HAM)[1]
        
        #GLU
        A_GLU = u_f.low_filter(StimL[i,GLU],tau,dt,A_prec_GLU)
        act_memory_GLU = np.vstack([act_memory_GLU,[A_GLU,current_t[i]]])
        lmtc_memory_GLU =np.vstack([lmtc_memory_GLU,[lmtu_GLU,current_t[i]]])
        lce_GLU = muscle.integrateur2000(lce_prec_GLU, dt, 5, lmtc_memory_GLU[:,0], act_memory_GLU[:,0], current_t[i-1], GLU,lmtc_memory_GLU[:,1])
        Fm_GLU = muscle.vce_compute(lce_GLU, lmtu_GLU, A_GLU, GLU)[1]
        
        #HFL
        A_HFL = u_f.low_filter(StimL[i,HFL],tau,dt,A_prec_HFL)
        act_memory_HFL = np.vstack([act_memory_HFL,[A_HFL,current_t[i]]])
        lmtc_memory_HFL =np.vstack([lmtc_memory_HFL,[lmtu_HFL,current_t[i]]])
        lce_HFL = muscle.integrateur2000(lce_prec_HFL, dt, 5, lmtc_memory_HFL[:,0], act_memory_HFL[:,0], current_t[i-1], HFL,lmtc_memory_HFL[:,1])
        Fm_HFL = muscle.vce_compute(lce_HFL, lmtu_HFL, A_HFL, HFL)[1]
    
    
    
    ### update memory Fm, lce
    if i == 0 :
        f_m_memory[VAS]= Fm_VAS
        f_m_memory[SOL]= Fm_SOL
        f_m_memory[GAS]= Fm_GAS
        f_m_memory[TA]= Fm_TA
        f_m_memory[HAM]= Fm_HAM
        f_m_memory[GLU]= Fm_GLU
        f_m_memory[HFL]= Fm_HFL
        
        lce_memory[VAS]= lce_VAS
        lce_memory[SOL]= lce_SOL
        lce_memory[GAS]= lce_GAS
        lce_memory[TA]= lce_TA
        lce_memory[HAM]= lce_HAM
        lce_memory[GLU]= lce_GLU
        lce_memory[HFL]= lce_HFL
    else :
        
        f_m_memory[-1,VAS]= Fm_VAS
        f_m_memory[-1,SOL]= Fm_SOL
        f_m_memory[-1,GAS]= Fm_GAS
        f_m_memory[-1,TA]= Fm_TA
        f_m_memory[-1,HAM]= Fm_HAM
        f_m_memory[-1,GLU]= Fm_GLU
        f_m_memory[-1,HFL]= Fm_HFL
        
        lce_memory[-1,VAS]= lce_VAS
        lce_memory[-1,SOL]= lce_SOL
        lce_memory[-1,GAS]= lce_GAS
        lce_memory[-1,TA]= lce_TA
        lce_memory[-1,HAM]= lce_HAM
        lce_memory[-1,GLU]= lce_GLU
        lce_memory[-1,HFL]= lce_HFL
    
    
    ### calcul des torques
    
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
    jl_ankle=muscle.joint_limits(ankle, ankleL_q[i], ankleL_qd[i])
    jl_knee=muscle.joint_limits(knee, kneeL_q[i], kneeL_qd[i])
    jl_hip=muscle.joint_limits(hip, hipL_q[i], hipL_qd[i])
    
    Torque_ankle[i] = jl_ankle + Torque_ankle_GAS + Torque_ankle_SOL - Torque_ankle_TA
    Torque_knee[i] = jl_knee + Torque_knee_VAS - Torque_knee_GAS - Torque_knee_HAM
    Torque_hip[i] = jl_hip + Torque_hip_GLU + Torque_hip_HAM - Torque_hip_HFL
    
    
    




###### Plot


# plt.scatter(t, test1-Ldx_control,s=1)
# plt.grid()
# plt.title("diff Ldx filtered")
# plt.show()

# plt.scatter(t, test2-Rdx_control,s=1)
# plt.grid()
# plt.title("diff Rdx filtered")
# plt.show()

# plt.plot(t,test1,label='Ldx python')
# plt.plot(t,Ldx_control,'--',label='ldx simulink')
# plt.plot(t,Ldx_nf_control,label='ldx no filtered simulink')
# plt.legend()
# plt.show()

# plt.scatter(t, test3-StanceL_control,s=2)
# plt.grid()
# plt.title("diff Stance L")
# plt.show()

# plt.scatter(t, test4-StanceR_control,s=1)
# plt.grid()
# plt.title("diff Stance R")
# plt.show()


# plt.scatter(t, StimL[:,VAS]-Stim_VAS_control,s=1)
# plt.grid()
# plt.title("diff Stim VAS")
# plt.show()

# plt.scatter(t, StimL[:,SOL]-Stim_SOL_control,s=1)
# plt.grid()
# plt.title("diff Stim SOL")
# plt.show()

# plt.scatter(t, StimL[:,GAS]-Stim_GAS_control,s=1)
# plt.grid()
# plt.title("diff Stim GAS")
# plt.show()

# plt.scatter(t, StimL[:,TA]-Stim_TA_control,s=1)
# plt.grid()
# plt.title("diff Stim TA")
# plt.show()

# plt.scatter(t, StimL[:,HAM]-Stim_HAM_control,s=1)
# plt.grid()
# plt.title("diff Stim HAM")
# plt.show()

# plt.scatter(t, StimL[:,GLU]-Stim_GLU_control,s=1)
# plt.grid()
# plt.title("diff Stim GLU")
# plt.show()

# plt.scatter(t, StimL[:,HFL]-Stim_HFL_control,s=1)
# plt.grid()
# plt.title("diff Stim HFL")
# plt.show()

plt.scatter(t, Torque_ankle-Torque_ankle_control,s=1)
plt.grid()
plt.title("diff ankle torque")
plt.show()

plt.scatter(t, Torque_knee-Torque_knee_control,s=1)
plt.grid()
plt.title("diff knee torque")
plt.show()

plt.scatter(t, Torque_hip-Torque_hip_control,s=1)
plt.grid()
plt.title("diff hip torque")
plt.show()





