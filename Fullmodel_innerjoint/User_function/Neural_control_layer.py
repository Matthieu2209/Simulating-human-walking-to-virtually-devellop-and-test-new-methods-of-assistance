#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 14:29:41 2022

@author: matthieuxaussems
"""
import numpy as np
import useful_functions as u_f
import MBsysPy

# 2) Fonction qui permet de savoir quelle jambe a le lead sur l'autre, c'est à dire que lorsque les 2 jambes sont en stance depuis plusieurs 
# itérations, c'est la jambe est rentrée en dernière en stance phase
#
# parametre in : StanceL_memory vecteur de boolean mémoire de l'etat de stance de la jambe gauche (chaque valeur du vecteur correspond à une 
#itération) et StanceR_memory vecteur de boolean mémoire de l'etat de stance de la jambe droite
# parametre out :  LonR (1 si L lead on R et 0 sinon) et RonL (1 si R lead on L 0 sinon)
#

countL = 0
countR = 0


def lead(Stance_memory, tsim):

    global countL
    global countR
    
    if tsim == 0:
        RonL = 0
        LonR = 0

    else:

        if Stance_memory[-1, 0]:

            countL += 1
        else:
            countL = 0

        if Stance_memory[-1, 1]:

            countR += 1
        else:
            countR = 0

        if Stance_memory[-1, 0] and Stance_memory[-1, 1]:

            if countL > countR:
                if countR == 1:
                    RonL = 0
                    LonR = 0
                else:
                    RonL = 1
                    LonR = 0
            if countL < countR:
                if countL == 1:
                    RonL = 0
                    LonR = 0
                else:
                    RonL = 0
                    LonR = 1
        else:
            RonL = 0
            LonR = 0

    return RonL,LonR

# 1) Fonction utile pour réaliser les lois de feedback des muscles en stance et en swing phase.
#
# parametre in : stance_or_swing(boolean true if stance, false if swing), Fm_memory c'est une matrice qui possède dans chaque ligne l'historique
# des Forces motrices à chaque pas de temps, current time, pas de temps qui a permis de générer la mémoire de la force motrice.
# parametre out : Vecteur Stim possédant l'update de la stimulation des muscles
#

PTO = 0
countLeft = 0
countRight = 0
memoryL = np.zeros(1)
memoryR = np.zeros(1)
current_PTO = 0
last_PTO_L = 0
last_PTO_R = 0
# count = 0
left_stim = np.array([0.01,0.01,0.01,0.01,0.01,0.01,0.01,0])
right_stim = np.array([0.01,0.01,0.01,0.01,0.01,0.01,0.01,0])
left_phase = np.zeros(2)
right_phase = np.zeros(2)
def Feedback(current_t, diff_t, Stance_memory, Fm_memory, ipsiDx_thigh, contraDx_thigh, lce_memory, theta_knee_memory,dtheta_knee_memory, theta_trunk_memory, dtheta_trunk_memory,leg):
    
    global left_stim
    global right_stim
    global left_phase
    global right_phase

    VAS = 0
    SOL = 1
    GAS = 2
    TA = 3
    HAM = 4
    GLU = 5
    HFL = 6

    n_muscle = 7

    # données nécessaires pour les lois de reflexes

    G_VAS = 2e-4  # gains for the VAS muscle normalised with the maximal force
    G_SOL = 1.2/4000  # gains for the SOL muscle normalised with the maximal force

    G_GAS = 1.1/1500
    G_TA = 1.1
    G_SOL_TA = 0.0001

    G_HAM = 2.166666666666667e-04
    G_GLU = 1/3000.
    G_HFL = 0.5
    G_HAM_HFL = 4
    G_delta_theta = 1.145915590261647

    # offset

    loff_TA = 0.72
    lopt_TA = 0.06

    loff_HAM = 0.85
    lopt_HAM = 0.10

    loff_HFL = 0.65
    lopt_HFL = 0.11

    # pre-stimulations

    So = 0.01
    So_VAS = 0.08  # VAS
    So_BAL = 0.05  # HAM,GLU,HFL en stance

    # parameters
    k_swing = 0.25
    k_p = 1.909859317102744
    k_d = 0.2
    phi_k_off = 2.967059728390360
    theta_ref = 0.104719755119660

    # delay

    t_l = current_t - 0.02  # (t- time delay) for long neural signal delay
    t_m = current_t - 0.01  # (t- time delay) for middle neural signal delay
    t_s = current_t - 0.005  # (t- time delay) for small neural signal delay

    index_t_l = round((current_t/diff_t)) - round((0.02/diff_t))
    index_t_m = round((current_t/diff_t)) - round((0.01/diff_t))
    index_t_s = round((current_t/diff_t)) - round((0.005/diff_t))
    index_current_t = round((current_t/diff_t))

    global countL
    global countR
  
    global current_PTO
    global last_PTO_L
    global last_PTO_R
    global count
    global PTO
    
    global memoryL
    global memoryR
    
    global countLeft
    global countRight
    

    Stim = np.zeros(n_muscle)
    
    PTO = 0

    if t_m <= 0:
        stance_or_swing = False
        leader = 0
    
    else:
        stance_or_swing = Stance_memory[index_current_t,leg]

    
    leader = lead(Stance_memory, current_t)[leg]
    # MBsysPy.set_output_value(float(lead(Stance_memory, current_t)[0]), 1, "RonL")
    # MBsysPy.set_output_value(theta_knee_memory - phi_k_off, 1, "delta_knee")
    # MBsysPy.set_output_value(dtheta_knee_memory, 1, "dtheta_knee")
    
    if index_current_t >= round((0.01/diff_t)):
        
        if (theta_knee_memory - phi_k_off ) > 0 and dtheta_knee_memory > 0:
            if leg == 0:
                
                memoryL = np.append(memoryL,theta_knee_memory - phi_k_off) 
            else: 
                memoryR =  np.append(memoryR,theta_knee_memory - phi_k_off) 
        else: 
            if leg == 0:
                
                memoryL = np.append(memoryL,0)
            else: 
                memoryR = np.append(memoryR,0)
        if leg == 0:
            knee_state = memoryL[index_t_m]
            # MBsysPy.set_output_value(knee_state, 1, "knee_state")
            # MBsysPy.set_output_value(theta_knee_memory - phi_k_off, 2, "knee_state")
            # MBsysPy.set_output_value(dtheta_knee_memory, 3, "knee_state") 
        else:
            knee_state = memoryR[index_t_m]
            # MBsysPy.set_output_value(knee_state, 4, "knee_state")
            # MBsysPy.set_output_value(theta_knee_memory - phi_k_off, 5, "knee_state")
            # MBsysPy.set_output_value(dtheta_knee_memory, 6, "knee_state")  
       
        
        
     
       
        
    elif index_current_t > 0:
        if leg == 0:
            memoryL  = np.append(memoryL,0)
        else: 
            memoryR = np.append(memoryR,0)
            
        knee_state = 0

    
    # STANCE

    if stance_or_swing:
        if leg == 0:
            countLeft = 0
        else:
            countRight = 0

        # delai court   okkkk
        if t_s <= 0:
            Stim[HAM] = So_BAL
            Stim[GLU] = So_BAL
            Stim[HFL] = So_BAL
           
            PTO = 0
        else:

            theta = theta_trunk_memory[index_t_s, 0]
            d_theta = dtheta_trunk_memory[index_t_s, 0]
            delta_theta = theta - theta_ref
            current_PTO = delta_theta
            PTO = delta_theta

            u1 =  (k_p * delta_theta + k_d* d_theta)

            u1_pos = max(0,u1) 
            u1_neg = min(u1,0)
          
            dx_thigh = ipsiDx_thigh[index_t_s,0]
            dx_thigh_buff_positif = max(0,dx_thigh)
            u2 = dx_thigh_buff_positif * 200

            #HAM
            pre_HAM = So_BAL + u1_pos
            Stim_Stance_HAM = pre_HAM * u2

            Stim[HAM] = max(0,min(Stim_Stance_HAM,1))
            #HFL
            pre_HFL = So_BAL - u1_neg
            Stim_Stance_HFL = pre_HFL * u2
            Stim[HFL] =max(0,min(Stim_Stance_HFL,1)) + leader *k_swing
            
            #GLU
            Stim[GLU] = 0.7 * Stim[HAM] - leader *0.7 * k_swing
            
        #delai moyen
        if t_m <= 0:
            Stim[VAS] = 0
        else: 
            Stim[VAS] = So_VAS + Fm_memory[index_t_m,VAS] * G_VAS - 2 * knee_state - leader *200 * max(0,contraDx_thigh[index_t_s,VAS])

        #delai long
        if t_l <= 0:
            
            Stim[GAS] = So
            Stim[SOL] = So
            Stim[TA] = So
        else:
            
            Stim[GAS] = So + Fm_memory[index_t_l,GAS] * G_GAS
            Stim[SOL] = So + Fm_memory[index_t_l,SOL] * G_SOL
            Stim[TA]  = So - Fm_memory[index_t_l,SOL] * G_SOL_TA + max(0,(lce_memory[index_t_l,TA]/lopt_TA - loff_TA)) * G_TA

    ###SWING
    else:
        
        #delai court
        if leg == 0:
            countLeft +=1
            count = countLeft
        else: 
            countRight +=1
            count =  countRight
            

        if t_s <= 0:
            Stim[HAM] = So
            Stim[GLU] = So
            Stim[HFL] = So
        else:
            
            Stim[HAM] = So + G_HAM * Fm_memory[index_t_s, HAM]
            Stim[GLU] = So + G_GLU * Fm_memory[index_t_s, GLU]
            
            if count == 1 :
                
                theta = theta_trunk_memory[index_t_s, 0]
                delta_theta = theta - theta_ref
                if leg == 0:
                    last_PTO_L = delta_theta
                else: 
                    last_PTO_R = delta_theta
            if leg == 0 :
                last_PTO = last_PTO_L
            else:
                last_PTO = last_PTO_R

            Stim[HFL] = So + G_delta_theta * last_PTO + G_HFL * max(0,(lce_memory[index_t_s,HFL]/lopt_HFL - loff_HFL)) - G_HAM_HFL * max(0,(lce_memory[index_t_s, HAM]/lopt_HAM - loff_HAM))
        
        
        #delai moyen
        Stim[VAS] = 0
        
        #delai long
        Stim[GAS] = So
        Stim[SOL] = So
        
        if t_l <= 0 : 
            Stim[TA] = So
        else: 
            
            Stim[TA] = So +  max(0,(lce_memory[index_t_l,TA]/lopt_TA - loff_TA)) * G_TA
    
        
    Stim[HAM] = max(0.01,min(Stim[HAM],1))
    Stim[GLU] = max(0.01,min(Stim[GLU],1))
    Stim[HFL] = max(0.01,min(Stim[HFL],1))
    Stim[SOL] = max(0.01,min(Stim[SOL],1))
    Stim[GAS] = max(0.01,min(Stim[GAS],1))
    Stim[TA] = max(0.01,min(Stim[TA],1))
    Stim[VAS] = max(0.01,min(Stim[VAS],1))


    # VAS = 0
    # SOL = 1
    # GAS = 2
    # TA = 3
    # HAM = 4
    # GLU = 5
    # HFL = 6
    # if leg == 0:
        # L = np.hstack([Stim,current_t])
        # left_stim = np.vstack([left_stim,[L]])
        # # left_stim = np.vstack([left_stim,[Stim]])
        # phase = np.hstack([stance_or_swing,current_t])
        # left_phase = np.vstack([left_phase,[phase]])
        # MBsysPy.set_output_value(float(stance_or_swing), 1, "stance")
        # MBsysPy.set_output_value(Stim[VAS], 1, "stim_left")
        # MBsysPy.set_output_value(Stim[SOL], 2, "stim_left")
        # MBsysPy.set_output_value(Stim[GAS], 3, "stim_left")
        # MBsysPy.set_output_value(Stim[TA], 4, "stim_left")
        # MBsysPy.set_output_value(Stim[HAM], 5, "stim_left")
        # MBsysPy.set_output_value(Stim[GLU], 6, "stim_left")
        # MBsysPy.set_output_value(Stim[HFL], 7, "stim_left")

    # if leg == 1:
        # R = np.hstack([Stim,current_t])
        # right_stim = np.vstack([right_stim,[R]])
        # # right_stim = np.vstack([right_stim,[Stim]])
        # phase = np.hstack([stance_or_swing,current_t])
        # right_phase = np.vstack([right_phase,[phase]])

        # MBsysPy.set_output_value(float(stance_or_swing), 2, "stance")
        # MBsysPy.set_output_value(Stim[VAS], 1, "stim_right")
        # MBsysPy.set_output_value(Stim[SOL], 2, "stim_right")
        # MBsysPy.set_output_value(Stim[GAS], 3, "stim_right")
        # MBsysPy.set_output_value(Stim[TA], 4, "stim_right")
        # MBsysPy.set_output_value(Stim[HAM], 5, "stim_right")
        # MBsysPy.set_output_value(Stim[GLU], 6, "stim_right")
        # MBsysPy.set_output_value(Stim[HFL], 7, "stim_right")
    return Stim

    


    
    
    
    
