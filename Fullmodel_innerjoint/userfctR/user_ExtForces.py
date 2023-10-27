# -*- coding: utf-8 -*-
"""Module for defining user function required to compute external forces."""
# Author: Robotran Team
# (c) Universite catholique de Louvain, 2019


import math
import numpy as np
import MBsysPy
# Useful data for the contact force model

import sys

import os
import time
# Get the directory where your script is located
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0,  os.path.join(parent_dir, "User_function"))
sys.path.insert(1,  os.path.join(parent_dir, "userfctR"))
import gait_graph

ground_limit = 0
vmax = 0.03
kz = 78480
kx = 7848
must = 0.9
musl = 0.8
v_limit = 0.01


Q_ballL=1
Qn_ballL=0
Stiction_test_ballL=0
Stiction_prec_test_ballL=0 
prec_slide_test_ballL=0
prec_stick_test_ballL=0

Q_heelL=1
Qn_heelL=0
Stiction_test_heelL=0
Stiction_prec_test_heelL=0 
prec_slide_test_heelL=0
prec_stick_test_heelL=0

x0_ballL = 0
x0_heelL = 0

Q_ballR=1
Qn_ballR=0
Stiction_test_ballR=0
Stiction_prec_test_ballR=0 
prec_slide_test_ballR=0
prec_stick_test_ballR=0

Q_heelR=1
Qn_heelR=0
Stiction_test_heelR=0
Stiction_prec_test_heelR=0 
prec_slide_test_heelR=0
prec_stick_test_heelR=0

x0_ballR = 0
x0_heelR = 0



flag_graph = True

# Useful function for the contact force model

def flip_flop_SR(S, R, Q, Qn):
    if S and not R: # Set
        Q = 1
        Qn = 0
    elif R and not S: # Reset
        Q = 0
        Qn = 1
    elif not R and not S: # Hold
        Q=Q
        Qn=Qn
    else: # Invalid input
        print("Invalid input: S and R cannot be both 1")
    return Q, Qn



def user_ExtForces(PxF, RxF, VxF, OMxF, AxF, OMPxF, mbs_data, tsim, ixF):
    

    """Compute an user-specified external force.

    Parameters
    ----------
    PxF : numpy.ndarray
        Position vector (index starting at 1) of the force sensor expressed in
        the inertial frame: PxF[1:4] = [P_x, P_y, P_z]
    RxF : numpy.ndarray
        Rotation matrix (index starting at 1) from the inertial frame to the
        force sensor frame: Frame_sensor = RxF[1:4,1:4] * Frame_inertial
    VxF : numpy.ndarray
        Velocity vector (index starting at 1) of the force sensor expressed in
        the inertial frame: VxF[1:4] = [V_x, V_y, V_z]
    OMxF : numpy.ndarray
        Angular velocity vector (index starting at 1) of the force sensor
        expressed in the inertial frame: OMxF[1:4] = [OM_x, OM_y, OM_z]
    AxF : numpy.ndarray
        Acceleration vector (index starting at 1) of the force sensor expressed
        in the inertial frame: AxF[1:4] = [A_x, A_y, A_z]
    OMPxF : numpy.ndarray
        Angular acceleration vector (index starting at 1) of the force sensor
        expressed in the inertial frame: OMPxF[1:4] = [OMP_x, OMP_y, OMP_z]
    mbs_data : MBsysPy.MbsData
        The multibody system associated to this computation.
    tsim : float
        The current time of the simulation.
    ixF : int
        The ID identifying the computed force sensor.

    Notes
    -----
    For 1D numpy.ndarray with index starting at 1, the first index (array[0])
    must not be modified. The first index to be filled is array[1].

    For 2D numpy.ndarray with index starting at 1, the first row (mat[0, :]) and
    line (mat[:,0]) must not be modified. The subarray to be filled is mat[1:, 1:].

    Returns
    -------
    Swr : numpy.ndarray
        An array of length 10 equal to [0., Fx, Fy, Fz, Mx, My, Mz, dxF].
        F_# are the forces components expressed in inertial frame.
        M_# are the torques components expressed in inertial frame.
        dxF is an array of length 3 containing the component of the forces/torque
        application point expressed in the BODY-FIXED frame.
        This array is a specific line of MbsData.SWr.
    """
    
    #global variable
    
    global Q_ballL
    global Qn_ballL
    global Stiction_test_ballL
    global Stiction_prec_test_ballL
    global x0_ballL
    global prec_slide_test_ballL
    global prec_stick_test_ballL
    
    global Q_heelL
    global Qn_heelL
    global Stiction_test_heelL
    global Stiction_prec_test_heelL
    global x0_heelL
    global prec_slide_test_heelL
    global prec_stick_test_heelL
    
    global Q_ballR
    global Qn_ballR
    global Stiction_test_ballR
    global Stiction_prec_test_ballR
    global x0_ballR
    global prec_slide_test_ballR
    global prec_stick_test_ballR
    
    global Q_heelR
    global Qn_heelR
    global Stiction_test_heelR
    global Stiction_prec_test_heelR
    global x0_heelR
    global prec_slide_test_heelR
    global prec_stick_test_heelR
    
    #Resetting variables

    Fx = 0.0
    Fy = 0.0
    Fz = 0.0
    Mx = 0.0
    My = 0.0
    Mz = 0.0
    idpt = mbs_data.xfidpt[ixF]
    dxF = mbs_data.dpt[1:, idpt]
    
    #Initialization of external force sensors:
    
    Force_BallL = mbs_data.extforce_id["Force_BallL"]
    Force_HeelL = mbs_data.extforce_id["Force_HeelL"]
    Force_BallR = mbs_data.extforce_id["Force_BallR"]
    Force_HeelR = mbs_data.extforce_id["Force_HeelR"]
    
    
    #Contact model for external force on BALL L
    
    if ixF == Force_BallL:

        MBsysPy.set_output_value(PxF[3], 4, "pos_Z")
        MBsysPy.set_output_value(VxF[1], 4, "velocity_X")
        MBsysPy.set_output_value(VxF[3], 4, "velocity_Z")
                
        #Contact model
        Fx, Fz, Q_ballL, Qn_ballL, Stiction_test_ballL , Stiction_prec_test_ballL ,  x0_ballL , prec_slide_test_ballL, prec_stick_test_ballL = GRF(PxF, RxF, VxF, OMxF, AxF, OMPxF, mbs_data, tsim, ixF, Q_ballL, Qn_ballL, Stiction_test_ballL , Stiction_prec_test_ballL ,  x0_ballL ,prec_slide_test_ballL, prec_stick_test_ballL)
        
          

    #Contact model for external force on BALL R
    
    if ixF == Force_BallR:
        
        MBsysPy.set_output_value(PxF[3], 3, "pos_Z")
        MBsysPy.set_output_value(VxF[1], 3, "velocity_X")
        MBsysPy.set_output_value(VxF[3], 3, "velocity_Z")
        gait_graph.collect_Px(PxF,tsim)

        #Contact model
        Fx, Fz, Q_ballR, Qn_ballR, Stiction_test_ballR, Stiction_prec_test_ballR ,  x0_ballR , prec_slide_test_ballR, prec_stick_test_ballR = GRF(PxF, RxF, VxF, OMxF, AxF, OMPxF, mbs_data, tsim, ixF, Q_ballR, Qn_ballR, Stiction_test_ballR , Stiction_prec_test_ballR ,  x0_ballR , prec_slide_test_ballR, prec_stick_test_ballR)

    #Contact model for external force on HEEL L

    if ixF == Force_HeelL:
        
        MBsysPy.set_output_value(PxF[3], 2, "pos_Z")
        MBsysPy.set_output_value(VxF[1], 2, "velocity_X")
        MBsysPy.set_output_value(VxF[3], 2, "velocity_Z")
        
        #Contact model
        Fx, Fz, Q_heelL, Qn_heelL, Stiction_test_heelL, Stiction_prec_test_heelL ,  x0_heelL , prec_slide_test_heelL, prec_stick_test_heelL  = GRF(PxF, RxF, VxF, OMxF, AxF, OMPxF, mbs_data, tsim, ixF, Q_heelL, Qn_heelL, Stiction_test_heelL , Stiction_prec_test_heelL ,  x0_heelL , prec_slide_test_heelL , prec_stick_test_heelL)

    
    #Contact model for external force on HEEL R

    if ixF == Force_HeelR:
        
        MBsysPy.set_output_value(PxF[3], 1, "pos_Z")
        MBsysPy.set_output_value(VxF[1], 1, "velocity_X")
        MBsysPy.set_output_value(VxF[3], 1, "velocity_Z")
        #Contact model
        Fx, Fz, Q_heelR, Qn_heelR, Stiction_test_heelR, Stiction_prec_test_heelR ,  x0_heelR , prec_slide_test_heelR , prec_stick_test_heelR  = GRF(PxF, RxF, VxF, OMxF, AxF, OMPxF, mbs_data, tsim, ixF, Q_heelR, Qn_heelR, Stiction_test_heelR , Stiction_prec_test_heelR ,  x0_heelR , prec_slide_test_heelR , prec_stick_test_heelR)



    # Concatenating force, torque and force application point to returned array.
    # This must not be modified.
    Swr = mbs_data.SWr[ixF]
    Swr[1:] = [Fx, Fy, Fz, Mx, My, Mz, dxF[0], dxF[1], dxF[2]]
    
    
    Fx_HeelR = mbs_data.SWr[Force_HeelR][1]
    Fz_HeelR = mbs_data.SWr[Force_HeelR][3]
    
    Fx_HeelL = mbs_data.SWr[Force_HeelL][1]
    Fz_HeelL = mbs_data.SWr[Force_HeelL][3]
    
    Fx_BallR = mbs_data.SWr[Force_BallR][1]
    Fz_BallR = mbs_data.SWr[Force_BallR][3]
    
    Fx_BallL = mbs_data.SWr[Force_BallL][1]
    Fz_BallL = mbs_data.SWr[Force_BallL][3]
    
    
    MBsysPy.set_output_value(Fx_HeelR, 1, "external_force_X")
    MBsysPy.set_output_value(Fx_HeelL, 2, "external_force_X")
    MBsysPy.set_output_value(Fx_BallR, 3, "external_force_X")
    MBsysPy.set_output_value(Fx_BallL, 4, "external_force_X")


    MBsysPy.set_output_value(Fz_HeelR, 1, "external_force_Z")
    MBsysPy.set_output_value(Fz_HeelL, 2, "external_force_Z")
    MBsysPy.set_output_value(Fz_BallR, 3, "external_force_Z")
    MBsysPy.set_output_value(Fz_BallL, 4, "external_force_Z")
    
    gait_graph.collect_ext(mbs_data,ixF,Fx,PxF[1],VxF[1],tsim)      

    return Swr

def GRF(PxF, RxF, VxF, OMxF, AxF, OMPxF, mbs_data, tsim, ixF, Q, Qn, Stiction_test , Stiction_prec_test,  x0 , prec_slide_test , prec_stick_test):

    #Contact model
    d_z = PxF[3] -ground_limit
    delta_x = 0
    delta_vx = 0
    Force_slide = 0
    Force_stick = 0
    Fx_mod = 0
    Fz = 0
    slide_test =0
    stick_test = 0
    
    #Contact
    if d_z >=0 :
        
        #Vertical ground force model
            
        if VxF[3]/vmax >= -1:
            Fz = (kz * -d_z *(1+VxF[3]/vmax))
        else :
            Fz = 0
        
        #update stiction mode
        old_Stiction = Stiction_test
        Stiction_test =Stiction_prec_test
        

    
        #Static friction model
    
        if Stiction_test == 1:
            
            #resetting parameters
            Force_slide = 0
            stick_test = 0
            
            #if(ixF==3):
                #print(PxF[1],x0,tsim)

            #force compute
            delta_x = PxF[1] - x0
            delta_vx = VxF[1] / vmax
            
            Fx_mod = -kx * delta_x * (1 + math.copysign(1, delta_x) * delta_vx)
            
            if Fz < 0 :
                Force_stick = Fx_mod  #Force stick
            else:
                Force_stick = 0
            
            #determination of friction mode
            if abs(Force_stick)-abs(Fz*must) > 0:
                slide_test = 1
            else :
                slide_test = 0
                

        
        
        #Kinetic friction model
        else :
            #update x0 (PosFP)
            
            if Stiction_test == 1:
                if old_Stiction == 0:
                    x0=0#PxF[3]
            else:
                x0=0#PxF[3] 

            #resetting parameters
            Force_stick = 0
            slide_test =0
            
            #force compute
            if VxF[1] != 0 :
                Fx_mod = math.copysign(1, VxF[1]) * musl * Fz
            else :
                Fx_mod = 0
        
            Force_slide = Fx_mod
            
            #determination of friction mode
            if abs(VxF[1]) - v_limit <= 0:
                stick_test = 1
            else :
                stick_test= 0
                

                
            
        #Horizontal ground force
        Fx = (Force_slide + Force_stick )

        
        #stiction mode
        Q, Qn = flip_flop_SR(stick_test,slide_test, Q, Qn)
        Stiction_prec_test=Q
        
    
    #no contact   
    else:
        
        Fz =0
        Fx =0
        
        Q=1
        Qn=0
        Stiction_prec_test=0
    
    
    
    if(ixF==1):
        gait_graph.collect_stiction(mbs_data, Stiction_test, Force_slide, delta_x, delta_vx,Fx_mod, Force_stick, d_z, Fz, slide_test, stick_test, tsim)

    
    return Fx, Fz, Q, Qn, Stiction_test , Stiction_prec_test,  x0 , prec_slide_test , prec_stick_test
      



# (c) Universite catholique de Louvain, 2020
import sys

import os
import time
# Get the directory where your script is located
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0,  os.path.join(parent_dir, "User_function"))
sys.path.insert(1,  os.path.join(parent_dir, "userfctR"))
sys.path.insert(2,  os.path.join(parent_dir, "workR"))


import TestworkR


if __name__ == "__main__":
    TestworkR.runtest(250e-7,0.9)