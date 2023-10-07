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
import os
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0,  os.path.join(parent_dir, "User_function"))
sys.path.insert(1,  os.path.join(parent_dir, "userfctR"))


# Change the current working directory to the script's directory
os.chdir(os.path.dirname(os.path.abspath(__file__))) # the script is running from it's parent directory ?


# ===========================================================================
# Project loading
# =============================================================================
import pandas as pd



def runtest(dt0,tf):
    mbs_data = Robotran.MbsData('../dataR/Fullmodel_innerjoint.mbs',)
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
    #mbs_dirdyn.set_options(dt0=250e-7, tf=1.5, save2file=1) #  9560.745 seconds
    #mbs_dirdyn.set_options(dt0=250e-7, tf=0.005, save2file=1) # 5.033s
    #mbs_dirdyn.set_options(dt0=250e-7, tf=0.01, save2file=1) # 7.333s
    #mbs_dirdyn.set_options(dt0=250e-7, tf=0.1, save2file=1) # 96
    
    mbs_dirdyn.set_options(dt0=dt0, tf=tf, save2file=1) # 96

    results = mbs_dirdyn.run()


    # Replace 'file1.anim' with the actual file path
    dirpath=str(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    df1 = pd.read_csv(dirpath+'/animationR/'+'dirdyn_q.anim', delimiter='\t')
    df2 = pd.read_csv(dirpath+'/animationR/'+str(tf)+'dirdyn_q.anim', delimiter='\t')

    are_equal = df1.equals(df2)

    if are_equal:
        print("\n\nThe two .anim files are the same for duration : "+str(tf))
    else:
        print("\n\nThe two .anim files are different for duration : "+str(tf))
        return

print("ish")
        
for tf in [0.002]:#,0.01,0.1]:
    print('starting tf' +str(tf))
    runtest(250e-7,tf)