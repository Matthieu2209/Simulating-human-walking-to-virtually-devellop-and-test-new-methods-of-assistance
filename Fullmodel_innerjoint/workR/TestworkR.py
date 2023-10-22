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
import time
import shutil



def runtest(dt0,tf,c=False):
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
    #mbs_dirdyn.set_options(dt0=250e-7, tf=0.5, save2file=1) #1124.726s
    #mbs_dirdyn.set_options(dt0=250e-7, tf=0.005, save2file=1) #5s
    #mbs_dirdyn.set_options(dt0=250e-7, tf=0.05, save2file=1) #40s
    #mbs_dirdyn.set_options(dt0=250e-7, tf=0.1, save2file=1) #95s
    #mbs_dirdyn.set_options(dt0=250e-7, tf=1, save2file=1) #4385s
    #mbs_dirdyn.set_options(dt0=250e-7, tf=0.0001, save2file=1) #1s
    #mbs_dirdyn.set_options(dt0=250e-7, tf=0.00001, save2file=1) #0.9s
    #mbs_dirdyn.set_options(dt0=250e-7, tf=0.000001, save2file=1) #0.9s
    
    from datetime import datetime

    now = datetime.now()

    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    
    mbs_dirdyn.set_options(dt0=dt0, tf=tf, save2file=1) # 96
    start_time = time.time()

    results = mbs_dirdyn.run()
    
    elapsed_time = time.time() - start_time


    if(c):
        import os
        # Replace 'file1.anim' with the actual file path
        dirpath=str(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

        df1 = pd.read_csv(dirpath+'/animationR/'+'dirdyn_q.anim', delimiter='\t')
        df2 = pd.read_csv(dirpath+'/animationR/'+str(tf)+'dirdyn_q.anim', delimiter='\t')

        are_equal = df1.equals(df2)

        if are_equal:
            print("\n\nCONGRATS : The two .anim files are the same for duration : "+str(tf))
        else:
            print("\n\nThe two .anim files are different for duration : "+str(tf))
            return
    


    # Calculate the elapsed time in seconds
    elapsed_time = time.time() - start_time

    # Convert the elapsed time to minutes and round to two decimal places
    elapsed_time_minutes = round(elapsed_time / 60, 2)

    # Print the elapsed time in minutes
    print(f"Time taken to run the line: {elapsed_time_minutes:.2f} minutes")



    import os
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    now = str(datetime.now())[:19]
    now = now.replace(":","_")

    tf=tf
    dt0=dt0
    src_dir=parent_dir+"/animationR/dirdyn_q.anim"
    dst_dir=parent_dir+"/animationR/archive/tf:"+str(tf)+"dt0"+str(dt0)+"rt"+str(elapsed_time_minutes)+".anim"

    shutil.copy(src_dir,dst_dir)

#0.02: 17s

if __name__ == "__main__":
    
    for tf in [0.5]:#,0.01,0.1]:
        print('starting tf' +str(tf))
        runtest(250e-6,tf)
    

""" 

import winsound    
winsound.Beep(1440, 200)    


import numpy as np

arr= np.zeros(2)

for i in range (100):
    arr=np.vstack([arr,[i,i]])[-20:]

print(arr)


 import numpy as np
print(np.arange(100)[:200])

dt=round(250e-7,10)
print(dt)

relevant_timesteps=int(0.21/dt)
import numpy as np

a=np.arange(100000)
print(a[-relevant_timesteps:])


counter=0
    global counter
    counter+=1
    if(counter%10000==0):
        print(time[-1],t)
    
                global elapsed_time
            start=time.time()
                        elapsed_time += time.time()-start
            if(tsim>0.095):
                print("Elapsed time"+str(elapsed_time))
 """
 