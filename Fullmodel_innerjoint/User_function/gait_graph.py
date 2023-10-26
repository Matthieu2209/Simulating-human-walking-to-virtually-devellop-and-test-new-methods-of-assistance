
    
import os
import sys 
import matplotlib.pyplot as plt
import numpy as np


import os
import time
# Get the directory where your script is located
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0,  os.path.join(parent_dir, "User_function"))
sys.path.insert(1,  os.path.join(parent_dir, "userfctR"))
sys.path.insert(2,  os.path.join(parent_dir, "workR"))
import TestworkR



t=np.zeros((0))
dt=0
BallL=np.zeros((0))
BallR=np.zeros((0))
HeelL=np.zeros((0))
HeelR=np.zeros((0))

BallR_stiction=np.zeros(0)

#HeelL_Px=np.zeros((0,2))

def initiate():
    global dt, t, BallL, BallR, HeelL, HeelR, BallR_stiction
    dt,tsim=np.load("paramaters.npy")
    t = np.arange(0, tsim + dt, dt)  # Create the time array directly
    l = len(t)  # Calculate the length based on t
    BallL=(np.zeros((3,l)))
    BallR=(np.zeros((3,l)))
    HeelL=np.zeros((3,l))
    HeelR=np.zeros((3,l))
    BallR_stiction=np.zeros((6,l))
    
    


flag_initiated=False

def collect_ext(mbs_data,type,fx,px,vx,tsim):
    global flag_initiated
    if (flag_initiated==False):
       initiate()
       flag_initiated=True
    
    global t, dt,  BallL , BallR , HeelL , HeelR
    ti= int(tsim/dt)
    if(type==mbs_data.extforce_id["Force_BallL"]):
        BallL[:,ti]=[fx,px,vx]
        
    if(type==mbs_data.extforce_id["Force_BallR"]):
        BallR[:,ti]=[fx,px,vx]
        
    if(type==mbs_data.extforce_id["Force_HeelL"]):
        HeelL[:,ti]=[fx,px,vx]
        
    if(type==mbs_data.extforce_id["Force_HeelR"]):
        HeelR[:,ti]=[fx,px,vx]    

def collect_stiction(mbs_data,Stiction_test_ballR,Force_slide_ballR,delta_x,delta_vx,Fx_mod, Force_stick_ballR,tsim):
    global flag_initiated
    if (flag_initiated==False):
       initiate()
       flag_initiated=True
       
    global t, dt, BallR_stiction
    ti= int(tsim/dt)
    
    BallR_stiction[:,ti]=[Stiction_test_ballR,Force_slide_ballR,delta_x,delta_vx,Fx_mod, Force_stick_ballR]

    
        


        
    
   

def show_ext():
    global t, dt,  BallL, BallR, HeelL, HeelR, BallR_stiction
    
    
    #BallL
    id="plot_archive/BallL_Fx"
    plt.plot(t, BallL[0], label="BallL FX")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close()
    
    id="plot_archive/BallL_Px"
    plt.plot(t, BallL[1], label="BallL PX")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close()
    
    id="plot_archive/BallL_Vx"
    plt.plot(t, BallL[2], label="BallL VX")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close() 
    
    
    #BallR
    id="plot_archive/BallR_Fx"
    plt.plot(t, BallR[0], label="BallR FX")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close()
    
    id="plot_archive/BallR_Px"
    plt.plot(t, BallR[1], label="BallR PX")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close()
    
    id="plot_archive/BallR_Vx"
    plt.plot(t, BallR[2], label="BallR VX")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close() 
    
    
    #HeelL
    id="plot_archive/HeelL_Fx"
    plt.plot(t, HeelL[0], label="HeelL FX")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close()
    
    id="plot_archive/HeelL_Px"
    plt.plot(t, HeelL[1], label="HeelL PX")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close()
    
    id="plot_archive/HeelL_Vx"
    plt.plot(t, HeelL[2], label="HeelL VX")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close() 
    
        
    #HeelR
    id="plot_archive/HeelR_Fx"
    plt.plot(t, HeelR[0], label="HeelR FX")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close()
    
    id="plot_archive/HeelR_Px"
    plt.plot(t, HeelR[1], label="HeelR PX")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close()
    
    id="plot_archive/HeelR_Vx"
    plt.plot(t, HeelR[2], label="HeelR VX")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close() 
    
    
    #mixed
    id="plot_archive/Fx"
    plt.plot(t, BallL[0], label="BallL FX")
    plt.plot(t, BallR[0], label="BallR FX")
    plt.plot(t, HeelL[0], label="HeelL FX")
    plt.plot(t, HeelR[0], label="HeelR FX")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close()
    
    
    id="plot_archive/Px"
    plt.plot(t, BallL[1], label="BallL pX")
    plt.plot(t, BallR[1], label="BallR PX")
    plt.plot(t, HeelL[1], label="HeelL PX")
    plt.plot(t, HeelR[1], label="HeelR PX")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close()
    
    id="plot_archive/Vx"
    plt.plot(t, BallL[2], label="BallL VX")
    plt.plot(t, BallR[2], label="BallR VX")
    plt.plot(t, HeelL[2], label="HeelL VX")
    plt.plot(t, HeelR[2], label="HeelR VX")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close()
    
    id="plot_archive/Stiction_test_ballR"
    plt.plot(t, BallR_stiction[0], label="Stiction_test_ballR")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close() 
    
    id="plot_archive/Force_slide_ballR"
    plt.plot(t, BallR_stiction[1], label="Force_slide_ballR")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close() 

    id="plot_archive/delta_x"
    plt.plot(t, BallR_stiction[2], label="delta_x")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close() 
    
    id="plot_archive/delta_vx"
    plt.plot(t, BallR_stiction[3], label="delta_vx")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close() 
    
    id="plot_archive/Fx_mod"
    plt.plot(t, BallR_stiction[4], label="Fx_mod")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close() 
    
    id="plot_archive/Force_stick_ballR"
    plt.plot(t, BallR_stiction[5], label="Force_stick_ballR")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close() 
    


if __name__ == "__main__":
    TestworkR.runtest(250e-7,0.6,c=False)

    