
    
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






flag_initiated=False
t=np.zeros((0))
dt=0

PxA=0
PzA = 0
VxA = 0
VzA = 0
StickingA = 0
SlidingA = 0
StictionA = 0
PosFPA = 0
dxA = 0
dvxA = 0
FxA = 0
FzA = 0
Fx_slidingA = 0
Fx_stickingA = 0
test_slideA= 0
test_stickA= 0

def initiate():
    global t, dt, PxA, PzA, VxA, VzA, StickingA, SlidingA, StictionA, PosFPA, dxA, dvxA, FxA, FzA, Fx_slidingA, Fx_stickingA, test_slideA, test_stickA
    dt,tsim=np.load("paramaters.npy")
    t = np.arange(0, tsim + dt, dt)  # Create the time array directly
    l = len(t)  # Calculate the length based on t
    PxA=np.zeros((5,l))
    PzA = np.zeros((5,l))
    VxA = np.zeros((5,l))
    VzA = np.zeros((5,l))
    StickingA = np.zeros((5,l))
    SlidingA = np.zeros((5,l))
    StictionA = np.zeros((5,l))
    PosFPA = np.zeros((5,l))
    dxA = np.zeros((5,l))
    dvxA = np.zeros((5,l))
    FxA = np.zeros((5,l))
    FzA = np.zeros((5,l))
    Fx_slidingA = np.zeros((5,l))
    Fx_stickingA = np.zeros((5,l))
    test_slideA = np.zeros((5,l))
    test_stickA = np.zeros((5,l))

    
def collect(ixF,Px,Pz,Vx,Vz,Sticking,Sliding,Stiction,PosFP,dx,dvx,Fx,Fz,Fx_sliding,Fx_sticking,test_slide,test_stick,tsim):
    global flag_initiated, t, dt, PxA, PzA, VxA, VzA, StickingA, SlidingA, StictionA, PosFPA, dxA, dvxA, FxA, FzA, Fx_slidingA, Fx_stickingA, test_slideA, test_stickA

    if (flag_initiated==False):
       initiate()
       flag_initiated=True
    ti= int(tsim/dt)
       
    PxA[ixF,ti]=Px
    PzA[ixF, ti] = Pz
    VxA[ixF, ti] = Vx
    VzA[ixF, ti] = Vz
    StickingA[ixF, ti] = Sticking
    SlidingA[ixF, ti] = Sliding
    StictionA[ixF, ti] = Stiction
    PosFPA[ixF, ti] = PosFP
    dxA[ixF, ti] = dx
    dvxA[ixF, ti] = dvx
    FxA[ixF, ti] = Fx
    FzA[ixF, ti] = Fz
    Fx_slidingA[ixF, ti] = Fx_sliding
    Fx_stickingA[ixF, ti] = Fx_sticking
    test_slideA[ixF, ti] = test_slide
    test_stickA[ixF, ti] = test_stick
    
    
    
def show_ext():
    global t, dt, PxA, PzA, VxA, VzA, StickingA, SlidingA, StictionA, PosFPA, dxA, dvxA, FxA, FzA, Fx_slidingA, Fx_stickingA ,test_slideA, test_stickA

    #name = ["temp", "HeelR", "HeelL", "BallR", "BallL"]
    name = ["temp", "BallR", "HeelL", "HeelR", "BallL"]
    
    variables = [PxA, PzA, VxA, VzA, StickingA, SlidingA, StictionA, PosFPA, dxA, dvxA, FxA, FzA, Fx_slidingA, Fx_stickingA, test_slideA, test_stickA]
    variable_names = ["Px", "Pz", "Vx", "Vz", "Sticking", "Sliding", "Stiction", "PosFP", "dx", "dvx", "Fx", "Fz", "Fx_sliding", "Fx_sticking", "test_slide", "test_stick"]

    for i in range(1, 5):
        """ 
        id = "plot/" + name[i]
        plt.plot(variables[i], label=name[i])
        plt.legend()
        plt.title(id)
        plt.savefig(id)
        plt.close()
        print(id) 
        """

        for j in range(len(variables)):
            id = "plot/"+ name[i] + "/" + name[i] + " " + variable_names[j]
            plt.figure(figsize=(12, 6))
            plt.plot(t, variables[j][i], label=variable_names[j])
            plt.legend()
            plt.grid()
            plt.title(id)
            plt.savefig(id)
            plt.close()
            
            np.save("numpy_archive/"+variable_names[j],variables[j])



def replay_show_ext():
    dt,tsim=np.load("paramaters.npy")
    t = np.arange(0, tsim + dt, dt)  # Create the time array directly
    
    variable_names = ["Px", "Pz", "Vx", "Vz", "Sticking", "Sliding", "Stiction", "PosFP", "dx", "dvx", "Fx", "Fz", "Fx_sliding", "Fx_sticking", "test_slide", "test_stick"]
    variables = [ np.load(os.getcwd()+"/numpy_archive/"+var+".npy") for var in variable_names]
    #name = ["temp", "HeelR", "HeelL", "BallR", "BallL"]
    name = ["temp", "BallR", "HeelL", "HeelR", "BallL"]

    for i in range(1, 5):
        """ 
        id = "plot/" + name[i]
        plt.plot(variables[i], label=name[i])
        plt.legend()
        plt.title(id)
        plt.savefig(id)
        plt.close()
        print(id) 
        """

        for j in range(len(variables)):
            id = "plot/"+ name[i] + "/" + name[i] + " " + variable_names[j]
            plt.figure(figsize=(12, 6))
            plt.scatter(t, variables[j][i], label=variable_names[j])
            delta=0.00035
            plt.axvline(x=delta+0.163742, color='red', linestyle='--')
            plt.axvline(x=delta+0.16445, color='red', linestyle='--')
            plt.axvline(x=delta+0.16585, color='red', linestyle='--')
            plt.xlim([0.1630,0.1690])
            #plt.xlim([0.13,0.23])
            #plt.ylim([-0.001,0.001])
            plt.grid()
            plt.legend()
            plt.title(id)
            plt.savefig(id)
            plt.close()
            
            if(variable_names[j]=="temp"):
                id = "plot/"+ name[i] + "/" + name[i] + " " + variable_names[j]+"zommed"
                plt.plot(t, variables[j][i], label=variable_names[j])
                plt.xlim([0.15,0.25])
                plt.ylim([-0.1,0.1])
                plt.grid()
                plt.legend()
                plt.title(id)
                plt.savefig(id)
                plt.close()
                
            
            
    
    
    
    
if __name__ == "__main__":
    replay_show_ext()
    #TestworkR.runtest(250e-7,0.3,c=False)    




















""" 


BallL=np.zeros((0))
BallR=np.zeros((0))
HeelL=np.zeros((0))
HeelR=np.zeros((0))

BallR_stiction=np.zeros(0)
Px=np.zeros(0)

#HeelL_Px=np.zeros((0,2))



def initiate_temp():
    global dt, t, BallL, BallR, HeelL, HeelR, BallR_stiction, Px
    dt,tsim=np.load("paramaters.npy")
    t = np.arange(0, tsim + dt, dt)  # Create the time array directly
    l = len(t)  # Calculate the length based on t
    BallL=(np.zeros((3,l)))
    BallR=(np.zeros((3,l)))
    HeelL=np.zeros((3,l))
    HeelR=np.zeros((3,l))
    BallR_stiction=np.zeros((11,l))
    Px=np.zeros((4,l))

    
    


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

def collect_stiction(mbs_data,Stiction_test_ballR,Force_slide_ballR,delta_x,delta_vx,Fx_mod, Force_stick_ballR,d_z, Fz,  slide_test, stick_test,delta_test, tsim):
    global flag_initiated
    if (flag_initiated==False):
       initiate()
       flag_initiated=True
       
    global t, dt, BallR_stiction
    ti= int(tsim/dt)
    
    BallR_stiction[:,ti]=[Stiction_test_ballR, Force_slide_ballR,delta_x, delta_vx,Fx_mod, Force_stick_ballR, d_z, Fz, slide_test, stick_test,delta_test]

    
      
def collect_Px(PxF,tsim):
    global flag_initiated
    if (flag_initiated==False):
       initiate()
       flag_initiated=True
       
    global t, dt, Px
    ti= int(tsim/dt)
    
    Px[:,ti]=PxF
    






   

def show_ext():
    global t, dt,  BallL, BallR, HeelL, HeelR, BallR_stiction, Px
    
    
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

    id="plot_archive/delta_x BallR"
    plt.plot(t, BallR_stiction[2], label="delta_x")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close() 
    
    id="plot_archive/delta_vx BallR"
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
    
    id="plot_archive/BallR D_Z"
    plt.plot(t, BallR_stiction[6], label="BallR D_Z")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close() 
    
    id="plot_archive/BallR Fz"
    plt.plot(t, BallR_stiction[7], label="BallR Fz")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close() 
  
    id="plot_archive/BallR slide test"
    plt.plot(t, BallR_stiction[8], label="BallR slide test")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close()   
    
    
    id="plot_archive/BallR stick test"
    plt.plot(t, BallR_stiction[9], label="BallR stick test")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close()   
    
    id="plot_archive/detla x_test"
    plt.plot(t, BallR_stiction[10], label="detla x_test")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close()  
    
    
    id="plot_archive/stiction_ballr"
    plt.figure(figsize=(12, 6))
    plt.plot(t, BallR_stiction[0], label="Stiction_test_ballR")
    plt.plot(t, BallR_stiction[1], label="Force_slide_ballR")
    plt.plot(t, BallR_stiction[2], label="delta_x")
    plt.plot(t, BallR_stiction[3], label="delta_vx")
    plt.plot(t, BallR_stiction[4], label="Fx_mod")
    plt.plot(t, BallR_stiction[5], label="Force_stick_ballR")
    plt.plot(t, BallR_stiction[6], label="BallR D_Z")
    #plt.xlim([0.2510,0.2525])
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close() 
    
    
    
    id="plot_archive/Px"
    plt.plot(t, Px[1], label="Px 1")
    plt.plot(t, Px[3], label="Px 3")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
    plt.close() 
    
    
    np.save("numpy_archive/stiction_BallR", BallR_stiction)
    np.save("numpy_archive/BallR", BallR)
    np.save("numpy_archive/BallL", BallL)
    np.save("numpy_archive/HeelL", HeelL)
    np.save("numpy_archive/HeelR", HeelR)


"""

    
""" initiate()
    BallR_stiction = np.load(os.getcwd()+"/numpy_archive/stiction_BallR.npy")
    print(BallR_stiction)
    
    id="plot_archive/Force_stick_ballR"
    print(t.shape,BallR_stiction[5].shape)
    
    plt.plot(t[10:], BallR_stiction[5,10:], label="Force_stick_ballR")
    plt.legend()
    plt.title(id)
    plt.savefig(id)
plt.close()  
"""

     