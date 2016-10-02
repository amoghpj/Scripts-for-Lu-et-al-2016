from scipy.integrate import ode
import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
import random
from tqdm import tqdm
#Define System:
D=3.8*10**(-3)
#Figure 4c looks at variation in parameter Omega, the phase of the system
O_range=np.arange(-4,4,0.05)

def f(t, z, p):
    K,F,D,O=p
    s=0.5*((K*z+F)-z**2*np.conj(K*z+F))-(D+1j*O)*z
    return s
l=0
#Setup simulation to identify steady state for each value of K.
t0=0
t1 = 1000
dt = 0.1
fp=[]
deltap=[np.pi/4,np.pi/2,np.pi*3/4,-np.pi,-np.pi/4,-np.pi/2,-np.pi*3/4]
names=['3W','6W','9W','12E/W','3E','6E','9E']
z0=0
tmax=1500
timescale=np.pi*8

for i in tqdm(range(0,len(deltap))):

    # Initial simulation
    #print p
    

    Y_VAL=[]

    for O in (O_range): 
        p=[4.5*D,3.5*D,D,O*D] # Notice that I've made the coefficient of D in the first term a variable.
        r = ode(f).set_integrator('zvode', method='bdf',with_jacobian=False)
        r.set_f_params(p)
        sol=[]
        r.set_initial_value(z0, t0)
        while r.successful() and r.t < t1:
            r.integrate(r.t+dt)
            sol=np.append(sol,r.y)
        ss=sol[-1]

        time=[]
        z_tau=[]
        z_tau=ss*(np.cos(deltap[i])+1j*np.sin(deltap[i])) # Setting initial condition by calculating z_st
        r = ode(f).set_integrator('zvode', method='bdf',with_jacobian=False)
        r.set_f_params(p)
        sol=[]
        r.set_initial_value(z_tau, t0)
        
        while r.successful() and r.t < tmax:
            r.integrate(r.t+dt)
            sol=np.append(sol,r.y)
            time=np.append(time,r.t+dt)
        #print sol
        y_data=[]

        for m in range(0,len(sol)):
            if np.abs(sol[m]-ss)>0.2:
            #y_data=np.append(y_data,np.abs(sol[i]-fp))
                y_data=np.append(y_data,time[m])
        
        Y_VAL=np.append(Y_VAL,y_data[-1])
    #TOTAL_Y_VAL=np.append(TOTAL_Y_VAL,Y_VAL)
    if deltap[i]<0:
        plt.plot(O_range,Y_VAL/timescale,label=names[i])
    else:
        plt.plot(O_range,Y_VAL/timescale,linestyle='dashed',label=names[i])

plt.ylabel('days')
plt.xlabel('O/Delta')
plt.xlim([-4,4])
plt.ylim([0,30])
plt.xticks(np.arange(-2,3,2))
plt.yticks(np.arange(0,35,5))
plt.legend()
#plt.show()
plt.savefig('/home/amogh/Documents/VTech/sem3GBCB/ModAndSim/Fig4_c_lu_oscillator.pdf')
