from scipy.integrate import ode
import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
import random
from tqdm import tqdm

D=3.8*10**(-3)
# Parameter set for Type I dynamics
p=[4.5*D,3.5*D,D,1.4*D]
# Arbitrary initial condition to obtain steady state
z0=-1+0.5j
# system definition
def f(t, z, p):
    K,F,D,O=p
    s=0.5*((K*z+F)-z**2*np.conj(K*z+F))-(D+1j*O)*z
    return s

t0=0
t1 = 1000
dt = 0.1
fp=[]

# initial run to get fixed point

for i in tqdm(range(0,len(p))):
    #axarr[i].plot(circ.real,circ.imag,'r',linestyle='dashed')
    r = ode(f).set_integrator('zvode', method='bdf',with_jacobian=False)
    r.set_f_params(p)
#    print ("entering parameter set %i" %(i+1))
#    print ("Simulating initial condition %i" %(j+1))
    sol=[]
    r.set_initial_value(z0, t0)
    while r.successful() and r.t < t1:
        r.integrate(r.t+dt)
        sol=np.append(sol,r.y)
    fp=sol[-1]


dp=[np.pi/4,np.pi/2,np.pi*3/4,-np.pi,-np.pi/4,-np.pi/2,-np.pi*3/4]
names=['3W','6W','9W','12E/W','3E','6E','9E']
dt=0.1
tmax=350
timescale=np.pi*8

# Code to simulate till steady state for each initial condition, i.e. for varying time zone differences.
for i in range(0,len(dp)):
    time=[]
 
    z_tau=fp*(np.cos(dp[i])+1j*np.sin(dp[i])) # Setting initial condition
    r = ode(f).set_integrator('zvode', method='bdf',with_jacobian=False)
    r.set_f_params(p)
#    print ("entering parameter set %i" %(i+1))
#    print ("Simulating initial condition %i" %(j+1))
    sol=[]
    r.set_initial_value(z_tau, t0)
    while r.successful() and r.t < tmax:
        r.integrate(r.t+dt)
        sol=np.append(sol,r.y)
        time=np.append(time,r.t+dt)
    x_data=[]
    y_data=[]
    for j in range(0,len(sol)):
        if np.abs(sol[j]-fp)>0.2:
            y_data=np.append(y_data,np.abs(sol[j]-fp))
            x_data=np.append(x_data,time[j])
    if dp[i]<0:
        plt.plot(x_data/(timescale),y_data,label=names[i])
    else:
        plt.plot(x_data/(timescale),y_data,linestyle='dashed',label=names[i])


###########################################################
# some shoddy code to generate a straight line. ignore this
d=len(time)-len(sol)
thresh_T=[]
thresh=[]
for i in range(0,len(sol)):
    thresh=np.append(thresh,[0.2])
    thresh_T=np.append(thresh_T,time[d+i]/timescale)
plt.plot(thresh_T,thresh,'r',linestyle='dashed')
###########################################################
plt.ylim([0,2])
plt.legend()
plt.xlabel('days')
plt.ylabel('|z-zst|')
plt.show()
#plt.savefig('Fig3_lu_oscillator.pdf')
