from scipy.integrate import ode
import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
import random
from tqdm import tqdm

major_tics=np.arange(-1,1,1)
minor_tics=np.arange(-1,-1,0.5)
rcParams['figure.figsize'] = 4, 12
t0 =  0
def f(t, z, p):
    K,F,D,O=p
    s=0.5*((K*z+F)-z**2*np.conj(K*z+F))-(D+1j*O)*z
    return s

D=3.8*10**(-3)

p=[[4.5*D,3.5*D,D,1.4*D],[4.5*D,0.65*D,D,1.4*D],[10.0*D,3.5*D,D,1.4*D]]
z0=[[-1+0.5j],[-1-0.1j],[-1-0.5j]]
#z0=[[-1+0.5j,-1-0.1j,-1-0.5j,0.5+1j,0.9-1j,1+0.4j],[-1+0.5j,-1-0.1j,-1-0.5j,-0.1-0.1j,-0.11-0.12j,0.65+1j,1+0.45j,0.9+1j],[-1+0.5j,-1-0.1j,-1-0.25j,-1-0.5j,0.75-1j,-0.435-0.13j,-0.13j-0.41,-0.445-0.2j,-0.45-0.15j,-0.5-0.3j,-0.5-0.2j]]

t1 = 10000
dt = 0.1
fp=[]
plt.xlim=[-1.1,1.1]
plt.ylim=[-1.1,1.1]

print ("Computing Circle")
theta=np.arange(-np.pi,np.pi,np.pi/100)
circ=np.exp(1j*theta)

f1, axarr = plt.subplots(3)


for i in range(0,3):
    axarr[i].set_xlabel('Re(z)')
    axarr[i].set_ylabel('Im(z)')
    axarr[i].set_xticks(major_tics)
    axarr[i].set_xticks(minor_tics, minor=True)
    axarr[i].set_yticks(major_tics)
    axarr[i].set_yticks(minor_tics, minor=True)
    axarr[i].grid(which='both')                                                            

#fp stores the fixed points which are stable for plot a and c.
for i in tqdm(range(0,len(p))):
    axarr[i].plot(circ.real,circ.imag,'r',linestyle='dashed')
    r = ode(f).set_integrator('zvode', method='bdf',with_jacobian=False)
    r.set_f_params(p[i])
    print ("entering parameter set %i" %(i+1))
    for j in tqdm(range(0,len(z0[i]))):
        print ("Simulating initial condition %i" %(j+1))
        sol=[]
        r.set_initial_value(z0[i][j], t0)
        while r.successful() and r.t < t1:
            r.integrate(r.t+dt)
            sol=np.append(sol,r.y)
        a=np.ceil(10*random.random())/10
        b=np.ceil(10*random.random())/10
        c=np.ceil(10*random.random())/10
        axarr[i].plot(sol.real,sol.imag,color=(a,b,c))
        
    fp=np.append(fp,sol[-1])

axarr[0].plot(fp.real[0],fp.imag[0],'ko')
axarr[2].plot(fp.real[2],fp.imag[2],'ko')
#axarr[2].plot(-0.43,-0.2,'wo') # generates the point close to unstable fixed point for plot c, not sure what the real value is.

plt.show()

#plt.savefig("Fig1_c_lu_oscillator.pdf",bbox_inches='tight')
