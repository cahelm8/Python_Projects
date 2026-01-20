# Cody Helm


import numpy as np
import math
from pylab import *

X10 = 0.5; X0 = 0.5; Y0 = 0.5; LMT = 1.0; gamma0 = 0.0001; theta = 0
STIM0 = 0.0001;
Tact= 0.015; Tdeact = 0.05;
t0 = 0;
v = 0.185;
Ftoe = 0.33; ktoe = 3.0 ;e0 = 0.04;
Eptoe = 0.609*e0 ;klin = 1.712/e0;
kPE, ePE = 5.0, 0.6
Vmax,FMFmax, Af = 10.0, 1.4, 0.25

ss, dur = 0.00001, 1.00 # step size (s), duration (s)
n = int(dur / ss) # number of steps

T, STIM, GAMMA, = np.zeros(n), np.zeros(n), np.zeros(n)
Xs, Ys, = np.zeros(n), np.zeros(n)
FT, alpha = np.zeros(n), np.zeros(n)

for i in arange (0,n,1):
    t1 = t0 + ss
        
    if t0 < 0.1:
        STIM1 = STIM0
    elif t0 < 0.6:
        STIM1 = 1
    else:
        STIM1 = STIM0
             
    if STIM1 > gamma0:
        T = Tact*(0.5 + 1.5*gamma0)
    else:
        T = Tdeact/(0.5 + 1.5*gamma0)
        
    gamma1 = ss*((STIM1-gamma0) / T) + gamma0
        
    alpha0 = math.exp(-(X10/X0-1)**2/v)
    
    Y1 = (LMT - X10*math.cos(theta))
    
    if (Y1-Y0)/Y0 <= Eptoe:
        FT1 = (Ftoe/(math.exp(ktoe)-1))*(math.exp(ktoe*((Y1-Y0)/Y0)/Eptoe)-1)
    elif (Y1-Y0)/Y0 > Eptoe:
        FT1 = (klin)*(((Y1-Y0)/Y0) - Eptoe) + Ftoe
        
    FPE = (math.exp(kPE*(X10/X0 - 1)/ePE) - 1) / (math.exp(kPE) - 1)
    
    FMF = (FT1 - FPE*math.cos(theta)) / math.cos(theta)
    
    if FMF > gamma0*alpha0:
        b = ((2 + 2/Af)*(gamma0*alpha0*FMFmax - FMF)) / (FMFmax - 1)
    elif FMF <= gamma0*alpha0:
        b = gamma0*alpha0 + (FMF/Af)
                 
    X_dot = ss*((0.25 + 0.75*gamma0)*Vmax*((FMF-gamma0*alpha0)/b))*X0
    
    X1 = (X_dot) + X10
    
    GAMMA[i] = gamma0
    gamma0 = gamma1
    STIM[i] = STIM1
    alpha[i] = alpha0
    Xs[i] = X10
    X10 = X1
    Ys[i] = Y1
    FT[i] = FT1
    t0 = t1
    
T = arange(ss,dur,ss)

plt.figure(1)
plot(T,STIM,label = 'STIM', linestyle = ':', color = 'lightsteelblue')
plot(T,GAMMA,label = 'Activation Dynamics (\u03B3)',linestyle = '--',color= 'cornflowerblue')
plot(T,FT,label = 'Tendon Force',linestyle = '--', color = 'teal')
plot(T,Xs,label = 'Muscle Fiber Length',linestyle = '-', color = 'magenta')
plot(T,Ys,label = 'Tendon Length',linestyle = '-', color = 'red')
ylabel('Normalized Variables')
xlabel('Time (sec)')
title('Hill-Model')
xticks(arange(0,1.2,0.2))
yticks(arange(0,1.75,0.25))
legend(loc='upper left', shadow=True, prop={'size':7})
