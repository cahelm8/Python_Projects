# Cody Helm

import numpy as np
import math
from pylab import *

## Activation Dynamics
ss, dur = 0.001, 1.00 # step size (s), duration (s)
n = int(dur / ss) # number of steps
Tact, Tdeact = 0.015, 0.05
STIM0 = 0.0001
gamma0, t0 = 0, 0
T, STIM, GAMMA = np.zeros(n), np.zeros(n), np.zeros(n)

for i in arange (0,n,1):
    t1 = t0 + ss
    T = t0
    
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
    
    GAMMA[i] = ss*((STIM1-gamma0) / T) + gamma0
    STIM[i] = STIM1
    gamma0 = GAMMA[i]
    t0 = t1
    
T = arange(0.0,dur,ss)

plt.figure(1)
plot(T,GAMMA,label = '\u03B3', color= 'orange')
plot(T,STIM,label = 'STIM', color = 'blue')
xlabel('Time (sec)')
ylabel('Muscle Activity (unitless)')
title('Muscle Activation Dynamics')
xticks(arange(0,1.25,0.25))
yticks(arange(0,1.25,0.25))
legend()

## Active Force Length Curve
v = 0.185 # width parameter
ss, dur = 0.001, 2.0
n = int(dur/ss)
X = arange(0,dur,ss)
X0 = 1
alpha = np.zeros(n)

for i in arange(0,n,1):

    alpha[i] = math.exp(-(X[i]/X0-1)**2/v)

figure(2)
plot(X,alpha, color = 'blue')
xlabel('Normalized Muscle Fiber Length (X/X0)')
ylabel('Normalized Force')
title('Active Force-Length Curve')
xticks(arange(0,2.50,0.50))
yticks(arange(0,1.75,0.25))

## Passive Force-Length Curve
kPE, ePE = 5.0, 0.6
FPE = np.zeros(n)

for i in arange(0,n,1):

    FPE[i] = (math.exp(kPE*(X[i]/X0 - 1)/ePE) - 1) / (math.exp(kPE) - 1)
    
FPE = 1.5/max(FPE)*FPE
X = X - 0.5

figure(3)
plot(X,FPE, color = 'blue')
xlabel('Normalized Muscle Fiber Length (X/X0)')
ylabel('Normalized Force')
title('Passive Force-Length Curve')
xticks(arange(0,2.50,0.50))
yticks(arange(0,1.75,0.25))

## Force (x-axis)-Velocity Curve
Vmax = 10.0 # maximum shortening velcoty
gamma,alpha,FMFmax,Af = 1.0, 1.0, 1.4, 0.25
ss, dur = 0.001, 1.40
n = int(dur/ss) + 1
FMF = arange(0,dur,ss)
X_dot, b, line = np.zeros(n), np.zeros(n), np.zeros(n)

for i in arange(0,n,1):
    if FMF[i] > gamma*alpha:
        b[i] = ((2 + 2/Af)*(gamma*alpha*FMFmax - FMF[i])) / (FMFmax - 1)
    elif FMF[i] <= gamma*alpha:
        b[i] = gamma*alpha + (FMF[i]/Af)
                 
    X_dot[i] = (0.25 + 0.75*gamma)*Vmax*((FMF[i]-gamma*alpha)/b[i])

figure(4)
plot(FMF[0:1385],X_dot[0:1385], color = 'blue')
plot(FMF[0:1385],line[0:1385], color = 'grey')
xlabel('Normalized Force')
ylabel("Muscle Fiber Velocity (X'/X0)")
title('Force-Velocity Curve')
xticks(arange(0,1.75,0.25))
yticks(linspace(-10,10,5))

## Tendon Strain-Force Relationship
Ftoe,ktoe,e0 = 0.33,3.0,0.04
Eptoe,klin = 0.609*e0,1.712/e0
ss, dur = 0.001, 0.05
n = int(dur/ss)
Y = arange(0,dur,ss)
FT = np.zeros(n)

for i in arange(0,n,1):
    if Y[i] <= Eptoe:
        FT[i] = (Ftoe/(math.exp(ktoe)-1))*(math.exp(ktoe*(Y[i])/Eptoe)-1)
    elif Y[i] > Eptoe:
        FT[i] = (klin)*(Y[i] - Eptoe) + Ftoe
                 
figure(5)
plot(Y,FT, color = 'blue')
xlabel('Normalized Strain ((y-y0)/y0)')
ylabel('Normalized Force')
title('Tendon Force-Strain')
xticks(arange(0,0.06,0.01))
yticks(arange(0,1.75,0.25))

## Tendon Compliance Curve
ss, dur = 0.04, 2
n = int(dur/ss)
FT = arange(0,dur,ss)
Y = np.zeros(n)

for i in arange(0,n,1):
        Y[i] = 0.46962 / (57.835 * FT[i] + 1)
                 
figure(6)
plot(FT,Y, color = 'blue')
xlabel('Tendon Force')
ylabel('dY/d(FT)')
title('Compliance Function')
xticks(arange(0,2.5,0.5))
yticks(arange(0,0.75,0.25))

## Velocity (x-axis)-Force Curve
Vmax = 10.0 # maximum shortening velcoty
gamma,alpha,FMFmax,Af = 1.0, 1.0, 1.4, 0.25
ss, dur = 0.001, 1.4
n = int(dur/ss) + 1
FMF, line = np.zeros(n), np.zeros(n)

for i in arange(0,n,1):       
    FMF[i] = X_dot[i]*(1/((0.25 + 0.75*gamma)*Vmax))*(b[i]) + gamma*alpha
     

figure(7)
plot(X_dot[0:1385],FMF[0:1385], color = 'blue')
plot(line[0:1385],FMF[0:1385], color = 'grey')
xlabel("Muscle Fiber Velocity (X'/X0)")
ylabel('Normalized Force')
title('Force-Velocity Curve')
xticks(linspace(-10,10,5))
yticks(arange(0,1.75,0.25))




