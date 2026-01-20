# Cody Helm

# Distribution-Moment Approximation (Zahalak 1981)

import numpy, math, matplotlib
from pylab import *

# Shortening Velocity
h, s = 0.0005, 0.8
n = int(s / h)
f1, g1, g2, g3 = 43.4, 10.0, 209.0, 0.0# fill in values
Q0, Q1, Q2 =  0.684, 0.438, 0.322# fill in values
Q0IC = Q0
Q1IC = Q1
Q2IC = Q2

STIFFNESS, FORCE, ENERGY = np.zeros(n), np.zeros(n), np.zeros(n)
TIME = np.zeros(n)
VEL = np.zeros(n)
B0, B1, B2 = f1 / 2.0, f1 / 3.0, f1 / 4.0 
t0 = 0

for i in arange(0,n,1):
    ########################################## initial conditions  ########################################
    tn = t0 + i*h
    v0 = -10 #fill in value
    p, q = Q1 / Q0, (Q2 / Q0 - (Q1 / Q0)**2)**0.5
    ta, tb = ((1 - p) / q), -p / q
    Pa, Pb = 0.5 * (1 + math.erf((ta) / (2**0.5))), 0.5 * (1 + math.erf((tb) / (2**0.5)))
    J0a = Pa
    J0b = Pb
    J1a = p * Pa - q * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5
    J1b = p * Pb - q * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5
    J2a = p**2 * Pa - 2 * p * q * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5 + q**2 * (Pa - ta * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5)
    J2b = p**2 * Pb - 2 * p * q * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5 + q**2 * (Pb - tb * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5)
    J3a = p**3 * Pa - 3 * p**2 * q * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5 + 3 * p * q**2 * (Pa - ta * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5)\
        - q**3 * (2 + ta**2) * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5
    J3b = p**3 * Pb - 3 * p**2 * q * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5 + 3 * p * q**2 * (Pb - tb * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5)\
        - q**3 * (2 + tb**2) * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5
    p0 = Q0 * (g2 * J0b + (f1 + g1) * (J1a - J1b) + g1 * (p - J1a) + g3 * (p - J1a - 1 + J0a))
    p1 = Q0 * (g2 * J1b + (f1 + g1) * (J2a - J2b) + g1 * (p**2 + q**2 - J2a) + g3 * (p**2 + q**2 - J2a - p + J1a))
    p2 = Q0 * (g2 * J2b + (f1 + g1) * (J3a - J3b) + g1 * (p**3 + 3 * p * q**2 - J3a) + g3 * (p**3 + 3 * p * q**2 - J3a - (p**2 + q**2) + J2a))
    Q01 = Q0 + (B0 - p0) * h
    Q11 = Q1 + (B1 - p1 - (v0) * Q0) * h
    Q21 = Q2 + (B2 - p2 - 2 * (v0) * Q1) * h
    STIFFNESS[i] = (Q0/Q0IC)
    FORCE[i] = (Q1/Q1IC) #the denominator the same value as the initial condition of Q1
    ENERGY[i] = (Q2/Q2IC)
    TIME[i] = tn
    VEL[i] = v0
    Q0, Q1, Q2 = Q01, Q11, Q21
    #print(tn)

figure(1)
plot(TIME,STIFFNESS,linestyle = '--', color = 'lime',label = 'Stiffness (Q0)')
plot(TIME,FORCE,linestyle = '--', color = 'teal',label = 'Force (Q1)')
plot(TIME,ENERGY,linestyle = '--', color = 'royalblue',label = 'Energy (Q2)')
ylabel('Normalized Force')
xlabel('Time(sec)')
title('DMA - Model (Lengthening)')
xticks(arange(0,1.0,0.2))
yticks(arange(0,6.5,0.5))
legend()

# Lengthening Velocity
h, s = 0.0005, 0.8
n = int(s / h)
f1, g1, g2, g3 = 43.4, 10.0, 209.0, 0.0# fill in values
Q0, Q1, Q2 =  0.684, 0.438, 0.322# fill in values
Q1IC = Q1
Q1IC = Q1
Q2IC = Q2

STIFFNESS, FORCE, ENERGY = np.zeros(n), np.zeros(n), np.zeros(n)
TIME = np.zeros(n)
VEL = np.zeros(n)
B0, B1, B2 = f1 / 2.0, f1 / 3.0, f1 / 4.0 
t0 = 0

for i in arange(0,n,1):
    ########################################## initial conditions  ########################################
    tn = t0 + i*h
    v0 = 10 #fill in value
    p, q = Q1 / Q0, (Q2 / Q0 - (Q1 / Q0)**2)**0.5
    ta, tb = ((1 - p) / q), -p / q
    Pa, Pb = 0.5 * (1 + math.erf((ta) / (2**0.5))), 0.5 * (1 + math.erf((tb) / (2**0.5)))
    J0a = Pa
    J0b = Pb
    J1a = p * Pa - q * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5
    J1b = p * Pb - q * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5
    J2a = p**2 * Pa - 2 * p * q * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5 + q**2 * (Pa - ta * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5)
    J2b = p**2 * Pb - 2 * p * q * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5 + q**2 * (Pb - tb * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5)
    J3a = p**3 * Pa - 3 * p**2 * q * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5 + 3 * p * q**2 * (Pa - ta * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5)\
        - q**3 * (2 + ta**2) * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5
    J3b = p**3 * Pb - 3 * p**2 * q * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5 + 3 * p * q**2 * (Pb - tb * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5)\
        - q**3 * (2 + tb**2) * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5
    p0 = Q0 * (g2 * J0b + (f1 + g1) * (J1a - J1b) + g1 * (p - J1a) + g3 * (p - J1a - 1 + J0a))
    p1 = Q0 * (g2 * J1b + (f1 + g1) * (J2a - J2b) + g1 * (p**2 + q**2 - J2a) + g3 * (p**2 + q**2 - J2a - p + J1a))
    p2 = Q0 * (g2 * J2b + (f1 + g1) * (J3a - J3b) + g1 * (p**3 + 3 * p * q**2 - J3a) + g3 * (p**3 + 3 * p * q**2 - J3a - (p**2 + q**2) + J2a))
    Q01 = Q0 + (B0 - p0) * h
    Q11 = Q1 + (B1 - p1 - (v0) * Q0) * h
    Q21 = Q2 + (B2 - p2 - 2 * (v0) * Q1) * h
    STIFFNESS[i] = (Q0/Q0IC)
    FORCE[i] = (Q1/Q1IC) #the denominator the same value as the initial condition of Q1
    ENERGY[i] = (Q2/Q2IC)
    TIME[i] = tn
    VEL[i] = v0
    Q0, Q1, Q2 = Q01, Q11, Q21
    #print(tn)

figure(2)
plot(TIME,STIFFNESS,linestyle = '--', color = 'lime',label = 'Stiffness (Q0)')
plot(TIME,FORCE,linestyle = '--', color = 'teal',label = 'Force (Q1)')
plot(TIME,ENERGY,linestyle = '--', color = 'royalblue',label = 'Energy (Q2)')
ylabel('Normalized Force')
xlabel('Time(sec)')
title('DMA - Model (Shortening)')
xticks(arange(0,1.0,0.2))
yticks(arange(0,3.5,0.5))
legend()

# Oscillating Muscle (T(0.05))
h, s = 0.0005, 0.8
n = int(s / h)
f1, g1, g2, g3 = 1.0, 10.0, 210.0, 100.0# fill in values
Q0, Q1, Q2 =  0.0712, 0.0436, 0.0295# fill in values
Q1IC = Q1
Q1IC = Q1
Q2IC = Q2

STIFFNESS, FORCE, ENERGY = np.zeros(n), np.zeros(n), np.zeros(n)
TIME = np.zeros(n)
VEL = np.zeros(n)
B0, B1, B2 = f1 / 2.0, f1 / 3.0, f1 / 4.0 
t0 = 0
T = [0.05]

for i in arange(0,n,1):
    ########################################## initial conditions  ########################################
    tn = t0 + i*h
    v0 = -25*sin(2*pi*tn/T) #fill in value
    p, q = Q1 / Q0, (Q2 / Q0 - (Q1 / Q0)**2)**0.5
    ta, tb = ((1 - p) / q), -p / q
    Pa, Pb = 0.5 * (1 + math.erf((ta) / (2**0.5))), 0.5 * (1 + math.erf((tb) / (2**0.5)))
    J0a = Pa
    J0b = Pb
    J1a = p * Pa - q * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5
    J1b = p * Pb - q * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5
    J2a = p**2 * Pa - 2 * p * q * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5 + q**2 * (Pa - ta * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5)
    J2b = p**2 * Pb - 2 * p * q * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5 + q**2 * (Pb - tb * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5)
    J3a = p**3 * Pa - 3 * p**2 * q * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5 + 3 * p * q**2 * (Pa - ta * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5)\
        - q**3 * (2 + ta**2) * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5
    J3b = p**3 * Pb - 3 * p**2 * q * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5 + 3 * p * q**2 * (Pb - tb * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5)\
        - q**3 * (2 + tb**2) * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5
    p0 = Q0 * (g2 * J0b + (f1 + g1) * (J1a - J1b) + g1 * (p - J1a) + g3 * (p - J1a - 1 + J0a))
    p1 = Q0 * (g2 * J1b + (f1 + g1) * (J2a - J2b) + g1 * (p**2 + q**2 - J2a) + g3 * (p**2 + q**2 - J2a - p + J1a))
    p2 = Q0 * (g2 * J2b + (f1 + g1) * (J3a - J3b) + g1 * (p**3 + 3 * p * q**2 - J3a) + g3 * (p**3 + 3 * p * q**2 - J3a - (p**2 + q**2) + J2a))
    Q01 = Q0 + (B0 - p0) * h
    Q11 = Q1 + (B1 - p1 - (v0) * Q0) * h
    Q21 = Q2 + (B2 - p2 - 2 * (v0) * Q1) * h
    STIFFNESS[i] = (Q0/Q0IC)
    FORCE[i] = (Q1/Q1IC) #the denominator the same value as the initial condition of Q1
    ENERGY[i] = (Q2/Q2IC)
    TIME[i] = tn
    VEL[i] = v0
    Q0, Q1, Q2 = Q01, Q11, Q21
    #print(tn)

figure(3)
plot(TIME,STIFFNESS,linestyle = '--', color = 'lime',label = 'Stiffness (Q0)')
plot(TIME,FORCE,linestyle = '--', color = 'teal',label = 'Force (Q1)')
plot(TIME,ENERGY,linestyle = '--', color = 'royalblue',label = 'Energy (Q2)')
ylabel('Normalized Force')
xlabel('Time(sec)')
title('DMA - Model (T=0.05)')
xticks(arange(0,1.0,0.2))
yticks(arange(0,3.5,0.5))
legend()

# Oscillating Muscle (T(0.1))
h, s = 0.0005, 0.8
n = int(s / h)
f1, g1, g2, g3 = 1.0, 10.0, 210.0, 100.0# fill in values
Q0, Q1, Q2 =  0.0712, 0.0436, 0.0295# fill in values
Q1IC = Q1
Q1IC = Q1
Q2IC = Q2

STIFFNESS, FORCE, ENERGY = np.zeros(n), np.zeros(n), np.zeros(n)
TIME = np.zeros(n)
VEL = np.zeros(n)
B0, B1, B2 = f1 / 2.0, f1 / 3.0, f1 / 4.0 
t0 = 0
T = [0.1]

for i in arange(0,n,1):
    ########################################## initial conditions  ########################################
    tn = t0 + i*h
    v0 = -25*sin(2*pi*tn/T) #fill in value
    p, q = Q1 / Q0, (Q2 / Q0 - (Q1 / Q0)**2)**0.5
    ta, tb = ((1 - p) / q), -p / q
    Pa, Pb = 0.5 * (1 + math.erf((ta) / (2**0.5))), 0.5 * (1 + math.erf((tb) / (2**0.5)))
    J0a = Pa
    J0b = Pb
    J1a = p * Pa - q * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5
    J1b = p * Pb - q * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5
    J2a = p**2 * Pa - 2 * p * q * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5 + q**2 * (Pa - ta * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5)
    J2b = p**2 * Pb - 2 * p * q * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5 + q**2 * (Pb - tb * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5)
    J3a = p**3 * Pa - 3 * p**2 * q * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5 + 3 * p * q**2 * (Pa - ta * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5)\
        - q**3 * (2 + ta**2) * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5
    J3b = p**3 * Pb - 3 * p**2 * q * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5 + 3 * p * q**2 * (Pb - tb * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5)\
        - q**3 * (2 + tb**2) * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5
    p0 = Q0 * (g2 * J0b + (f1 + g1) * (J1a - J1b) + g1 * (p - J1a) + g3 * (p - J1a - 1 + J0a))
    p1 = Q0 * (g2 * J1b + (f1 + g1) * (J2a - J2b) + g1 * (p**2 + q**2 - J2a) + g3 * (p**2 + q**2 - J2a - p + J1a))
    p2 = Q0 * (g2 * J2b + (f1 + g1) * (J3a - J3b) + g1 * (p**3 + 3 * p * q**2 - J3a) + g3 * (p**3 + 3 * p * q**2 - J3a - (p**2 + q**2) + J2a))
    Q01 = Q0 + (B0 - p0) * h
    Q11 = Q1 + (B1 - p1 - (v0) * Q0) * h
    Q21 = Q2 + (B2 - p2 - 2 * (v0) * Q1) * h
    STIFFNESS[i] = (Q0/Q0IC)
    FORCE[i] = (Q1/Q1IC) #the denominator the same value as the initial condition of Q1
    ENERGY[i] = (Q2/Q2IC)
    TIME[i] = tn
    VEL[i] = v0
    Q0, Q1, Q2 = Q01, Q11, Q21
    #print(tn)

figure(4)
plot(TIME,STIFFNESS,linestyle = '--', color = 'lime',label = 'Stiffness (Q0)')
plot(TIME,FORCE,linestyle = '--', color = 'teal',label = 'Force (Q1)')
plot(TIME,ENERGY,linestyle = '--', color = 'royalblue',label = 'Energy (Q2)')
ylabel('Normalized Force')
xlabel('Time(sec)')
title('DMA - Model (T=0.1)')
xticks(arange(0,1.0,0.2))
yticks(arange(0,3.5,0.5))
legend()


# Oscillating Muscle (T(0.2))
h, s = 0.0005, 0.8
n = int(s / h)
f1, g1, g2, g3 = 1.0, 10.0, 210.0, 100.0# fill in values
Q0, Q1, Q2 =  0.0712, 0.0436, 0.0295# fill in values
Q1IC = Q1
Q1IC = Q1
Q2IC = Q2

STIFFNESS, FORCE, ENERGY = np.zeros(n), np.zeros(n), np.zeros(n)
TIME = np.zeros(n)
VEL = np.zeros(n)
B0, B1, B2 = f1 / 2.0, f1 / 3.0, f1 / 4.0 
t0 = 0
T = [0.2]

for i in arange(0,n,1):
    ########################################## initial conditions  ########################################
    tn = t0 + i*h
    v0 = -25*sin(2*pi*tn/T) #fill in value
    p, q = Q1 / Q0, (Q2 / Q0 - (Q1 / Q0)**2)**0.5
    ta, tb = ((1 - p) / q), -p / q
    Pa, Pb = 0.5 * (1 + math.erf((ta) / (2**0.5))), 0.5 * (1 + math.erf((tb) / (2**0.5)))
    J0a = Pa
    J0b = Pb
    J1a = p * Pa - q * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5
    J1b = p * Pb - q * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5
    J2a = p**2 * Pa - 2 * p * q * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5 + q**2 * (Pa - ta * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5)
    J2b = p**2 * Pb - 2 * p * q * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5 + q**2 * (Pb - tb * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5)
    J3a = p**3 * Pa - 3 * p**2 * q * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5 + 3 * p * q**2 * (Pa - ta * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5)\
        - q**3 * (2 + ta**2) * math.exp(-ta**2 / 2.0) / (2 * math.pi)**0.5
    J3b = p**3 * Pb - 3 * p**2 * q * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5 + 3 * p * q**2 * (Pb - tb * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5)\
        - q**3 * (2 + tb**2) * math.exp(-tb**2 / 2.0) / (2 * math.pi)**0.5
    p0 = Q0 * (g2 * J0b + (f1 + g1) * (J1a - J1b) + g1 * (p - J1a) + g3 * (p - J1a - 1 + J0a))
    p1 = Q0 * (g2 * J1b + (f1 + g1) * (J2a - J2b) + g1 * (p**2 + q**2 - J2a) + g3 * (p**2 + q**2 - J2a - p + J1a))
    p2 = Q0 * (g2 * J2b + (f1 + g1) * (J3a - J3b) + g1 * (p**3 + 3 * p * q**2 - J3a) + g3 * (p**3 + 3 * p * q**2 - J3a - (p**2 + q**2) + J2a))
    Q01 = Q0 + (B0 - p0) * h
    Q11 = Q1 + (B1 - p1 - (v0) * Q0) * h
    Q21 = Q2 + (B2 - p2 - 2 * (v0) * Q1) * h
    STIFFNESS[i] = (Q0/Q0IC)
    FORCE[i] = (Q1/Q1IC) #the denominator the same value as the initial condition of Q1
    ENERGY[i] = (Q2/Q2IC)
    TIME[i] = tn
    VEL[i] = v0
    Q0, Q1, Q2 = Q01, Q11, Q21
    #print(tn)

figure(5)
plot(TIME,STIFFNESS,linestyle = '--', color = 'lime',label = 'Stiffness (Q0)')
plot(TIME,FORCE,linestyle = '--', color = 'teal',label = 'Force (Q1)')
plot(TIME,ENERGY,linestyle = '--', color = 'royalblue',label = 'Energy (Q2)')
ylabel('Normalized Force')
xlabel('Time(sec)')
title('DMA - Model (T=0.2)')
xticks(arange(0,1.0,0.2))
yticks(arange(0,3.5,0.5))
legend()
