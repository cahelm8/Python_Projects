# Cody Helm

# Forward Dynamics Calculations

import numpy as np
import math
from pylab import *

# Integration using Euler
dur, h, m, r, g, I = 10.0, 0.0001, 1.65, 0.5, 9.81, 0.025
theta0 = 90 * math.pi / 180 # radians
theta_dot0 = 0 # rad/sec
t0 = 0 # intial time
Q = 0 # input torque
n = int(dur / h)

T, thetaN, theta_dotN = np.zeros(n), np.zeros(n), np.zeros(n)
thetaN1, theta_dotN1 = np.zeros(n), np.zeros(n)

for i in range(n):
    t1 = t0 + h
    theta1 = theta0 + h*(theta_dot0)
    theta_dot1 = theta_dot0 + h*((Q-m*g*r*math.sin(theta0))/(m*(r**2)+I))
    T[i], thetaN[i], theta_dotN[i] = t0, theta0, theta_dot0
    thetaN1[i],theta_dotN1[i]= theta1, theta_dot1
    theta0, theta_dot0, t0 = theta1, theta_dot1, t1
      
thetaN = thetaN * 180 / math.pi
theta_dotN = theta_dotN * 180 / math.pi

figure(1)
plot(T,thetaN, linestyle = '-', linewidth = 2.0, color = 'orange', label = '\u03B8(deg)')
plot(T,theta_dotN, linestyle = '-', linewidth = 2.0, color = 'cyan',label = "\u03B8'(deg/s)")
title('1-Link Arm (Euler)')
xlabel('Time (s)')
ylabel('States')
xticks(arange(0,12,2))
yticks(arange(-400,500,100))
legend(bbox_to_anchor=(1.05, 1))
show()


# Integration using odeint
from scipy.integrate import odeint

# function defining differential equation
def LinkArm(state, t):
    theta, theta_dot = state # unpack the state vector
    m, r, g, I, Q = 1.65, 0.5, 9.81, 0.025, 0
    theta_ddot = (Q-m*g*r*math.sin(theta))/(m*r**2+I) # compute acceleration xdd
    # return the two state derivatives
    return [theta_dot, theta_ddot]

# initial conditions
theta0 = 90 * math.pi / 180 # radians
theta_dot0 = 0 # rad/sec
state0 = np.array([theta0, theta_dot0])

# time vector
tstart, tend, timestep = 0.0, 10.0, 0.001
t = arange(tstart, tend, timestep)

state = odeint(LinkArm, state0, t)
thetaDeg = state[:,0] * 180 / math.pi
theta_dotDeg = state[:,1]* 180 / math.pi

figure(2)
plot(t,thetaDeg, linestyle = '-', linewidth = 2.0, color = 'orange', label = '\u03B8(deg)')
plot(t,theta_dotDeg, linestyle = '-', linewidth = 2.0, color = 'cyan', label = "\u03B8'(deg/s)")
title('1-Link Arm (ODEINT)')
xlabel('Time (s)')
ylabel('States')
xticks(arange(0,12,2))
yticks(arange(-400,500,100))
legend(bbox_to_anchor=(1.05, 1))
show()


# Energy
k = len(state[:,0])
thetaRad,theta_dotRad = np.zeros(k), np.zeros(k)
U,Trot,Tlin,Etot = np.zeros(k),np.zeros(k),np.zeros(k),np.zeros(k)
for i in range(k):
    thetaRad[i] = state[i,0]
    theta_dotRad[i] = state[i,1]
    U[i] = m*g*r*(1-math.cos(thetaRad[i]))
    Trot[i] = 1/2*I*theta_dotRad[i]**2
    Tlin[i] = 1/2*m*(r**2)*theta_dotRad[i]**2
    Etot[i] = Trot[i] + Tlin[i] + U[i]

figure(3)
plot(t,Trot, label = 'Trot')
plot(t,Tlin, label = 'Tlin')
plot(t,U, label = 'U')
plot(t,Etot, label = 'Etotal')
title('Kinetic and Potential Energy')
xlabel('Time (s)')
ylabel('Energy (J)')
legend(bbox_to_anchor=(1.05, 1))
show()


from scipy.integrate import odeint


def arm_2dof(state, t):
    # unpack the state vector
    theta1, theta2 = state[0], state[1]
    theta1dot, theta2dot = state[2], state[3]
    m1, m2, I1, I2, l1, l2 = 2.1, 1.65, 0.025, 0.075, 0.3384, 0.4554
    r1, r2, g, Q1, Q2 = 0.1692, 0.2277, 9.81, 0, 0
    M11 = I1 + I2 + m1*(r1**2) + m2*(l1**2 + r2**2 + 2*l1*r2*math.cos(theta2))
    M12 = I2 + m2*(r2**2 + l1*r2*math.cos(theta2))
    M21 = M12
    M22 = I2 + m2*(r2**2)
    M = np.array([[M11, M12], [M21, M22]])
    C1 = -l1*m2*r2*math.sin(theta2)*(theta2dot**2) - 2*l1*m2*r2*math.sin(theta2)*theta1dot*theta2dot
    C2 = l1*m2*r2*math.sin(theta2)*(theta1dot**2)
    C = np.array([[C1], [C2]])
    G1 = g*math.sin(theta1)*(m2*l1+m1*r1) + g*m2*r2*math.sin(theta1+theta2)
    G2 = g*m2*r2*math.sin(theta1 + theta2)
    G = np.array([[G1], [G2]])
    Q = np.array([[Q1], [Q2]])
    Minver = np.linalg.pinv(M)
    thetaddot = np.dot(Minver, (Q - C - G))
    theta1ddot = thetaddot[0,0]
    theta2ddot = thetaddot[1,0]
    # return the state derivatives
    return [theta1dot, theta2dot, theta1ddot, theta2ddot]

# initial conditions
theta1_0, theta2_0 = 180.0 * math.pi / 180, 1.0 * math.pi / 180
theta1dot_0, theta2dot_0 = 0.0 * math.pi / 180.0, 0.0 * math.pi / 180

state0 = np.array([theta1_0, theta2_0, theta1dot_0, theta2dot_0])

# time vector
tstart, tend, timestep = 0.0, 10.0, 0.001
t = arange(tstart, tend, timestep)

# differential equations
state = odeint(arm_2dof, state0, t)

theta1Deg = state[:,0] * 180 / math.pi
theta2Deg = state[:,1] * 180 / math.pi
theta1dotDeg = state[:,2] * 180 / math.pi
theta2dotDeg = state[:,3] * 180 / math.pi

plot(t,theta1Deg, label = '\u03B81(deg)')
plot(t,theta2Deg, label = '\u03B82(deg)')
plot(t,theta1dotDeg, label = "\u03B81'(deg/s)")
plot(t,theta2dotDeg, label = "\u03B82'(deg/s)")
title('2-Link Arm')
xlabel('Time (s)')
ylabel('States')
legend(bbox_to_anchor=(1.30, 1))
show()

# Energy
m1, m2, I1, I2, l1, l2 = 2.1, 1.65, 0.025, 0.075, 0.3384, 0.4554
r1, r2, g, Q1, Q2 = 0.1692, 0.2277, 9.81, 10, 10
n = len(state[:,0])
Ttotal, Utotal, Etotal = np.zeros(n), np.zeros(n), np.zeros(n)

for i in range(n):
    theta = state[i,0]
    alpha = state[i,1]
    thetad = state[i,2]
    alphad = state[i,3]

    x1 = r1 * math.sin(theta)
    y1 = -r1 * math.cos(theta)
    x2 = l1 * math.sin(theta) + r2 * math.sin(theta + alpha)
    y2 = -l1 * math.cos(theta) -r2 * math.cos(theta + alpha)
    
    x1d = r1*math.cos(theta)*thetad
    y1d = r1*math.sin(theta)*thetad
    x2d = l1*math.cos(theta)*thetad + r2*(alphad + thetad)*math.cos(theta + alpha)
    y2d = l1*math.sin(theta)*thetad + r2*(alphad + thetad)*math.sin(theta + alpha)

    Tlin1 = 1/2. * m1 * (x1d ** 2 + y1d ** 2)
    Tlin2 = 1/2. * m2 * (x2d ** 2 + y2d ** 2)
    Trot1 = 1/2. * I1 * (thetad) ** 2
    Trot2 = 1/2. * I2 * (thetad + alphad) ** 2
    Ttotal[i] = Tlin1 + Tlin2 + Trot1 + Trot2
    U1 = r1 * m1 * g * (1-math.cos(theta))
    U2 = l1 * m2 * g * (1-math.cos(theta)) + r2 * m2 * g * (1-math.cos(theta + alpha))
    Utotal[i] = U1 + U2
    Etotal[i] = Ttotal[i] + Utotal[i]
    

h = 0.001
axiscolor = 'k'
titlefsize = 32
axisfsize = 28
scalefsize = 24
fig, ax = plt.subplots(figsize=(8,8))
axis([0,n * h, 0 ,30])
plt.title('Kinetic and Potential Energy',fontsize = titlefsize)
plt.xlabel('Time (s)',fontsize = scalefsize)
plt.ylabel('Energy (J)',fontsize = scalefsize)
d1 = errorbar(t,Ttotal, linestyle = '-', linewidth = 2.0, color = 'blue')
d2 = errorbar(t,Utotal, linestyle = '-', linewidth = 2.0, color = 'green')
d3 = errorbar(t,Etotal, linestyle = '-', linewidth = 2.0, color = 'red')
lg = legend([d1,d2,d3], ['Lin','Potential','Total'])
legend(bbox_to_anchor=(1.05, 1))
show()




