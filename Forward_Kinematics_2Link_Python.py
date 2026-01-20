# Cody Helm
# Forward kinematics for a 2link arm

import numpy as np
import math

# Given
l1 = 0.34;l2 = 0.46
# Hand Position
Hx = 0.36;Hy = 0.65
H = np.array([[Hx],[Hy]])
Hdot = np.array([[-3.89],[1.30]])
Hddot = np.array([[-7.79],[-26.18]])

# Calculating Thetas
c = math.sqrt(Hx**2 + Hy**2)
beta = math.acos((c**2 - l1**2 - l2**2) / (-2*l1*l2)) * 180 / math.pi
Theta2 = 180 - beta
alpha = math.acos((l2**2 - l1**2 - c**2) / (-2*l1*c)) * 180 / math.pi
gamma = math.atan2(Hy,Hx)  * 180 / math.pi
Theta1 = gamma - alpha
Theta = np.array([[Theta1],[Theta2]]) * math.pi / 180

# Joint Angle
Angle = Theta * 180 / math.pi
print('Joint Angle = ',Angle,'deg')

# Jacobian
J11 = -l1 * math.sin(Theta[0,0]) - l2 * math.sin(Theta[0,0] + Theta[1,0])
J12 = -l2 * math.sin(Theta[0,0] + Theta[1,0])
J21 = l1 * math.cos(Theta[0,0]) + l2 * math.cos(Theta[0,0] + Theta[1,0])
J22 = l2 * math.cos(Theta[0,0] + Theta[1,0])
J = np.array([[J11,J12],[J21,J22]])

# Angular Velocity
Jinv = np.linalg.pinv(J)
AngVel = np.dot(Jinv,Hdot)
print('Angular Velocity = ',AngVel,'rad/s')

# Jacobian Dot
Jdot11 = -l1 * math.cos(Theta[0,0]) * AngVel[0,0] - l2 * math.cos(Theta[0,0] + Theta[1,0]) * (AngVel[0,0] + AngVel[1,0])
Jdot12 = -l2 * math.cos(Theta[0,0] + Theta[1,0]) * (AngVel[0,0] + AngVel[1,0])
Jdot21 = -l1 * math.sin(Theta[0,0]) * AngVel[0,0] - l2 * math.sin(Theta[0,0] + Theta[1,0]) * (AngVel[0,0] + AngVel[1,0])
Jdot22 = -l2 * math.sin(Theta[0,0] + Theta[1,0]) * (AngVel[0,0] + AngVel[1,0])
Jdot = np.array([[Jdot11,Jdot12],[Jdot21,Jdot22]])

# Angular Acceleration
AngAcc = np.dot(Jinv,(Hddot - np.dot(Jdot,AngVel)))
print('Angular Acceleration = ',AngAcc,'rad/s^2') 



