# Cody Helm
# Forward kinematics for a 1 link arm

import numpy as np
import math

# Define the variables and constant values
I = 0.46 # meters
Theta = 60.0 * math.pi / 180. # convert from degrees to radians
Theta_dot = np.array([[2.0]]) # rad/s
Theta_ddot = np.array([[1.0]]) # rad/s^2

## Problem 1
# Hand Position
hx = I * math.cos(Theta)
hy = I * math.sin(Theta)
H = np.array([[hx],[hy]])
print('Hand Position = ',H,'m')

# %%
# Jacobian (how does hand position change with a change in joint angle)
J11 = -I * math.sin(Theta) # first derivative of the hx
J21 = I * math.cos(Theta) # first derivative of the hy
J = np.array([[J11],[J21]])

# Hand Velocity (change in H velocity with a change in angle times the velocity)
Hdot = np.dot(J,Theta_dot) # np.dot does matrix multiplication
print('Hand Velocity = ',Hdot,'m/s')




# Jacobian Dot
Jdot11 = -I * math.cos(Theta) * Theta_dot[0,0] # [0,0] extract element from 1x1 array, first derivative of jacobian
Jdot21 = -I * math.sin(Theta) * Theta_dot[0,0] # first derivative of the jacobian
Jdot = np.array([[Jdot11],[Jdot21]])

# Hand Acceleration
Hddot = np.dot(Jdot,Theta_dot) + np.dot(J,Theta_ddot) # 
print('Hand Acceleration = ',Hddot,'m/s^2')



# %%
## Problem 2 - going from hand position to joint angle
# Joint Angle
Angle = np.array([[math.atan2(H[1,0],H[0,0]) * 180 / math.pi]])
print('Joint Angle = ',Angle,'deg')

# Angular Velocity
Jinv = np.linalg.pinv(J)
AngVel = np.dot(Jinv,Hdot)
print('Angular Velocity = ',AngVel,'rad/s') 

# Angular Acceleration
AngAcc = np.dot(Jinv,(Hddot - np.dot(Jdot,AngVel)))
print('Angular Acceleration = ',AngAcc,'rad/s^2') 



