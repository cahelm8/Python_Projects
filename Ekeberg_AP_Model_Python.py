import numpy, math, matplotlib
from pylab import *

ss, dur = 0.00001, 0.2 # step size (s), duration (s)
n = int(dur / ss) # number of steps
E_soma0, m0, h0, n0, t0 = -70*10**(-3), 0, 1, 0, 0 # Initial Conditions
E_SOMA, M, H, N, T = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
I_ext = 1.0*10**(-10) #I_ext = 100.0e-10

for i in arange (0,n,1):
    t1 = t0 + ss
    T[i] = t0
    
    # sodium activation channel
    # define alpha_m term
    A_alpha_m, B_alpha_m, C_alpha_m = 2.0*10**5, -4.0*10**-2, 1.0*10**(-3)
    alpha_m = (A_alpha_m*(E_soma0-B_alpha_m)) / (1-exp((B_alpha_m-E_soma0)/C_alpha_m))
    # define beta_m term
    A_beta_m, B_beta_m, C_beta_m = 6.0*10**(4), -4.9*10**(-2), 2.0*10**(-2)
    beta_m = (A_beta_m*(B_beta_m-E_soma0)) / (1-exp((E_soma0-B_beta_m)/C_beta_m))
    # define Na activation channel - differential equation here
    M[i] = m0
    m1 = m0 + ss*(alpha_m*(1-m0)-beta_m*m0)
    
    # sodium inactivation channel
    # define alpha_h term
    A_alpha_h, B_alpha_h, C_alpha_h = 8.0*10**(4), -4.0*10**(-2), 1.0*10**(-3)
    alpha_h = (A_alpha_h*(B_alpha_h-E_soma0)) / (1-exp((E_soma0-B_alpha_h)/C_alpha_h))
    # define beta_h term
    A_beta_h, B_beta_h, C_beta_h = 4.0*10**(2), -3.6*10**(-2), 2.0*10**(-3)
    beta_h = (A_beta_h) / (1+exp((B_beta_h-E_soma0)/C_beta_h))
    # define Na inactivation channel - differential equation here
    H[i] = h0 
    h1 = h0 + ss*(alpha_h*(1-h0)-beta_h*h0)
    
    # define NA channel current
    E_Na, G_Na = 5.0*10**(-2), 1.0*10**(-6)
    I_Na = (E_Na - E_soma0)*G_Na*(m0**3)*h0
    
    # potassium channel
    # define alpha_n term
    A_alpha_n, B_alpha_n, C_alpha_n = 2.0*10**(4), -3.1*10**(-2), 8.0*10**(-4)
    alpha_n = (A_alpha_n*(E_soma0-B_alpha_n)) / (1-exp((B_alpha_n-E_soma0)/C_alpha_n))
    # define beta_n term
    A_beta_n, B_beta_n, C_beta_n = 5.0*10**(3), -2.8*10**(-2), 4.0*10**(-4)
    beta_n = (A_beta_n*(B_beta_n-E_soma0)) / (1-exp((E_soma0-B_beta_n)/C_beta_n))
    # define K channel - differential equation
    N[i] = n0
    n1 = n0 + ss*(alpha_n*(1-n0) - beta_n*n0)
    
    # define K channel current
    E_K, G_K = -9.0*10**(-2), 2.0*10**(-7)
    I_K = (E_K - E_soma0)*G_K*(n0**4)
       
    # sum current from all channels (e.g., I_Na + I_K)
    # I = I_Na + I_K + I_Ca + I_K_Ca + I_ext
    I = I_Na + I_K + I_ext
    
    # leak current
    E_leak, G_m, C_m = -7.0*10**(-2), 3.0*10**(-9), 3.0*10**(-11)
    # define (E_soma) memberane current - differential equation
    E_SOMA[i] = E_soma0
    E_soma1 = E_soma0 + ss*(((E_leak - E_soma0)*G_m + I) / C_m)
    
    # update differential equation terms from n+1 to n
    t0, m0, n0, h0, E_soma0 = t1, m1, n1, h1, E_soma1


plt.figure(1)
plot(T,M,label = "Nam")
xlabel("Time (sec)")
ylabel("Na Activation Channels")

plt.figure(2)
plot(T,H,label = "Nah")
xlabel("Time (sec)")
ylabel("Na Inactivation Channels")

plt.figure(3)
plot(T,N,label = "K")
xlabel("Time (sec)")
ylabel("Potassium Channels")

plt.figure(4)
plot(T,E_SOMA,label ="ESOMA")
xlabel("Time (sec)")
ylabel("Membrane Potential (V)")


plt.figure(5)
plot(M,H)
xlabel('Na Activation')
ylabel('Na Inactivation')
plt.figure(6)
plot(M,N)
xlabel('Na Activation')
ylabel('Potassium')
plt.figure(7)
plot(M,E_SOMA)
xlabel('Na Activation')
ylabel('Membrane Potential (V)')
plt.figure(8)
plot(H,N)
xlabel('Na Inactivation')
ylabel('Potassium')
plt.figure(9)
plot(H,E_SOMA)
xlabel('Na Inactivation')
ylabel('Membrane Potential (V)')
plt.figure(10)
plot(N,E_SOMA)
xlabel('Potassium')
ylabel('Membrane Potential (V)')




from scipy.integrate import odeint
def ActionPotential(state,t):
    E_soma0, m0, h0, n0 = state
    # sodium activation channel
    # define alpha_m term
    A_alpha_m, B_alpha_m, C_alpha_m = 2.0*10**5, -4.0*10**-2, 1.0*10**(-3)
    alpha_m = (A_alpha_m*(E_soma0-B_alpha_m)) / (1-exp((B_alpha_m-E_soma0)/C_alpha_m))
    # define beta_m term
    A_beta_m, B_beta_m, C_beta_m = 6.0*10**(4), -4.9*10**(-2), 2.0*10**(-2)
    beta_m = (A_beta_m*(B_beta_m-E_soma0)) / (1-exp((E_soma0-B_beta_m)/C_beta_m))
    # define Na activation channel - differential equation here
    m1 = (alpha_m*(1-m0)-beta_m*m0)
    # sodium inactivation channel
    # define alpha_h term
    A_alpha_h, B_alpha_h, C_alpha_h = 8.0*10**(4), -4.0*10**(-2), 1.0*10**(-3)
    alpha_h = (A_alpha_h*(B_alpha_h-E_soma0)) / (1-exp((E_soma0-B_alpha_h)/C_alpha_h))
    # define beta_h term
    A_beta_h, B_beta_h, C_beta_h = 4.0*10**(2), -3.6*10**(-2), 2.0*10**(-3)
    beta_h = (A_beta_h) / (1+exp((B_beta_h-E_soma0)/C_beta_h))
    # define Na inactivation channel - differential equation here
    h1 = (alpha_h*(1-h0)-beta_h*h0)
    
    # define NA channel current
    E_Na, G_Na = 5.0*10**(-2), 1.0*10**(-6)
    I_Na = (E_Na - E_soma0)*G_Na*(m0**3)*h0
    
    # potassium channel
    # define alpha_n term
    A_alpha_n, B_alpha_n, C_alpha_n = 2.0*10**(4), -3.1*10**(-2), 8.0*10**(-4)
    alpha_n = (A_alpha_n*(E_soma0-B_alpha_n)) / (1-exp((B_alpha_n-E_soma0)/C_alpha_n))
    # define beta_n term
    A_beta_n, B_beta_n, C_beta_n = 5.0*10**(3), -2.8*10**(-2), 4.0*10**(-4)
    beta_n = (A_beta_n*(B_beta_n-E_soma0)) / (1-exp((E_soma0-B_beta_n)/C_beta_n))
    # define K channel - differential equation
    n1 = (alpha_n*(1-n0) - beta_n*n0)
    
    # define K channel current
    E_K, G_K = -9.0*10**(-2), 2.0*10**(-7)
    I_K = (E_K - E_soma0)*G_K*(n0**4)
       
    # sum current from all channels (e.g., I_Na + I_K)
    # I = I_Na + I_K + I_Ca + I_K_Ca + I_ext
    I = I_Na + I_K + I_ext
    
    # leak current
    E_leak, G_m, C_m = -7.0*10**(-2), 3.0*10**(-9), 3.0*10**(-11)
    # define (E_soma) memberane current - differential equation
    E_soma1 = ((E_leak - E_soma0)*G_m + I) / C_m
    
    
    return [E_soma1,m1,h1,n1]
    
# initial conditions
state0 = np.array([-70*10**(-3), 0.0, 1.0, 0.0])

# time vector
tstart, tend, timestep = 0.0, dur, ss
t = arange(tstart, tend, timestep)

state = odeint(ActionPotential, state0, t)

plt.figure(11)
plot(t,state[:,0])
xlabel('Time (sec)')
ylabel('Membrane Potential (V)')
legend()
 
    
    
    
    
