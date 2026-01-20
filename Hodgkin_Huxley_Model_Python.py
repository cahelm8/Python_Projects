import numpy, math, matplotlib
from pylab import *

ss, dur = 0.00001, 0.2 # step size (s), duration (s)
n = int(dur / ss) # number of steps
E_soma0, m0, h0, n0, q0, CaAP0, t0 = -70*10**(-3), 0, 1, 0, 0, 0, 0 # Initial Conditions
E_SOMA, M, H, N, Q, CAAP, T = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
I_ext = 2*10**(-9) #I_ext = 100.0e-10

for i in arange (0,n,1):
    t1 = t0 + ss
    T[i] = t0
    # sodium activation channel
    # define alpha_m term
    A_alpha_m, B_alpha_m, C_alpha_m = 2.0*10**(+5), -4.0*10**(-2), 1.0*10**(-3)
    alpha_m = (A_alpha_m*(E_soma0-B_alpha_m)) / (1-exp((B_alpha_m-E_soma0)/C_alpha_m))
    # define beta_m term
    A_beta_m, B_beta_m, C_beta_m = 6.0*10**(+4), -4.9*10**(-2), 2.0*10**(-2)
    beta_m = (A_beta_m*(B_beta_m-E_soma0)) / (1-exp((E_soma0-B_beta_m)/C_beta_m))
    # define Na activation channel - differential equation here
    M[i] = m0
    m1 = m0 + ss*(alpha_m*(1-m0)-beta_m*m0)

    
    # sodium inactivation channel
    # define alpha_h term
    A_alpha_h, B_alpha_h, C_alpha_h = 8.0*10**(+4), -4.0*10**(-2), 1.0*10**(-3)
    alpha_h = (A_alpha_h*(B_alpha_h-E_soma0)) / (1-exp((E_soma0-B_alpha_h)/C_alpha_h))
    # define beta_h term
    A_beta_h, B_beta_h, C_beta_h = 4.0*10**(+2), -3.6*10**(-2), 2.0*10**(-3)
    beta_h = (A_beta_h) / (1+exp((B_beta_h-E_soma0)/C_beta_h))
    # define Na inactivation channel - differential equation here
    H[i] = h0 
    h1 = h0 + ss*(alpha_h*(1-h0)-beta_h*h0)
    
    # define NA channel current
    E_Na, G_Na = 5.0*10**(-2), 1.0*10**(-6)
    I_Na = (E_Na - E_soma0)*G_Na*(m0**3)*h0
    
    # potassium channel
    # define alpha_n term
    A_alpha_n, B_alpha_n, C_alpha_n = 2.0*10**(+4), -3.1*10**(-2), 8.0*10**(-4)
    alpha_n = (A_alpha_n*(E_soma0-B_alpha_n)) / (1-math.exp((B_alpha_n-E_soma0)/C_alpha_n))
    # define beta_n term
    A_beta_n, B_beta_n, C_beta_n = 5.0*10**(+3), -2.8*10**(-2), 4.0*10**(-4)
    beta_n = (A_beta_n*(B_beta_n-E_soma0)) / (1-exp((E_soma0-B_beta_n)/C_beta_n))
    # define K channel - differential equation
    N[i] = n0
    n1 = n0 + ss*(alpha_n*(1-n0) - beta_n*n0)
    
    # define K channel current
    E_K, G_K = -9.0*10**(-2), 2.0*10**(-7)
    I_K = (E_K - E_soma0)*G_K*(n0**4)
    
    # calcium channel
    # define alpha_q term
    A_alpha_q, B_alpha_q, C_alpha_q = 0.08*10**(+6), -10*10**(-3), 11*10**(-3)
    alpha_q = (A_alpha_q*(E_soma0-B_alpha_q)) / (1-exp((B_alpha_q-E_soma0)/C_alpha_q))
    # define beta_q term
    A_beta_q, B_beta_q, C_beta_q = 0.001*10**(+6), -10*10**(-3), 0.5*10**(-3)
    beta_q = (A_beta_q*(B_beta_q-E_soma0)) / (1-exp((E_soma0-B_beta_q)/C_beta_q))
    # define Ca channel - differential equation
    Q[i] = q0
    q1 = q0 + ss*(alpha_q*(1-q0) - beta_q*q0)
    
    # define Ca channel current
    E_Ca, G_Ca = 150*10**(-3), 1.0*10**(-8)
    I_Ca = (E_Ca - E_soma0)*G_Ca*q0**5
    
    # define intracellular Ca - differential equation
    rho_ap, delta_ap = 4.0*10**(+3), 30.0
    CAAP[i] = CaAP0
    CaAP1 = CaAP0 + ss*((E_Ca - E_soma0)*rho_ap*q0**5 - delta_ap*CaAP0)
    
    # define Ca dependent Potassium channel current
    G_K_Ca = 0.01*10**(-6)
    I_K_Ca = (E_K - E_soma0)*G_K_Ca*CaAP0
    
    # sum current from all channels (e.g., I_Na + I_K)
    I = I_Na + I_K + I_Ca + I_K_Ca + I_ext

    
    # leak current
    E_leak, G_m, C_m = -7.0*10**(-2), 3.0*10**(-9), 3.0*10**(-11)
    # define (E_soma) memberane current - differential equation
    E_SOMA[i] = E_soma0
    E_soma1 = E_soma0 + ss*(((E_leak - E_soma0)*G_m + I) / C_m)
    
    # update differential equation terms from n+1 to n
    t0, m0, n0, h0, q0, CaAP0, E_soma0 = t1, m1, n1, h1, q1, CaAP1, E_soma1



plt.figure(1)
plot(T,M,label = "Nam")
xlabel("Time(sec)")
ylabel("Na Activation")

plt.figure(2)
plot(T,H,label = "Nah")
xlabel("Time(sec)")
ylabel("Na Inactivation")

plt.figure(3)
plot(T,N,label = "K")
xlabel("Time(sec)")
ylabel("Potassium")

plt.figure(4)
plot(T,Q,label = "Q")
xlabel("Time(sec)")
ylabel("Calcium Channel")

plt.figure(5)
plot(T,CAAP,label ="CaAP")
xlabel("Time(sec)")
ylabel("CaAP")

plt.figure(6)
plot(T,E_SOMA,label ="ESOMA")
xlabel("Time(sec)")
ylabel("Membrane Potential (V)")


# Find time of each membrane potential peak
from scipy.signal import find_peaks
peaks,_ = find_peaks(E_SOMA, height = 0)

r = T[peaks]
print(r)





    
    
    
    
    