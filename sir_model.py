import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numba
from scipy import integrate

N = pow(10,9) #Total number of individuals, N
I0, R0 = 10000 , pow(10,9)*0.6 #Initial number of infected and recovered individuals
S0 = N - I0 - R0 #Susceptible individuals to infection initially is deduced
beta, gamma = 5*1/14, 1/14 #Contact rate and mean recovery rate
tmax = 365 #A grid of time points (in days)
Nt = 100
t = np.linspace(0, tmax, Nt+1)
X0 = S0, I0, R0 #Initial conditions vector
def derivative(X, t):
    S, I, R = X
    dotS = -beta * S * I / N
    dotI = beta * S * I / N - gamma * I
    dotR = gamma * I
    return np.array([dotS, dotI, dotR])
res = integrate.odeint(derivative, X0, t)
S, I, R = res.T
Seuil = 1 - 1 / (beta/gamma)
plt.figure()
plt.grid()
plt.title("odeint method")
plt.plot(t, S, 'orange', label='Susceptible')
plt.plot(t, I, 'r', label='Infected')
plt.plot(t, R, 'g', label='Recovered with immunity')
plt.xlabel('Time t, [days]')
plt.ylabel('Numbers of individuals')
plt.ylim([0,N])
plt.legend()

plt.show()