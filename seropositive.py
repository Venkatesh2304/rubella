import numpy as np
from scipy.optimize import minimize
import pandas as pd 
import matplotlib.pyplot as plt
""" 
   Uses Vellore 1999 seropositive data , fitted a exponential decay curve 
   S(x) = 1 - a*e^(-b*x) , where x is age and S(x) is propotion of suspectible population 
   with log likelihood as the cost function . 
   Returns :: a1,b1 from age <= 13 and a2,b2 > 13 
"""
_df = pd.read_excel("data/seropositive_vellore.xlsx")

def exponential(x, a, b):
    if type(x) == list : x = np.array(x)
    return a*np.exp(-b*x)

def log_likelihood(params, x):
    a, b = params  
    y_pred = exponential(x, a, b)
    return -np.sum( w1*np.log(y_pred) + w2*np.log(1-y_pred) )

initial_guess = np.random.uniform(0.1,0.9,2)
sub_model = [[13,39],[0,13]]
Y,Y_PRED = np.array([]) , np.array([])
params = []
for [n1,n2] in sub_model :
    df = _df.iloc[n1:n2]
    w1 = np.array(df["total negative"])
    w2 = np.array(df["total test"])
    x = np.array(df["Age"])
    result = minimize(log_likelihood, initial_guess, args=(x,) , bounds=((0+pow(10,-5),1-pow(10,-5)),(0+pow(10,-5),None)))
    a, b = result.x
    params.append((a,b))
    Y=np.append(1 - w1/(w1+w2),Y)
    Y_PRED=np.append(1 - exponential(x,a,b),Y_PRED)

[[a1,b1],[a2,b2]] = params 
plt.plot(Y)
plt.plot(Y_PRED)
plt.savefig("observations/seropostive_vellore.png")
print("Model is a*exp(-b*x)")
for i in range(len(sub_model)): 
    print(f"For Age between {sub_model[i]} ,  a = {params[i][0]:.5f}, b = {params[i][1]:.5f}")


