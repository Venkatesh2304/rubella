import pandas as pd 
import numpy as np 
from numpy import linalg as LA
from transform import transform

cm = pd.read_excel("data/extended_polymod_fitted.xlsx").to_numpy()
#cm = pd.read_excel("data/polymod_fitted.xlsx").to_numpy()
intial_pl  = pd.read_excel("data/population_2001.xlsx")
agroups = list(range(0,81,5)) + [100]
pl = transform(intial_pl.POPULATION.to_numpy(),agroups).round()

import seropositive
# seropositive : a*exp(-b*x) , here a1,b1 for < 13 yrs & a2,b2 >= 13 
ss_ratio = [ 0.2 ] +  list(seropositive.exponential(range(1,13),seropositive.a1,seropositive.b1)) + list(
                           seropositive.exponential(range(13,50),seropositive.a2,seropositive.b2))


def R0_finder(factor) :
    n = pl.sum()
    y = np.eye(50)
    l = np.eye(100)
    for i in range(50) : 
        for j in range(50) : 
            y[i,j] = (1 - np.exp(-factor*cm[i,j]*1/n))*pl[i] #pl[j]
            l[i,j] = factor*cm[i,j]*pl[i]/pl[j]

    #linear estimation 
    r = (y*pl/n).sum()
    #largest positive eigen value of the next generation matrix 
    #(e,v) = LA.eig(factor*cm)
    (e,v) = LA.eig(l)
    print( r , e[0])
    return r,0

factor =  pow(10,0)*2.5
e,v = R0_finder(factor)
print("R0 fitted in contact matrix :: ",e)
cm = cm*factor