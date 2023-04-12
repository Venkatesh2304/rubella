import pandas as pd 
import numpy as np 
from transform import transform

fr = [0,11.3,113.6,139.6,84.4,35.6,11.7,4.7]
fr_ages = [0,15,20,25,30,35,40,45,50]
fr_rates =  transform(fr,fr_ages,no_scale=True) / 1000 

#intial conditions 
import seropositive 
intial_pl  = pd.read_excel("data/population_2001.xlsx")
agroups = list(range(0,81,5)) + [100]
intial_pl = transform(intial_pl.POPULATION.to_numpy(),agroups).round()

# seropositive : a*exp(-b*x) , here a1,b1 for < 13 yrs & a2,b2 >= 13 
ss_ratio = [ 0.2 ] +  list(seropositive.exponential(range(1,13),seropositive.a1,seropositive.b1)) + list(
                           seropositive.exponential(range(13,50),seropositive.a2,seropositive.b2))
intial_ss = (np.array(ss_ratio) * intial_pl)
intial_m = np.zeros(50)
intial_m[0] = intial_pl[0]*0.8 
intial_i =  np.concatenate((np.ones(15)*30,np.ones(35)*5))*pow(10,3)*5/( 30*15 + 5*35) #np.random.choice(5,50,p=[0.3,0.3,0.15,0.15,0.10])
intial_im = intial_pl - intial_m - intial_ss  - intial_i
intial_im[intial_im < 0 ] = 0  
states = [None]*(50*104)
states[0] = np.stack((intial_m,intial_ss,intial_i,intial_im),axis=-1).flatten()
import cm_rescaler
import mortality 

#returns a 200x200 transition matrix 
def transition_matrix(d,v,infected,cm,sr,br,immigrations,total_population,t) : 
    t1 = []
    cm = cm*( 1 + 0.15*np.cos(2*np.pi*t) )
    for i in range(50) : 
        foi = 1 - np.exp((-cm[i,:] * np.nan_to_num(infected**0.97)).sum()/ total_population)
        a_i = [
             [ 1-d , 0 , 0 , 0 ] , 
             [ d, (1-foi)*(1-v[i]) , 0 , 0] , 
             [ 0, foi*(1-v[i]) ,0, 0],
             [0, v[i] , 1 ,1 ] 
            ]
        a_i = np.array(a_i)
        if np.isnan(a_i).any() :
           print(infected)
           input(a_i)
        t1.append(a_i)
    return t1 


R0=5
t=14
year=2000
years = list(range(2000,2040))
cm = cm_rescaler.cm

for i,year in enumerate(years) :
  for idx in range(i*24,(i+1)*24) : #range(i*104,(i+1)*105) :   
    state = states[idx]
    pl_from_state = lambda state : state.reshape(50,4).sum(axis=-1)
    sr = mortality.death_rate( pl_from_state(state)  , year)
    sr = sr*t/365
    sr = 1 - sr 
    d = 1/((9*30)/t) #probability of losing maternal immunity 
    infected = state.reshape(50,4)[:,2]
    v = np.zeros(50) #vaccination for that state  

    br = mortality.births( state.sum() , year  )*t/365
    if idx%24 == 0 : print( year ,  (br*365)/ (((1-sr)*pl_from_state(state)).sum()*365)    )
    immigrations = np.zeros(50)
    total_population = state.sum()
    
    a = transition_matrix(d,v,infected,cm,sr,br,immigrations,total_population,idx%24)
    next_state = []
    for i in range(50) :
        next_state_i = (a[i]*sr[i]*( 1 - t/365 )).dot(state[4*i:4*(i+1)])
        if i!=0 : next_state_i += (a[i-1]*sr[i-1]*t/365).dot(state[4*(i-1):4*i])
        else : next_state_i[0] += br 
        next_state += list(next_state_i)

    next_state = np.array(next_state)    
    states[idx+1] = next_state

states = [ i for i in states if i is not None ]

#sus1  = states[0].reshape(50,4)[:,1]/states[0].reshape(50,4)[:,1].sum()
#sus2  = states[-1].reshape(50,4)[:,1]/states[-1].sum()


inf = [ state.reshape(50,4)[:,2].sum()/(state.reshape(50,4).sum())*pow(10,9) for state in states if state is not None ]
sus = [ 1 - (state.reshape(50,4)[:,1].sum()/(state.reshape(50,4).sum())) for state in states if state is not None  ]
plot_ages = [1,3,7,10,15,30,45]
for age in plot_ages : 
    _sus = [ 1 - (state.reshape(50,4)[age,1]/(state.reshape(50,4)[age,:].sum())) for state in states if state is not None  ]
    #plt.plot(_sus,label = str(age))

_sus = [ state.reshape(50,4)[16:40,2].sum() for state in states if state is not None  ]
_sus1 = [ state.reshape(50,4)[0:10,2].sum() for state in states if state is not None  ]

crs = [ (state.reshape(50,4)[:,2]*fr_rates*0.65*14/52).sum() for state in states if state is not None ]
crs = [ sum(crs[i:i+24]) for i in range(0,len(crs)-24,24) ]

from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

x = [ i for i in years ][10:]
y = crs[10:] 
X_Y_Spline = make_interp_spline(x, y)
X_ = np.linspace(min(x), max(x), 500)
Y_ = X_Y_Spline(X_)
plt.xlim(left = 2010, right=2040)
plt.plot(X_  ,Y_)
plt.show()


for i in range(0,50,10): 
    _sus = [ state.reshape(50,4)[i:i+10,1].sum() for state in states if state is not None  ]
    _pop = [ state.reshape(50,4)[i:i+10,:].sum() for state in states if state is not None  ]
    #plt.plot( np.array(_sus),label=f"{i} - {i+10}")
    #plt.plot( np.array(_sus)/np.array(_pop),label=f"{i} - {i+10}")

#print( _sus )
mat = [ state.reshape(50,4)[:,0].sum() for state in states if state is not None ]
im = [ state.reshape(50,4)[:,3].sum() for state in states if state is not None  ]
#plt.plot(_sus1)
#plt.plot(inf)
#plt.plot(sus)
#plt.plot(mat)
#plt.plot(im)
















 


