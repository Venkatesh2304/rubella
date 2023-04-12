import pandas as pd 
import numpy as np 
from transform import transform

ratios = np.array([13.6,3.1,2.3,30,50.9])
agroups = [0,1] +  [ i*5 for i in range(1,12) ] 
mr = pd.read_excel("data/mortality.xlsx")
cbr = pd.read_excel("data/CBR.xlsx")
dr = [0.04014,0.00768,0.00379,0.00315,0.00484,0.00698,0.00777,0.00975,0.01341,0.01834,0.02623,0.04254,0.06241,0.09228,0.13633,0.20462,0.29153,0.44991]
dr =  dr[:len(agroups)-1]
def death_rate(population_age_wise,year) :
    #total_deaths = round( population_age_wise.sum() * mr[mr.year == year].iloc[0].MR / 1000 )
    #print(f"Total deaths in {year} :: {total_deaths}")
    #deaths = (ratios/100)*total_deaths
    #sr = transform(deaths,agroups)/population_age_wise
    #print( transform(dr,agroups) )
    #print( population_age_wise )
    #print( transform(dr,agroups) )
    #input()
    #print( ".." , population_age_wise.sum() , (population_age_wise*transform(dr,agroups)).sum())
    return transform(dr,agroups)

def births(population,year) : 
    if type(population) == np.array : population = population.sum()
    #print( population , cbr[cbr.year == year].iloc[0].CBR , population * cbr[cbr.year == year].iloc[0].CBR / 1000)
    return round( population * cbr[cbr.year == year].iloc[0].CBR / 1000 )

"""
Test :: 
x = np.random.randint(10,100,50)
print( f"Total population :: {x.sum()}" ) 
print( x )
print( survival_rates(x,2001) )
print( births(1000,2001) )
"""

    




