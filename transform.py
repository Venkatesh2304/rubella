#returns the linearly transformed data for all age groups 
import numpy as np
def transform(data,ages1,no_scale=False) : 
    if type(ages1) == list : ages1 = np.array(ages1)
    lb = ages1[:-1]
    ub = ages1[1:] 
    diff = ub-lb 
    if no_scale : temp = data 
    else : temp = data/diff
    data1 = [] 
    for i in range(len(data)): 
        data1 += [temp[i]]*diff[i]
    return np.array(data1[:50])

