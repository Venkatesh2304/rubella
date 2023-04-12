import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
"""  
  Data Source :: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697#pcbi-1005697-g003
  5  year age strata's from 0 - 80 years (16 categories)
  Extended into 80 age groups (1-80) using a linear fit 
"""

df = pd.read_excel("data/extended_polymod.xlsx")

f = interpolate.interp2d(np.arange(2.5,82.5,5), np.arange(2.5,82.5,5), df.to_numpy() , kind='cubic')
r = np.matrix(np.eye(81))
import matplotlib.pyplot as plt
for i in np.arange(0,81) :
    for j in np.arange(0,81) : 
         r[i,j] = f(i,j)
r = r[:50,:50]
pd.DataFrame(r).to_excel("data/extended_polymod_fitted.xlsx",index=False)
#plt.imshow(r, origin="lower")
#plt.savefig("observations/extended_polymod_fitted.png")

df = pd.read_excel("data/polymod.xlsx")
f = interpolate.interp2d(np.arange(2.5,72.5,5), np.arange(2.5,72.5,5), df.to_numpy() , kind='cubic')
r = np.matrix(np.eye(71))
import matplotlib.pyplot as plt
for i in np.arange(0,71) :
    for j in np.arange(0,71) : 
         r[i,j] = f(i,j)
r = r[:50,:50]
pd.DataFrame(r).to_excel("data/polymod_fitted.xlsx",index=False)
plt.imshow(r, origin="lower")
plt.savefig("observations/polymod_fitted.png")


plt.colorbar()
#plt.show()


