import matplotlib.pyplot as plt
import pandas as pd
#from scipy.optimize import curve_fit
msd = pd.read_table('msd.txt', header=None, sep='\t')
#print (curve_fit(msd[0], msd[1]))
plt.plot(msd[0], msd[1])
plt.show()
