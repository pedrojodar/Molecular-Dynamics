import pandas as pd
import matplotlib.pyplot as plt
res = pd.read_table('position.txt', header=None, sep='\t')
plt.scatter(res[0], res[1])
res = pd.read_table('position_ac.txt', header=None, sep='\t')
plt.scatter(res[0], res[1])
plt.show()