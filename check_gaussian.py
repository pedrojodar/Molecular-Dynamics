import matplotlib.pyplot as plt
import pandas as pd
res = pd.read_table("gaussian_test.txt", header=None)
print (len(res))
print (res.mean())
print (res.std())
plt.hist(res)
plt.show()