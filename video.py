from matplotlib import animation
import matplotlib.pyplot as plt
import pandas as pd
import sys as sys
res=[]
fig, ax = plt.subplots()
for i in range (int(sys.argv[1])):
	res.append(pd.read_table('video2/'+str(i)+'.txt',header=None, sep='\t'))
len(res)
def update (i):
	str_val=str(i)+'.txt'
	ax.clear()
	ax.scatter(res[i][0], res[i][1],c='b')
update(0)
anim = animation.FuncAnimation(fig, update, frames = int(sys.argv[1]))
anim.save('linechart2.gif')