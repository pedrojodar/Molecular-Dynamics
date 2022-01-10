from matplotlib import animation
import matplotlib.pyplot as plt
import pandas as pd
import sys as sys
res=[]
res_ac=[]

fig, ax = plt.subplots()
for i in range (int(sys.argv[1])):
    res_ac.append(pd.read_table('video2/'+str(i)+'ac.txt',header=None, sep='\t'))
    res.append(pd.read_table('video2/'+str(i)+'.txt',header=None, sep='\t'))
    
len(res)
def update (i):
    str_val=str(i)+'.txt'
    ax.clear()
    ax.scatter(res_ac[i][0], res_ac[i][1],c='purple')
    ax.scatter(res[i][0], res[i][1],c='cyan')

update(0)
anim = animation.FuncAnimation(fig, update, frames = int(sys.argv[1]))
anim.save('linechart2_ac.gif')