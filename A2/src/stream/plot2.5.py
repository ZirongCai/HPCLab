import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas

matplotlib.use("Agg")

thread=[]
copy=[]
scale=[]
add=[]
triad=[]
file = open("output_2.5_size40M.txt", "r")
file = file.readlines()[1:]
for line in file:
    data=line.rstrip().split()
    thread.append(int(data[0]))
    copy.append(float(data[1]))
    scale.append(float(data[2]))
    add.append(float(data[3]))
    triad.append(float(data[4]))


plt.plot(thread,copy,marker='.', label='Copy')
plt.plot(thread,scale,marker='.', label='Scale')
plt.plot(thread,add,marker='.', label='Add')
plt.plot(thread,triad,marker='.', label='Triad')
plt.xlabel('nThreads')
#plt.xscale('log', base=2)
plt.ylabel('Best rate (MB/s)')
plt.grid()
plt.legend()
#plt.title('Initial performance on laptop')
#plt.show();
plt.savefig("../../figs/2.5_StreamScale.png",dpi=300,bbox_inches="tight")
#plt.close()

