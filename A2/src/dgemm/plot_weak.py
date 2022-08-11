import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas

matplotlib.use("Agg")

thread=[]
size=[]
time=[]

file = open("output_weak.txt", "r")
file = file.readlines()[1:]
for line in file:
    data=line.rstrip().split()
    thread.append(int(data[0]))
    size.append(float(data[1]))
    time.append(float(data[2]))

tmp = time[0]
for i in range(len(time)):
    time[i] /= tmp 


plt.plot(thread,time,marker='.')
plt.xlabel('nThreads')
plt.xscale('log', base=2)
plt.ylabel('Normalized Time')
plt.grid()
plt.legend()
#plt.title('Initial performance on laptop')
#plt.show();
plt.savefig("../../figs/4.3_WeakScale.png",dpi=300,bbox_inches="tight")
#plt.close()

