import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas

matplotlib.use("Agg")

thread=[]
copy=[]
add=[]
file = open("output_2.6.txt", "r")
file = file.readlines() [1:]
for line in file:
    data=line.rstrip().split()
    thread.append(int(data[0]))
    copy.append(float(data[1]))
    add.append(float(data[3]))


plt.plot(thread,copy,marker='.', label='Copy')
plt.plot(thread,add,marker='.', label='Add')


file = open("output_2.6_noft.txt", "r")
file = file.readlines()[1:]
copy=[]
add=[]
for line in file:
    data=line.rstrip().split()
    copy.append(float(data[1]))
    add.append(float(data[3]))


plt.plot(thread,copy,marker='.', label='Copy-noFT')
plt.plot(thread,add,marker='.', label='Add-noFT')




plt.xlabel('nThreads')
#plt.xscale('log', base=2)
plt.ylabel('Best rate (MB/s)')
plt.grid()
plt.legend()
#plt.title('Initial performance on laptop')
#plt.show();
plt.savefig("../../figs/2.6_StreamCompare.png",dpi=300,bbox_inches="tight")
#plt.close()

