import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas

matplotlib.use("Agg")

array=[]
copy=[]
scale=[]
add=[]
triad=[]
file = open("output_scanning_2.3.txt", "r")
file = file.readlines()[1:]
for line in file:
    data=line.rstrip().split()
    array.append(int(data[0]))
    copy.append(float(data[1]))
    scale.append(float(data[2]))
    add.append(float(data[3]))
    triad.append(float(data[4]))


plt.plot(array,copy,marker='.', label='Copy')
plt.plot(array,scale,marker='.', label='Scale')
plt.plot(array,add,marker='.', label='Add')
plt.plot(array,triad,marker='.', label='Triad')
plt.xlabel('ArraySize')
plt.xscale('log', base=2)
plt.ylabel('Best rate (MB/s)')
plt.grid()
plt.legend()
#plt.title('Initial performance on laptop')
#plt.show();
plt.savefig("../../figs/2.4_StreamCaches.png",dpi=300,bbox_inches="tight")
#plt.close()

