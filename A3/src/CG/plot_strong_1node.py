import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas

matplotlib.use("Agg")

size=[]
time=[]
nmpi=[]
file = open("result_strong1nodes.txt", "r")
file = file.readlines()[1:]
for line in file:
    data=line.rstrip().split()
    nmpi.append(data[0])
    size.append(float(data[0]))
    time.append(float(data[5]))

#--- TIME PLOT -----------------------------------------------------
time = np.array(time)
speedup = time[0]/time
eff  = speedup/(np.array(size))
#plt.plot(size, time, marker='.')
#
#plt.xlabel('10 * Gridsize', fontsize=14)
#plt.xscale('log', base=2)
#plt.yscale('log')
#plt.ylabel('Normalized Time', fontsize=14)
#plt.xlim(2**(-4), 2**(-14))
#plt.grid()
#plt.xticks(size, my_xticks)
#plt.savefig("../../figs/5.gridsizes_1.png",dpi=300,bbox_inches="tight")

# create figure and axis objects with subplots()
fig,ax = plt.subplots()
# make a plot
ax.bar(nmpi, speedup, color="blue")
#plt.xscale('log', base=2)
#plt.xlim(2**(-4), 2**(-14))
#plt.xticks(size, my_xticks)
plt.grid(axis='y')
#plt.ylim(0,10.5)
# set x-axis label
ax.set_xlabel("Number of MPI processes",fontsize=14)
# set y-axis label
ax.tick_params(axis='y', colors='blue')
ax.set_ylabel("Speedup",color="blue",fontsize=14)

# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
# make a plot with different y-axis using second axis object
ax2.plot(nmpi, eff ,color="red",marker="o")
ax2.set_ylabel("Efficiency",color="red",fontsize=14)
ax2.tick_params(axis='y', colors='red')
#plt.ylim(0,1.05)
#plt.xticks(size, my_xticks)
#plt.grid()
# save the plot as a file
#plt.savefig("../../figs/5.strong_1.png",dpi=300,bbox_inches="tight")
plt.savefig("../../figs/5.strong_4.png",dpi=300,bbox_inches="tight")

