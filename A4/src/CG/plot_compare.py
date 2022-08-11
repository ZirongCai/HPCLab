import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas

matplotlib.use("Agg")

size=[]
nmpi=[]
time1=[]
time2=[]
file = open("output_purempi.txt", "r")
file = file.readlines()[1:]
for line in file:
    data=line.rstrip().split()
    nmpi.append(data[0])
    size.append(float(data[0]))
    time1.append(float(data[6]))

file2 = open("output_hybrid.txt", "r")
file2 = file2.readlines()[1:]
for line in file2:
    data=line.rstrip().split()
    time2.append(float(data[6]))

#--- TIME PLOT -----------------------------------------------------
time1 = np.array(time1)
time2 = np.array(time2)

plt.xlabel('nCores', fontsize=14)
#plt.xscale('log', base=2)
#plt.yscale('log')
plt.ylabel('RunTime (s)', fontsize=14)
plt.grid()
plt.plot(size, time1, marker="o", label="PureMPI")
plt.plot(size, time2, marker="o", label="Hybrid, nOMP=7")
plt.legend()
plt.savefig("../../figs/3.hybrid_vs_mpi.png",dpi=300,bbox_inches="tight")

## create figure and axis objects with subplots()
#fig,ax = plt.subplots()
## make a plot
#ax.bar(nmpi, speedup, color="blue")
##plt.xscale('log', base=2)
##plt.xlim(2**(-4), 2**(-14))
##plt.xticks(size, my_xticks)
#plt.grid(axis='y')
##plt.ylim(0,10.5)
## set x-axis label
#ax.set_xlabel("Number of MPI processes",fontsize=14)
## set y-axis label
#ax.tick_params(axis='y', colors='blue')
#ax.set_ylabel("Speedup",color="blue",fontsize=14)
#
## twin object for two different y-axis on the sample plot
#ax2=ax.twinx()
## make a plot with different y-axis using second axis object
#ax2.plot(nmpi, eff ,color="red",marker="o")
#ax2.set_ylabel("Efficiency",color="red",fontsize=14)
#ax2.tick_params(axis='y', colors='red')
##plt.ylim(0,1.05)
##plt.xticks(size, my_xticks)
##plt.grid()
## save the plot as a file
##plt.savefig("../../figs/5.strong_1.png",dpi=300,bbox_inches="tight")

