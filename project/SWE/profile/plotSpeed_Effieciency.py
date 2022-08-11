import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
import sys

numThread = np.array([1,2,4])


with open("./profile/outputs/StrongScale_8192.txt","r") as out:
	lines = out.readlines()

cpuTime_8192 = []
commTime_8192 = []
my_xticks=[]

for line in lines:
	cpuTime_8192.append(float(line.split()[4]))
	commTime_8192.append(float(line.split()[10]))

out.close()

for i in range(len(cpuTime_8192)):
    my_xticks.append(str(2**i))


speedupCPU = cpuTime_8192[0] / np.array(cpuTime_8192) 
effCPU  = speedupCPU/numThread

# create figure and axis objects with subplots()
fig,ax = plt.subplots()
# make a plot
ax.bar(my_xticks, speedupCPU, color="blue")
plt.grid(axis='y')
#plt.ylim(0,10.5)
# set x-axis label
ax.set_xlabel("Number of Nodes",fontsize=14)
# set y-axis label
ax.tick_params(axis='y', colors='blue')
ax.set_ylabel("Normalized Speedup",color="blue",fontsize=14)

# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
# make a plot with different y-axis using second axis object
ax2.plot(my_xticks, effCPU ,color="red",marker="o")
ax2.set_ylabel("Efficiency",color="red",fontsize=14)
ax2.tick_params(axis='y', colors='red')
plt.ylim(0.8,1.05)

plt.savefig("./profile/SpeedUp_Eff_CPU.png",dpi=300,bbox_inches="tight")



#++++++++++++++++++++ figure 2+++++++++++++++++++++++++++++++
plot2 = plt.figure(2)

speedupCOMM = commTime_8192[0] / np.array(commTime_8192) 
effCOMM  = speedupCOMM/numThread

# create figure and axis objects with subplots()
fig,ax = plt.subplots()
# make a plot
ax.bar(my_xticks, speedupCOMM, color="blue")
plt.grid(axis='y')
#plt.ylim(0,10.5)
# set x-axis label
ax.set_xlabel("Number of Nodes",fontsize=14)
# set y-axis label
ax.tick_params(axis='y', colors='blue')
ax.set_ylabel("Normalized Speedup",color="blue",fontsize=14)

# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
# make a plot with different y-axis using second axis object
ax2.plot(my_xticks, effCOMM ,color="red",marker="o")
ax2.set_ylabel("Efficiency",color="red",fontsize=14)
ax2.tick_params(axis='y', colors='red')
#plt.ylim(0.8,1.05)

plt.savefig("./profile/SpeedUp_Eff_COMM.png",dpi=300,bbox_inches="tight")