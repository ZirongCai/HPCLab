import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

numThread = np.array([1,2,4,7,8,14,16,28,56])

#--------------------1000---------------------------------
#--------------------balanced-----------------------------
with open("../../result/calculatePi/scalingbalanced_1000.txt","r") as out:
	lines = out.readlines()

serial_balanced_1000 = []
reduction_balanced_1000 = []
for line in lines[1:]:
	serial_balanced_1000.append(float(line.split()[2]))
	reduction_balanced_1000.append(float(line.split()[4]))

speedup_balanced_1000 = np.divide(np.array(serial_balanced_1000),np.array(reduction_balanced_1000))
out.close()

#--------------------compact-----------------------------
with open("../../result/calculatePi/scalingcompact_1000.txt","r") as out:
	lines = out.readlines()

serial_compact_1000 = []
reduction_compact_1000 = []
for line in lines[1:]:
	serial_compact_1000.append(float(line.split()[2]))
	reduction_compact_1000.append(float(line.split()[4]))

speedup_compact_1000 = np.divide(np.array(serial_compact_1000),np.array(reduction_compact_1000))
out.close()

#--------------------scatter-----------------------------
with open("../../result/calculatePi/scalingscatter_1000.txt","r") as out:
	lines = out.readlines()

serial_scatter_1000 = []
reduction_scatter_1000 = []
for line in lines[1:]:
	serial_scatter_1000.append(float(line.split()[2]))
	reduction_scatter_1000.append(float(line.split()[4]))

speedup_scatter_1000 = np.divide(np.array(serial_scatter_1000),np.array(reduction_scatter_1000))
out.close()




#--------------------4000---------------------------------
#--------------------balanced-----------------------------
with open("../../result/calculatePi/scalingbalanced_4000.txt","r") as out:
	lines = out.readlines()

serial_balanced_4000 = []
reduction_balanced_4000 = []
for line in lines[1:]:
	serial_balanced_4000.append(float(line.split()[2]))
	reduction_balanced_4000.append(float(line.split()[4]))

speedup_balanced_4000 = np.divide(np.array(serial_balanced_4000),np.array(reduction_balanced_4000))
out.close()

#--------------------compact-----------------------------
with open("../../result/calculatePi/scalingcompact_4000.txt","r") as out:
	lines = out.readlines()

serial_compact_4000 = []
reduction_compact_4000 = []
for line in lines[1:]:
	serial_compact_4000.append(float(line.split()[2]))
	reduction_compact_4000.append(float(line.split()[4]))

speedup_compact_4000 = np.divide(np.array(serial_compact_4000),np.array(reduction_compact_4000))
out.close()

#--------------------scatter-----------------------------
with open("../../result/calculatePi/scalingscatter_4000.txt","r") as out:
	lines = out.readlines()

serial_scatter_4000 = []
reduction_scatter_4000 = []
for line in lines[1:]:
	serial_scatter_4000.append(float(line.split()[2]))
	reduction_scatter_4000.append(float(line.split()[4]))

speedup_scatter_4000 = np.divide(np.array(serial_scatter_4000),np.array(reduction_scatter_4000))
out.close()

#--------------------8000-----------------------------
#--------------------balanced-----------------------------
with open("../../result/calculatePi/scalingbalanced_8000.txt","r") as out:
	lines = out.readlines()

serial_balanced_8000 = []
reduction_balanced_8000 = []
for line in lines[1:]:
	serial_balanced_8000.append(float(line.split()[2]))
	reduction_balanced_8000.append(float(line.split()[4]))

speedup_balanced_8000 = np.divide(np.array(serial_balanced_8000),np.array(reduction_balanced_8000))
out.close()

#--------------------compact-----------------------------
with open("../../result/calculatePi/scalingcompact_8000.txt","r") as out:
	lines = out.readlines()

serial_compact_8000 = []
reduction_compact_8000 = []
for line in lines[1:]:
	serial_compact_8000.append(float(line.split()[2]))
	reduction_compact_8000.append(float(line.split()[4]))

speedup_compact_8000 = np.divide(np.array(serial_compact_8000),np.array(reduction_compact_8000))
out.close()

#--------------------scatter-----------------------------
with open("../../result/calculatePi/scalingscatter_8000.txt","r") as out:
	lines = out.readlines()

serial_scatter_8000 = []
reduction_scatter_8000 = []
for line in lines[1:]:
	serial_scatter_8000.append(float(line.split()[2]))
	reduction_scatter_8000.append(float(line.split()[4]))

speedup_scatter_8000 = np.divide(np.array(serial_scatter_8000),np.array(reduction_scatter_8000))
out.close()


#--------------------16000-----------------------------
#--------------------balanced-----------------------------
with open("../../result/calculatePi/scalingbalanced_16000.txt","r") as out:
	lines = out.readlines()

serial_balanced_16000 = []
reduction_balanced_16000 = []
for line in lines[1:]:
	serial_balanced_16000.append(float(line.split()[2]))
	reduction_balanced_16000.append(float(line.split()[4]))

speedup_balanced_16000 = np.divide(np.array(serial_balanced_16000),np.array(reduction_balanced_16000))
out.close()

#--------------------compact-----------------------------
with open("../../result/calculatePi/scalingcompact_16000.txt","r") as out:
	lines = out.readlines()

serial_compact_16000 = []
reduction_compact_16000 = []
for line in lines[1:]:
	serial_compact_16000.append(float(line.split()[2]))
	reduction_compact_16000.append(float(line.split()[4]))

speedup_compact_16000 = np.divide(np.array(serial_compact_16000),np.array(reduction_compact_16000))
out.close()

#--------------------scatter-----------------------------
with open("../../result/calculatePi/scalingscatter_16000.txt","r") as out:
	lines = out.readlines()

serial_scatter_16000 = []
reduction_scatter_16000 = []
for line in lines[1:]:
	serial_scatter_16000.append(float(line.split()[2]))
	reduction_scatter_16000.append(float(line.split()[4]))

speedup_scatter_16000 = np.divide(np.array(serial_scatter_16000),np.array(reduction_scatter_16000))
out.close()

#--------------------28000-----------------------------
#--------------------balanced-----------------------------
with open("../../result/calculatePi/scalingbalanced_28000.txt","r") as out:
	lines = out.readlines()

serial_balanced_28000 = []
reduction_balanced_28000 = []
for line in lines[1:]:
	serial_balanced_28000.append(float(line.split()[2]))
	reduction_balanced_28000.append(float(line.split()[4]))

speedup_balanced_28000 = np.divide(np.array(serial_balanced_28000),np.array(reduction_balanced_28000))
out.close()

#--------------------compact-----------------------------
with open("../../result/calculatePi/scalingcompact_28000.txt","r") as out:
	lines = out.readlines()

serial_compact_28000 = []
reduction_compact_28000 = []
for line in lines[1:]:
	serial_compact_28000.append(float(line.split()[2]))
	reduction_compact_28000.append(float(line.split()[4]))

speedup_compact_28000 = np.divide(np.array(serial_compact_28000),np.array(reduction_compact_28000))
out.close()

#--------------------scatter-----------------------------
with open("../../result/calculatePi/scalingscatter_28000.txt","r") as out:
	lines = out.readlines()

serial_scatter_28000 = []
reduction_scatter_28000 = []
for line in lines[1:]:
	serial_scatter_28000.append(float(line.split()[2]))
	reduction_scatter_28000.append(float(line.split()[4]))

speedup_scatter_28000 = np.divide(np.array(serial_scatter_28000),np.array(reduction_scatter_28000))
out.close()

#--------------------56000-----------------------------
#--------------------balanced-----------------------------
with open("../../result/calculatePi/scalingbalanced_56000.txt","r") as out:
	lines = out.readlines()

serial_balanced_56000 = []
reduction_balanced_56000 = []
for line in lines[1:]:
	serial_balanced_56000.append(float(line.split()[2]))
	reduction_balanced_56000.append(float(line.split()[4]))

speedup_balanced_56000 = np.divide(np.array(serial_balanced_56000),np.array(reduction_balanced_56000))
out.close()

#--------------------compact-----------------------------
with open("../../result/calculatePi/scalingcompact_56000.txt","r") as out:
	lines = out.readlines()

serial_compact_56000 = []
reduction_compact_56000 = []
for line in lines[1:]:
	serial_compact_56000.append(float(line.split()[2]))
	reduction_compact_56000.append(float(line.split()[4]))

speedup_compact_56000 = np.divide(np.array(serial_compact_56000),np.array(reduction_compact_56000))
out.close()

#--------------------scatter-----------------------------
with open("../../result/calculatePi/scalingscatter_56000.txt","r") as out:
	lines = out.readlines()

serial_scatter_56000 = []
reduction_scatter_56000 = []
for line in lines[1:]:
	serial_scatter_56000.append(float(line.split()[2]))
	reduction_scatter_56000.append(float(line.split()[4]))

speedup_scatter_56000 = np.divide(np.array(serial_scatter_56000),np.array(reduction_scatter_56000))
out.close()









fig, axes = plt.subplots(nrows=3, ncols=2) # two axes on figure
fig.suptitle('Comparison of Different Affinity with Different DataSize')


axes[0,0].plot(numThread, speedup_balanced_1000, color='Blue', marker='.',)
axes[0,0].plot(numThread, speedup_compact_1000, color='Red', marker='.',)
axes[0,0].plot(numThread, speedup_scatter_1000, color='Green', marker='.', )
axes[0,0].set_ylabel('SpeedUp',fontsize=8)
axes[0,0].get_xaxis().set_visible(False)

#plt.show()
axes[0,1].plot(numThread, speedup_balanced_4000, color='Blue', marker='.',label='balanced')
axes[0,1].plot(numThread, speedup_compact_4000, color='Red', marker='.',label='compact')
axes[0,1].plot(numThread, speedup_scatter_4000, color='Green', marker='.', label='scatter')
axes[0,1].get_xaxis().set_visible(False)
axes[0,1].legend(loc='upper right',prop={'size': 4})


axes[1,0].plot(numThread, speedup_balanced_8000, color='Blue', marker='.',)
axes[1,0].plot(numThread, speedup_compact_8000, color='Red', marker='.',)
axes[1,0].plot(numThread, speedup_scatter_8000, color='Green', marker='.', )
axes[1,0].set_ylabel('SpeedUp',fontsize=8)

axes[1,1].plot(numThread, speedup_balanced_16000, color='Blue', marker='.',)
axes[1,1].plot(numThread, speedup_compact_16000, color='Red', marker='.',)
axes[1,1].plot(numThread, speedup_scatter_16000, color='Green', marker='.', )

axes[2,0].plot(numThread, speedup_balanced_28000, color='Blue', marker='.',)
axes[2,0].plot(numThread, speedup_compact_28000, color='Red', marker='.',)
axes[2,0].plot(numThread, speedup_scatter_28000, color='Green', marker='.', )
axes[2,0].set_ylabel('SpeedUp',fontsize=8)
axes[2,0].set_xlabel('NumThread',fontsize=8)

axes[2,1].plot(numThread, speedup_balanced_56000, color='Blue', marker='.',)
axes[2,1].plot(numThread, speedup_compact_56000, color='Red', marker='.',)
axes[2,1].plot(numThread, speedup_scatter_56000, color='Green', marker='.', )
axes[2,1].set_xlabel('NumThread',fontsize=8)

#plt.show();
fig.savefig("../../figs/Compare_Aff_AllinOne.png",dpi=300,box_inches="tight")