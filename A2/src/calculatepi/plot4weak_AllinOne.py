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

serial_balanced_1000.append(float(lines[1].split()[3]))
reduction_balanced_1000.append(float(lines[1].split()[4]))
speedup_balanced_1000 = np.divide(np.array(serial_balanced_1000),np.array(reduction_balanced_1000))
out.close()



#--------------------2000---------------------------------
#--------------------balanced-----------------------------
with open("../../result/calculatePi/scalingbalanced_2000.txt","r") as out:
	lines = out.readlines()

serial_balanced_2000 = []
reduction_balanced_2000 = []
serial_balanced_2000.append(float(lines[2].split()[3]))
reduction_balanced_2000.append(float(lines[2].split()[4]))

speedup_balanced_2000 = np.divide(np.array(serial_balanced_2000),np.array(reduction_balanced_2000))
out.close()


#--------------------4000---------------------------------
#--------------------balanced-----------------------------
with open("../../result/calculatePi/scalingbalanced_4000.txt","r") as out:
	lines = out.readlines()

serial_balanced_4000 = []
reduction_balanced_4000 = []

serial_balanced_4000.append(float(lines[3].split()[3]))
reduction_balanced_4000.append(float(lines[3].split()[4]))

speedup_balanced_4000 = np.divide(np.array(serial_balanced_4000),np.array(reduction_balanced_4000))
out.close()




#--------------------7000-----------------------------
#--------------------balanced-----------------------------
with open("../../result/calculatePi/scalingbalanced_7000.txt","r") as out:
	lines = out.readlines()

serial_balanced_7000 = []
reduction_balanced_7000 = []
serial_balanced_7000.append(float(lines[4].split()[3]))
reduction_balanced_7000.append(float(lines[4].split()[4]))

speedup_balanced_7000 = np.divide(np.array(serial_balanced_7000),np.array(reduction_balanced_7000))
out.close()




#--------------------8000-----------------------------
#--------------------balanced-----------------------------
with open("../../result/calculatePi/scalingbalanced_8000.txt","r") as out:
	lines = out.readlines()

serial_balanced_8000 = []
reduction_balanced_8000 = []
serial_balanced_8000.append(float(lines[5].split()[3]))
reduction_balanced_8000.append(float(lines[5].split()[4]))

speedup_balanced_8000 = np.divide(np.array(serial_balanced_8000),np.array(reduction_balanced_8000))
out.close()


#--------------------14000-----------------------------
#--------------------balanced-----------------------------
with open("../../result/calculatePi/scalingbalanced_14000.txt","r") as out:
	lines = out.readlines()

serial_balanced_14000 = []
reduction_balanced_14000 = []
serial_balanced_14000.append(float(lines[6].split()[3]))
reduction_balanced_14000.append(float(lines[6].split()[4]))

speedup_balanced_14000 = np.divide(np.array(serial_balanced_14000),np.array(reduction_balanced_14000))
out.close()



#--------------------16000-----------------------------
#--------------------balanced-----------------------------
with open("../../result/calculatePi/scalingbalanced_16000.txt","r") as out:
	lines = out.readlines()

serial_balanced_16000 = []
reduction_balanced_16000 = []
serial_balanced_16000.append(float(lines[7].split()[3]))
reduction_balanced_16000.append(float(lines[7].split()[4]))

speedup_balanced_16000 = np.divide(np.array(serial_balanced_16000),np.array(reduction_balanced_16000))
out.close()



#--------------------28000-----------------------------
#--------------------balanced-----------------------------
with open("../../result/calculatePi/scalingbalanced_28000.txt","r") as out:
	lines = out.readlines()

serial_balanced_28000 = []
reduction_balanced_28000 = []
serial_balanced_28000.append(float(lines[8].split()[3]))
reduction_balanced_28000.append(float(lines[8].split()[4]))

speedup_balanced_28000 = np.divide(np.array(serial_balanced_28000),np.array(reduction_balanced_28000))
out.close()


#--------------------56000-----------------------------
#--------------------balanced-----------------------------
with open("../../result/calculatePi/scalingbalanced_56000.txt","r") as out:
	lines = out.readlines()

serial_balanced_56000 = []
reduction_balanced_56000 = []
serial_balanced_56000.append(float(lines[9].split()[3]))
reduction_balanced_56000.append(float(lines[9].split()[4]))

speedup_balanced_56000 = np.divide(np.array(serial_balanced_56000),np.array(reduction_balanced_56000))
out.close()






fig1, axe1 = plt.subplots() # two axes on figure

xdata = ('(1,1000)', '(2,2000)', '(4,4000)', '(8,8000)', '(16,16000)')
ydata = []
ydata.append(reduction_balanced_1000)
ydata.append(reduction_balanced_2000)
ydata.append(reduction_balanced_4000)
ydata.append(reduction_balanced_8000)
ydata.append(reduction_balanced_16000)
ydata = np.multiply(ydata,1000000)

ydata_critical = []
ydata_critical.append(serial_balanced_1000)
ydata_critical.append(serial_balanced_2000)
ydata_critical.append(serial_balanced_4000)
ydata_critical.append(serial_balanced_8000)
ydata_critical.append(serial_balanced_16000)
ydata_critical = np.multiply(ydata_critical,1000)

axe1.set_ylabel('RunTime($\mu s$)',fontsize=8)
axe1.set_xlabel('(NumThread, Data Size)',fontsize=8)

axe1.plot(xdata,ydata, color='Blue', marker='.')


fig3, axe3 = plt.subplots() # two axes on figure


axe3.set_ylabel('RunTime(ms)',fontsize=8)
axe3.set_xlabel('(NumThread, Data Size)',fontsize=8)
axe3.plot(xdata,ydata_critical, color='Green', marker='.')

fig1.savefig("../../figs/WeakScaling1_AllinOne.png",dpi=300,box_inches="tight")


# fig2, axe2 = plt.subplots() # two axes on figure

# xdata = ('(7,7000)', '(14,14000)', '(28,28000)', '(56,56000)')
# ydata = []
# ydata.append(reduction_balanced_7000)
# ydata.append(reduction_balanced_14000)
# ydata.append(reduction_balanced_28000)
# ydata.append(reduction_balanced_56000)
# ydata = np.multiply(ydata,1000000)

# axe2.set_ylabel('RunTime($\mu s$)',fontsize=8)
# axe2.set_xlabel('(NumThread, Data Size)',fontsize=8)


# axe2.plot(xdata,ydata, marker='.')

fig3.savefig("../../figs/WeakScaling2_AllinOne.png",dpi=300,box_inches="tight")