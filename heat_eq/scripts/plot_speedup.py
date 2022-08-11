import numpy as np
import matplotlib.pyplot as plt

# import all the data
from data import *
from data_blocking import *

topos = [
    "1 x 48", "2 x 24", "3 x 16", "4 x 12", "6 x 8", "8 x 6", "12 x 4",
    "16 x 3", "24 x 2", "48 x 1"
]

print(f"1 node  best topology based on average {topos[best_1node_index]}")
print(f"2 nodes best topology based on average {topos[best_2node_index]}")
print(f"3 nodes best topology based on average {topos[best_3node_index]}")
print(f"4 nodes best topology based on average {topos[best_4node_index]}")

#
# Speed up comparison of influence_topology
#
timeSq = np.array(
    [0.089263, 0.363984, 0.824760, 1.470021, 2.271087 , 3.266086, 4.527948, 5.833436, 7.363150, 9.066547, 10.942300])


fig, ax1 = plt.subplots()
ax1.set_title("Comparison with sequential verison")
ax1.set_xlabel("Resolution")
ax1.set_ylabel("Speedup")

for i, t in enumerate(time_1node):
    _ = ax1.plot(res, timeSq/t, label=f'Topo {topos[i]}', linestyle='-', markevery=1)

plt.legend(loc='upper right')

plt.savefig('influence_topology_speedup.png', format='png', dpi=300)

