import numpy as np
import matplotlib.pyplot as plt

# import all the data
from data import *
from data_blocking import *
from data_hybrid import *

topos = [
    "1 x 48", "2 x 24", "3 x 16", "4 x 12", "6 x 8", "8 x 6", "12 x 4",
    "16 x 3", "24 x 2", "48 x 1"
]

print(f"1 node  best topology based on average {topos[best_1node_index]}")
print(f"2 nodes best topology based on average {topos[best_2node_index]}")
print(f"3 nodes best topology based on average {topos[best_3node_index]}")
print(f"4 nodes best topology based on average {topos[best_4node_index]}")
 
print(f"1 node  best hybrid topology based on average {omp_1node_best_topo[0]} with {omp_1node_best_topo[1]}")
print(f"2 nodes best hybrid topology based on average {omp_2node_best_topo[0]} with {omp_2node_best_topo[1]}")
print(f"3 nodes best hybrid topology based on average {omp_3node_best_topo[0]} with {omp_3node_best_topo[1]}")

#
# Speedup vs sequential version
#
#  fig, ax1 = plt.subplots()
#  ax1.set_title("Speedup compared to sequential")
#  ax1.set_xlabel("Number of Nodes")
#  ax1.set_ylabel("Speedup")
#
#  nodes = np.array([1, 2, 3, 4])
#
#  perf = np.array([best_1node, best_2node, best_3node, best_4node])
#
#  plt.xticks(np.arange(min(nodes), max(nodes) + 1, 1.0))
#
#  NUM_COLORS = perf[0].size
#  LINE_STYLES = ['solid', 'dashed', 'dashdot', 'dotted']
#  NUM_STYLES = len(LINE_STYLES)
#
#  cm = plt.get_cmap('gist_rainbow')
#
#  ax1.set_prop_cycle('color',
#                     [cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])
#
#  for i, row in enumerate(perf.transpose()[:6]):
#      lines = ax1.plot(nodes,
#                       sequential_time[i] / row,
#                       label=f'Res {(i+1) * 1000}',
#                       linestyle='-',
#                       markevery=1)
#      lines[0].set_color(cm(i // NUM_STYLES * float(NUM_STYLES) / NUM_COLORS))
#      lines[0].set_linestyle(LINE_STYLES[i % NUM_STYLES])
#
#  plt.legend(loc='upper left')
#  plt.savefig('seq_speedup_first5.png', format='png', dpi=300)
#
#  fig, ax1 = plt.subplots()
#  ax1.set_title("Speedup compared to sequential")
#  ax1.set_xlabel("Number of Nodes")
#  ax1.set_ylabel("Speedup")
#
#  nodes = np.array([1, 2, 3, 4])
#
#  perf = np.array([best_1node, best_2node, best_3node, best_4node])
#
#  plt.xticks(np.arange(min(nodes), max(nodes) + 1, 1.0))
#
#  NUM_COLORS = perf[0].size
#  LINE_STYLES = ['solid', 'dashed', 'dashdot', 'dotted']
#  NUM_STYLES = len(LINE_STYLES)
#
#  cm = plt.get_cmap('gist_rainbow')
#
#  ax1.set_prop_cycle('color',
#                     [cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])
#
#  for i, row in enumerate(perf.transpose()[5:]):
#      lines = ax1.plot(nodes,
#                       sequential_time[i] / row,
#                       label=f'Res {(i+6) * 1000}',
#                       linestyle='-',
#                       markevery=1)
#      lines[0].set_color(cm(i // NUM_STYLES * float(NUM_STYLES) / NUM_COLORS))
#      lines[0].set_linestyle(LINE_STYLES[i % NUM_STYLES])
#
#  plt.legend(loc='upper left')
#  plt.savefig('seq_speedup_last6.png', format='png', dpi=300)

#
# How does the topology of the cartesian grid influence the throughput
#
#  fig, ax1 = plt.subplots()
#  ax1.set_title("Test")
#  ax1.set_xlabel("Resolution")
#  ax1.set_ylabel("GFlop/s")
#
#  for i, p in enumerate(perf_1node):
#      _ = ax1.plot(res, p, label=f'Topo {topos[i]}', linestyle='-', markevery=1)
#
#  plt.legend(loc='upper right')
#
#  plt.savefig('influence_topology.png', format='png', dpi=300)

#
# How does our code scale with more nodes?
#
#  fig, ax1 = plt.subplots()
#  ax1.set_title("Scaling per Node")
#  ax1.set_xlabel("Number of Nodes")
#  ax1.set_ylabel("Scaling")
#
#  nodes = np.array([1, 2, 3, 4])
#
#  perf = np.array([best_1node, best_2node, best_3node, best_4node])
#
#  ax1.axhline(2, color="grey", ls="dashed")
#  ax1.axhline(3, color="lightgrey", ls="dashed")
#  ax1.axhline(4, color="silver", ls="dashed")
#
#  ax1.axvline(2, color="grey", ls="dashed")
#  ax1.axvline(3, color="lightgrey", ls="dashed")
#  ax1.axvline(4, color="silver", ls="dashed")
#
#  plt.xticks(np.arange(min(nodes), max(nodes) + 1, 1.0))
#
#  NUM_COLORS = perf[0].size
#  LINE_STYLES = ['solid', 'dashed', 'dashdot', 'dotted']
#  NUM_STYLES = len(LINE_STYLES)
#
#  cm = plt.get_cmap('gist_rainbow')
#
#  ax1.set_prop_cycle('color',
#                     [cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])
#
#  for i, row in enumerate(perf.transpose()):
#      lines = ax1.plot(nodes,
#                       row[0] / row,
#                       label=f'Res {(i+1) * 1000}',
#                       linestyle='-',
#                       markevery=1)
#      lines[0].set_color(cm(i // NUM_STYLES * float(NUM_STYLES) / NUM_COLORS))
#      lines[0].set_linestyle(LINE_STYLES[i % NUM_STYLES])
#
#  plt.legend(loc='upper left')
#  plt.savefig('node_scaling.png', format='png', dpi=300)
#
#  #
#  # Plot to compare non blocking and blocking communication
#  #
#
#  fig, ax1 = plt.subplots()
#  ax1.set_title("Blocking vs non-blocking communication")
#  ax1.set_xlabel("Resolution")
#  ax1.set_ylabel("GFlop/s")
#
#  _ = ax1.plot(res,
#               perf_1node[best_1node_index],
#               label=f'1 Node nb',
#               color='red',
#               linestyle='-',
#               markevery=1)
#  _ = ax1.plot(res,
#               perf_2node[best_2node_index],
#               label=f'2 Node nb',
#               color='blue',
#               linestyle='-',
#               markevery=1)
#  _ = ax1.plot(res,
#               perf_3node[best_3node_index],
#               label=f'3 Node nb',
#               color='green',
#               linestyle='-',
#               markevery=1)
#  _ = ax1.plot(res,
#               perf_4node[best_4node_index],
#               label=f'4 Node nb',
#               color='black',
#               linestyle='-',
#               markevery=1)
#
#  _ = ax1.plot(res,
#               bbest_1node,
#               label=f'1 Node b',
#               color='red',
#               linestyle='--',
#               markevery=1)
#  _ = ax1.plot(res,
#               bbest_2node,
#               label=f'2 Node b',
#               color='blue',
#               linestyle='--',
#               markevery=1)
#  _ = ax1.plot(res,
#               bbest_3node,
#               label=f'3 Node b',
#               color='green',
#               linestyle='--',
#               markevery=1)
#  _ = ax1.plot(res,
#               bbest_4node,
#               label=f'4 Node b',
#               color='black',
#               linestyle='--',
#               markevery=1)
#
#  plt.legend(loc='upper right')
#  plt.savefig('communication.png', format='png', dpi=300)

fig, ax1 = plt.subplots()
ax1.set_title("Hybrid topologies")
ax1.set_xlabel("Resolution")
ax1.set_ylabel("GFlop/s")

NUM_COLORS = omp_1node_perf[0].size
LINE_STYLES = ['solid', 'dashed', 'dashdot', 'dotted']
NUM_STYLES = len(LINE_STYLES)

cm = plt.get_cmap('gist_rainbow')

ax1.set_prop_cycle('color',
                   [cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])

a = np.mean(omp_1node_perf, axis=1)
ind = np.argpartition(a, -8)[-8:]
sorted_ind = ind[np.argsort(a[ind])]
topos = np.flip(omp_1node_topo[sorted_ind], axis=0)
perf = np.flip(omp_1node_perf[sorted_ind], axis=0)

for i, row in enumerate(perf):
    lines = ax1.plot(res,
                     row,
                     label=f'Topo {topos[i][0]}, {topos[i][1]} OMP threads',
                     linestyle='-',
                     markevery=1)

plt.legend(loc='upper right')
plt.savefig('hybrid_best1node.png', format='png', dpi=300)
 
fig, ax1 = plt.subplots()
ax1.set_title("Scaling per Node")
ax1.set_xlabel("Number of Nodes")
ax1.set_ylabel("Scaling")

nodes = np.array([1, 2, 3])

perf = np.array([omp_1node_best_perf, omp_2node_best_perf, omp_3node_best_perf])

ax1.axhline(2, color="grey", ls="dashed")
ax1.axhline(3, color="lightgrey", ls="dashed")
ax1.axhline(4, color="silver", ls="dashed")

ax1.axvline(2, color="grey", ls="dashed")
ax1.axvline(3, color="lightgrey", ls="dashed")
ax1.axvline(4, color="silver", ls="dashed")

plt.xticks(np.arange(min(nodes), max(nodes) + 1, 1.0))

NUM_COLORS = perf[0].size
LINE_STYLES = ['solid', 'dashed', 'dashdot', 'dotted']
NUM_STYLES = len(LINE_STYLES)

cm = plt.get_cmap('gist_rainbow')

ax1.set_prop_cycle('color',
                   [cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])

for i, row in enumerate(perf.transpose()):
    lines = ax1.plot(nodes,
                     row / row[0],
                     label=f'Res {(i+1) * 1000}',
                     linestyle='-',
                     markevery=1)
    lines[0].set_color(cm(i // NUM_STYLES * float(NUM_STYLES) / NUM_COLORS))
    lines[0].set_linestyle(LINE_STYLES[i % NUM_STYLES])

plt.legend(loc='upper left')
plt.savefig('hybrid_scaling.png', format='png', dpi=300)

plt.show()
