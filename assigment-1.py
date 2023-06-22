import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='Times')
plt.rc('mathtext', fontset='cm')
#------------------------------------------------------------------------------------#
# ASSIGMENT 1 (http://konect.cc/networks/arenas-meta/)
#------------------------------------------------------------------------------------#
# Get the number of nodes
N = 0
with open('edge_list.dat') as f:
       for line in f:
              N = max(N, int(line.split()[0]), int(line.split()[1]))
#------------------------------------------------------------------------------------#
# Get the degree and neighbors of each node
D = np.zeros(N, dtype=int)
V = [[] for i in range(N)]
with open('edge_list.dat') as f:
       for line in f:
              i = int(line.split()[0])
              j = int(line.split()[1])
              if i != j:
                     if j not in V[i-1]:
                            D[i-1] += 1
                            V[i-1].append(j)
                     if i not in V[j-1]:
                            D[j-1] += 1
                            V[j-1].append(i)
# Sort the neighbors of each node
for i in range(N):
       V[i].sort()
#------------------------------------------------------------------------------------#
# Print the number of nodes, links, average degree and list of degrees
print('Number of nodes: ', N)
print('Number of links: ', int(sum(D)/2))
k_avg = np.mean(D)
print('Average degree: ', k_avg)
#print('List of degrees: ', D)
#------------------------------------------------------------------------------------#
# Write in memory number of nodes, links, average degree, list of degrees and neighbors
with open('info.dat', 'w') as f:
       f.write('number of nodes, links, average degree' + '\n')
       f.write(f'{N} {int(sum(D)/2)} {k_avg}')
       f.write('\n')
       f.write('list of degrees' + '\n')
       f.write(' '.join(map(str, D)))
       f.write('\n')
       f.write('neighbors' + '\n')
       for i in range(N):
              f.write(' '.join(map(str, V[i])) + '\n')
#------------------------------------------------------------------------------------#