import numpy as np
import matplotlib.pyplot as plt
import random
plt.rc('font', family='Times')
plt.rc('mathtext', fontset='cm')
#------------------------------------------------------------------------------------#
# ASSIGMENT 5
#------------------------------------------------------------------------------------#
# Read in memory number of nodes, links, average degree, list of degrees and neighbors
with open('info.dat', 'r') as f:
       f.readline()
       N, L, k_avg = map(float, f.readline().split())
       N, L = int(N), int(L)
       f.readline()
       D = list(map(int, f.readline().split()))
       f.readline()
       V = []
       for i in range(N):
              V.append(list(map(int, f.readline().split())))
#------------------------------------------------------------------------------------#
# Given a list of degrees, create a function that returns the adjacency matrix of a
# configuration model network with the given degree sequence.
def CM_djacency_matrix(degrees):
       # Create a list of nodes with degrees equal to the given list of degrees
       nodes = []
       for i in range(len(degrees)):
              for j in range(degrees[i]):
                     nodes.append(i+1)
       # Shuffle the list of nodes
       random.shuffle(nodes)
       # Create a list of edges
       edges = []
       for i in range(0, len(nodes), 2):
              edges.append([nodes[i], nodes[i+1]])
       # Remove self and multiple loops
       for i in range(len(edges)):
              if edges[i][0] == edges[i][1]:
                     edges[i][1] = random.choice(nodes)
       # Create the adjacency matrix
       adjacency_matrix = np.zeros((len(degrees), len(degrees)), dtype=int)
       for i in range(len(edges)):
              adjacency_matrix[edges[i][0]-1, edges[i][1]-1] = 1
              adjacency_matrix[edges[i][1]-1, edges[i][0]-1] = 1
       return adjacency_matrix
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
# Given a list of degrees, create a function that returns the adjacency matrix of a
# configuration model network with the given degree sequence.
def CM_djacency_matrix_1(degrees):
       # Create a list of nodes with degrees equal to the given list of degrees
       nodes = []
       V = []
       for i in range(len(degrees)):
              for j in range(degrees[i]):
                     nodes.append(i+1)
                     V = [[] for i in range(len(degrees))]
       # Fill the list of neighbors randomly without self and multiple loops
       for i in range(len(degrees)):
              while len(V[i]) < degrees[i]:
                     to_connect = random.choice(nodes)
                     if to_connect != i+1 and to_connect not in V[i]:
                            V[i].append(to_connect)
                            V[to_connect-1].append(i+1)
       # Create the adjacency matrix
       adjacency_matrix = np.zeros((len(degrees), len(degrees)), dtype=int)
       for i in range(len(V)):
              for j in range(len(V[i])):
                     adjacency_matrix[i, V[i][j]-1] = 1
       return adjacency_matrix
#------------------------------------------------------------------------------------#
A = CM_djacency_matrix(D)
#------------------------------------------------------------------------------------#
V = []
for i in range(N):
       V.append([])
       for j in range(N):
              if A[i][j] == 1:
                     V[i].append(j+1)
#------------------------------------------------------------------------------------#
# Calculate the number of nodes with degree k
N_k = np.zeros(N, dtype=int)       
for i in range(N):
       N_k[D[i]] += 1
#------------------------------------------------------------------------------------#
# Calculate the probability of each degree value
P = np.zeros(N, dtype=float)
for i in range(N):
       P[D[i]] += 1
P /= N
#------------------------------------------------------------------------------------#
# Calculate the average nearest neighbor degree
k_nn = np.zeros(N, dtype=float)
for i in range(N):
       for j in range(len(V[i])):
              k_nn[D[i]] += D[V[i][j]-1]/(D[i]*N_k[D[i]])
#------------------------------------------------------------------------------------#
# Plot the average nearest neighbor degree
fig = plt.figure(figsize=(3.2,3))
ax = plt.subplot()
plt.subplots_adjust(wspace=0.5,hspace=0.5,left=0.2,top=0.9,right=0.95,bottom=0.15)
ax.tick_params(direction='in', top=True, right=True, which='both')
ax.plot(range(0,N), k_nn/k_avg, marker='o', markersize=1, linewidth=.4, linestyle='-', color='crimson')
ax.set_ylabel('$\\overline{k}_{nn}(k)/\\langle k \\rangle$')
ax.set_xlabel('$k$')
ax.set_xlim(1, 1e2)
#ax.set_ylim(0.9,10)
ax.set_xscale('log')
ax.set_yscale('log')
plt.yticks([0.1,1,10,100])
plt.savefig('fig4.pdf')
#------------------------------------------------------------------------------------#
# Calculate the clustering coefficient for each node
c_avg = np.zeros(N, dtype=float)
for i in range(N):
       if D[i] < 2:
              c_avg[D[i]] += 0
       else:
              num_triangles = 0
              for j in range(len(V[i])):
                     for k in range(j+1, len(V[i])):
                            if V[i][k] in V[V[i][j]-1]:
                                   num_triangles += 1
              c_avg[D[i]] += 2*num_triangles / (D[i]*(D[i]-1)*N_k[D[i]])
c = sum(c_avg*P)
#------------------------------------------------------------------------------------#
# Plot the clustering coefficient
plt.figure(figsize=(3.2,3))
plt.tick_params(direction='in', top=True, right=True, which='both')
plt.plot(range(len(c_avg)), c_avg, marker='o', markersize=1, linewidth=.4, linestyle='-', color='forestgreen')
plt.subplots_adjust(wspace=0.5,hspace=0.5,left=0.2,top=0.9,right=0.95,bottom=0.15)
plt.xlabel('$k$')
plt.ylabel('$\\overline{c}(k)$')
plt.xlim(1, 2e2)
#plt.ylim()
plt.xscale('log')
plt.yscale('log')
plt.yticks([0.0001,0.001,0.01,0.1,1])
plt.title('$\\overline{c} = $'+f'{c:.5f}', fontsize=12)
plt.savefig('fig5.pdf')
#------------------------------------------------------------------------------------#
import time
# get the start time
st = time.time()
# Do the same for the configuration model but average over 100 realizations
R = 100
k_nn = np.zeros(N, dtype=float)
c_avg = np.zeros(N, dtype=float)
for r in range(R):
       A = CM_djacency_matrix(D)
       V = []
       for i in range(N):
              V.append([])
              for j in range(N):
                     if A[i][j] == 1:
                            V[i].append(j+1)
       for i in range(N):
              for j in range(len(V[i])):
                     k_nn[D[i]] += D[V[i][j]-1]/(D[i]*N_k[D[i]])
       for i in range(N):
              if D[i] < 2:
                     c_avg[D[i]] += 0
              else:
                     num_triangles = 0
                     for j in range(len(V[i])):
                            for k in range(j+1, len(V[i])):
                                   if V[i][k] in V[V[i][j]-1]:
                                          num_triangles += 1
                     c_avg[D[i]] += 2*num_triangles / (D[i]*(D[i]-1)*N_k[D[i]])
       if r==0:
              # get the end time
              et = time.time()
              print(f'Estimated time = {time.strftime("%H:%M:%S", time.gmtime((et-st)*R))}')
k_nn = k_nn/R
c_avg = c_avg/R
c = sum(c_avg*P)
#------------------------------------------------------------------------------------#
fig = plt.figure(figsize=(3.2,3))
ax = plt.subplot()
plt.subplots_adjust(wspace=0.5,hspace=0.5,left=0.2,top=0.9,right=0.95,bottom=0.15)
ax.tick_params(direction='in', top=True, right=True, which='both')
ax.plot(range(0,N), k_nn/k_avg, marker='o', markersize=1, linewidth=.4, linestyle='-', color='crimson')
ax.set_ylabel('$\\overline{k}_{nn}(k)/\\langle k \\rangle$')
ax.set_xlabel('$k$')
ax.set_xlim(1, 1e2)
#ax.set_ylim(0.9,10)
ax.set_xscale('log')
ax.set_yscale('log')
plt.yticks([0.1,1,10,100])
plt.savefig('fig6.pdf')
#------------------------------------------------------------------------------------#
plt.figure(figsize=(3.2,3))
plt.tick_params(direction='in', top=True, right=True, which='both')
plt.plot(range(len(c_avg)), c_avg, marker='o', markersize=1, linewidth=.4, linestyle='-', color='forestgreen')
plt.subplots_adjust(wspace=0.5,hspace=0.5,left=0.2,top=0.9,right=0.95,bottom=0.15)
plt.xlabel('$k$')
plt.ylabel('$\\overline{c}(k)$')
plt.xlim(1, 2e2)
#plt.ylim()
plt.xscale('log')
plt.yscale('log')
plt.yticks([0.0001,0.001,0.01,0.1,1])
plt.title('$\\overline{c} = $'+f'{c:.5f}', fontsize=12)
plt.savefig('fig7.pdf')

