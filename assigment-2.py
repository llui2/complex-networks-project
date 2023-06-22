import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='Times')
plt.rc('mathtext', fontset='cm')
#------------------------------------------------------------------------------------#
# ASSIGMENT 2
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
# Calculate the cumulative and complementary cumulative distribution function
CDF = np.zeros(N, dtype=float)
CCDF = np.zeros(N, dtype=float)
for i in range(len(D)):
       CDF[i] = sum(P[:i+1])
       CCDF[i] = sum(P[i:])
#------------------------------------------------------------------------------------#
# Plot the degree distribution, cumulative and complementary cumulative distribution
fig = plt.figure(figsize=(3.4,3))
ax = plt.subplot()
plt.subplots_adjust(wspace=0.5,hspace=0.5,left=0.15,top=0.9,right=0.9,bottom=0.15)
ax.tick_params(direction='in', top=True, right=True, which='both')
ax.plot(range(0,N), CDF, label='$P(k_i<k)$', marker='o', markersize=1, linewidth=.4, linestyle='-', color='blue')
ax.plot(range(0,N), CCDF, label='$P_c(k)$', marker='o', markersize=1, linewidth=.4, linestyle='-', color='red')
ax.plot(range(0,N), P, label='$P(k)$', marker='o', markersize=1, linewidth=.4, linestyle='-', color='black')
ax.set_xlabel('$k$')
ax.set_xlim(1, 1e3)
#plt.ylim()
ax.set_xscale('log')
ax.set_yscale('log')
plt.legend(loc=1,frameon=False,bbox_to_anchor=(1.0, 0.9))
plt.savefig('fig1.pdf')
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
plt.savefig('fig2.pdf')
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
fig = plt.figure(figsize=(3.2,3))
plt.tick_params(direction='in', top=True, right=True, which='both')
plt.plot(range(len(c_avg)), c_avg, marker='o', markersize=1, linewidth=.4, linestyle='-', color='forestgreen')
plt.subplots_adjust(wspace=0.5,hspace=0.5,left=0.2,top=0.9,right=0.95,bottom=0.15)
plt.xlabel('$k$')
plt.ylabel('$\\overline{c}(k)$')
plt.xlim(1, 2e2)
plt.ylim(0.01,1)
plt.yticks([0.01,0.1,1])
plt.xscale('log')
plt.yscale('log')
plt.title('$\\overline{c} = $'+f'{c:.5f}', fontsize=12)
plt.savefig('fig3.pdf')