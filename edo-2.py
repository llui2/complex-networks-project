# NUMERICAL SOLUTION OF THE SUSCEPTIBLE-INFECTED-SUSCEPTIBLE (SIS) MODEL OF THE HETEROGENEOUS MEAN-FIELD APPROXIMATION
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
print('EDO-2')
#------------------------------------------------------------------------------------#
# Parameters
delta = 1 # Recovery rate
#lambda_list = np.linspace(0.1, 0.9, 5) # Infection rates to be tested
lambda_list =[0.1,0.2,0.4,0.9]
#------------------------------------------------------------------------------------#
# Read the network from a file edge_list.dat
edge_list = np.loadtxt('edge_list.dat', dtype=int)
# Number of nodes
N = np.max(edge_list)
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
# Create the adjacency matrix
A = np.zeros((N,N), dtype=int)
for i in range(N):
       for j in V[i]:
              A[i][j-1] = 1
#------------------------------------------------------------------------------------#
# Create the initial condition
N_infected = 0.05*N # Number of initially infected nodes
p_0 = np.zeros(N, dtype=float)
for i in range(int(N_infected)):
       p_0[i] = 1
np.random.shuffle(p_0)
#------------------------------------------------------------------------------------#
def sis_model(p, t, delta, lambda_, adjacency_matrix):
       n = len(p)
       dp_dt = np.zeros(n)

       for i in range(n):
              sum_neighbors = np.sum(adjacency_matrix[i] * p)
              dp_dt[i] = -delta * p[i] + lambda_ * sum_neighbors * (1 - p[i])

       return dp_dt     

def solve_sis_model(p0, t, delta, lambda_, adjacency_matrix):
       p_solution = odeint(sis_model, p0, t, args=(delta, lambda_, adjacency_matrix))
       p_avg = np.mean(p_solution, axis=1)
       return p_avg
#------------------------------------------------------------------------------------#
# Solve the SIS model and obtain the average p(t)

# Plot the results
plt.rc('font', family='Times', size=12)
plt.rc('mathtext', fontset='cm')

plt.figure(figsize=(7,3))
ax = plt.subplot()
plt.tick_params(direction='in', top=True, right=True, which='both')
plt.subplots_adjust(wspace=0.5,hspace=0.5,left=0.15,top=0.9,right=0.9,bottom=0.15)

import time
# get the start time
st = time.time()

colors = ['#DB2D43','#007FFF','#006B3C','#FF8C00','#9400D3']
for a,lambda_ in enumerate(lambda_list):

       data = np.loadtxt('results/lambda_'+f'{lambda_:.2}'.replace('.', '')+'.dat', dtype=float)
       t = data[:,0]
       nt = data[:,1]
       plt.plot(t,nt, marker=' ', markersize=1, linewidth=1, linestyle='-', color=colors[a], label='$\\lambda = {}$'.format(lambda_))

       p_avg = solve_sis_model(p_0, t, delta, lambda_, A)
       plt.plot(t, p_avg, marker=' ', linewidth=1, linestyle='--', alpha=0.8, color=colors[a])

       if lambda_ == lambda_list[0]:
              # estimate time in HH:MM:SS
              et = time.time()
              print(f'Estimated time = {time.strftime("%H:%M:%S", time.gmtime((et-st)*len(lambda_list)))}')


plt.xlabel('$t$')
plt.ylabel('$\\rho(t)$')
plt.xlim(0,10)
plt.ylim(0,1)

plt.yticks([0,0.2,0.4,0.6,0.8,1],['0','0.2','0.4','0.6','0.8','1'])

labels = [label.get_text() for i, label in enumerate(ax.xaxis.get_majorticklabels())]
labels[0] = '0'
ax.set_xticklabels(labels)

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.savefig('fig13.pdf')