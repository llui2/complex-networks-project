# STOCHASTIC SIMULATION OF THE SUSCEPTIBLE-INFECTED-SUSCEPTIBLE (SIS) MODEL USING THE GILLESPIE ALGORITHM
import numpy as np
import matplotlib.pyplot as plt
print('SIMULATION-1')
#------------------------------------------------------------------------------------#
# Parameters
delta = 1 # Recovery rate
#lambda_list = np.linspace(0.1, 0.9, 4) # Infection rates to be tested
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
# Initial condition
N_infected = 0.05*N # Number of initially infected nodes
n_0 = np.zeros(N, dtype=int)
for i in range(int(N_infected)):
       n_0[i] = 1
       np.random.shuffle(n_0)
#------------------------------------------------------------------------------------#
# Simulation time
T = 20000
#------------------------------------------------------------------------------------#
import time
# get the start time
st = time.time()
# Create directory to store the results
import os
os.system('mkdir results/')
# Dynamics
for lambda_ in lambda_list:
       #----------------------------------------------------------#
       # Initial Infected nodes list and Initial Active links list
       N_0 = []
       E_0 = []
       for i in range(N):
              if n_0[i] == 1:
                     N_0.append(i+1)
              for j in V[i]:
                     if n_0[i] == 1 and n_0[j-1] == 0:
                            E_0.append([i+1,j])
       #----------------------------------------------------------#
       # Suport list for active links list update
       support_list_0 = [[] for i in range(N)]
       for i in range(N):
              for label,j in enumerate(V[i]):
                     support_list_0[i].append(0)
                     for index,k in enumerate(E_0):
                            if i+1 == k[0] and j == k[1] or i+1 == k[1] and j == k[0]:
                                   support_list_0[i][label] = index+1
       #----------------------------------------------------------#
       # Write the time series of the number of infected nodes
       f = open('results/lambda1_'+f'{lambda_:.2}'.replace('.', '')+'.dat', 'w')
       f.write('# t <n_t>\n')
       f.write(f'{0} {np.mean(n_0)}'+'\n')
       #----------------------------------------------------------#
       # Physical time
       tau = 0
       #----------------------------------------------------------#
       # Initial condition
       n_t = n_0.copy()
       N_t = N_0.copy()
       E_t = E_0.copy()
       support_list = support_list_0.copy()
       #----------------------------------------------------------#
       for t in range(T):

              # End if there are no infected nodes or active links
              if len(N_t) == 0 or len(E_t) == 0:
                     break

              # Compute the infection and recovery probabilities
              prob_infection = len(E_t)*lambda_/(lambda_*len(E_t)+delta*len(N_t))
              prob_recovery = len(N_t)*delta/(lambda_*len(E_t)+delta*len(N_t))
              # Pick one of the possible events with the corresponding probability
              event = np.random.choice(['infection', 'recovery'], p=[prob_infection, prob_recovery])

              # Infection event
              if event == 'infection' and len(E_t) > 0:
#------------------------------------------------------------------------------------------------#
                     # Pick one of the active links at random
                     link = np.random.choice(range(len(E_t)))
                     # Pick the susceptible node of the link
                     if n_t[E_t[link][0]-1] == 0:
                            node_0 = E_t[link][0]
                            node_1 = E_t[link][1]
                     else:
                            node_0 = E_t[link][1]
                            node_1 = E_t[link][0]
                     # Update the state of the susceptible node
                     n_t[node_0-1] = 1
                     # Add the susceptible node to the list of infected nodes
                     N_t.append(node_0)
                     # Remove the link from the list of active links
                     if link != len(E_t)-1:
                            E_t[link] = E_t[-1]
                            support_list[node_0-1][V[node_0-1].index(node_1)] = 0
                            support_list[node_1-1][V[node_1-1].index(node_0)] = 0
                            support_list[E_t[-1][0]-1][V[E_t[-1][0]-1].index(E_t[-1][1])] = link+1
                            support_list[E_t[-1][1]-1][V[E_t[-1][1]-1].index(E_t[-1][0])] = link+1
                            E_t.pop()
                     else:
                            E_t.pop()
                            support_list[node_0-1][V[node_0-1].index(node_1)] = 0
                            support_list[node_1-1][V[node_1-1].index(node_0)] = 0
                     # Add new active links and remove the new inactive links
                     for j in V[node_0-1]:
                            if n_t[j-1] == 0:
                                   E_t.append([node_0,j])
                                   support_list[node_0-1][V[node_0-1].index(j)] = len(E_t)
                                   support_list[j-1][V[j-1].index(node_0)] = len(E_t)
                            elif n_t[j-1] == 1 and j != node_1:
                                   link_position = support_list[node_0-1][V[node_0-1].index(j)]
                                   if link_position != len(E_t):
                                          E_t[link_position-1] = E_t[-1]
                                          support_list[node_0-1][V[node_0-1].index(j)] = 0
                                          support_list[j-1][V[j-1].index(node_0)] = 0
                                          support_list[E_t[-1][0]-1][V[E_t[-1][0]-1].index(E_t[-1][1])] = link_position
                                          support_list[E_t[-1][1]-1][V[E_t[-1][1]-1].index(E_t[-1][0])] = link_position
                                          E_t.pop()
                                   else:
                                          E_t.pop()
                                          support_list[node_0-1][V[node_0-1].index(j)] = 0
                                          support_list[j-1][V[j-1].index(node_0)] = 0
#------------------------------------------------------------------------------------------------#
              # Recovery event
              elif event == 'recovery' and len(N_t) > 0:
                     # Pick one of the infected nodes with the corresponding probability
                     node_position = np.random.choice(len(N_t))
                     node = N_t[node_position]
                     # Pick the infected node
                     N_t[node_position] = N_t[-1]
                     # Remove the last infected node
                     N_t.pop()
                     # Update the state of the infected node
                     n_t[node-1] = 0
                     # Add new active links and remove the new inactive links
                     for j in V[node-1]:
                            if n_t[j-1] == 1:
                                   E_t.append([node,j])
                                   support_list[node-1][V[node-1].index(j)] = len(E_t)
                                   support_list[j-1][V[j-1].index(node)] = len(E_t)
                            else:
                                   link_position = support_list[node-1][V[node-1].index(j)]
                                   if link_position != len(E_t):
                                          E_t[link_position-1] = E_t[-1]
                                          support_list[node-1][V[node-1].index(j)] = 0
                                          support_list[j-1][V[j-1].index(node)] = 0
                                          support_list[E_t[-1][0]-1][V[E_t[-1][0]-1].index(E_t[-1][1])] = link_position
                                          support_list[E_t[-1][1]-1][V[E_t[-1][1]-1].index(E_t[-1][0])] = link_position
                                          E_t.pop()
                                   else:
                                          E_t.pop()
                                          support_list[node-1][V[node-1].index(j)] = 0
                                          support_list[j-1][V[j-1].index(node)] = 0
#------------------------------------------------------------------------------------------------#
              # Update physical time with a Poisson process
              tau += np.random.exponential(1/(lambda_*len(E_t)+delta*len(N_t)))
              f.write(f'{tau} {np.mean(n_t)}'+'\n')

              if lambda_ == lambda_list[0] and t == T-1:
                     # estimate time in HH:MM:SS
                     et = time.time()
                     print(f'Estimated time = {time.strftime("%H:%M:%S", time.gmtime((et-st)*len(lambda_list)))}')
# estimate time in HH:MM:SS
cput = time.time()
print(f'CPU time = {time.strftime("%H:%M:%S", time.gmtime((cput-st)))}')     
#----------------------------------------------------------------------------------------#

nt_avg_list = []

import matplotlib.pyplot as plt
import numpy as np

plt.rc('font', family='Times', size=12)
plt.rc('mathtext', fontset='cm')

plt.figure(figsize=(7,3))
ax = plt.subplot()
plt.tick_params(direction='in', top=True, right=True, which='both')
plt.subplots_adjust(wspace=0.5,hspace=0.5,left=0.14,top=0.9,right=0.9,bottom=0.15)

colors = ['#DB2D43','#007FFF','#006B3C','#FF8C00','#9400D3']
for a,lambda_ in enumerate(lambda_list):
       data = np.loadtxt('results/lambda1_'+f'{lambda_:.2}'.replace('.', '')+'.dat', dtype=float)
       t = data[:,0]
       nt = data[:,1]

       nt_avg = np.mean(nt[int(len(nt)*0.3):])
       nt_avg_list.append(nt_avg)

       t1 = np.linspace(0,100,10)

       plt.plot(t1,nt_avg*np.ones(len(t1)), marker=' ', markersize=1, linewidth=1, linestyle='--', color=colors[a])
       plt.plot(t,nt, marker=' ', markersize=1, linewidth=1, linestyle='-', label=f'$\\lambda = $ {lambda_:.2}', color=colors[a])

plt.xlabel('$t$')
plt.ylabel('$\\rho_\mathrm{stoc}$')
plt.xlim(0,25)
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

plt.savefig('fig8.pdf')
#----------------------------------------------------------------------------------------#

plt.rc('font', family='Times', size=10)

plt.figure(figsize=(3.2,3))
ax = plt.subplot()
plt.tick_params(direction='in', top=True, right=True, which='both')

plt.plot(lambda_list,nt_avg_list, marker=' ', markersize=1, linewidth=1, linestyle='-', color='black')

plt.subplots_adjust(wspace=0.5,hspace=0.5,left=0.2,top=0.9,right=0.95,bottom=0.15)
plt.xlabel('$\\lambda$')
plt.ylabel('$\\rho^S(\\lambda)$')
plt.xlim(0,1)
plt.ylim(0,1)

plt.yticks([0,0.2,0.4,0.6,0.8,1],['0','0.2','0.4','0.6','0.8','1'])

labels = [label.get_text() for i, label in enumerate(ax.xaxis.get_majorticklabels())]
labels[0] = '0'
ax.set_xticklabels(labels)

plt.savefig('fig9.pdf')
#----------------------------------------------------------------------------------------#