# Extracting epsilon's

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

def RK4(nu,epsilon):                             # RK4 method with x,y,z being intial conditions, h is the incriments, and N is the number that you want x to go to
    
    # Initialise variables:
    i=0
    tau = 0
    theta = 0.03
    z = 0
    h = 0.05
    N = 500
    Y = np.zeros(int(N*(1/h))+1)
    X = np.zeros(int(N*(1/h))+1)
    def f(tau,theta,z):                                       # z=dy/dx
        return z    
    def g(tau,theta,z):                                       # dz/dx = d^2y/dx^2=...(solved for dz/dx)
        ## Constants:
        a = 0.01
        l = 1
        #nu = 1+(a/l)*(0.5)
        #epsilon = nu*(a/l)
    
        ## Equation: 
        dz = - nu*theta - epsilon*np.cos(2*tau)*theta
        return dz    

    ## Solving the equation:
    while tau<N:                                      # N is the nnumber of the x-value that you wish to go to
        
        ## Determining k1,k2,k3,k4:
    
        k1 = h * f(tau,theta,z)
    
        l1 = h * g(tau,theta,z)
    
        k2 = h * f(tau + h/2, theta + k1/2, z + l1/2)
    
        l2 = h * g(tau + h/2, theta + k1/2, z + l1/2)
    
        k3 = h * f(tau + h/2, theta + k2/2, z + l2/2)
    
        l3 = h * g(tau + h/2, theta + k2/2, z + l2/2)
    
        k4 = h * f(tau+h, theta+k3, z+l3)
    
        l4 = h * g(tau+h, theta+k3, z+l3)
        
        
        ## Plugging into equation:
        
        theta = theta + 1/6 * (k1 + 2*k2 + 2*k3 + k4)       # This equation outputs a point Y(x+h) on the "graph" and updates the y-value
        
        Y[i] = theta                                    # Stores the point y in the array Y to be plotted later
        
        z = z + 1/6 * (l1 + 2*l2 + 2*l3 + l4)       # This equation updates the z-value
        
        X[i] = tau
        
        i+=1                                        # Updates the index of the array Y and X
        
        tau+=h                                        # x+h incriments
       
    #if np.amax(Y) > 1.0505864368267241:        #0.47051971068636994:     #0.1422079355928062:                                 # Plots the solution curve
        
    return Y 

a=0.01
l=1
eps = np.linspace( 0.01, 1.5 ,200 )
nuA = np.linspace( 0.75,1.25,50 )
StorEps = []
StorNu = []
count =0
for i in nuA:
    
    count+=1
    print(count)
    for j in eps:
        
        #RK4(i,j)
        
        if np.amax(RK4(i,j)) > 0.04:
            StorEps.append(j)
            StorNu.append(i)
            break

#nuA2 = np.linspace( 0.5,0,100 )
#nuA1 = np.linspace( 0,1.5,100 )
epsil1 = abs(2 * nuA - 2)

plt.plot(StorNu,StorEps,label = 'Numerical Solution for $\epsilon_1(\\nu)$')
plt.plot( nuA, epsil1,label = 'Analytical Solution for $\epsilon_1(\\nu)$' )
plt.xlabel('$\\nu$')
plt.ylabel('$\epsilon_1(\\nu)$')
plt.legend()
plt.grid()
for k in range(50):
    print(StorNu[k],StorEps[k])