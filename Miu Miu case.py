# 
#
#
# Fourth Order Runge-Kutta method for Python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


def f(tau,theta,z):                                             # z=dy/dx
    return z

def g(tau,theta,z):                                             # dz/dx = d^2y/dx^2=...(solved for dz/dx)
    ## Constants:
    a = 0.01
    l = 1
    nu = 4 + ((0.04)**2)*(5/48)-0.00001
    epsilon = 0.04 #nu*(a/l)
    
    ## Equation: 
    dz = - nu*theta - epsilon*np.cos(2*tau)*theta
    return dz

def RK4(tau,theta,z,h,N):                                       # RK4 method with x,y,z being intial conditions, h is the incriments, and N is the number that you want x to go to
    
    # Initialise variables:
    i = 0
    Y = np.zeros(int(N*(1/h))+1)
    X = np.zeros(int(N*(1/h))+1)
    

    ## Solving the equation:
    while tau<N:                                                # N is the nnumber of the x-value that you wish to go to
        
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
        
        theta = theta + 1/6 * (k1 + 2*k2 + 2*k3 + k4)           # This equation outputs a point Y(x+h) on the "graph" and updates the y-value
        
        Y[i] = theta                                            # Stores the point y in the array Y to be plotted later
        
        z = z + 1/6 * (l1 + 2*l2 + 2*l3 + l4)                   # This equation updates the z-value
        
        X[i] = tau
        
        i+=1                                                    # Updates the index of the array Y and X
        
        tau+=h                                                  # x+h incriments
        
    plt.plot(X,Y)                                               # Plots the solution curve
    plt.xlabel('$\\tau$')
    plt.ylabel('$\\theta$')
    plt.grid()
    print(np.amax(Y))



## Plotting ------------------------------------------------------------------------------------------------------------------------------------------

RK4(0,0.03,0,0.01,50000)