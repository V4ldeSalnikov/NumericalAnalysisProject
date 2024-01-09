#!/usr/bin/env python
# coding: utf-8

# In[9]:


from numpy import linalg as LA
import numpy as np
def ImpEulerStep(f,y0,t0,h,tol):

    # Initialize
    t1 = t0 + h
    
    # Initial estimate
    f0 = f(t0,y0)
    Y1 = y0 + h*f0
    
    errNorm = LA.norm(Y1 - y0)
    
    # Fixed point iteration
    i = 0;
    while errNorm > h*tol:
        Y1old = Y1
        Y1 = y0 + h*f(t1,Y1) # y0 constant
        errNorm = LA.norm(Y1 - Y1old)
        i = i+1
        if(i>=1e5):
            print('Fixed point iteration terminated after %d steps with error %g\n',i,errNorm)
            break
    
    # Define solution at next timestep
    y1 = Y1
    
    return y1


# In[10]:


def ImpEulerSolve(f,y0,T,h,tol):

    # Initialize variables
    N = int(np.ceil((T[1]-T[0])/h))# number of steps
    h = (T[1]-T[0])/N # adjust step size to fit interval length
    d = len(y0) # dimension of solution
    t = np.linspace(T[0],T[1],N+1) # time grid
    y = np.zeros((d,N+1)) # solution
    y[:,0] = y0 # set initial value
    
    # Compute solution
    for j in range(0,N):
        y[:,j+1] = ImpEulerStep(f,y[:,j],t[j],h,tol)
      
    return np.array([t,y])


# In[14]:


def f(t,y):
    return np.array([y[1],-9.81/1.0 * np.sin(y[0])])

# Define problem and discretization parameters
y0 = np.array([np.pi/2,0])
T = np.array([0,5])
h = 0.1
tol = 1e-7

# Serie 8, Aufgabe 3, Part a)

print(ImpEulerSolve(f,y0,T,h,tol))


# In[23]:


import matplotlib.pyplot as plt
import pandas as pd

plt.plot(t[:],yImp[0,:],'r^-',label="Euler Implicit")


# In[17]:



plt.plot(yImp[0,:], yImp[1,:],'r^-',label="Euler Implicit")


# In[ ]:




