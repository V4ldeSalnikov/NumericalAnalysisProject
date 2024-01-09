#!/usr/bin/env python
# coding: utf-8

# In[34]:


import numpy as np
import pandas as pd
def f(y):
    g = 9.81
    L = 1.0
    theta, v = y
    return np.array([v, -g / L * np.sin(theta)])

def runge_kutta_4(y0, t0, tf, h):
    t = np.arange(t0, tf, h)
    np.append(y,y0)
    for i in range(1, len(t)):
        k1 = h * f(y[i-1])
        k2 = h * f(y[i-1] + 0.5 * k1)
        k3 = h * f(y[i-1] + 0.5 * k2)
        k4 = h * f(y[i-1] + k3)
        np.append(y, y[i-1] + (k1 + 2*k2 + 2*k3 + k4) / 6)

    return y

g = 9.81  
L = 1.0 
y0 = np.array([np.pi/2, 0])
t0 = 0    
tf = 10 
h = 0.1  
y = runge_kutta_4(y0, t0, tf, h)
df = pd.DataFrame(y, columns =['theta', 'omega']) 


# In[35]:


df


# In[33]:


df.plot("omega", "theta")


# In[23]:


df['time'] = pd.Series(pd.np.arange(0, 100, 0.1))


# In[24]:


df.plot("time", "theta")


# In[27]:


df.plot("time", "omega")


# In[ ]:




