#!/usr/bin/env python
# coding: utf-8

# # Initial conditions
import numpy as np
theta_0 = np.pi /2  # Initial angle
omega_0 = 0.0 # Initial angular velocity
h = 0.1 #Step size
y = np.zeros((1000, 2))
y[0][0] = theta_0
y[0][1] = omega_0


# 1st step of Integration using Implicit Euler method
def Jacobian_Euler(y):
    g = 9.81
    h = 0.1
    l = 1.0
    return np.array([
    [-1, h],
    [-h * (g / l) * np.cos(y[0]), -1]])

def negative(y):
    x = []
    for j in y:
        if j != 0:
            x.append(-j)
        else :
            x.append(0)
    return x
            
def F_Euler(y, guess):

    g = 9.81
    h = 0.1
    l = 1.0
    
    theta_0 = y[0][0]
    omega_0 = y[0][1]
    
    theta_1, omega_1 = guess[0], guess[1]
    eq1 = theta_0 + h * omega_1 - theta_1
    eq2 = omega_0 - h * (g / l) * np.sin(theta_1) - omega_1
    return [eq1, eq2]
        
def Newton_system_Euler(y,guess): 
    F_value = F_Euler(y,guess)
    eps = 1e-7
    F_norm = np.linalg.norm(F_value, ord=2)  # l2 norm of vector
    iteration_counter = 0
    while abs(F_norm) > eps and iteration_counter < 100:
        
        delta = np.linalg.solve(Jacobian_Euler(guess), negative(F_value))
        guess = guess + delta
        F_value = F_Euler(y, guess)
        F_norm = np.linalg.norm(F_value, ord=2)
        iteration_counter += 1
       
    # Here, either a solution is found, or too many iterations
    if abs(F_norm) > eps:
        iteration_counter = -1
    theta_1, omega_1 = guess[0], guess[1]
    return theta_1, omega_1
theta_1, omega_1 = Newton_system_Euler(y,[theta_0 + h * omega_0,omega_0]) #Solving first system of equations
y[1][0] = theta_1
y[1][1] = omega_1
print(y[1][0],y[1][1])


# Further Integration Using BDF-2 Method
def F_BDF_2(i,y,guess): #calculating value of function
    
    g = 9.81
    h = 0.1
    l = 1.0
    theta_np1, omega_np1 = guess[0], guess[1]
    #Newton_system_ImplicitEuler(y)
    f1 = -theta_np1 + (4/3)*y[i][0] - (1/3)*y[i-1][0] + (2/3)*h * omega_np1
    f2 = -omega_np1 + (4/3)*y[i][1] - (1/3)*y[i-1][1] - (2/3)*h * (g / l) * np.sin(theta_np1)
    return [f1, f2]

def Jacobian_BDF_2(y):
    g = 9.81
    h = 0.1
    l = 1.0
    return np.array([[-1, (2/3)*h],
                      [-(2/3)*h*(g/l)*np.cos(y[0]), -1]])

def Newton_system_BDF(y,guess):
    eps = 1e-6
    F_value = F_BDF_2(i,y,guess)
    F_norm = np.linalg.norm(F_value, ord=2)  # l2 norm of vector
    iteration_counter = 0  
    while abs(F_norm) > eps and iteration_counter < 1000:
        delta = np.linalg.solve(Jacobian_BDF_2(guess),negative(F_value))
        guess = guess + delta
        F_value = F_BDF_2(i,y,guess)
        F_norm = np.linalg.norm(F_value, ord=2)
        iteration_counter += 1

    # Here, either a solution is found, or too many iterations
    if abs(F_norm) > eps:
        iteration_counter = -1
    theta_np1, omega_np1 = guess[0], guess[1]
    return theta_np1, omega_np1

g = 9.81
h = 0.1
l = 1.0
for i in range(1,999):
    theta, omega = Newton_system_BDF(y, [y[i][0]+ h * y[i][1],y[i][1] - h*g/l*sin(y[i][0])])
    y[i+1][0] = theta
    y[i+1][1] = omega

import pandas as pd #Creating data frame for results of simulation without obstacle
df = pd.DataFrame(y, columns =['theta', 'omega']) 

df.head()
df['time'] = pd.Series(pd.np.arange(0, 1000, 0.1)) #initializing column with time values
df.iloc[1:100].plot(x='time', y='omega', figsize=(10, 5)) 
plt.plot(df["theta"].iloc[0:500], df["omega"].iloc[0:500], 'r^-')
df.iloc[1:100].plot(x='time', y=['theta', 'omega'])
df.iloc[900:1000].plot(x='theta', y = "omega")


#Lagrange polynomial

import numpy as np
from numpy import (array, transpose, searchsorted, atleast_1d, atleast_2d,
                   ravel, poly1d, asarray, intp)

def lagrange(x, w):
    p = poly1d(0.0)
    for j in range(len(x)):
        pt = poly1d(w[j])
        for k in range(len(x)):
            if k == j:
                continue
            fac = x[j]-x[k]
            pt *= np.poly1d([1.0, -x[k]])/fac
        p += pt
    return p

t = np.arange(0, 1000, 0.1)
array_lagrange_theta = [0]
array_lagrange_omega = [0]


#Creating interpolation polynomials for every pair of 3 solution points (theta,omega)
from numpy.polynomial.polynomial import Polynomial
for i in range(0,999):
    if i == 999: 
        array_lagrange_theta.append(0)
        array_lagrange_omega.append(0)
    else :
        x_current = [t[i-1], t[i], t[i+1]]
        theta_current = [y[i-1][0], y[i][0], y[i+1][0]]
        omega_current = [y[i-1][1], y[i][1], y[i+1][1]]
        array_lagrange_theta.append(Polynomial(lagrange(x_current, theta_current)).coef)
        array_lagrange_omega.append(Polynomial(lagrange(x_current, omega_current)).coef)

#Creating dataframe with polynomials for theta and omega
df["poly of theta"] = array_lagrange_theta
df["poly of omega"] = array_lagrange_omega

def is_solution_check(coef,t_n,t_np1):
    a = np.polyval(coef, t_n) + np.pi/3
    b = np.polyval(coef, t_np1) + np.pi/3
    if (a*b > 0):
        return False
    else :
        return True


# # Finding root of lagrange polynomial


def is_solution_check(coef,t_n,t_np1):
    a = np.polyval(coef, t_n) + np.pi/3
    b = np.polyval(coef, t_np1) + np.pi/3
    if (a*b > 0):
        return False
    else :
        return True


# For every time t[i] and every polynomial theta[i] checking if solution exist
for i in range (1, 1000):
    if(is_solution_check(array_lagrange_theta[i], t[i],t[i+1])):
        print(array_lagrange_theta[i])
        print(t[i])
        print(i)
        break


def func(t):
    return np.polyval(array_lagrange_theta[9], t) + np.pi/3

def bisection_method(func, low, high, tolerance=1e-6, max_iterations=1000):

    for _ in range(max_iterations):
        midpoint = (low + high) / 2
        f_mid = func(midpoint)
        if abs(f_mid) < tolerance:
            return midpoint
        elif func(low) * f_mid < 0:
            high = midpoint
        else:
            low = midpoint

    return (low + high) / 2


# Definition of interval of solution
tn = 0.8  
tn_plus_1 = 0.9

# Find the root using the bisection method
tobst = bisection_method(func, tn, tn_plus_1)

if tobst is not None:
    print("The solution tobst is:", tobst)
else:
    print("No solution found in the interval.")


np.polyval(array_lagrange_theta[8], 0.877) + np.pi/3

from scipy.optimize import brentq
def func(t):
    return np.polyval(array_lagrange_theta[8], t) + np.pi/3

t_obst = brentq(func, 0.8, 0.9)
print(t_obst)


theta_obstacle = -np.pi/3
omega_obstacle = np.polyval(array_lagrange_omega[8],tobst)
print(omega_obstacle)


# Restarting the integration
theta_new_0 = -np.pi /3  # Initial angle
omega_new_0 = -omega_obstacle # Initial angular velocity
h = 0.1
y_new = np.zeros((1000, 2))
y_new[0][0] = theta_new_0
y_new[0][1] = omega_new_0

def Jacobian_Euler(y_new):
    g = 9.81
    h = 0.1
    l = 1.0
    return np.array([
    [-1, h],
    [-h * (g / l) * np.cos(y_new[0]), -1]])

def negative(y):
    x = []
    for j in y:
        if j != 0:
            x.append(-j)
        else :
            x.append(0)
    return x
            
def F_Euler_new(y_new, guess):

    g = 9.81
    h = 0.1
    l = 1.0
    
    theta_0 = y_new[0][0]
    omega_0 = y_new[0][1]
    
    theta_1, omega_1 = guess[0], guess[1]
    eq1 = theta_0 + h * omega_1 - theta_1
    eq2 = omega_0 - h * (g / l) * np.sin(theta_1) - omega_1
    return [eq1, eq2]


def Newton_system_Euler(y_new,guess): 
    F_value = F_Euler_new(y_new,guess)
    eps = 1e-7
    F_norm = np.linalg.norm(F_value, ord=2)  # l2 norm of vector
    iteration_counter = 0
    while abs(F_norm) > eps and iteration_counter < 100:
        
        delta = np.linalg.solve(Jacobian_Euler(guess), negative(F_value))
        guess = guess + delta
        F_value = F_Euler(y_new, guess)
        F_norm = np.linalg.norm(F_value, ord=2)
        iteration_counter += 1
       

    # Here, either a solution is found, or too many iterations
    if abs(F_norm) > eps:
        iteration_counter = -1
    theta_1, omega_1 = guess[0], guess[1]
    return theta_1, omega_1

theta_new_1, omega_new_1 = Newton_system_Euler(y_new,[theta_new_0 + h * omega_new_0,omega_new_0])
y_new[1][0] = theta_new_1
y_new[1][1] = omega_new_1
print(y_new[1][0],y_new[1][1])


# # Restarting integration usign BDF-2

def F_BDF_2_new(i,y_new,guess):
    
    g = 9.81
    h = 0.1
    l = 1.0
    theta_np1, omega_np1 = guess[0], guess[1]
    #Newton_system_ImplicitEuler(y)
    f1 = -theta_np1 + (4/3)*y_new[i][0] - (1/3)*y_new[i-1][0] + (2/3)*h * omega_np1
    f2 = -omega_np1 + (4/3)*y_new[i][1] - (1/3)*y_new[i-1][1] - (2/3)*h * (g / l) * np.sin(theta_np1)
    return [f1, f2]

def Jacobian_BDF_2(y_new):
    g = 9.81
    h = 0.1
    l = 1.0
    return np.array([[-1, (2/3)*h],
                      [-(2/3)*h*(g/l)*np.cos(y_new[0]), -1]])

def Newton_system_BDF(y_new,guess):
    eps = 1e-6
    F_value = F_BDF_2_new(i,y_new,guess)
    F_norm = np.linalg.norm(F_value, ord=2)  # l2 norm of vector
    iteration_counter = 0  
    while abs(F_norm) > eps and iteration_counter < 1000:
        delta = np.linalg.solve(Jacobian_BDF_2(guess),negative(F_value))
        guess = guess + delta
        F_value = F_BDF_2_new(i,y_new,guess)
        F_norm = np.linalg.norm(F_value, ord=2)
        iteration_counter += 1

    # Here, either a solution is found, or too many iterations
    if abs(F_norm) > eps:
        iteration_counter = -1
    theta_np1, omega_np1 = guess[0], guess[1]
    return theta_np1, omega_np1

#Solving Bdf-2 using Newton
for i in range(1,1000):
    theta, omega = Newton_system_BDF(y_new, [y_new[i-1][0]+ h * y_new[i-1][1],y_new[i-1][1]])
    y_new[i+1][0] = theta
    y_new[i+1][1] = omega
    print(y_new[i+1][0],y_new[i+1][1])

df_new['time'] = pd.Series(pd.np.arange(t_obst, 992, 0.1))
df['time'] = pd.Series(pd.np.arange(0, 1000, 0.1))
#Plot after collision
plt.figure(figsize=(15,12))
plt.plot(df_new["theta_new"], df_new["omega_new"], 'r^-')

df_new.rename(columns = {'omega_new':'omega'}, inplace = True)
df_new

from pylab import *

df_collision = df.iloc[0:9]
concatenated_df = pd.concat([df_collision, df_new])

#Plots
concatenated_df.reset_index(drop=True, inplace=True)
concatenated_df.iloc[0:75].plot(x = "theta", y = "omega")

ax = concatenated_df.iloc[0:75].plot(x = "theta", y = "omega", label='with collision')
df.iloc[0:75].plot(x = "theta", y = "omega", label='without collision', ax = ax)

