import numpy as np 
from math import *
import matplotlib.pyplot as plt
from scipy.sparse import diags

#initial condition 1
def init_condition_1(x):
    u_n = []
    for temp_x in x:
        if temp_x <= 0.2:
            u_n.append(1)
        else:
            u_n.append(0)
    return u_n

def init_condition_2(x):
    u_n = []
    for temp_x in x:
        if temp_x < 0.05:
            u_n.append(0)
        elif temp_x < 0.35:
            u_n.append(sin(4*np.pi*((temp_x-0.05)/0.3)))
        else:
            u_n.append(0)
    return u_n

def init_condition_3(x):
    u_n = []
    for temp_x in x:
        if temp_x < 0.05:
            u_n.append(0)
        elif temp_x < 0.35:
            u_n.append(sin(8*np.pi*((temp_x-0.05)/0.3)))
        else:
            u_n.append(0)
    return u_n

def init_condition_4(x):
    u_n = []
    for temp_x in x:
        if temp_x < 0.05:
            u_n.append(0)
        elif temp_x < 0.35:
            u_n.append(sin(12*np.pi*((temp_x-0.05)/0.3)))
        else:
            u_n.append(0)
    return u_n

def init_condition_5(x):
    u_n = []
    for temp_x in x:
        u_n.append(np.exp(-50.0*np.power(temp_x-0.2,2)/0.16))
    return u_n

def implicit_scheme(u_0, a, CFL, delta_x, t_final):
    diag_main=np.zeros((99))
    diag_upper=np.zeros((98,1))
    diag_lower=np.zeros((98,1))
    RHS=np.zeros((99,1))
    u_n = u_0.copy()
    u_np1 = u_0.copy()
    t=0.0
    delta_t = CFL * delta_x / a
    for i in range(len(u_0)-2):
        diag_main[i] = 1
    for i in range(len(u_0)-3):
        diag_upper[i] = CFL / 2
        diag_lower[i] = -CFL / 2

    RHS[0] = u_n[1] + ((CFL / 2) * u_n[0])
    RHS[98] = u_n[99] - ((CFL / 2) * u_n[100])
        
    for i in range(1,len(u_0)-2):
        RHS[i] = u_n[i + 1]

    while t <= t_final:
        diagonals = [np.squeeze(diag_upper) ,np.squeeze(diag_main), np.squeeze(diag_lower)]
        A=diags(diagonals,[1,0,-1]).toarray()
        U = np.dot(np.linalg.inv(A), RHS )
        RHS[0] = u_n[1] + ((CFL / 2) * u_n[0])
        RHS[98] = u_n[99] - ((CFL / 2) * u_n[100])
        t += delta_t 

    return U 


        

def initialization(x,ic):
    match ic:
        case 1:
            return init_condition_1(x)
        case 2:
            return init_condition_2(x)
        case 3:
            return init_condition_3(x)
        case 4:
            return init_condition_4(x)
        case 5:
            return init_condition_5(x)

t_final=0.35
l=1
n=101
delta_x=l/(n-1)
a=1.0
nu = float(input("\nEnter CFL Number:\n"))
ic = float(input("\nEnter the initial condition: [1,2,3,4,5] \n"))
x = np.round(np.arange(0,1+delta_x,delta_x), 2)
U_0 = initialization(x, ic)
implicit_scheme(U_0, a, nu, delta_x, t_final)
print(np.round(np.arange(0.01,0.99+delta_x,delta_x),2)) 
U_f = implicit_scheme(u_0=U_0, a=a, CFL=nu, delta_x=delta_x, t_final=t_final)
Fig1 = plt.figure()
plt.plot(x, U_0, '-')
plt.plot(np.round(np.arange(0.01,0.99+delta_x,delta_x),2), U_f, '--')
plt.xlabel('Length')
plt.ylabel('U')
plt.xlim(0,1)
plt.ylim(-1.5,1.5)
plt.grid()
plt.show()
