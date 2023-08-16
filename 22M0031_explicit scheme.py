import numpy as np 
from math import *
import matplotlib.pyplot as plt

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

def explicit_scheme(u_0, a, CFL, delta_x, t_final, scheme):
    u_n = u_0.copy()
    u_np1 = u_0.copy()
    t=0.0
    delta_t = CFL * delta_x / a
    while t <= t_final:
        for i in range(1,len(u_0)-1):
            if scheme == "FTFS":
                u_np1[i] = u_n[i] - (CFL * (u_n[i+1] - u_n[i]))
            elif scheme == "FTCS":
                u_np1[i] = u_n[i] - ((CFL/2) * (u_n[i+1] - u_n[i-1]))
            elif scheme == "FTBS":
                u_np1[i] = u_n[i] - (CFL * (u_n[i] - u_n[i-1]))
            elif scheme == "LW":
                u_np1[i] = u_n[i] - ((CFL/2) * (u_n[i+1] - u_n[i-1]))+((pow(CFL, 2)/2)*(u_n[i+1]-2*u_n[i]+u_n[i-1]))
            elif scheme == "BW":
               u_np1[i] = u_n[i] - ((CFL/2) * (3*u_n[i] - 4*u_n[i-1] + u_n[i-2])) + ((pow(CFL, 2)/2)*(u_n[i]-2*u_n[i-1]+u_n[i-2]))
            elif scheme == "FR":
                temp_LW= u_n[i] - ((CFL/2) * (u_n[i+1] - u_n[i-1]))+((pow(CFL, 2)/2)*(u_n[i+1]-2*u_n[i]+u_n[i-1]))
                temp_BW= u_n[i] - ((CFL/2) * (3*u_n[i] - 4*u_n[i-1]+u_n[i-2]))+((pow(CFL, 2)/2)*(u_n[i]-2*u_n[i-1]+u_n[i-2]))
                u_np1[i] = 0.5*(temp_LW+temp_BW)
        u_n = u_np1.copy()
        t += delta_t
    return u_n

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
#scheme = input("\nEnter Scheme FTFS/FTCS/FTBS\n")
x = np.round(np.arange(0,l+delta_x,delta_x),2)

U_0 = initialization(x, ic)
U_f_FTCS = explicit_scheme(u_0=U_0, a=a, CFL=nu, delta_x=delta_x, t_final=t_final, scheme="FTCS")
U_f_FTBS = explicit_scheme(u_0=U_0, a=a, CFL=nu, delta_x=delta_x, t_final=t_final, scheme="FTBS")
U_f_LW = explicit_scheme(u_0=U_0, a=a, CFL=nu, delta_x=delta_x, t_final=t_final, scheme="LW")
U_f_BW = explicit_scheme(u_0=U_0, a=a, CFL=nu, delta_x=delta_x, t_final=t_final, scheme="BW")
U_f_FR = explicit_scheme(u_0=U_0, a=a, CFL=nu, delta_x=delta_x, t_final=t_final, scheme="FR")
Fig1 = plt.figure()
plt.plot(x, U_0, '-')
plt.plot(x, U_f_FTCS, '-')
plt.plot(x, U_f_FTBS, '--')
plt.plot(x, U_f_LW, '--')
plt.plot(x, U_f_BW, '--')
plt.plot(x, U_f_FR, '--')
plt.legend(['U_0','FTBS', 'LW', 'BW', 'FR'])
plt.xlabel('Length')
plt.ylabel('U')
plt.xlim(0,1)
plt.ylim(-1.5,1.5)
plt.grid()
plt.show()




