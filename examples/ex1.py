'''
**************************************************
Mulitlayer Chemical Vapor Deposition Diffusion Analytical solution
Date:2024-11-14
Author:zhanggang
**************************************************
'''
from  diffusion.Diffusion import ML_CVD_Model
import numpy as np
import matplotlib.pyplot as plt

# 1 layer 
filepath = 'para_1layer_exp.yaml'
model = ML_CVD_Model()
model(filepath)
x = model.x_positon
c1 = model.c1_res
c2 = model.c2_res

data_exp = np.loadtxt('N_P_Constant.csv',delimiter=',',skiprows=2)

plt.scatter(data_exp[:,0],data_exp[:,1],s=1)
plt.scatter(data_exp[:,0],data_exp[:,2],s=1)
plt.plot(x,c1)
plt.plot(x,c2)
plt.ylim(5e9,3.0e16)
plt.yscale('log')
plt.show()

# 3 layer
filepath = 'para_3layer.yaml'
model(filepath)
x = model.x_positon
c1 = model.c1_res
c2 = model.c2_res

plt.plot(x,c1)
plt.plot(x,c2)
plt.ylim(5e9,3.0e16)
plt.yscale('log')
plt.show()