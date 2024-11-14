'''
**************************************************
Mulitlayer Chemical Vapor Deposition Diffusion Analytical solution
Date:2024-11-14
Author:zhanggang
**************************************************
'''
import numpy as np
import yaml
from pathlib import Path
from scipy.special import erf

class ML_CVD_Model:

    def __init__(self):
        self.c1_res = np.array([])
        self.c2_res = np.array([])
        self.x_positon = np.array([])
        self.D = []
    
    def __call__(self, filepath:str):
        file_y = Path(filepath)
        with open(file_y,'r',encoding='utf-8') as fr:
            para_data = yaml.load(fr,Loader=yaml.FullLoader)
        # set process and grid parameters
        layer = np.array(para_data['layer'],dtype=np.float64)
        dep_t = np.array(para_data['dep_t'],dtype=np.float64)
        c1 = np.array(para_data['c1'],dtype=np.float64)
        c2 = np.array(para_data['c2'],dtype=np.float64)

        d1_coff = para_data['d1_coff']
        d1_temp_ref = para_data['d1_temp_ref']
        d1_exp_c = para_data['d1_exp_c']
        d2_coff = para_data['d2_coff']
        d2_temp_ref = para_data['d2_temp_ref']
        d2_exp_c = para_data['d2_exp_c']

        step_temperature = np.array(para_data['step_temperature'],dtype=np.float64)
        step_time = np.array(para_data['step_time'],dtype=np.float64)

        layer = np.array([layer[:i].sum() for i in range(len(layer)+1)])
        layer = (layer-layer[1])

        def gen_mesh(layer:np.array):
            dx = 1e-4 #μm
            geo = layer-layer[1]
            mesh = []
            for i in range(1,len(geo)-1):
                x0 = (geo[i-1]+geo[i])/2
                xm = (geo[i]+geo[i+1])/2
                if(i==1):
                    x0 = geo[i-1]
                elif(i==(len(geo)-1)):
                    xm = geo[i+1]
                ele_num = int((xm-x0)/dx+1)
                mesh.append(np.linspace(x0,xm,ele_num,endpoint=True)[:-1])
            return mesh

        mesh = gen_mesh(layer)

        def d_int_process(dep_t,d1,d1_ref,d1c,d2,d2_ref,d2c):
            dt = 1e-3
            D = []
            global title
            for i in range(len(dep_t)-1):    
                dt_num = int((dep_t[-1]-dep_t[i])/dt+1)
                dt = (dep_t[-1]-dep_t[i])/(dt_num-1)
                time_int = np.linspace(dep_t[i],dep_t[-1],dt_num,endpoint=True)[:-1]
                d1_int = np.interp(time_int,step_time,step_temperature)
                d1_int = d1*np.exp(d1c*(d1_int-d1_ref))*dt
                d2_int = np.interp(time_int,step_time,step_temperature)
                d2_int = d2*np.exp(d2c*(d2_int-d2_ref))*dt
                D.append([d1_int.sum(),d2_int.sum()])
                
            return D

        D = d_int_process(dep_t,d1_coff,d1_temp_ref,d1_exp_c,d2_coff,d2_temp_ref,d2_exp_c)
        
        self.D = D
        c1_res = np.array([])
        c2_res = np.array([])
        x_positon = np.array([])

        for i in range(len(mesh)):
            c1_temp = 0.5*(c1[i]+c1[i+1])-0.5*(c1[i]-c1[i+1])*erf((mesh[i]-layer[i+1])/(2.0*D[i][0]**0.5))
            c1_res = np.hstack([c1_res,c1_temp])
            c2_temp = 0.5*(c2[i]+c2[i+1])-0.5*(c2[i]-c2[i+1])*erf((mesh[i]-layer[i+1])/(2.0*D[i][1]**0.5))
            c2_res = np.hstack([c2_res,c2_temp])
            
            x_positon = np.hstack([x_positon,mesh[i]])
        self.c1_res = c1_res
        self.c2_res = c2_res
        self.x_positon = x_positon


    