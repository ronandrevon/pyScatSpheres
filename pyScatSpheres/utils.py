import os,pandas as pd,numpy as np
from . import qdot_sphere_array as qsa #;imp.reload(qsa)
from utils import glob_colors as colors

def get_param_list(params):
    ns = np.unique([np.array(p).size for p in params if type(p) in [np.ndarray,list]])
    if not any(ns):return None
    ns=ns[0]
    for i,p in enumerate(params) :
        if not type(p) in [np.ndarray,list]: p = np.array([p]*ns)
        params[i] = np.array(p)
    return params

def solve_set(df_name,kas,kps,kdkas,Ns,new=0):
    cols = ['N','ka','kp','kd','nmax','ap','bp']
    cols_add = ['bp0','err_0','bpa','err_a','bp2','err_2']
    df = pd.DataFrame(columns=cols+cols_add)
    thetas = np.deg2rad(np.linspace(0,60,361))
    dt = thetas[1]-thetas[0]
    for i,(ka,kp,kdka,N) in enumerate(zip(kas,kps,kdkas,Ns)):
        nmax = max(int(np.ceil(1.3*ka)),int(ka)+4)
        print('%d/%d' %(i,kas.size) +
            colors.yellow+' nmax=%d, N=%d, ndofs=(nmax+1)*N*2=%d' %(nmax,N,2*nmax*N)+colors.black)#;print(kp,kd,nmax))
        kd = kdka*ka  #;print(N,ka,kd,kp,nmax)
        s = qsa.QdotSphereArray(N=N,ka=ka,kp=kp,kd=kd,nmax=nmax,solve=1,copt=N>1,v=0)
        df.loc[i,cols]=[N,ka,kp,kd,nmax,s.ap,s.bp] #;print(df.loc[i])
        #uncoupled
        ap,bp0 = s.get_cp0()#solve(copt=0,v=0)
        df.loc[i,'bp0'] = bp0
        #Secondary
        ap,bp2 = s.get_cp2()#solve(copt=0,v=0)
        # print(bp0.shape,bp2.shape)
        df.loc[i,'bp2'] = bp2
        #forward coupling
        ap,bpa = s.get_cpa()#solve(copt=0,v=0)
        df.loc[i,'bpa'] = bpa
    if os.path.exists(df_name) and not new:
        df0 = pd.read_pickle(df_name)
        df.append(df0)#;print(df.index)
        # print('updating df')
    df.to_pickle(df_name)
    # print(df.iloc[0],s.bp)
    print(colors.green+'DataFrame saved:\n'+ colors.yellow+df_name+colors.black)
    return df
