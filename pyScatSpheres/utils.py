import importlib as imp
import os,pandas as pd,numpy as np
from utils import glob_colors as colors
from . import qdot_sphere_array as qsa ;imp.reload(qsa)


def load_qdot(df,cond,opts='co'):
    '''opts:returns c(df_c) o(qdot object)'''
    if isinstance(cond,dict):
        cond = ' & '.join(['(%s==%.4f)' %(k,v) for k,v in cond.items()])#; print(cond)
    # print(df.loc[df.eval(cond)])
    c = df.loc[df.eval(cond)].iloc[0]

    ret=[]
    for i in opts:
        if   'c'==i:ret+=[c]
        elif 'o'==i:
            qdot = qsa.QdotSphereArray(N=c.N,ka=c.ka,kp=c.kp,kd=c.kd,nmax=c.nmax,solve=0)
            qdot.ap,qdot.bp= c.ap.copy(),c.bp.copy()
            ret+=[qdot]
    if len(opts)<2:ret=ret[0]
    return ret

def show_ff_multiscat(df,cond,npts=2001,**kwargs):

    c = df.loc[df.eval(cond)].iloc[0]
    qdot2   = qsa.QdotSphereArray(N=c.N,ka=c.ka,kp=c.kp,kd=c.kd,nmax=c.nmax,solve=0)

    theta_d = np.linspace(0,30,npts)
    theta   = np.deg2rad(theta_d)
    qdot2.bp = c.bp
    ff =  qdot2.get_ff(theta)

    # qdot2.bp = c.bp0
    # ff0 = qdot2.get_ff(theta)
    qdot2.bp = c.bp2
    ff2 = qdot2.get_ff(theta)
    qdot2.bp = c.bpa
    ffa = qdot2.get_ff(theta)
    ffmax = abs(ff).max()

    err = lambda ff0:ff0
    err = lambda ff0:np.log10(abs(abs(ff0)-abs(ff))/abs(ff))
    plts=[]
    # plts =[[theta_d,np.log10(abs(ff)) ,'b-' ,'exact']]
    plts+=[[theta_d,err(ffa),'g--','forward']]
    plts+=[[theta_d,err(ff2),'c--','secondary']]
    plts+=[[theta_d,err(ff0),'r--','uncoupled']]
    tle = r'$N=%d, ka=%d, k_p=%.4f, kd=%d$ ' %(c.N,c.ka,c.kp,c.kd)
    dsp.stddisp(plts,labs=[r"$\theta(^{\circ})$",r"$f(\theta)$"],#title=tle,
        lw=2,**kwargs)


def set_errors(df_name,err='bp'):
    '''err:
        - 'bp' : sum(bpa-bp)
        - 'ff' : sum(ff-cp)
        - 'fm' : error at the peaks
    '''
    df = pd.read_pickle(df_name)
    thetas_d = np.linspace(0,180,180*10+1)
    thetas = np.deg2rad(thetas_d)
    dt = thetas[1]-thetas[0]
    ns  = len([col for col in df.columns if 'bp_' in col])
    cols_bp  = ['bp_%d' %(i+1) for i in range(ns)]
    cols_err = ['err_%d' %(i+1) for i in range(ns)]
    for i,c in df.iterrows():
        s = qsa.QdotSphereArray(N=c.N,ka=c.ka,kp=c.kp,kd=c.kd,nmax=c.nmax-1,solve=0,copt=1,v=0)
        if err=='bp':
            # df.loc[i,'err_a'] = abs(c.bp-c.bpa).sum()
            bpm = abs(c.bp).sum()
            df.loc[i,'erra'] = abs(c.bp-c.bpa).sum()/bpm
            df.loc[i,'err0'] = abs(c.bp-c.bp0).sum()/bpm
            df.loc[i,'err2'] = abs(c.bp-c.bp2).sum()/bpm
            #bp_n : n x (N*nmax*2)
            bp_n = np.array([np.array(b,dtype=complex) for b in c[cols_bp]])
            bps = np.cumsum(bp_n,axis=0)
            # print(bps.shape)
            df.loc[i,cols_err] = abs(c.bp-bps).sum(axis=1)/bpm
        else:
            s.bp = c.bp
            ff = s.get_ff(thetas)
            ff2 = abs(ff)**2
            ffmax = abs(ff).max()

            idx = find_peaks(abs(ff))#;print(idx);
            idx = idx[0]
            theta_i = thetas[idx]
            ff_2i = ff2[idx]

            s.bp = c.bpa
            df.loc[i,'erra'] = (abs(abs(s.get_ff(theta_i))**2-ff_2i)/ff_2i).mean()
            # df.loc[i,'err_a'] = abs(s.get_ff(thetas)-ff).sum()*dt/ffmax
            s.bp = c.bp2
            df.loc[i,'err2'] = (abs(abs(s.get_ff(theta_i))**2-ff_2i)/ff_2i).mean()
            # df.loc[i,'err_2'] = abs(s.get_ff(thetas)-ff).sum()*dt/ffmax
            s.bp = c.bp0
            df.loc[i,'err0'] = (abs(abs(s.get_ff(theta_i))**2-ff_2i)/ff_2i).mean()
            # df.loc[i,'err_0'] = abs(s.get_ff(thetas)-ff).sum()*dt/ffmax

        # if df.loc[i,'err_0']<8e-3 and df.loc[i,'err_0']>3e-3:
        #     tle = 'ka=%d,kp=%.4f,kd=%d,nmax=%d,err=%.2E' %(s.ka,s.kp,s.kd,s.nmax,df.loc[i,'err_0'])
        #     plts=[[thetas_d,abs(ff),'b-','exact'],[thetas_d,abs(ffa),'r--','forward'],[thetas_d,abs(ff0),'g-.','uncoupled'],]
        #     dsp.stddisp(plts,title=tle,lw=2);plt.show()
    df.to_pickle(df_name)
    print(colors.green+'DataFrame updated with errors:\n'+ colors.yellow+df_name+colors.black)
    return df


def get_param_list(params):
    ns = np.unique([np.array(p).size for p in params if type(p) in [np.ndarray,list]])
    if not any(ns):return None
    ns=ns[0]
    for i,p in enumerate(params) :
        if not type(p) in [np.ndarray,list]: p = np.array([p]*ns)
        params[i] = np.array(p)
    return params

def solve_set(df_name,kas,kps,kdkas,Ns,new=0,n=2):
    cols = ['N','ka','kp','kd','nmax','ap','bp']
    cols_add = ['bpa','bp2','bp0']#,'err_0','err_a']
    bpn_cols = ['bp_%d' %(i+1) for i in range(n)]
    df = pd.DataFrame(columns=cols+cols_add+bpn_cols)
    # thetas = np.deg2rad(np.linspace(0,60,361))
    # dt = thetas[1]-thetas[0]
    for i,(ka,kp,kdka,N) in enumerate(zip(kas,kps,kdkas,Ns)):
        nmax = max(int(np.ceil(1.3*ka)),int(ka)+4)
        print('%d/%d' %(i,kas.size) +
            colors.yellow+' nmax=%d, N=%d, ndofs=(nmax+1)*N*2=%d' %(nmax,N,2*nmax*N)+colors.black)#;print(kp,kd,nmax))
        kd = kdka*ka  #;print(N,ka,kd,kp,nmax)
        s = qsa.QdotSphereArray(N=N,ka=ka,kp=kp,kd=kd,nmax=nmax,solve=1,copt=N>1,v=0)
        df.loc[i,cols]=[N,ka,kp,kd,nmax,s.ap,s.bp] #;print(df.loc[i])
        # #uncoupled
        # ap,bp0 = s.get_cp0()#solve(copt=0,v=0)
        #forward coupling
        # ap,bpa = s.get_cpa()#solve(copt=0,v=0)
        # ap,bpa = s.get_cpa()#solve(copt=0,v=0)
        df.loc[i,'bpa'] = s.get_cpa()[1]
        df.loc[i,'bp2'] = s.get_cp2()[1]
        df.loc[i,'bp0'] = s.get_cp0()[1]
        #multiple
        # apn,bpn = s.get_cpn(n=n)#solve(copt=0,v=0)
        # print(bp0.shape,bp2.shape)
        df.loc[i,bpn_cols] = s.get_cpn(n=n)[1]
    if os.path.exists(df_name) and not new:
        df0 = pd.read_pickle(df_name)
        df=df.append(df0)#;print(df.index)
        print('updating df')
    df.to_pickle(df_name)
    # print(df.iloc[0],s.bp)
    print(colors.green+'DataFrame saved:\n'+ colors.yellow+df_name+colors.black)
    return df
