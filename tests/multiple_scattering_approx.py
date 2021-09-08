from utils import *                                 ;imp.reload(dsp)
import scipy,copy,os
import scipy.special as spe
from scipy.signal import find_peaks
from scipy.integrate import trapz,quad
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)
from pyScatSpheres import utils as ut ;imp.reload(ut)
# from pyScatSpheres import spherical_utils as      #;imp.reload(spu)
plt.close('all')
path='../docs/figures/'
nameB=path+'qdotSphereArray'
name=nameB+'Napprox_'

opts = '' #2(show error for N=2), s(solve) S(sigma_tot)
# name=nameB+'2approx_'
# df_name = 'data/df_qdotSpheres_2.pkl'
df_name = 'data/df_qdotSpheresN2.pkl'


# def solve_set(df_name,kas,kps,kdkas,Ns):
#     cols = ['N','ka','kp','kd','nmax','ap','bp']
#     cols_add = ['bp0','err_0','bpa','err_a','bp2','err_2']
#     df = pd.DataFrame(columns=cols+cols_add)
#     thetas = np.deg2rad(np.linspace(0,60,361))
#     dt = thetas[1]-thetas[0]
#     for i,(ka,kp,kdka,N) in enumerate(zip(kas,kps,kdkas,Ns)):
#         nmax = max(int(np.ceil(1.3*ka)),int(ka)+4)
#         print('%d/%d' %(i,kas.size) +
#             colors.yellow+' nmax=%d, N=%d, ndofs=(nmax+1)*N*2=%d' %(nmax,N,2*nmax*N)+colors.black)#;print(kp,kd,nmax))
#         kd = kdka*ka
#         s = qsa.QdotSphereArray(N=N,ka=ka,kp=kp,kd=kd,nmax=nmax,solve=1,copt=1,v=0)
#         df.loc[i,cols]=[N,ka,kp,kd,nmax,s.ap,s.bp]
#         #uncoupled
#         ap,bp0 = s.get_cp0()#solve(copt=0,v=0)
#         df.loc[i,'bp0'] = bp0
#         #Secondary
#         ap,bp2 = s.get_cp2()#solve(copt=0,v=0)
#         # print(bp0.shape,bp2.shape)
#         df.loc[i,'bp2'] = bp2
#         #forward coupling
#         ap,bpa = s.get_cpa()#solve(copt=0,v=0)
#         df.loc[i,'bpa'] = bpa
#     if os.path.exists(df_name):
#         df0 = pd.read_pickle(df_name)
#         df.append(df0)#;print(df.index)
#         # print('updating df')
#     df.to_pickle(df_name)
#     print(colors.green+'DataFrame saved:\n'+ colors.yellow+df_name+colors.black)
#     return df

def set_errors(df_name):
    df = pd.read_pickle(df_name)
    thetas_d = np.linspace(0,180,180*10+1)
    thetas = np.deg2rad(thetas_d)
    dt = thetas[1]-thetas[0]
    for i,c in df.iterrows():
        s = qsa.QdotSphereArray(N=c.N,ka=c.ka,kp=c.kp,kd=c.kd,nmax=c.nmax-1,solve=0,copt=1,v=0)
        s.bp = c.bp
        ff = s.get_ff(thetas)
        ff2 = abs(ff)**2
        ffmax = abs(ff).max()


        idx = find_peaks(abs(ff))#;print(idx);
        idx = idx[0]
        theta_i = thetas[idx]
        ff_2i = ff2[idx]

        s.bp = c.bpa
        df.loc[i,'err_a'] = (abs(abs(s.get_ff(theta_i))**2-ff_2i)/ff_2i).mean()
        # df.loc[i,'err_a'] = abs(s.get_ff(thetas)-ff).sum()*dt/ffmax
        s.bp = c.bp2
        df.loc[i,'err_2'] = (abs(abs(s.get_ff(theta_i))**2-ff_2i)/ff_2i).mean()
        # df.loc[i,'err_2'] = abs(s.get_ff(thetas)-ff).sum()*dt/ffmax
        s.bp = c.bp0
        df.loc[i,'err_0'] = (abs(abs(s.get_ff(theta_i))**2-ff_2i)/ff_2i).mean()
        # df.loc[i,'err_0'] = abs(s.get_ff(thetas)-ff).sum()*dt/ffmax

        # if df.loc[i,'err_0']<8e-3 and df.loc[i,'err_0']>3e-3:
        #     tle = 'ka=%d,kp=%.4f,kd=%d,nmax=%d,err=%.2E' %(s.ka,s.kp,s.kd,s.nmax,df.loc[i,'err_0'])
        #     plts=[[thetas_d,abs(ff),'b-','exact'],[thetas_d,abs(ffa),'r--','forward'],[thetas_d,abs(ff0),'g-.','uncoupled'],]
        #     dsp.stddisp(plts,title=tle,lw=2);plt.show()
    df.to_pickle(df_name)
    print(colors.green+'DataFrame updated with errors:\n'+ colors.yellow+df_name+colors.black)
    return df

def show_ff(df,cond,npts=2001,**kwargs):

    c = df.loc[df.eval(cond)].iloc[0]
    qdot2   = qsa.QdotSphereArray(N=c.N,ka=c.ka,kp=c.kp,kd=c.kd,nmax=c.nmax,solve=0)

    theta_d = np.linspace(0,30,npts)
    theta   = np.deg2rad(theta_d)
    qdot2.bp = c.bp
    ff =  qdot2.get_ff(theta)

    qdot2.bp = c.bp0
    ff0 = qdot2.get_ff(theta)
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
    dsp.stddisp(plts,title=tle,labs=[r"$\theta(^{\circ})$",r"$f(\theta)$"],lw=2,**kwargs)

def load_qdot(df,cond):
    if isinstance(cond,dict):
        cond = ' & '.join(['(%s==%.4f)' %(k,v) for k,v in cond.items()])
        # print(cond)
    c = df.loc[df.eval(cond)].iloc[0]
    qdot = qsa.QdotSphereArray(N=c.N,ka=c.ka,kp=c.kp,kd=c.kd,nmax=c.nmax,solve=0)
    qdot.ap,qdot.bp= c.ap.copy(),c.bp.copy()
    return qdot



if 's' in opts:
    # kas0   = np.array([0.5,1,2,3,4,5,7,10,15,20,30])
    # kdkas0 = np.array([2.1,3,5,10,25,20])
    # kps0   = np.array([1.2,1.1,1.01,1.001,1.0001])
    # kas,kdkas,kps = np.meshgrid(kas0,kdkas0,kps0)
    # kas,kdkas,kps =kas.flatten(),kdkas.flatten(),kps.flatten()
    # Ns=np.array([2]*kas.size,dtype=int)

    kas0   = np.array([5,7,10,15])
    kps0   = np.array([1.1,1.01,1.001,1.0001])
    Ns0    = np.array([2,3,4,6,8,10,15,20,30]) #np.stack([np.arange(2,11,1),np.arange(15,26,5)])
    kas,Ns,kps = np.meshgrid(kas0,Ns0,kps0)
    kas,Ns,kps = kas.flatten(),Ns.flatten(),kps.flatten()
    kdkas = np.array([3]*kas.size)

    df = solve_set(df_name,kas,kps,kdkas,Ns)



df = pd.read_pickle(df_name)
# if 'e' in opts:df = set_errors(df_name)

if 'F' in opts:
    ka,kp=15,1.01
    npts = 3600
    dfc = df.loc[(df.ka==ka) &(df.kp==kp)].copy()
    Ns = np.array(dfc.N.values,dtype=float)
    # cs = dsp.getCs('jet',nkas)
    dfc['sigma']=0
    dfc['sigma1']=0
    for l,c in dfc.iterrows():
        qdot = load_qdot(dfc,{'ka':c.ka,'kp':c.kp,'N':c.N})
        qdot.bp=c.bpa
        dfc.loc[l,'sigma'] = qdot.get_s(npts=npts,norm=True)
        qdot.bp=c.bp0
        dfc.loc[l,'sigma1'] = qdot.get_s(npts=npts,norm=True)

    S =dfc.sigma.max()*10
    plts = [[Ns,1-dfc.sigma/S,'r-','coh']]
    plts+= [[Ns,dfc.sigma1/S,'g-','kin']]
    plts+= [[Ns,(dfc.sigma-dfc.sigma1)/S,'b-','dyn']]
    # plts+= [[Ns,sigma2,'b-','dyn']]
    dsp.stddisp(plts,labs=['$N$',r'$\sigma$'],lw=2,name=nameB+'Pi.svg',opt='p')

if 'S' in opts:
    kas,kps = [5,10,15],[1.1,1.01,1.001]
    nkas = df.ka.unique().size
    cs = dsp.getCs('jet',nkas)
    npts = 3600
    for i,kp in enumerate(kps):
        plts=[]
        for j,ka in enumerate(kas):
            dfc = df.loc[(df.ka==ka) &(df.kp==kp)].copy()
            Ns = np.array(dfc.N.values,dtype=float)
            dfc['sigma']=0
            for l,c in dfc.iterrows():
                qdot = load_qdot(dfc,{'ka':c.ka,'kp':c.kp,'N':c.N})
                dfc.loc[l,'sigma'] = qdot.get_s(npts=npts,norm=True)
            plts+=[[Ns,dfc.sigma,[cs[j],'-o'],r'$ka=%d$' %ka]]

            qdot = qsa.QdotSphereArray(N=1,ka=c.ka,kp=c.kp,kd=c.kd,nmax=c.nmax,solve=1,copt=0)
            sig0 = qdot.get_s(npts=npts,norm=True)
            plts += [[Ns,Ns**2*sig0,[cs[j],'--o'],''],[Ns,Ns*sig0,[cs[j],'-.o'],'']]
            legElt = {r'$N^2\sigma_0$':'k--o',r'$N\sigma_0$':'k-.o'}
        tle = r'$k_p=%.4f$' %kp
        dsp.stddisp(plts,labs=['$N$',r'$\sigma_{tot}$'],lw=2,title=tle,legElt=legElt,
            name=nameB+'_sigmaN%d.svg' %i,opt='sc')


# thetas = np.deg2rad(np.linspace(0,180,2001))
# ff = qdot.get_ff(thetas)
# # idp = np.hstack([ p*qdot.nmax+np.arange(4,qdot.nmax) for p in range(qdot.N)])
# idp = abs(qdot.bp)<1e-4
# qdot.bp[idp]=0#; qdot.show_ff()
# ff0 = qdot.get_ff(thetas)
# dsp.stddisp([[thetas,np.log10(abs(ff0-ff)),'b-']])

# qdot.test_convergence(nmaxs=[5,7,10])
################################################################################
#### N=N
################################################################################
if 'N' in opts:
    optsN = 'f' #e(errors), f(ff)
    kas0=np.array(df.ka.unique(),dtype=float)
    kps0=np.array(df.kp.unique(),dtype=float)
    Ns0 =np.array(df.N.unique(),dtype=float)
    argsM = {'pOpt':'im','cmap':'RdYlGn_r'}
    ikp = 3
    kp = kps0[ikp]
    # kderrs,errs = np.zeros(kas0.shape),np.zeros(kas0.shape)
    if 'e' in optsN:
        dfkp = df.loc[df.kp==kp]
        nkas = kas0.size
        errs_a = np.zeros((nkas,Ns0.size))
        errs_0 = np.zeros((nkas,Ns0.size))
        errs_2 = np.zeros((nkas,Ns0.size))
        for i,ka in enumerate(kas0):
            dfka = dfkp.loc[dfkp.ka==ka]
            errs_0[i,:] = np.log10(np.array(dfka.err_0.values,dtype=float))
            errs_2[i,:] = np.log10(np.array(dfka.err_2.values,dtype=float))
            errs_a[i,:] = np.log10(np.array(dfka.err_a.values,dtype=float))

        xylims=['y',min(errs_0.min(),errs_a.min()),max(errs_0.max(),errs_a.max())]
        cs = dsp.getCs('jet',nkas)
        args = {'labs':['$N$',r'$log_{10}(err)$'],'xylims':['y',-8,0],'lw':2,
            'opt':'ps'}
        plts0 = [[Ns0,errs_0[i,:],[cs[i],'--o'],''] for (i,ka) in enumerate(kas0) ]
        plts0+= [[Ns0,errs_2[i,:],[cs[i],'-.o'],''] for (i,ka) in enumerate(kas0) ]
        plts0+= [[Ns0,errs_a[i,:],[cs[i],'-o' ],''] for (i,ka) in enumerate(kas0) ]
        legElt = {'forward':'k-o','secondary':'k-.o','uncoupled':'k--o'}
        legElt.update({'$ka=%.1f$'%ka:[cs[i],'o'] for i,ka in enumerate(kas0)})
        dsp.stddisp(plts0,title='error uncoupled, $k_p=%.4f$' %kp,
            name=name+'err_kp%d.svg' %ikp,legElt=legElt,**args)
        # dsp.stddisp(pltsa,title='error forward coupling, $k_p=%.4f$' %kp,
        # name=name+'err_forward.svg',**args)
    if 'f' in optsN:
        ka0,N0 = 15,20
        show_ff(df,cond='(kp==%.4f) &(ka==%d) &(N==%d)' %(kp,ka0,N0),npts=2000,
            name=name+'ff_kp%d.svg' %ikp,xylims=['y',-10,0],opt='ps')

#############################################################################
#### N=2
#############################################################################
if '2' in opts:
    optsA='' #d(ka,kd) n(ka,kp) e(err vs param)
    df['kdka'] = df.kd/df.ka
    df = df.loc[df.kdka<100]
    kas0=np.array(df.ka.unique(),dtype=float)
    kps0=np.array(df.kp.unique(),dtype=float)
    kds0=np.array(df.kdka.unique(),dtype=float)
    argsM = {'pOpt':'im','cmap':'RdYlGn_r'}

    if 'd' in optsA:
        ikp = 3
        kp = kps0[ikp]
        # kderrs,errs = np.zeros(kas0.shape),np.zeros(kas0.shape)
        dfkp = df.loc[df.kp==kp]
        nkas = kas0.size
        errs_a = np.zeros((nkas,kds0.size))
        errs_0 = np.zeros((nkas,kds0.size))
        for i,ka in enumerate(kas0):
            dfka = dfkp.loc[dfkp.ka==ka]
            errs_a[i,:] = np.log10(np.array(dfka.err_a.values,dtype=float))
            errs_0[i,:] = np.log10(np.array(dfka.err_0.values,dtype=float))


        xylims=['y',min(errs_0.min(),errs_a.min()),max(errs_0.max(),errs_a.max())]
        # cs = dsp.getCs('Spectral',nkas)
        # args = {'labs':['$kd/ka$',r'$log_{10}(err)$'],'xylims':xylims,'lw':2}
        # plts0 = [[kds0,errs_0[i,:],[cs[i],'-o'],'$ka=%.1f$'%ka] for (i,ka) in enumerate(kas0) ]
        # pltsa = [[kds0,errs_a[i,:],[cs[i],'-o'],'$ka=%.1f$'%ka] for (i,ka) in enumerate(kas0) ]
        # dsp.stddisp(plts0,title='error uncoupled, $k_p=%.4f$' %kp       ,**args)
        # dsp.stddisp(pltsa,title='error forward coupling, $k_p=%.4f$' %kp,**args)

        # kds,kas =np.meshgrid(kas0,kds0)
        #### maps ka,kdka
        lkds0 = np.arange(kps0.size)#;fkp = interp1d(kps0,lkps0)
        lkas0 = np.arange(kas0.size)#;fka = interp1d(kas0,lkas0)
        argsM.update({'labs':['$kd/ka$','$ka$'],'caxis':[-8,0],
            'xyTicks':[lkds0,lkas0],'xyTickLabs':[np.array(kds0,dtype=str),np.array(kas0,dtype=str)],
            'opt':'ps'})
        dsp.stddisp(im=[lkds0,lkas0,errs_a],title='error forward   kp=%.4f' %kp,
            name=name+'kakd_forward%d.png' %ikp,**argsM)#caxis=[errs_a.min(),errs_a.max()]
        dsp.stddisp(im=[lkds0,lkas0,errs_0],title='error uncoupled kp=%.4f' %kp,
            name=name+'kakd_uncoupled%d.png' %ikp,**argsM)#caxis=[errs_0.min(),errs_0.max()]

    if 'c' in optsA:
        kd0=3
        dfkd = df.loc[df.kdka==kd0]
        errs_pa = np.zeros((nkas,kps0.size))
        errs_p0 = np.zeros((nkas,kps0.size))
        for i,ka in enumerate(kas0):
            dfkap = dfkd.loc[dfkd.ka==ka]
            errs_pa[i,:] = np.log10(np.array(dfkap.err_a.values,dtype=float))
            errs_p0[i,:] = np.log10(np.array(dfkap.err_0.values,dtype=float))

        ###points for carbone
        from scipy.interpolate import interp1d
        E = 200              #keV
        lam = cst.keV2lam(E) #A
        r,V = np.loadtxt('data/C_V.txt').T; V /=1e3
        fV = interp1d(r,V)
        kasC = np.array(kas0[kas0>=6],dtype=float)
        r0 = kasC/(2*np.pi/lam)
        kpsC = np.sqrt(1+fV(r0)/E)

        lkps0 = np.arange(kps0.size);fkp = interp1d(kps0,lkps0)
        lkas0 = np.arange(kas0.size);fka = interp1d(kas0,lkas0)
        plts = [fkp(kpsC),fka(kasC),'bo-','kp(C)@200keV']
        argsM.update({'labs':['$k_p$','$ka$'],'plots':plts,'ms':10,'caxis':[-5,0],
            'xyTicks':[lkps0,lkas0],'xyTickLabs':[np.array(kps0,dtype=str),np.array(kas0,dtype=str)],
            'opt':'ps'})
        dsp.stddisp(im=[lkps0,lkas0,errs_pa],title='error forward   kd/ka=%.4f' %kd0,
            name=name+'kakp_forward.png',**argsM)
        dsp.stddisp(im=[lkps0,lkas0,errs_p0],title='error uncoupled kd/ka=%.4f' %kd0,
            name=name+'kakp_uncoupled.png',**argsM)


    # idx = abs(ff)/ffmax>1e-3
    # err = abs((ff[idx]-ffu[idx])/ff[idx]).mean()
    # print(err)

    # if err<0.05:qdot2.show_f()

    # plts = [[kas0,kderrs,'b-o','$err=%1.E$' %err]]
    # txts = [[ka,kdka,'%.2E' %err,'b'] for ka,kdka,err in zip(kas0,kderrs,errs)]
    # dsp.stddisp(plts,texts=txts,labs=['$ka$','$kdka$'],fonts={'text':15})

    # s2s,s0s = [],[]
    # for i,N in enumerate(Ns):
    #     s2 = qsa.QdotSphereArray(N=N,ka=ka,kp=kp,kd=kd*ka,nmax=nmax);s2.solve(v=1)
    #     s0 = copy.copy(s2); s0.solve(copt=0)
    #     fig,ax = s0.show_ff(npts=2000,fopts='m',opt='',leg=0)
    #     s2.show_ff(npts=2000,ax=ax,fopts='m',legElt={'uncoupled':'k--'},name=name+'N%d_' %i,opt='p')
    #     s2s+=[s2]
    #     s0s+=[s0]


# if 'd' in opts:
#     N,ka,kp,nmax = 5,2,1.2,5
#     kds = [2.5,4,10]
#     for i,kd in enumerate(kds):
#         s2 = qsa.QdotSphereArray(N=N,ka=ka,kp=kp,kd=kd*ka,nmax=nmax,solve=True)
#         s0 = copy.copy(s2); s0.solve(copt=0)
#         fig,ax = s0.show_ff(fopts='m',opt='',leg=0)
#         s2.show_ff(ax=ax,fopts='m',legElt={'uncoupled':'k--'},name=name+'kd%d_' %i,opt='sc')
#
