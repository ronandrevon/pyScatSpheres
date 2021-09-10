from utils import * ;imp.reload(dsp)
from pyScatSpheres import utils as ut ;imp.reload(ut)
plt.close('all')
name='figures/Tmatrix_'
opt='ps'
opts='a' #b (fig10b scattering amplitudes) a(fig10a Probas)


theta = np.deg2rad(np.linspace(0,90,1000))#3601))
if 'a' in opts:
    df_name = 'data/df_qdot_ka7_kp2_Ns.pkl'
    df = pd.read_pickle(df_name)
    Nmax,n_max = df.shape[0],9#max scattering orders

    dfc  = df[:Nmax].copy()
    Ns   = np.array(dfc.N.values,dtype=float)[:Nmax]
    sig_cols = ['sigma%d' %i for i in range(1,n_max+1)]
    dfc[['sigma']+sig_cols]=0
    dt = theta[2]-theta[1]
    ka = dfc.ka[0]
    sig = lambda sigf:np.sum(sigf*np.sin(theta))*2*np.pi*dt*4/ka**2
    for l,c in dfc.iterrows():
        qdot = ut.load_qdot(dfc ,{'ka':c.ka,'kp':c.kp,'N':c.N},'o')

        ff  = qdot.get_ffn(theta,c.bpa)
        sigT = sig(abs(ff)**2)
        dfc.loc[l,'sigma']  = sigT
        ffm=np.zeros(theta.shape,dtype=complex)#(n_max,theta.size))
        for i in range(n_max,0,-1):
        # for i in range(1,n_max):
            ffi = qdot.get_ffn(theta,c['bp_%d' %i])
            sig_n = sig(abs(ffi)**2+2*np.real(ffi*np.conj(ffm)))
            dfc.loc[l,'sigma%d' %i] = sig_n
            ffm+=ffi
            # ff_n[i,:]=ffi

    S = dfc.sigma.max()*10
    sigma_a = (dfc.sigma/Ns).mean()
    le = S*c.kd/sigma_a
    ns = np.arange(Ns.max())*1
    Ppoisson = 1-np.exp(-ns*c.kd/le)
    dfc['sigma_dyn']=dfc[sig_cols[1:]].sum(axis=1)
    dfc['scat']=dfc[['sigma1','sigma_dyn']].sum(axis=1)

    plts,cs=[],dsp.getCs('viridis',n_max+1)
    # plts+= [[Ns,1-dfc.sigma/S,'r-','coh']]
    plts+= [[ns,Ppoisson,'k-.','$P_{scat}^{P}$']]
    # plts+= [[Ns,dfc.scat/S  ,'k--' ,'$P_{kin}+P_{dyn}$']]
    plts+= [[Ns,dfc.sigma/S ,'k-'  ,'$P_{scat}$']]
    plts+= [[Ns,dfc.sigma1/S,'c-o','$P_{kin}$']]
    plts+= [[Ns,dfc.sigma_dyn/S,'b-o','$P_{dyn}$']]
    plts+= [[Ns,dfc['sigma%d' %i]/S,[cs[i],'--o'],''] for i in range(1,n_max+1)]
    legElt={'$P_{n}$':'k--o'}
    dsp.stddisp(plts,labs=['$N$',r'$Fraction$'],lw=2,fonts='2',legElt=legElt,
        name=name+'Pn.eps',opt=opt)

if 'b' in opts:
    df_name = 'data/df_qdot_ka7_kp2_Ns.pkl'
    nmax = 3
    df = pd.read_pickle(df_name)
    c,qdot = ut.load_qdot(df ,{'ka':df.ka[0],'kp':df.kp[0],'N':15},'co')

    ff  = abs(qdot.get_ffn(theta,c.bpa))
    ffs=[abs(qdot.get_ffn(theta,c['bp_%d' %(i+1)])) for i in range(nmax)]

    theta_d = np.rad2deg(theta)
    plts,cs=[],dsp.getCs('viridis',nmax)
    plts+= [[theta_d, ff ,'k-',r'$\sum f_n$',1]]
    plts+= [[theta_d, ffs[i],[cs[i],'--'],'$f_%d$' %(i+1),2] for i in range(nmax)]
    dsp.stddisp(plts,labs=[r'$\theta(deg)$',r'$|f_n|$'],fonts='2',
        xylims=[0,91,0,ffs[1].max()],
        name=name+'fn.eps',opt=opt)





# if 'F' in opts:
#     kas,kps = [5,10,15],[1.1,1.01,1.001]
#     nkas = df.ka.unique().size
#     cs = dsp.getCs('jet',nkas)
#     npts = 3600
#     for i,kp in enumerate(kps):
#         plts=[]
#         for j,ka in enumerate(kas):
#             dfc = df.loc[(df.ka==ka) &(df.kp==kp)].copy()
#             Ns = np.array(dfc.N.values,dtype=float)
#             dfc['sigma']=0
#             for l,c in dfc.iterrows():
#                 qdot = ut.load_qdot(dfc,{'ka':c.ka,'kp':c.kp,'N':c.N},'o')
#                 dfc.loc[l,'sigma'] = qdot.get_s(npts=npts,norm=True)
#             plts+=[[Ns,dfc.sigma,[cs[j],'-o'],r'$ka=%d$' %ka]]
#
#             qdot = qsa.QdotSphereArray(N=1,ka=c.ka,kp=c.kp,kd=c.kd,nmax=c.nmax,solve=1,copt=0)
#             sig0 = qdot.get_s(npts=npts,norm=True)
#             plts += [[Ns,Ns**2*sig0,[cs[j],'--o'],''],[Ns,Ns*sig0,[cs[j],'-.o'],'']]
#             legElt = {r'$N^2\sigma_0$':'k--o',r'$N\sigma_0$':'k-.o'}
#         tle = r'$k_p=%.4f$' %kp
#         dsp.stddisp(plts,labs=['$N$',r'$\sigma_{tot}$'],lw=2,title=tle,legElt=legElt,
#             name=name+'_sigmaN%d.eps' %i,opt=opt)
