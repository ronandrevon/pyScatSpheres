from utils import * ;imp.reload(dsp)
from pyScatSpheres import utils as ut ;imp.reload(ut)
plt.close('all')
name='figures/Tmatrix_'
opt='ps'
opts='ab' #b (fig10b scattering amplitudes) a(fig10a Probas)

df_name = 'data/df_qdotSpheresN2.pkl'
df = pd.read_pickle(df_name)

ka,kp=7,1.01
dfc  = df.loc[(df.ka==ka) &(df.kp==kp)].copy()
n_max = 2#max scattering orders
theta = np.deg2rad(np.linspace(0,90,1000))#3601))
plts,cs=[],dsp.getCs('viridis',n_max+1)

if 'a' in opts:
    Ns   = np.array(dfc.N.values,dtype=float)
    dfc[['sigma']+['sigma%d' %i for i in range(1,n_max+1)]]=0
    dt = theta[2]-theta[1]
    sig = lambda sigf:np.sum(np.abs(sigf)*np.sin(theta))*2*np.pi*dt*4/ka**2
    for l,c in dfc.iterrows():
        qdot = ut.load_qdot(dfc ,{'ka':c.ka,'kp':c.kp,'N':c.N},'o')

        ff  = qdot.get_ffn(theta,c.bpa)
        ff1 = qdot.get_ffn(theta,c.bp0)
        ff2 = qdot.get_ffn(theta,c.bp2-c.bp0)

        sigT = sig(abs(ff)**2)
        sig1 = sig(abs(ff1)**2+2*np.real(ff1*np.conj(ff2)))
        sig2 = sig(abs(ff2)**2)
        sigN = sigT-sig1

        dfc.loc[l,'sigma']  = sigT
        dfc.loc[l,'sigma1'] = sig1
        dfc.loc[l,'sigma2'] = sig2
        dfc.loc[l,'sigmaN'] = sigN
        # dfc.loc[l,'sigmaN'] =
        # qdot.bp=c.bp2
        # dfc.loc[l,'sigma1'] = qdot.get_s(npts=npts,norm=True)

    S = dfc.sigma.max()*10
    # plts+= [[Ns,1-dfc.sigma/S,'r-','coh']]
    le = (dfc.sigma/Ns).mean()*S*c.kd
    ns = np.arange(1000)#Ns.max())*1
    Pscat = 1-np.exp(-ns*c.kd/le)
    # plts+= [[ns,Pscat,'k--','$P_{scat}^{P}$']]
    plts+= [[Ns,dfc.sigma/S ,'k-' ,'$P_{scat}$']]
    plts+= [[Ns,dfc.sigmaN/S,'b-o','$P_{dyn}$']]
    # plts+= [[Ns,dfc.sigma1/S,'g-o','$P_{kin}$']]
    plts+= [[Ns, dfc['sigma%d' %i]/S,[cs[i],'--o'],'$P_%d$' %i] for i in range(1,n_max+1)]
    # plts+= [[Ns,sigma2,'b-','dyn']]
    dsp.stddisp(plts,labs=['$N$',r'$Fraction$'],lw=2,fonts='2',
        name=name+'Pn.eps',opt=opt)

if 'b' in opts:

    c,qdot = ut.load_qdot(dfc ,{'ka':ka,'kp':kp,'N':15},'co')

    ff  = abs(qdot.get_ffn(theta,c.bpa))
    ff1 = abs(qdot.get_ffn(theta,c.bp0))
    ff2 = abs(qdot.get_ffn(theta,c.bp2-c.bp0))
    ffs=[ff,ff1,ff2]

    theta_d = np.rad2deg(theta)
    plts = [[theta_d, ff ,'k-',r'$\sum f_n$',1]]
    plts+= [[theta_d, ffs[i],[cs[i],'--'],'$f_%d$' %i,2] for i in range(1,n_max+1)]
    dsp.stddisp(plts,labs=[r'$\theta(deg)$',r'$|f_n|$'],fonts='2',
        xylims=[0,91,0,ff2.max()],
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
