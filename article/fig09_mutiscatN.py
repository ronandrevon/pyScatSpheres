from utils import *                                 ;imp.reload(dsp)
import scipy,copy,os,scipy.special as spe
from scipy.signal import find_peaks
from scipy.integrate import trapz,quad
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)
from pyScatSpheres import utils as ut ;imp.reload(ut)
plt.close('all')
name='figures/qdotSphereArray2approx_kakp_'
df_name = 'data/df_qdotSpheresN2.pkl'
opts = '' # s(solve) e(set_error)
opt = 'p'
optsN = 'a' #a(fig9a errors), f(ff)


df = pd.read_pickle(df_name)
if 'e' in opts:df = set_errors(df_name)

kas0=np.array([7,15])#np.array(df.ka.unique(),dtype=float)
kps0=np.array(df.kp.unique(),dtype=float)
Ns0 =np.array(df.N.unique(),dtype=float)
argsM = {'pOpt':'im','cmap':'RdYlGn_r'}
ikp = 3
kp = kps0[ikp]

# kderrs,errs = np.zeros(kas0.shape),np.zeros(kas0.shape)
if 'a' in optsN:
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
    # cs = dsp.getCs('jet',nkas)
    cs = ['g','b']
    args = {'labs':['$N$',r'$log_{10}(b_{p;l,m})$'],'xylims':['y',-10,0],'lw':2,
        'fonts':'2','opt':opt}
    plts0 = [[Ns0,errs_0[i,:],[cs[i],'--o'],''] for (i,ka) in enumerate(kas0) ]
    plts0+= [[Ns0,errs_2[i,:],[cs[i],'-.o'],''] for (i,ka) in enumerate(kas0) ]
    plts0+= [[Ns0,errs_a[i,:],[cs[i],'-o' ],''] for (i,ka) in enumerate(kas0) ]
    legElt = {'forward':'k-o','secondary':'k-.o','uncoupled':'k--o'}
    legElt.update({'$ka=%.1f$'%ka:[cs[i],'o'] for i,ka in enumerate(kas0)})
    dsp.stddisp(plts0,#title='error uncoupled, $k_p=%.4f$' %kp,
        name=name+'err_kp%d.eps' %ikp,legElt=legElt,**args)
    # dsp.stddisp(pltsa,title='error forward coupling, $k_p=%.4f$' %kp,
    # name=name+'err_forward.svg',**args)


if 'f' in optsN:
    ka0,N0 = 15,20
    show_ff(df,cond='(kp==%.4f) &(ka==%d) &(N==%d)' %(kp,ka0,N0),npts=2000,
        name=name+'ff_kp%d.svg' %ikp,xylims=['y',-10,0],opt='ps')
