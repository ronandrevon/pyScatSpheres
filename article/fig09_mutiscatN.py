from utils import *                                 ;imp.reload(dsp)
import scipy,copy,os,scipy.special as spe
from scipy.signal import find_peaks
from scipy.integrate import trapz,quad
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)
from pyScatSpheres import utils as ut ;imp.reload(ut)
plt.close('all')
name='figures/TmatrixApproxErr'
# df_name = 'data/df_qdotSpheresN2.pkl'
opts = 'be' #a(fig9a errors), f(ff) e(set_error)
opt = 'p'
cmap_ns = 'viridis'

args = {'labs':['$N$',r'$log_{10}(err(b_p))$'],'lw':2,
    'fonts':'2','opt':opt}

# kderrs,errs = np.zeros(kas0.shape),np.zeros(kas0.shape)
if 'a' in opts:
    ns = 5
    df_name='data/df_qdot_ka7_kp3_Ns.pkl'
    df = pd.read_pickle(df_name)
    if 'e' in opts:df = ut.set_errors(df_name)

    Ns0 =np.array(df.N.unique(),dtype=float)
    #### uncomment if kas,kps need be selected
    # kas0=np.array([7,15])#np.array(df.ka.unique(),dtype=float)
    # kps0=np.array(df.kp.unique(),dtype=float)
    # ikp = 3
    # kp = kps0[ikp]
    # dfkp = df.loc[df.kp==kp]
    # nkas = kas0.size
    # errs_a = np.zeros((nkas,Ns0.size))
    # errs_0 = np.zeros((nkas,Ns0.size))
    # errs_2 = np.zeros((nkas,Ns0.size))
    # for i,ka in enumerate(kas0):
    #     dfka = dfkp.loc[dfkp.ka==ka]
    #     errs_0[i,:] = np.log10(np.array(dfka.err_0.values,dtype=float))
    #     errs_2[i,:] = np.log10(np.array(dfka.err_2.values,dtype=float))
    #     errs_a[i,:] = np.log10(np.array(dfka.err_a.values,dtype=float))

    errsa = np.log10(np.array(df.erra.values,dtype=float))
    errs0 = np.log10(np.array(df.err0.values,dtype=float))
    errs2 = np.log10(np.array(df.err2.values,dtype=float))
    errs_n = np.log10([np.array(df['err_%d' %i].values,dtype=float) for i in range(1,ns+1)])
    # xylims=['y',min(errs_n.min(),errsa.min()),max(errs_n.max(),errsa.max())]
    # cs = ['g','b']

    plts0,cs = [],dsp.getCs(cmap_ns,ns)
    plts0+= [[Ns0,errsa,['k','s'],'forward']]
    plts0+= [[Ns0,errs0,[cs[0],'d'],'kinematic']]
    plts0+= [[Ns0,errs2,[cs[1],'^'],'secondary']]
    # plts0+= [[Ns0,errs_a[i,:],[cs[i],'-o' ],''] for (i,ka) in enumerate(kas0) ]
    # plts0 = [[Ns0,errs_0[i,:],[cs[i],'--o'],''] for (i,ka) in enumerate(kas0) ]
    # plts0+= [[Ns0,errs_2[i,:],[cs[i],'-.o'],''] for (i,ka) in enumerate(kas0) ]
    # legElt = {'forward':'k-o','secondary':'k-.o','uncoupled':'k--o'}
    # legElt.update({'$ka=%.1f$'%ka:[cs[i],'o'] for i,ka in enumerate(kas0)})
    fig,ax=dsp.stddisp(plts0,opt='',pargs={'mfc':'none'},ms=13,lw=2)
    legElt={'$n=%d$' %(i+1):[cs[i],'--o'] for i in range(ns)}
    leg = plt.legend(handles= dsp.get_legElt(legElt),fontsize='25',loc=4)#,bbox_to_anchor=out)
    ax.add_artist(leg)
    plts1= [[Ns0,errs_n[i,:],[cs[i],'--o'],''] for i in range(ns) ]
    dsp.stddisp(plts1,ax=ax,legLoc='lower left',xylims=[0,105,-6.5,0.5],#title='error uncoupled, $k_p=%.4f$' %kp,
        name=name+'.eps',**args)
    # dsp.stddisp(pltsa,title='error forward coupling, $k_p=%.4f$' %kp,
    # name=name+'err_forward.svg',**args)

if 'b' in opts:
    plts,cs=[],dsp.getCs('jet',4)
    df_name='data/df_qdot_ka7_N_100_np.pkl'
    df = pd.read_pickle(df_name)
    kps = np.array(df.kp.unique(),dtype=float)
    if 'e' in opts:df = ut.set_errors(df_name)
    ns = 10
    errs_n = np.log10([np.array(df['err_%d' %n].values,dtype=float) for n in range(1,ns+1)]).T
    errs_n = errs_n[:ns]
    plts += [[np.arange(ns)+1,errs_n[i],[cs[i-1],'-o'],'$n_p=%.4f$' %kp] for i,kp in enumerate(kps)]
    erra = np.log10(df.erra.values,dtype=float)
    # legElt = {'approx':'k-s'}
    # leg = plt.legend(handles= dsp.get_legElt(legElt),fontsize='25',loc=4)#,bbox_to_anchor=out)
    # ax.add_artist(leg)
    plts0 = [[ [ns]*kps.size,erra,'ks','forward']]
    fig,ax=dsp.stddisp(plts0,opt='',pargs={'mfc':'none'},ms=10,lw=2)
    args['labs'][0]='$n$'
    dsp.stddisp(plts,xylims=[0,ns+0.5,-6,0],ax=ax,legLoc='lower right',
        name=name+'_n.eps',**args)

if 'f' in opts:
    # ka0,N0 = 15,20
    qdot = ut.load_qdot(df,{'N':2},'o')
    qdot.show_f('t',npts=100,r=(0,10,-10,30),pOpt='e',fonts='2',tle='')
    # show_ff(df,cond='(kp==%.4f) &(ka==%d) &(N==%d)' %(kp,ka0,N0),npts=2000,
    #     name=name+'ff_kp%d.svg' %ikp,xylims=['y',-10,0],opt='ps')
