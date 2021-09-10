from utils import *                                 ;imp.reload(dsp)
import scipy,copy,os,scipy.special as spe
from scipy.signal import find_peaks
from scipy.integrate import trapz,quad
from scipy.interpolate import interp1d
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)
from pyScatSpheres import utils as ut ;imp.reload(ut)
plt.close('all')
name='figures/qdotSphereArray2approx_kakp_'
df_name = 'data/df_qdotSpheres_2.pkl'
opt = 'ps'
opts = 'd' #c(fig7.c bp) d(fig7.d ff)
#### load dataframe with precomputed solutions
df = pd.read_pickle(df_name)
ka,kp = 30,1.01
c,qdot2 = ut.load_qdot(df,{'ka':ka,'kp':kp,'kd':3*ka},'co')


if 'c' in opts:
    ### gather coeffs for all kas with values kp,kdka
    kps = df.kp.unique()
    df['kdka'] = df.kd/df.ka
    dfka = df.loc[df.eval('(kp==%.4f) & (kdka==%d)' %(kp,3))]
    bp,bpa,bp0 = dfka[['bp','bpa','bp0']].values.T
    ####select sphere 2 and adapt order (trail with zeros)
    nmaxs,kds =  np.array(dfka.nmax.unique(),dtype=int),np.array(dfka.kd.unique(),float)
    bp2m = lambda bp:np.array([np.hstack([b[nmax_i+1:],
        [0]*(nmaxs.max()-nmax_i)]) for b,nmax_i in zip(bp,nmaxs)]).T/np.exp(-1J*kds)
    s2,ska = np.array([0,1,2,3,8,13,18]),slice(0,-2,1)
    # bp1,bp1a,bp10 = bp1m(bp)[s1,ska],bp1m(bpa)[s1,ska],bp1m(bp0)[s1,ska]
    bp2,bp2a,bp20 = bp2m(bp)[s2,ska],bp2m(bpa)[s2,ska],bp2m(bp0)[s2,ska]

    # s1 = slice(0,c.nmax,1)
    # s2 = slice(c.nmax,None,1)
    # bp1,bp1a,bp10 = bp2m(bp)[s1,ska],bp2m(bpa)[s1,ska],bp2m(bp0)[s1,ska]
    # bp2,bp2a,bp20 = bp2m(bp)[s2,ska],bp2m(bpa)[s2,ska],bp2m(bp0)[s2,ska]

    nmax_s=s2.size
    plts,cs = [],dsp.getCs('jet',nmax_s )
    plts+= [[np.real(bp2a[i]),np.imag(bp2a[i]),[cs[i],'--s'],''] for i in range(nmax_s)]
    plts+= [[np.real(bp20[i]),np.imag(bp20[i]),[cs[i],'-.d'],''] for i in range(nmax_s)]
    plts0= [[np.real(bp2[i]) ,np.imag(bp2[i]) ,[cs[i],'-o' ],''] for i in range(nmax_s)]
    fig,ax = dsp.stddisp(plts0,pargs={'mfc':'none'},ms=13,lw=2,opt='')
    #separate legend
    legElt_l={'$l=%d$' %nmax:[c,'o'] for nmax,c in zip(s2,cs)}
    leg = plt.legend(handles= dsp.get_legElt(legElt_l),fontsize='25',loc=4)#,bbox_to_anchor=out)
    ax.add_artist(leg)
    ### final display
    legElt = {'exact':('k-o',{'mfc':'none','markersize':13}),'forward':'k--s','kin':'k-.d'}
    dsp.stddisp(plts,ax=ax,labs=[r"$\Re e(b_{2;l,0})$",r"$\Im m(b_{2;l,0})$"],#title=tle,
        fonts='2',legElt=legElt,legLoc='lower center',lw=2,ms=6,
        xyTicks=[np.arange(-1,2)]*2,xylims=[-1.7,2,-1.9,1.4],
        opt=opt,name=name+'bp.eps')

if 'd'in opts:
    bp,bp0,bpa = c.bp,c.bp0,c.bpa

    theta_d = np.linspace(0,90,3600)
    theta   = np.deg2rad(theta_d)
    qdot2.bp=c.bp
    ff =  qdot2.get_ff(theta)
    qdot2.bp=c.bp0
    ff0 = qdot2.get_ff(theta)
    qdot2.bp=c.bpa
    ffa = qdot2.get_ff(theta)
    ffmax = abs(ff).max()

    err = lambda ff0:abs(ff0)/ffmax
    # err = lambda ff0:np.log10(abs(abs(ff0)-abs(ff))/abs(ff))
    plts=[]
    plts =[[theta_d,err(ff) ,'b-' ,'exact']]
    plts+=[[theta_d,err(ffa),'g--','forward']]
    # plts+=[[theta_d,err(ff2),'c--','secondary']]
    plts+=[[theta_d,err(ff0),'r-.','kinematic']]
    tle = r'$N=%d, ka=%d, k_p=%.4f, kd=%d$ ' %(c.N,c.ka,c.kp,c.kd)
    dsp.stddisp(plts,labs=[r"$\theta(^{\circ})$",r"$|f(\theta)|$"],#title=tle,
        lw=2,xylims=[0,90,0,0.08],fonts='2',
        opt=opt,name=name+'ff.eps')
