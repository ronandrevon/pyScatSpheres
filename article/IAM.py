from utils import *                                ;imp.reload(dsp)
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)
from scipy.integrate import trapz,quad
from scipy.interpolate import interp1d
import scipy.special as spe
plt.close('all')
name='figures/qdotSphereSingle_shells_'
opt='p'
opts = 'V'   #V(potential),Q(Fq)


E = 200              #keV
lam = cst.keV2lam(E) #A
k0 = 2*np.pi/lam

if 'V' in opts:
    r,V = np.loadtxt('../tests/data/C_V.txt').T
    fV = interp1d(r,V/1e3)
    n_p = lambda V:(np.sqrt(1+V/E)-1)*1e3
    #### potential Carbone and fits
    r0 = np.arange(0.02,1,0.02)
    V0 = fV(r0)
    Vc = 0.06633476/r0-0.08779877
    Vy = 0.07683447*np.exp(-2.91720804*r0)/r0
    # plts=[]
    plts = [[r0*k0,n_p(V0),'b-o',r'$V(EJK)$']]
    plts+= [[r0*k0,n_p(Vy),'b--',r'$e^{br}/r$']]
    plts+= [[r0*k0,n_p(Vc),'r--',r'$a/r+b$']]

    #### shells
    # kas = np.arange(5,101,5)
    # ri = np.linspace(0.02,2,40)
    kas = r0*k0
    nps = n_p(V0)
    # csf = dsp.getCs('jet',kas.size)
    csf = [(0.5,0.5,1)]*kas.size
    # kass = np.hstack([0,kas])
    kass = np.hstack([0,(kas[1:]+kas[:-1])/2,kas[-1]])
    # kass = np.hstack([0,kas])
    patches = [dsp.Rectangle((kass[i],0),kass[i+1]-kass[i],npi,
        fill=1,color=csf[i],ec='b') for i,npi in enumerate(nps)]

    xmax,ymax = 101,5
    fig,ax=dsp.stddisp(plts,lw=2,patches=patches)
    kat,kpt = np.arange(0,101,20),np.arange(0,ymax+1)
    kp2np  = lambda n_p:E*((n_p*1e-3+1)**2-1)
    rat,Vpt = kat[1:]/k0, kp2np(kpt[1:])
    rat,Vpt = np.round(rat*1e2)/1e2,np.round(Vpt*1e2)/1e2
    rat_s,Vpt_s = np.array(rat,dtype=str),np.array(Vpt,dtype=str)
    fonts = dsp.get_font('2')[:4]
    dsp.addxAxis(ax=ax,plots=[],xLab=r'$r(\AA)$',c=(0.5,)*3,xTicks=rat,xTickLabs=rat_s,gridOn=0,fonts=fonts,xylims=['x',0,xmax/k0],)
    dsp.addyAxis(ax=ax,plots=[],yLab=r'$V(kV)$' ,c=(0.5,)*3,yTicks=Vpt,yTickLabs=Vpt_s,gridOn=0,fonts=fonts,xylims=['y',0,kp2np(ymax)],fig=fig)
    dsp.stddisp(ax=ax,xylims=[0,xmax,0,ymax],xyTicks=[kat,kpt],
        xyTickLabs=[np.array(kat,dtype=str),np.array(kpt,dtype=str)],
        fonts='2',labs=['$ka$',r'$n_p-1 (\times 1e^{-3})$'],axPos=[0.12,0.09,0.76,0.8],
        name=name+'pot1.eps',opt=opt)

if 'Q' in opts:
    theta_deg = np.linspace(1e-3,15,361)
    theta = np.deg2rad(theta_deg)
    # s1 = [qsa.QdotSphereArray(N=1,ka=ka,kp=kp,kd=0,nmax=80,solve=1,copt=0) for ka in kas]
    # f1 = np.array([s.get_ff(theta) for s in s1])
    # cs = dsp.getCs(cm,kas.size)
    # ct,st = np.cos(theta),np.sin(theta)
    ka_max = kas.max()
    q = 2*np.sin(theta/2)

    eps = np.hstack([eps,0])
    fi = np.array([ (eps[i]-eps[i+1])*ka**3*(np.sin(ka*q)/(ka*q)-np.cos(ka*q))/(ka*q)**2 for i,ka in enumerate(kas)])
    # fm = fi.sum(axis=0)
    # fmax = abs(fm).max()
    # f=abs(fm)*2.5/fmax
    # plts = [[K*q,f,'k','total($ka_max=%d$)' ]]

    xylims = [-0.1,10,-0.1,2.7]
    nmaxs = np.hstack([np.arange(5,npts-1,10),[npts-1]])
    fim = np.array([abs(np.array(fi[:nmax+1,:]).sum(axis=0)) for nmax in nmaxs])
    cs = dsp.getCs('Spectral',npts)
    plts= [[K*q,fim[i]*2.5/abs(fim[i]).max(),[cs[nmax],'-'],'$ka_{max}=%d$' %kas[nmax] ] for i,nmax in enumerate(nmaxs)]
    # plts += [[k0*q,fi[i],[csf[i],'--'],''] for i,ka in enumerate(kas) ]
    labs = [r'$q(\AA^{-1})$',r'$|f(q)|$']
    # plts += [[theta_deg,fi[i],[csf[i],'--'],'$ka=%d$' %ka ] for i,ka in enumerate(kas) ]
    dsp.stddisp(plts,lw=2,labs=labs,#inset=inset,iOpt='GtX',
        xylims=xylims,fonts='2',#legElt=legElt,title = r'$\epsilon=%.2f$, $n_{ref}=%.3f$' %(eps,kp),
        name=name+'fka1.eps',opt=opt)
