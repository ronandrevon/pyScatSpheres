from utils import *                                ;imp.reload(dsp)
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)
from scipy.integrate import trapz,quad
from scipy.interpolate import interp1d
import scipy.special as spe
plt.close('all')
name='figures/qdotSphereSingle_shells_'
opt='ps'
opts = 'Q'   #V(potential),Q(Fq)

E = 200              #keV
lam = cst.keV2lam(E) #A
k0 = 2*np.pi/lam

r,V = np.loadtxt('../tests/data/C_V.txt').T
r,V = np.hstack([0,r]),np.hstack([1e4,V])
fV = interp1d(r,V/1e3,kind='cubic')
n_p = lambda V:(np.sqrt(1+V/E)-1)*1e3
ka2eps=lambda ka:fV(ka/k0)/E

if 'V' in opts:

    #### potential Carbone and fits
    r0 = np.arange(0.02,1,0.02)
    V0 = fV(r0)
    Vc = 0.06633476/r0-0.08779877
    Vy = 0.07683447*np.exp(-2.91720804*r0)/r0
    rs=np.linspace(0.02,1,1000)
    plts=[]
    # plts += [[r0*k0,n_p(V0),'bo',r'']]
    plts += [[rs*k0,n_p(fV(rs)),'b-',r'$V(IAM)$']]
    # plts+= [[r0*k0,n_p(Vy),'b--',r'$e^{br}/r$']]
    # plts+= [[r0*k0,n_p(Vc),'r--',r'$a/r+b$']]

    #### shells
    kas = np.arange(5,106,10)#r0*k0
    nps = n_p(fV(kas/k0))
    csf = [(0.5,0.5,1)]*kas.size
    kass = np.hstack([0,(kas[1:]+kas[:-1])/2,kas[-1]])
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
    #### Born multi shell approx
    kas = np.arange(5,401,5)
    theta_deg = np.linspace(1e-3,15,361)
    theta = np.deg2rad(theta_deg)
    ka_max = kas.max()
    q = 2*np.sin(theta/2)

    eps = ka2eps(kas)
    eps = np.hstack([eps,0])
    fi = np.array([
        (eps[i]-eps[i+1])*ka**3/(ka*q)**2*(np.sin(ka*q)/(ka*q)-np.cos(ka*q))
            for i,ka in enumerate(kas)])

    ka_maxs = [50,100,200,400]
    npts = len(ka_maxs)
    nmaxs = np.array([ abs(kas-ka).argmin() for ka in ka_maxs])
    fim = np.array([abs(np.array(fi[:nmax+1,:]).sum(axis=0)) for nmax in nmaxs])
    cs = dsp.getCs('Spectral',npts)
    plts= [[k0/(2*np.pi)*q,fim[i]*2.5/abs(fim[i]).max(),[cs[i],'-'],
        '$ka_{max}=%d$' %kas[nmax]] for i,nmax in enumerate(nmaxs)]

    xylims = [-0.1,10,-0.1,2.7]
    labs = [r'$q(\AA^{-1})$',r'$|f(q)|$']
    dsp.stddisp(plts,lw=2,labs=labs,#inset=inset,iOpt='GtX',
        xylims=xylims,fonts='2',#legElt=legElt,title = r'$\epsilon=%.2f$, $n_{ref}=%.3f$' %(eps,kp),
        name=name+'fka1.eps',opt=opt)
