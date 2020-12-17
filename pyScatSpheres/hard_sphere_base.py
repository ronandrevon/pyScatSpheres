# import importlib as imp
import scipy.special as spe
import numpy as np
from . import displayStandards as dsp #;imp.reload(dsp)
from . import spherical_utils as spu #;imp.reload(spu)

class HardSphereArrayBase():
    def __init__(self,N=1,ka=1,kd=2,kp=np.inf,k=1,nmax=7,Cp=None,solve=False,copt=1):
        '''equally disitributed linear array of hard spheres
        - N,ka,kp,kd : nb spheres,wave number inside spheres, normalized radius, normalized distances
        - k : incident wave number
        - nmax : max included order expansion
        '''
        self.k,self.N,self.nmax = k,N,nmax+1
        self.ka,self.kp,self.kd,self.kd_p = ka,kp,kd,kd*np.arange(N)
        self.d_p  = self.kd_p/k
        if isinstance(Cp,np.ndarray):self.set_Cp(Cp)
        if solve:self.solve(copt=copt)

    def get_s(self,npts=361):
        theta = np.linspace(0,np.pi,npts)
        dt = theta[1]-theta[0]
        f = self.get_ff(theta)
        s = np.sum(np.abs(f)**2*np.sin(theta))*2*np.pi*dt
        return s

    ###############################################################
    #  display
    ###############################################################
    def show_f(self,cart_opt=True,npts=200,r=None,opts='tT',idp=None,fz=np.abs,name='',def_args=True,**kwargs):
        '''
        opts : i(incident), s(scattered), t(total), P(sphere idp),T(All), G(Gradient)
        '''
        ka,kd,N = self.ka,self.kd,self.N
        if not isinstance(r,tuple) or not isinstance(r,list) or not isinstance(r,np.ndarray):
            r=(1e-3,4*ka+N*kd,-2*ka,N*kd+2*ka)
        r,theta,y,z = spu.polar_mesh2D(cart_opt,npts,r)
        k,N,nmax = self.k,self.N,self.nmax
        ap,dp = self.ka/k,self.kd/k

        if idp:self._check_idp(idp=None)
        args = {}
        if def_args:
            args = {'lw':2,'labs':['$y$','$z$'],'imOpt':'c','axPos':'V','fonts':{'title':25},}


        t = np.linspace(-np.pi/2,np.pi/2,100)
        ct,st = np.cos(t), np.sin(t)
        plts = [ [ap*ct, dp*p+ap*st,'k-',''] for p in range(N)]
        fstr =  r'$%s \psi(r,\theta)$' %(['',r'\partial_r']['G' in opts])
        name+=['f','df']['G' in opts]
        if 'P' in opts:
            for p in idp:
                print("...compute amplitude at p=%d..." %p)
                fi,fs = self.compute_f(r,theta,0,ftype='a',Gopt='G' in opts,idp=p)
                if 'i' in opts:
                    f = fz(fi)
                    dsp.stddisp(plts,im=[y,z,f],
                        title = r"Incident %s at sphere %d, $N_{max}=%d$" %(fstr,p,nmax),
                        name=name+'i_sphere%d.png' %p,**args,**kwargs)
                if 's' in opts:
                    f = fz(fs)
                    dsp.stddisp(plts,im=[y,z,f],
                        title = r"Scattered %s at sphere %d, $N_{max}=%d$" %(fstr,p,nmax),
                        name=name+'s_sphere%d.png' %p,**args,**kwargs)
                if 't' in opts:
                    f = fz(fi+fs)
                    dsp.stddisp(plts,im=[y,z,f],
                        title = r"Total %s at sphere %d, $N_{max}=%d$" %(fstr,p,nmax),
                        name=name+'t_sphere%d.png' %p,**args,**kwargs)
        if 'T' in opts:
            print("...compute total amplitude ...")
            fi,fs = self.compute_f(r,theta,0,ftype='a',Gopt='G' in opts)
            params = self._params()
            if 'i' in opts:
                f = fz(fi)
                dsp.stddisp(plts,im=[y,z,f],
                    title = r"Incident %s, %s " %(fstr,params),
                    name=name+'i.png',**args,**kwargs)
            if 's' in opts:
                f = fz(fs)
                dsp.stddisp(plts,im=[y,z,f],
                    title = r"Scattered %s, %s" %(fstr,params),
                    name=name+'s.png',**args,**kwargs)
            if 't' in opts:
                f = fz(fi+fs)
                dsp.stddisp(plts,im=[y,z,f],
                    title = r"Total %s, %s" %(fstr,params),
                    name=name+'t.png',**args,**kwargs)
        return f.max(),f.min()

    def show_ff(self,npts=361,fopts='m',lw=2,leg=1,name='',title=None,**kwargs):
        ''' display far field scattering amplitude
        - fopts : r(real), i(imag), m(mag), a(angle), 2(mag2)
        '''
        theta  = np.linspace(0,np.pi,npts)
        theta_d = np.rad2deg(theta)
        ff = self.get_ff(theta)


        plts,ls = [], ['--','-'][leg]
        if 'r' in fopts: plts+= [[theta_d,np.real(ff)    ,['c'          ,ls],['',r'$\Re$' ][leg] ]]
        if 'i' in fopts: plts+= [[theta_d,np.imag(ff)    ,['r'          ,ls],['',r'$\Im$' ][leg] ]]
        if 'a' in fopts: plts+= [[theta_d,np.angle(ff)   ,['m'          ,ls],['',r'$\phi$'][leg] ]]
        if 'm' in fopts: plts+= [[theta_d,np.abs(ff)     ,['b'          ,ls],['',r'$||$'  ][leg] ]]
        if '2' in fopts: plts+= [[theta_d,np.abs(ff)**2  ,[[0.25,0.75,1],ls],['',r'$||^2$'][leg] ]]
        if not isinstance(title,str):
            title='Scattering amplitude for %s' %self._params()
        return dsp.stddisp(plts,labs=[r"$\theta(^{\circ})$",r"$f(\theta)$"],lw=lw,
            title=title, xylims=[0,180,0,np.abs(ff).max()],
            name=name+'fka.svg',**kwargs)

    def _params(self):
        ss = '$N=%d$, $n_{ref}=%g$, $ka=%.1f$, $kd=%.1f$' %(self.N,self.kp, self.ka,self.kd)
        return ss
    ############################################################################
    # misc
    ############################################################################
    def _check_idp(self,idp=None):
        if isinstance(idp,int) :
            if idp>=0 and idp<self.N:
                idp=[idp]
            else:
                print(colors.red+'warn:sphere index not within 0<=idp<N'+colors.black)
                idp=np.arange(self.N)
        else:
            idp=np.arange(self.N)
        return idp

    def test_convergence(self,nmaxs,npts=361,name='', **kwargs):
        ns = np.array(nmaxs).size
        theta = np.linspace(0,np.pi,npts)
        f = np.zeros((ns,npts),dtype=complex)
        for i,nmax in enumerate(nmaxs):
            print('solving nmax=%d' %nmax)
            self.solve(nmax=nmax,copt=1,v=0)
            f[i] = self.get_ff(theta)
        theta_d = np.rad2deg(theta)

        #### display
        cs   = dsp.getCs('Spectral',ns)
        nref = ['%g' %self.kp,r'\infty'][self.kp==np.inf]
        tle  = 'Convergence test $ka=%.1f$, $kd=%.1f$, $n_{ref}=%s$' %(self.ka,self.kd,nref)
        labs =[r"$\theta(^{\circ})$",r"$|f(\theta)|$"]
        plts = [[theta_d,np.abs(f[i]), cs[i] ,'$n_{max}=%d$' %nmax] for i,nmax in enumerate(nmaxs)]
        dsp.stddisp(plts,labs=labs,lw=2,
            title=tle,xylims=['x',0,180],
            name=name+'fconv.svg',**kwargs)
