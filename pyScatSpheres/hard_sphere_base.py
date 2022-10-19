# import importlib as imp
import scipy.special as spe
import numpy as np
from . import displayStandards as dsp #;imp.reload(dsp)
from . import spherical_utils as spu #;imp.reload(spu)

class HardSphereArrayBase():
    def __init__(self,N=1,ka=1,kd=2,kp=np.inf,k=1,nmax=7,Cp=None,solve=True,copt=1,v=1):
        '''equally disitributed linear array of hard spheres
        - N,ka,kp,kd : nb spheres,wave number inside spheres, normalized radius, normalized distances
        - k : incident wave number
        - nmax : max included order expansion
        '''
        self.k,self.N,self.nmax = k,N,nmax+1
        self.ka,self.kp,self.kd,self.kd_p = ka,kp,kd,kd*np.arange(N)
        self.d_p  = self.kd_p/k
        if isinstance(Cp,np.ndarray):self.set_Cp(Cp)
        if solve:self.solve(copt=copt,v=v)

    def get_s(self,npts=3601,norm=False):
        theta = np.linspace(0,np.pi,npts)
        dt = theta[1]-theta[0]
        f = self.get_ff(theta)
        s = np.sum(np.abs(f)**2*np.sin(theta))*2*np.pi*dt
        if norm:
            s *= 4/self.ka**2
        return s

    ###############################################################
    #  display
    ###############################################################
    def show_f(self,cart_opt=True,npts=200,r=None,opts='tT',idp=None,fz=np.abs,
        name='',def_args=True,short_title=False,v=0,**kwargs):
        '''
        opts : i(incident), s(scattered), t(total), P(sphere idp),T(All), G(Gradient), F(flux)
        '''
        ka,kd,N = self.ka,self.kd,self.N
        if not isinstance(r,tuple) and not isinstance(r,list) and not isinstance(r,np.ndarray):
            r=(1e-3,4*ka+N*kd,-2*ka,N*kd+2*ka)
        # print(r)
        r,theta,y,z = spu.polar_mesh2D(cart_opt,npts,r)
        k,N,nmax = self.k,self.N,self.nmax
        ap,dp = self.ka/k,self.kd/k

        self._check_idp(idp);
        args = {}
        if def_args:
            args = {'lw':2,'labs':['$y$','$z$'],'imOpt':'c','axPos':'V','fonts':{'title':25},}

        if not ('T'  in opts or 'P' in opts) : opts+='T'

        t = np.linspace(-np.pi/2,np.pi/2,100)
        ct,st = np.cos(t), np.sin(t)
        plts = [ [ap*ct, dp*p+ap*st,'k-',''] for p in range(N)]
        fstr =  r'$%s \psi(r,\theta)$' %(['',r'\partial_r']['G' in opts])
        Gopt =  ''.join([c for c in opts if c in 'GF'] )
        name+=['f','df']['G' in opts]
        if 'P' in opts:
            p=idp
            if v:print("...compute field at sphere p=%d..." %p)
            fi,fs = self.compute_f(r,theta,0,ftype='a',Gopt=Gopt,idp=p)
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
            if v:print("...compute field from all spheres ...")
            fi,fs = self.compute_f(r,theta,0,ftype='a',Gopt=Gopt)
            params = self._params()
            if short_title:params,fstr = '',[r'$\Psi$',r'$\partial_r\Psi$']['G' in opts]
            # if ''
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

    def show_ff(self,npts=361,fopts='m',lw=2,leg=1,name='',title=None,xylims=[],**kwargs):
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
        # if not xylims:
        #     m,M = 0,0
        #     if 'r' in fopts:m,M = min(m,np.abs(ff).min()
        #     if 'i' in fopts:
        #     if 'a' in fopts:
        #     if 'm' in fopts:
        #     if '2' in fopts:
        #     xylims=[0,180,m,M]
        return dsp.stddisp(plts,labs=[r"$\theta(^{\circ})$",r"$f(\theta)$"],lw=lw,
            title=title, xylims=xylims,
            name=name+'fka.svg',**kwargs)

    def _params(self):
        ss = '$N=%d$, $n_{ref}=%g$, $ka=%.1f$, $kd=%.1f$' %(self.N,self.kp, self.ka,self.kd)
        return ss
    ############################################################################
    # misc
    ############################################################################
    def _check_idp(self,idp=None):
        if isinstance(idp,int) :
            if not (idp>=0 and idp<self.N):
                raise Exception(colors.red+'sphere index not within 0<=idp<N'+colors.black)

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



    def test_convergence_error(self,nu_max,npts=360,**kwargs):
        nmaxs=np.arange(2,nu_max+1,2)
        ns = np.array(nmaxs).size
        # theta = np.linspace(0,np.pi,npts)

        eps=1e-16
        phi=np.pi/2

        #### sphere 1
        r1_p=np.array(self.ka+eps)
        r1_m=np.array(self.ka-eps)
        theta0=np.linspace(0,np.pi ,npts)
        theta1 = np.hstack([theta0,np.flipud(theta0)])
        # phi=np.array([np.pi/2]*npts2+[-np.pi/2]*npts2)

        x,y,z = spu.sphere2cart(r1_p,theta1,phi=phi)
        r_p1_p,theta_p1_p,phi_p1_p = spu.cart2sphere(x,y,z+self.d_p[0])

        x,y,z = spu.sphere2cart(r1_m,theta1,phi=phi)
        r_p1_m,theta_p1_m,phi_p1_m = spu.cart2sphere(x,y,z+self.d_p[0])

        #### sphere 2
        r2_p=np.array(self.ka+eps)
        r2_m=np.array(self.ka-eps)
        x,y,z = spu.sphere2cart(r2_p,theta1,phi=phi)
        r_p2_p,theta_p2_p,phi_p2_p = spu.cart2sphere(x,y,z+self.d_p[1])
        x,y,z = spu.sphere2cart(r2_m,theta1,phi=phi)
        r_p2_m,theta_p2_m,phi_p2_m = spu.cart2sphere(x,y,z+self.d_p[1])

        #concat
        r_tot_p     = np.hstack((r_p1_p,r_p2_p))
        theta_tot_p = np.hstack((theta_p1_p,theta_p2_p))
        phi_tot_p   = np.hstack((phi_p1_p,phi_p2_p))
        r_tot_m     = np.hstack((r_p1_m,r_p2_m))
        theta_tot_m = np.hstack((theta_p1_m,theta_p2_m))
        phi_tot_m   = np.hstack((phi_p1_m,phi_p2_m))

        err=np.zeros(nmaxs.shape)
        for i,nmax in enumerate(nmaxs):
            print('solving nmax=%d' %nmax)
            self.solve(nmax=nmax,copt=1,v=0)
            f_p=self.compute_f(r_tot_p,theta_tot_p,phi_tot_p,ftype='t')
            f_m=self.compute_f(r_tot_m,theta_tot_m,phi_tot_m,ftype='t')
            f_tot=abs(f_p-f_m)
            err[i]=np.sum(f_tot)
        # print(err)
        #### display
        labs =[r"$\nu_{max}$",r"$\log_{10}(|err|)$"]
        plts=[nmaxs,np.log10(abs(err)),'b-o','']
        dsp.stddisp(plts,labs=labs,lw=2,**kwargs)
