from utils import *   ;imp.reload(dsp)
from matplotlib.patches import Circle
plt.close('all')
name='figures/Tmatrix.png'

txts,plts=[],[]
cs={1:[(0.75,)*3,0.6],6:[(0.1,)*3,1],7:['b',1.1],8:['r',1.2]}

#### Sphere locations, Circles and n_p
spheres = {
    'p':[(3,6), 6],'q':[(5,2), 7],
    '2':[(8.5,4.5), 6],'N':[(8.5,8), 8],'1':[(7.25,0.75),1],}
spheres = {k:v+cs[v[1]] for k,v in spheres.items()}
pps = [Circle(yz,radius=r,color=c,alpha=0.5) for (yz,Z,c,r) in spheres.values()]
txts+= [[yz[0]+r/3,yz[1]-r/3,'$n_%s$' %i,'k'] for i,(yz,Z,c,r) in spheres.items()]

####radius arrows
theta = 1.15*3*np.pi/4#+np.pi/2
rm=0.92
ct,st = np.cos(theta),np.sin(theta)
arrows = [ [yz[0],yz[1],rm*r*ct,rm*r*st,'k'] for i,(yz,Z,c,r) in spheres.items()]
txts+= [[yz[0]+0.7*r*ct,yz[1]+0.7*r*st+0.05,'$a_%s$' %i,'k'] for i,(yz,Z,c,r) in spheres.items()]

#### triedres
c = (0.35,)*3
u,(y,z) = 1,(0,0)
arrows+= [ [y,z,u*i,u*j,c] for i,j in zip([0,9,-0.6],[9,0,-0.5])]
txts+=[[u*i+0.1,u*j-0.1,s,c] for i,j,s in zip([0,9,-0.6],[9,0,-0.5],['$z$','$y$','$x$']) ]
u,(y,z) = 1.5,spheres['p'][0]
arrows+= [[y,z,u*i,u*j,c] for i,j in zip([0,1,-0.6],[1,0,-0.5])]
u,(y,z) = 1.5,spheres['q'][0]
arrows+= [[y,z,u*i,u*j,c] for i,j in zip([0,1,-0.6],[1,0,-0.5])]

#### distances
yP,zP = [6.25,7.5]
yp,zp = spheres['p'][0]
yq,zq = spheres['q'][0]
plts = [[[0,yp],[0,zp],'k-','']]
plts+= [[[0,yq],[0,zq],'k-','']]
plts+= [[[0,yP],[0,zP],'k-','']]
plts+= [[[yp,yq],[zp,zq],'k-','']]
plts+= [[[yp,yP],[zp,zP],'k-','']]
plts+= [[[yq,yP],[zq,zP],'k-','']]
txtsB  = [[yP,zP,'$P$','k']]
txtsB += [[0.8*yP/2+0.2,0.8*zP/2  ,r'$\vec{r}$'     ,'k']]
txtsB += [[1.3*yp/2+0.1,1.3*zp/2  ,r'$\vec{d_p}$'   ,'k']]
txtsB += [[1.3*yq/2,1.3*zq/2-0.25 ,r'$\vec{d_q}$'   ,'k']]
txtsB += [[(yq+yp)/2+0.1,(zp+zq)/2,r'$\vec{d_{pq}}$','k']]
txtsB += [[(yp+yP)/2,(zp+zP)/2+0.3,r'$\vec{r_p}$'   ,'k']]
txtsB += [[(yq+yP)/2+0.1,(zq+zP)/2,r'$\vec{r_q}$'   ,'k']]
#### arcs theta
theta_p = lambda yp,zp:np.arccos(zp/np.sqrt(yp**2+zp**2))
t,Rp = np.linspace(0,theta_p(yp,zp),20),5
plts+=[[Rp*np.sin(t),Rp*np.cos(t),'k-','',2]]
txts+= [[0.5,Rp+0.1,r'$\Theta_{p}$','k']]
t,R  = np.linspace(0,theta_p(yP,zP),20),4
plts+=[[R*np.sin(t),R*np.cos(t),'k-','',2]]
txts+= [[0.5,R+0.1,r'$\theta$','k']]
t,Rq = np.linspace(0,theta_p(yq,zq),20),3
plts+=[[Rq*np.sin(t),Rq*np.cos(t),'k-','',2]]
txts+= [[0.5,Rq+0.1,r'$\Theta_{q}$','k']]
t,R  = np.linspace(0,theta_p(yP-yp,zP-zp),20),1.3
plts+=[[yp+R*np.sin(t),zp+R*np.cos(t),c,'',2]]
txts+= [[yp+0.5,zp+R+0.1,r'$\theta_p$',c]]

#### incident wave
## theta_i
theta_i = np.deg2rad(50)
t,Ra  = np.linspace(0,theta_i,20),2
plts+=[[Ra*np.sin(t),Ra*np.cos(t),'r-','',2]]
txts+= [[0.5,Ra+0.1,r'$\theta_i$','r']]
ya,za = 1.1*Ra*np.sin(theta_i),1.1*Ra*np.cos(theta_i)
arrows+= [[0,0,ya,za,'r']]
txts+= [[ya+0.2,za-0.3,r'$\vec{k_0}$','r']]

##### phi,phi,Phi_q
##phi
Rb1,Rb2 = 0.9,0.25
phi= np.deg2rad(-16)
t  = np.linspace(0,phi,20)
xP = -0.5
plts+=[[Rb1*yP*np.cos(t),Rb2*yP*np.sin(t),'k-','',2]]
plts+=[[[yP,yP],[xP,zP],'k--','',1]]
plts+=[[[0,yP],[0,xP],'k--','',1]]
txts+= [[yP-0.5,xP+0.15,r'$\phi$','k']]
##Phi_q
Phi_q= np.deg2rad(-30)
t  = np.linspace(0,Phi_q,20)
xq = -0.8
plts+=[[Rb1*yq*np.cos(t),Rb2*yq*np.sin(t),'k-','',2]]
plts+=[[[yq,yq],[xq,zq+0.1],'k--','',1]]
plts+=[[[0,yq],[0,xq],'k--','',1]]
txts+= [[yq-0.5,xq+0.15,r'$\Phi_q$','k']]
##phi_i
phi_i= np.deg2rad(-56.5)
t = np.linspace(0,phi_i,20)
xa,ya = -0.75,ya+0.1
plts+=[[Rb1*ya*np.cos(t),Rb2*ya*np.sin(t),'r-','',2]]
plts+=[[[ya,ya],[xa,za+0.1],'r--','',1]]
plts+=[[[0,ya],[0,xa],'r--','',1]]
txts+= [[ya-0.35,xa+0.25,r'$\phi_i$','r']]


fig,ax = dsp.stddisp(texts=txtsB,targs={'ha':'left','va':'center','weight':'bold'},opt='p')
dsp.stddisp(plts,ax=ax,arrows=arrows,texts=txts,patches=pps,targs={'ha':'left'},
    xylims=[-1,10,-1,10],pOpt='ptX',axPos=[0,0,1,1],fonts={'text':25},
    name=name,opt='ps')
