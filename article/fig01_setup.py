from utils import *   ;imp.reload(dsp)
from matplotlib.patches import Circle
plt.close('all')
name='figures/Tmatrix.png'
opt='ps'

txts,plts=[],[]
cs={1:[(0.75,)*3,0.8],6:[(0.1,)*3,1],7:['b',1.1],8:['r',1.2]}

#### Sphere locations, radius, n_p
spheres = {
    'p':[(3.0,6.0), 6,0.5],'q':[(5.0,2.0), 7,0.5],
    '2':[(8.5,4.5), 6,0.2],'N':[(8.5,8.0), 8,0.2],'1':[(7.25,0.75),1,0.2],}
spheres = {k:v+cs[v[1]] for k,v in spheres.items()}
pps = [Circle(yz,radius=r,color=c,alpha=a) for (yz,Z,a,c,r) in spheres.values()]
txts+= [[yz[0]+r/3,yz[1]-r/3,'$n_%s$' %i,(1-2*a,)*3] for i,(yz,Z,a,c,r) in spheres.items()]

####radius arrows
theta = 1.15*3*np.pi/4#+np.pi/2
rm=0.91
ct,st = np.cos(theta),np.sin(theta)
arrows = [ [yz[0],yz[1],rm*r*ct,rm*r*st,(1-2*a,)*3] for i,(yz,Z,a,c,r) in spheres.items()]
txts+= [[yz[0]+0.7*r*ct,yz[1]+0.7*r*st+0.05,'$a_%s$' %i,(1-2*a,)*3] for i,(yz,Z,a,c,r) in spheres.items()]

#### triedres
c,Nf = (0.35,)*3,3
u,(y,z) = 1,(0,0)
arrows+= [ [y,z,u*i,u*j,c] for i,j in zip([0,9,-0.6*Nf],[9,0,-0.5*Nf])]
txts+=[[u*i+0.1,u*j-0.1,s,c] for i,j,s in zip([0,9,-0.6*Nf],[9,0,-0.5*Nf],['$z$','$y$','$x$']) ]
txts[-1][0]-=0.25;txts[-1][1]+=0.3
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
txtsB  = [[yP+0.1,zP,'$P$','k']]
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
t,Ra  = np.linspace(0,theta_i,20),1.55
plts+=[[Ra*np.sin(t),Ra*np.cos(t),'r-','',2]]
txts+= [[0.5,Ra+0.1,r'$\theta_i$','r']]
ya,za = 1.1*Ra*np.sin(theta_i),1.1*Ra*np.cos(theta_i)
arrows+= [[0,0,ya,za,'r']]
txts+= [[ya+0.2,za-0.3,r'$\vec{k_0}$','r']]

##### phi,phi,Phi_q
##phi
ta = 0.5/0.6
Rab=4.5#Rb1,Rb2 = 1.2,0.25
phi= np.deg2rad(75)
t  = np.linspace(0,phi,20)
xzP,xP = -1.4,-0.5
plts+=[[xzP/ta-Rab*xzP*np.sin(t),xzP*np.cos(t),'k-','',2]]
# plts+=[[Rb1*yP*np.cos(t),Rb2*yP*np.sin(t),'k-','',2]]
plts+=[[[yP,yP],[xP,zP],'k--','',1]]
plts+=[[[0,yP],[0,xP],'k--','',1]]
txts+= [[yP/2+0.1,xzP+0.1,r'$\phi$','k']]
##Phi_q
Phi_q= np.deg2rad(65)
t  = np.linspace(0,Phi_q,20)
xzq,xq = -0.9,-0.7
plts+=[[xzq/ta-Rab*xzq*np.sin(t),xzq*np.cos(t),'k-','',2]]
# plts+=[[Rb1*yq*np.cos(t),Rb2*yq*np.sin(t),'k-','',2]]
plts+=[[[yq,yq],[xq,zq+0.1],'k--','',1]]
plts+=[[[0,yq],[0,xq],'k--','',1]]
txts+= [[1.5,xzq-0.1,r'$\Phi_q$','k']]
##phi_i
phi_i= np.deg2rad(42)
t = np.linspace(0,phi_i,20)
xza,xa,ya = -0.4,-0.55,ya+0.1
plts+=[[xza/ta-Rab*xza*np.sin(t),xza*np.cos(t),'r-','',2]]
# plts+=[[Rb1*ya*np.cos(t),Rb2*ya*np.sin(t),'r-','',2]]
plts+=[[[ya,ya],[xa,za+0.1],'r--','',1]]
plts+=[[[0,ya],[0,xa],'r--','',1]]
txts+= [[0,-0.75,r'$\phi_i$','r']]


fig,ax = dsp.stddisp(texts=txtsB,targs={'ha':'left','va':'center','weight':'bold'},opt='p')
dsp.stddisp(plts,ax=ax,arrows=arrows,texts=txts,patches=pps,targs={'ha':'left'},
    xylims=[-1,10,-1,10],pOpt='peX',axPos=[0,0,1,1],fonts={'text':25},
    name=name,opt=opt)
