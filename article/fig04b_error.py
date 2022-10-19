from utils import*
import pyScatSpheres.spherical_utils as spu         ;imp.reload(spu)
import utils.displayStandards as dsp                ;imp.reload(dsp)
from pyScatSpheres import qdot_sphere_array as qsa  ;imp.reload(qsa)
path='../docs/figures/'
plt.close('all')

q=qsa.QdotSphereArray(N=2,kp=1.1,kd=2,ka=0.5)

q.test_convergence_error(16,fonts={'lab':30,'tick':20})
