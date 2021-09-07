import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt
from math import factorial
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import imageio
images = []
images_contour = []

beta_low=0.6
beta_high=0.85
nbetas=64
rho=1.
e=0.
h=1
nu=0.39
eta=0
nMAX=30
H_init=0

L=16

MaxP=129

if( (nu).is_integer()): nu=int(nu)
if( (e).is_integer()): e=int(e)
if( (rho).is_integer()): rho=int(rho)

BASEDIR=("/Users/ilaria/Desktop/Multicomponent_VillainModel/Output_Villain_2C/e_%s/nu_%s/h_%s" %(e, nu, h))
folder_in=("%s/L%d_rho%s_eta%s_e%s_h%s_nu%s_bmin%s_bmax%s_nMAX%s_init%s"  %(BASEDIR, L, rho, eta, e,  h, nu, beta_low, beta_high, nMAX, H_init))

beta=np.zeros((nbetas))

for b in range(nbetas):
	beta[b]=beta_low +b*(beta_high -beta_low)/(nbetas-1)

	filename=("%s/beta_%s/Villain_potential.txt" %(folder_in,b))

	x, y, potential1=np.loadtxt(filename, usecols=(0,1,2), unpack=True)
	data=np.loadtxt(filename, usecols=(0,1,2))
	
	potential=np.reshape(potential1, (MaxP, MaxP))
	
	fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
	fig.suptitle(r"$\nu=%s; \,\, \rho=%s; \,\, \beta=%s$" %(nu, rho, beta[b]))
	
	axs = plt.axes(projection='3d')
	
	x_surf = np.outer(np.arange(-64,65,1), np.ones(MaxP))
	print(x_surf, np.shape(x_surf))
	y_surf = x_surf.copy().T # transpose
	axs.set_zlim(-0.2, 3.5)
	axs.set_xlabel(r'$u_1$')
	axs.set_ylabel(r'$u_2$')
	axs.set_zlabel(r'$V[u_1, u_2; \beta]$')
	surf = axs.plot_surface(x_surf,y_surf, potential, cmap="viridis", linewidth=0, antialiased=False)
	
	fig.subplots_adjust(left=0, bottom=0, right=0.9, top=0.9)
	
	fig.tight_layout()
	figname=("%s/Villain_Potential_beta%s.png" %(BASEDIR, b))
	fig.savefig("%s/Villain_Potential_beta%s.png" %(BASEDIR, b))
	plt.close()
	#plt.show()

	fig1, axs1 = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
	axs1.set_title(r"$\nu=%s; \,\, \rho=%s; \,\, \beta=%s$" %(nu, rho, beta[b]))
	axs1.axis('equal')
	axs1.set(xlim=(-(MaxP-1)/2, (MaxP-1)/2), ylim=(-(MaxP-1)/2, (MaxP-1)/2))
	counts = axs1.contour(np.arange(-64,65,1),np.arange(-64,65,1) , potential)
	axs1.clabel(counts, inline=1, fontsize=10)
	axs1.axvline(x=0, ls="--", c="gray")
	axs1.axhline(y=0, ls="--", c="gray")

	fig1.subplots_adjust(left=0, bottom=0, right=0.9, top=0.8)
	
	fig1.tight_layout()
	figname_contour=("%s/Villain_Potential_contour_beta%s.png" %(BASEDIR, b))
	fig1.savefig("%s/Villain_Potential_contour_beta%s.png" %(BASEDIR, b))
	plt.close()

	images.append(imageio.imread(figname))
	images_contour.append(imageio.imread(figname_contour))


imageio.mimsave('%s/VillainPotential.gif' %(BASEDIR), images)
imageio.mimsave('%s/VillainPotential_contour.gif' %(BASEDIR), images_contour)

	#plt.show()
