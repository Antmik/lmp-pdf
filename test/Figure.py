import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

pgf_with_rc_fonts = {"pgf.texsystem": "pdflatex"}
matplotlib.rcParams.update(pgf_with_rc_fonts)

from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit


def find_delta(fun_vec,pos_vec):
        roots=[]
        roots_idx=[]
        for i in range(len(fun_vec)-1):
                if (  fun_vec[i] == 0 ):
                        roots.append( pos_vec[i] )
                if (  (fun_vec[i]<0 and fun_vec[i+1]>0 ) or ( fun_vec[i]>0 and fun_vec[i+1]<0 ) ): 
                        roots.append( (pos_vec[i] +pos_vec[i+1])/2 )
                        roots_idx.append(i)
                        
        if(len(roots)==2):
                free_space= roots[1]-roots[0]
        elif(len(roots)==4):
                free_space= roots[1]-roots[0] + roots[3]-roots[2]
        elif(len(roots)==6):
                free_space= roots[1]-roots[0] + roots[3]-roots[2] + roots[5]-roots[4]
                
        return max(roots)-min(roots), int(max(roots)), max(roots_idx), free_space
####delta=find_delta([-1,0,1,0,-2,-3,7],[0,1,2,3,4,5,6])
##print(delta)

        

figWidth = 8.6*0.39
figHeigh = 6.88*0.39

n_bin=200
n_bin_start=37
n_bin_end=n_bin-130#21
#n_wall_start=30
#
#n_bin_center= n_bin_start+ int((n_bin_end - n_bin_start)/2)
#
#n_approx_start=34+6
#n_approx_end=n_bin-44-5
#
#n_wall_end=n_bin_end#n_bin_start+70#n_bin_start+20
#
#h=12.5
#shear_rate=np.array([2,6,10])/h
#
#T=np.array([1.0,2.0,3.0,4.0])
#epsilon=np.array([0.6,1.0,1.4])

################################################################################################################
##Omega22=np.zeros(4)
##Omega23=np.zeros(4)
##Omega24=np.zeros(4)
##
##B22=[ 2.3508044 , 0.50110649, -0.47193769,  0.15806367, -2.6367184*10**(-2), 1.8120118*10**(-3)]
##B23=[-1.8569443 , 0.96985775, -0.39888526, 0.090063692, -1.09181991*10**(-2), 0.56646797*10**(-3)]
##B24=[-0.67406115, 0.42671907, -0.10177069, 0.00061857136, 0.31225358*10**(-2), -0.35206051*10**(-3)]
##
##C22=[1.6330213  , -0.69795156, 0.16096572, -2.2109440*10**(-2), 1.7031434*10**(-3), -0.56699986*10**(-4)]
##C23=[-1.4586197 , 0.52947262, -0.11946363, 1.6264589*10**(-2), -1.2354315*10**(-3), 0.40366357*10**(-4)]
##C24=[-0.62774499, 0.20700644, -0.047601690, 0.67153792*10**(-2), -0.52706167*10**(-3), 0.17705708*10**(-4) ]
##
##for i_T in range(4):
##        Omega22[i_T]=-0.92032979 
##        Omega23[i_T]= 2.5955799
##        Omega24[i_T]= 1.6042745
##        for i in range(6):
##                Omega22[i_T] += B22[i] / (T[i_T]**(i+1)) + C22[i] * (np.log(T[i_T]))**(i+1)
##                Omega23[i_T] += B23[i] / (T[i_T]**(i+1)) + C23[i] * (np.log(T[i_T]))**(i+1)
##                Omega24[i_T] += B24[i] / (T[i_T]**(i+1)) + C24[i] * (np.log(T[i_T]))**(i+1)
##
##                
##                
##f11 = 4 * Omega22
##f12 = 7 * Omega22 - 8* Omega23
##f22 = 301/12 * Omega22 -28* Omega23 + 20 * Omega24
##
##f=1+ f12*f12 / (f11*f22 - f12*f12 )
##
##eta0=5/16/np.sqrt(np.pi)* np.sqrt(T) * f / Omega22 #vector of eta0
##
##print(eta0)
##############################################################################################################



bin_vec=np.zeros((3,4,n_bin))
density=np.zeros((3,4,n_bin))
viscosity=np.zeros((3,4,n_bin))
velocity_grad=np.zeros((3,4,n_bin))
velocity=np.zeros((3,4,n_bin))
stressK12=np.zeros((3,4,n_bin))
stressV12=np.zeros((3,4,n_bin))
temperature=np.zeros((3,4,n_bin))

#bin_vec[0,0,:]=np.genfromtxt("epsilon06/T10/bin_vec.txt",skip_header=2, dtype=None)
#bin_vec[0,1,:]=np.genfromtxt("epsilon06/T20/bin_vec.txt",skip_header=2, dtype=None)
#bin_vec[0,2,:]=np.genfromtxt("epsilon06/T30/bin_vec.txt",skip_header=2, dtype=None)
#bin_vec[0,3,:]=np.genfromtxt("epsilon06/T40/bin_vec.txt",skip_header=2, dtype=None)
#
#bin_vec[1,0,:]=np.genfromtxt("epsilon10/T10/bin_vec.txt",skip_header=2, dtype=None)
#bin_vec[1,1,:]=np.genfromtxt("epsilon10/T20/bin_vec.txt",skip_header=2, dtype=None)
#bin_vec[1,2,:]=np.genfromtxt("epsilon10/T30/bin_vec.txt",skip_header=2, dtype=None)
#bin_vec[1,3,:]=np.genfromtxt("epsilon10/T40/bin_vec.txt",skip_header=2, dtype=None)
#
#bin_vec[2,0,:]=np.genfromtxt("epsilon14/T10/bin_vec.txt",skip_header=2, dtype=None)
#bin_vec[2,1,:]=np.genfromtxt("epsilon14/T20/bin_vec.txt",skip_header=2, dtype=None)
#bin_vec[2,2,:]=np.genfromtxt("epsilon14/T30/bin_vec.txt",skip_header=2, dtype=None)
#bin_vec[2,3,:]=np.genfromtxt("epsilon14/T40/bin_vec.txt",skip_header=2, dtype=None)

bin_vec=np.genfromtxt("bin_vec.txt",skip_header=2, dtype=None)

pdf_bin_vec=np.genfromtxt("pdf_bin_vec.txt",skip_header=2, dtype=None)

pdf_vy=np.genfromtxt("pdf_vy.txt",skip_header=2, dtype=None)



#######################################################################################################

color_vec=['tab:red','tab:blue', 'tab:green', 'tab:orange', 'tab:brown', 'tab:purple', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
line_style=['-','-','-','-']
marker_style=['x','o','+','d']




######################################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex= True)
plt.rc('font', family='serif')

ax = fig.gca(projection='3d')

X, Y = np.meshgrid(pdf_bin_vec, bin_vec[n_bin_start:n_bin_end])
Z=pdf_vy[n_bin_start:n_bin_end,:]

print(X.shape)
print(Y.shape)
print(Z.shape)
ax.plot_surface(X,Y, Z,cmap=cm.coolwarm, label='parametric curve')
#    plt.plot(pdf_bin_vec,pdf_vy[i_bin,:],lw=1.,ms=4,mew=.1,mfc='none', label=r"$U_{eff} \ (\epsilon=0.6,T=4.0) $")


plt.ylabel(r'$y-y_{w}$', fontsize=12)
plt.ylabel(r'$ pdf $', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
#plt.legend(loc='upper right', fontsize=7, edgecolor='k', frameon=False)
plt.locator_params(nbins=6, axis='x')
plt.locator_params(nbins=5, axis='y')
#ax.set_ylim(4, 16)
#ax.set_xlim(-5, 5)

ax.view_init(elev=0., azim=90)

#plt.ylim(-2,5)
fig.tight_layout()
fig.savefig('pdf.pdf',dpi=fig.dpi)

#plt.show()

##########################################################################################
#delta vs Wa
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex= True)
plt.rc('font', family='serif')

#plt.plot(Wall[:,0],delta[:,0] ,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
#plt.plot(Wall[:,1],delta[:,1] , lw=1.,ls=line_style[1],marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
#plt.plot(Wall[:,2],delta[:,2] , lw=1.,ls=line_style[2],marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
#plt.plot(Wall[:,3],delta[:,3], lw=1.,ls=line_style[3],marker=marker_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")

plt.plot(Wall[:,0],delta[:,0] ,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
plt.plot(Wall[:,1],delta[:,1] , lw=1.,ls=line_style[1],marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
plt.plot(Wall[:,2],delta[:,2] , lw=1.,ls=line_style[2],marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
plt.plot(Wall[:,3],delta[:,3], lw=1.,ls=line_style[3],marker=marker_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")

plt.plot(Wall[:,0],delta_10_4[:,0] ,lw=1.,ls='--',marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0])
plt.plot(Wall[:,1],delta_10_4[:,1] , lw=1.,ls='--',marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1])
plt.plot(Wall[:,2],delta_10_4[:,2] , lw=1.,ls='--',marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2])
plt.plot(Wall[:,3],delta_10_4[:,3], lw=1.,ls='--',marker=marker_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3])

plt.xlabel(r'$Wa$', fontsize=12)
plt.ylabel(r'$ \delta $', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='upper left', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)

plt.xlim(xmin=0)
plt.locator_params(nbins=6, axis='x') 
plt.locator_params(nbins=5, axis='y')
fig.tight_layout()
fig.savefig('delta.pdf',dpi=fig.dpi)
#plt.show()


##################################################################################################
##########################################################################################
#density-bulk vs Wa
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(Wall[:,0],rho_bulk[:,0] ,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
plt.plot(Wall[:,1],rho_bulk[:,1] , lw=1.,ls=line_style[1],marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
plt.plot(Wall[:,2],rho_bulk[:,2] , lw=1.,ls=line_style[2],marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
plt.plot(Wall[:,3],rho_bulk[:,3], lw=1.,ls=line_style[3],marker=marker_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")

plt.xlabel(r'$Wa$', fontsize=12)
plt.ylabel(r'$ \rho_b $', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='upper left', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
plt.xlim(xmin=0)
plt.locator_params(nbins=7, axis='x')
plt.locator_params(nbins=5, axis='y')
plt.xlim(0,1.5)
plt.ylim(0.81,0.85)
fig.tight_layout()
fig.savefig('density-bulk.pdf',dpi=fig.dpi)
#plt.show()

##################################################################################################
##########################################################################################
#density-Confined vs Wa
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(Wall[:,0]*T[0]**0.5,rho_confined[:,0] ,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
plt.plot(Wall[:,1]*T[1]**0.5,rho_confined[:,1] , lw=1.,ls=line_style[1],marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
plt.plot(Wall[:,2]*T[2]**0.5,rho_confined[:,2] , lw=1.,ls=line_style[2],marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
plt.plot(Wall[:,3]*T[3]**0.5,rho_confined[:,3], lw=1.,ls=line_style[3],marker=marker_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")

plt.plot(Wall[:,0]*T[0]**0.5,rho_free[:,0] ,lw=1.,ls='--',marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0])
plt.plot(Wall[:,1]*T[1]**0.5,rho_free[:,1] , lw=1.,ls='--',marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1])
plt.plot(Wall[:,2]*T[2]**0.5,rho_free[:,2] , lw=1.,ls='--',marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2])
plt.plot(Wall[:,3]*T[3]**0.5,rho_free[:,3] , lw=1.,ls='--',marker=marker_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3])

plt.xlabel(r'$Wa \ \sqrt{T}$', fontsize=12)
plt.ylabel(r'$ \rho_c / \eta_b , \rho_f / \eta_b $', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='upper left', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
plt.xlim(xmin=0)
plt.locator_params(nbins=7, axis='x')
plt.locator_params(nbins=5, axis='y')
plt.xlim(0,1.5)
plt.ylim(0.8,1.1)
fig.tight_layout()
fig.savefig('density-confined.pdf',dpi=fig.dpi)
#plt.show()

##################################################################################################
##########################################################################################
#density-Free vs Wa
#########################################################################################
#fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#
#plt.plot(epsilon,rho_free[:,0] ,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
#plt.plot(epsilon,rho_free[:,1] , lw=1.,ls=line_style[1],marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
#plt.plot(epsilon,rho_free[:,2] , lw=1.,ls=line_style[2],marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
#plt.plot(epsilon,rho_free[:,3] , lw=1.,ls=line_style[3],marker=marker_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")
#
#plt.xlabel(r'$\epsilon_(f,w)$', fontsize=12)
#plt.ylabel(r'$ \rho_f / \eta_b $', fontsize=12)
#plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
#plt.legend(loc='upper left', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
#plt.locator_params(nbins=6, axis='x')
#plt.locator_params(nbins=5, axis='y')
##plt.ylim(0.8,1.2)
#fig.tight_layout()
#fig.savefig('rho-free.pdf',dpi=fig.dpi)
##plt.show()

##################################################################################################
##########################################################################################
#viscosity-bulk vs Wa
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(Wall[:,0],eta_bulk[:,0] ,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
plt.plot(Wall[:,1],eta_bulk[:,1] , lw=1.,ls=line_style[1],marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
plt.plot(Wall[:,2],eta_bulk[:,2] , lw=1.,ls=line_style[2],marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
plt.plot(Wall[:,3],eta_bulk[:,3], lw=1.,ls=line_style[3],marker=marker_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")

plt.xlabel(r'$Wa$', fontsize=12)
plt.ylabel(r'$\eta_b$', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='upper left', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
plt.xlim(xmin=0)
plt.locator_params(nbins=7, axis='x')
plt.locator_params(nbins=5, axis='y')
#plt.xlim(0,1.5)
plt.ylim(1.8,3.0)
fig.tight_layout()
fig.savefig('viscosity-bulk.pdf',dpi=fig.dpi)
#plt.show()


##################################################################################################
##########################################################################################
#viscosity-Confined vs Wa
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot((epsilon-1)*T[0]**(-1),eta_confined[:,0] ,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
plt.plot((epsilon-1)*T[1]**(-1),eta_confined[:,1] , lw=1.,ls=line_style[1],marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
plt.plot((epsilon-1)*T[2]**(-1),eta_confined[:,2] , lw=1.,ls=line_style[2],marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
plt.plot((epsilon-1)*T[3]**(-1),eta_confined[:,3], lw=1.,ls=line_style[3],marker=marker_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")

plt.plot((epsilon-1)*T[0]**(-1),eta_free[:,0] ,lw=1.,ls='--',marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0])
plt.plot((epsilon-1)*T[0]**(-1),eta_free[:,1] , lw=1.,ls='--',marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1])
plt.plot((epsilon-1)*T[0]**(-1),eta_free[:,2] , lw=1.,ls='--',marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2])
plt.plot((epsilon-1)*T[0]**(-1),eta_free[:,3], lw=1.,ls='--',marker=marker_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3])

plt.xlabel(r'$ \left( \epsilon_{w,f} -1 \right)/T $', fontsize=12)
plt.ylabel(r'$ \eta_c / \eta_b, \ \eta_f / \eta_b $', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='upper left', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
plt.xlim(-0.5,0.5)
plt.locator_params(nbins=6, axis='x')
plt.locator_params(nbins=5, axis='y')
fig.tight_layout()
plt.ylim(0.7,1.5)
fig.savefig('viscosity-confined.pdf',dpi=fig.dpi)
#plt.show()

##################################################################################################
##########################################################################################
#viscosity-Free vs Wa
#########################################################################################
#fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#
#plt.plot(epsilon,eta_free[:,0] ,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
#plt.plot(epsilon,eta_free[:,1] , lw=1.,ls=line_style[1],marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
#plt.plot(epsilon,eta_free[:,2] , lw=1.,ls=line_style[2],marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
#plt.plot(epsilon,eta_free[:,3] , lw=1.,ls=line_style[3],marker=marker_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")
#
#plt.xlabel(r'$\epsilon_(f,w)$', fontsize=12)
#plt.ylabel(r'$ \eta_f / \eta_b $', fontsize=12)
#plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
#plt.legend(loc='upper left', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
#plt.locator_params(nbins=6, axis='x')
#plt.locator_params(nbins=5, axis='y')
#plt.ylim(0.8,1.2)
#fig.tight_layout()
#fig.savefig('viscosity-free.pdf',dpi=fig.dpi)
##plt.show()


##################################################################################################
##########################################################################################
#viscosity-approx
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#plt.plot(bin_vec[0,1,n_bin_start-1:n_bin_end],viscosity[1,0,n_bin_start-1:n_bin_end] ,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
#plt.plot(bin_vec[0,1,n_bin_start-1:n_bin_end],viscosity[2,0,n_bin_start-1:n_bin_end] , lw=1.,ls=line_style[1],marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")

#plt.plot(bin_vec[2,1,n_bin_start-1:n_bin_end],viscosity[2,0,n_bin_start-1:n_bin_end] , lw=1.,ls=line_style[1],marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")


#plt.plot(bin_vec[0,1,n_bin_start:n_bin_end],viscosity[1,2,n_bin_start:n_bin_end] , lw=1.,ls=line_style[2],marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
#plt.plot(bin_vec[0,1,n_bin_start:n_bin_end],viscosity[1,3,n_bin_start:n_bin_end] , lw=1.,ls=line_style[3],marker=marker_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")

#plt.plot(bin_vec[0,1,n_bin_start:n_bin_end],eta_approx[1,0,n_bin_start:n_bin_end] ,lw=1.,ls=line_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
#plt.plot(bin_vec[0,1,n_bin_start:n_bin_end],eta_approx[2,0,n_bin_start:n_bin_end] , lw=1.,ls=line_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")

#plt.plot(bin_vec[0,1,n_bin_start:n_bin_end],eta_approx[1,2,n_bin_start:n_bin_end] , lw=1.,ls=line_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
#plt.plot(bin_vec[0,1,n_bin_start:n_bin_end],eta_approx[1,3,n_bin_start:n_bin_end] , lw=1.,ls=line_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")


plt.xlabel(r'$Wa$', fontsize=12)
plt.ylabel(r'$  \frac{ \left \Vert \eta(y)-\hat{\eta}(y) \right \Vert }{\eta^l} $', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='upper left', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
plt.locator_params(nbins=6, axis='x') 
plt.locator_params(nbins=5, axis='y') 
fig.tight_layout()
fig.savefig('viscosity-approx.pdf',dpi=fig.dpi)
#plt.show()

##################################################################################################
##########################################################################################
#viscosity-approx
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.semilogy(bin_vec[0,1,n_bin_start:n_bin_end],viscosity[0,0,n_bin_start:n_bin_end]/eta_bulk[0,0] ,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
plt.semilogy(bin_vec[0,1,n_bin_start:n_bin_end],viscosity[1,0,n_bin_start:n_bin_end]/eta_bulk[1,0] , lw=1.,ls=line_style[1],marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
plt.semilogy(bin_vec[0,1,n_bin_start:n_bin_end],viscosity[2,0,n_bin_start:n_bin_end]/eta_bulk[2,0] , lw=1.,ls=line_style[2],marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")

#plt.semilogy(bin_vec[0,1,n_bin_start:n_bin_end],eta_approx[1,0,n_bin_start:n_bin_end] ,lw=1.,ls=line_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
#plt.semilogy(bin_vec[0,1,n_bin_start:n_bin_end],eta_approx[2,0,n_bin_start:n_bin_end] , lw=1.,ls=line_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
#plt.plot(bin_vec[0,1,n_bin_start:n_bin_end],eta_approx[1,2,n_bin_start:n_bin_end] , lw=1.,ls=line_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
#plt.plot(bin_vec[0,1,n_bin_start:n_bin_end],eta_approx[1,3,n_bin_start:n_bin_end] , lw=1.,ls=line_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")


plt.xlabel(r'$Wa$', fontsize=12)
plt.ylabel(r'$  \frac{1}{N} \frac{ \left \Vert \eta(y)-\hat{\eta}(y) \right \Vert }{\eta^l} $', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='upper left', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
#plt.locator_params(nbins=6, axis='x') 
#plt.locator_params(nbins=5, axis='y') 
fig.tight_layout()
fig.savefig('viscosity-approxLog.pdf',dpi=fig.dpi)
#plt.show()

##################################################################################################
##########################################################################################
#viscosity-approx
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(bin_vec[0,1,n_bin_start:n_bin_end],viscosity[2,0,n_bin_start:n_bin_end]-eta_approx[2,0,n_bin_start:n_bin_end] ,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
plt.plot(bin_vec[0,1,n_bin_start:n_bin_end],viscosity[2,2,n_bin_start:n_bin_end] -eta_approx[2,2,n_bin_start:n_bin_end], lw=1.,ls=line_style[1],marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
#plt.plot(bin_vec[0,1,n_bin_start:n_bin_end],viscosity[1,2,n_bin_start:n_bin_end] , lw=1.,ls=line_style[2],marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
#plt.plot(bin_vec[0,1,n_bin_start:n_bin_end],viscosity[1,3,n_bin_start:n_bin_end] , lw=1.,ls=line_style[3],marker=marker_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")

plt.plot(bin_vec[0,1,n_bin_start:n_bin_end],viscosity[2,0,n_bin_start:n_bin_end]-eta_approx_corrected[2,0,n_bin_start:n_bin_end] ,lw=1.,ls='--',marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
plt.plot(bin_vec[0,1,n_bin_start:n_bin_end],viscosity[2,2,n_bin_start:n_bin_end] -eta_approx_corrected[2,2,n_bin_start:n_bin_end], lw=1.,ls='--',marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")

plt.xlabel(r'$Wa$', fontsize=12)
plt.ylabel(r'$  \eta(y)-\hat{\eta}(y) $', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='upper left', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
plt.locator_params(nbins=6, axis='x') 
plt.locator_params(nbins=5, axis='y') 
fig.tight_layout()
fig.savefig('viscosity-diff.pdf',dpi=fig.dpi)
#plt.show()


##################################################################################################
##########################################################################################
#viscosity-density StdErr vs Wa
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(Wall[:,0],StdErr[:,0] ,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
plt.plot(Wall[:,1],StdErr[:,1] , lw=1.,ls=line_style[1],marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
plt.plot(Wall[:,2],StdErr[:,2] , lw=1.,ls=line_style[2],marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
plt.plot(Wall[:,3],StdErr[:,3] , lw=1.,ls=line_style[3],marker=marker_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")

#plt.plot(Wall[:,0],StdErr_corrected[:,0] ,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
#plt.plot(Wall[:,1],StdErr_corrected[:,1] , lw=1.,ls=line_style[1],marker=marker_style[1],ms=4,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
#plt.plot(Wall[:,2],StdErr_corrected[:,2] , lw=1.,ls=line_style[2],marker=marker_style[2],ms=4,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
#plt.plot(Wall[:,3],StdErr_corrected[:,3] , lw=1.,ls=line_style[3],marker=marker_style[3],ms=4,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")


plt.xlabel(r'$Wa$', fontsize=12)
plt.ylabel(r'$  \frac{ \left \Vert \eta(y)-\hat{\eta}(y) \right \Vert }{\bar{\eta}} $', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='upper left', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
plt.xlim(xmin=0)
plt.locator_params(nbins=6, axis='x') 
plt.locator_params(nbins=5, axis='y') 
fig.tight_layout()
fig.savefig('viscosity-densityStdErr.pdf',dpi=fig.dpi)
#plt.show()




##########################################################################################
#viscosity-density
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot( density[1,0,n_approx_start:n_approx_end]/rho_bulk[1,0] , eta_approx[1,0,n_approx_start:n_approx_end]/eta_bulk[1,0] ,lw=1.,ls=line_style[0], color=color_vec[0])
plt.plot( density[1,1,n_approx_start:n_approx_end]/rho_bulk[1,1] , eta_approx[1,1,n_approx_start:n_approx_end]/eta_bulk[1,1] ,lw=1.,ls=line_style[1], color=color_vec[1])
plt.plot( density[1,2,n_approx_start:n_approx_end]/rho_bulk[1,2] , eta_approx[1,2,n_approx_start:n_approx_end]/eta_bulk[1,2] ,lw=1.,ls=line_style[2], color=color_vec[2])
plt.plot( density[1,3,n_approx_start:n_approx_end]/rho_bulk[1,3] , eta_approx[1,3,n_approx_start:n_approx_end]/eta_bulk[1,3] ,lw=1.,ls=line_style[3], color=color_vec[3])

plt.scatter(density[1,0,n_approx_start:n_approx_end]/rho_bulk[1,0] ,viscosity[1,0,n_approx_start:n_approx_end]/eta_bulk[1,0] , color=color_vec[0],s=5.,marker=marker_style[0], label=r"$T=1.0 $")
plt.scatter(density[1,1,n_approx_start:n_approx_end]/rho_bulk[1,1] ,viscosity[1,1,n_approx_start:n_approx_end]/eta_bulk[1,1] , color=color_vec[1],s=5.,marker=marker_style[1], label=r"$T=2.0 $")
plt.scatter(density[1,2,n_approx_start:n_approx_end]/rho_bulk[1,2] ,viscosity[1,2,n_approx_start:n_approx_end]/eta_bulk[1,2] , color=color_vec[2],s=5.,marker=marker_style[2],label=r"$T=3.0 $")
plt.scatter(density[1,3,n_approx_start:n_approx_end]/rho_bulk[1,3] ,viscosity[1,3,n_approx_start:n_approx_end]/eta_bulk[1,3] , color=color_vec[3],s=5.,marker=marker_style[3], label=r"$T=4.0 $")


plt.xlabel(r'$\rho(y)/ \rho_b$', fontsize=12)
plt.ylabel(r'$\eta(y)/ \eta_b$', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='upper left', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
plt.locator_params(nbins=6, axis='x') 
plt.locator_params(nbins=5, axis='y') 
fig.tight_layout()
fig.savefig('viscosity-density.pdf',dpi=fig.dpi)
#plt.show()



##################################################################################################
##########################################################################################
#viscosity-density StdErr vs T
#########################################################################################
##fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
##plt.rc('text', usetex=True)
##plt.rc('font', family='serif')
##
##plt.plot(T,Omega22,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
##plt.plot(T,Omega23,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
##plt.plot(T,Omega24,lw=1.,ls=line_style[0],marker=marker_style[0],ms=4,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
##
##plt.xlabel(r'$Wa$', fontsize=12)
##plt.ylabel(r'$  \frac{1}{N} \frac{ \left \Vert \eta(y)-\hat{\eta}(y) \right \Vert }{\eta^l} $', fontsize=12)
##plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
##plt.legend(loc='upper left', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
##plt.locator_params(nbins=6, axis='x') 
##plt.locator_params(nbins=5, axis='y') 
##fig.tight_layout()
##fig.savefig('Omega.pdf',dpi=fig.dpi)
###plt.show()

##########################################################################################
#density-T
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(bin_vec[0,0,n_bin_start:n_bin_end]-bin_vec[0,0,n_bin_start],density[0,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[0],marker=marker_style[0],ms=2,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $" %shear_rate[0] )
plt.plot(bin_vec[0,1,n_bin_start:n_bin_end]-bin_vec[0,1,n_bin_start],density[0,1,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[1],marker=marker_style[1],ms=2,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $" %shear_rate[1])
plt.plot(bin_vec[0,2,n_bin_start:n_bin_end]-bin_vec[0,2,n_bin_start],density[0,2,n_bin_start:n_bin_end],lw=.5,ls=line_style[2],marker=marker_style[2],ms=2,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0  $" %shear_rate[2])
plt.plot(bin_vec[0,3,n_bin_start:n_bin_end]-bin_vec[0,3,n_bin_start],density[0,3,n_bin_start:n_bin_end],lw=.5,ls=line_style[3],marker=marker_style[3],ms=2,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0  $" %shear_rate[2])

plt.xlabel(r'$y-y_{w}$', fontsize=12)
plt.ylabel(r'$\rho(y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='upper center', fontsize=7, edgecolor='k', frameon=False)
plt.locator_params(nbins=6, axis='x') 
plt.locator_params(nbins=4, axis='y') 
fig.tight_layout()
fig.savefig('densityT.pdf',dpi=fig.dpi)
#plt.show()

##########################################################################################
#density-epsilon
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(bin_vec[0,0,n_bin_start:n_bin_end]-bin_vec[0,0,n_bin_start],density[0,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[0],marker=marker_style[0],ms=2,mew=.1,mfc='none', color=color_vec[0], label=r"$\epsilon_{f,w}=0.6 $" %shear_rate[0] )
plt.plot(bin_vec[1,0,n_bin_start:n_bin_end]-bin_vec[1,0,n_bin_start],density[1,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[1],marker=marker_style[1],ms=2,mew=.1,mfc='none', color=color_vec[1], label=r"$\epsilon_{f,w}=1.0  $" %shear_rate[1])
plt.plot(bin_vec[2,0,n_bin_start:n_bin_end]-bin_vec[2,0,n_bin_start],density[2,0,n_bin_start:n_bin_end],lw=.5,ls=line_style[2],marker=marker_style[2],ms=2,mew=.1,mfc='none', color=color_vec[2], label=r"$\epsilon_{f,w}=1.4 $" %shear_rate[2])

plt.xlabel(r'$y-y_{w}$', fontsize=12)
plt.ylabel(r'$\rho(y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='upper center', fontsize=7, edgecolor='k', frameon=False)
plt.locator_params(nbins=6, axis='x') 
plt.locator_params(nbins=4, axis='y') 
fig.tight_layout()
fig.savefig('densityEpsilon.pdf',dpi=fig.dpi)
#plt.show()

##########################################################################################
#velocity T
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(bin_vec[0,0,n_bin_start:n_bin_end]-bin_vec[0,0,n_bin_start],velocity[0,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[0],marker=marker_style[0],ms=2,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0$" %shear_rate[0] )
plt.plot(bin_vec[0,1,n_bin_start:n_bin_end]-bin_vec[0,1,n_bin_start],velocity[0,1,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[1],marker=marker_style[1],ms=2,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0$" %shear_rate[1])
plt.plot(bin_vec[0,2,n_bin_start:n_bin_end]-bin_vec[0,2,n_bin_start],velocity[0,2,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[2],marker=marker_style[2],ms=2,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0$" %shear_rate[2])
plt.plot(bin_vec[0,3,n_bin_start:n_bin_end]-bin_vec[0,3,n_bin_start],velocity[0,3,n_bin_start:n_bin_end],lw=.5,ls=line_style[3],marker=marker_style[3],ms=2,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0$" %shear_rate[0])

plt.xlabel(r'$y-y_{w}$', fontsize=12)
plt.ylabel(r'$ v(y) $', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='lower center', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
plt.locator_params(nbins=6, axis='x')
plt.locator_params(nbins=4, axis='y')
#plt.ylim(-1,0 )
fig.tight_layout()
fig.savefig('velocityT.pdf',dpi=fig.dpi)
#plt.show()


##########################################################################################
#velocity Epsilon
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(bin_vec[0,0,n_bin_start:n_bin_end]-bin_vec[0,0,n_bin_start],velocity[0,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[0],marker=marker_style[0],ms=2,mew=.1,mfc='none', color=color_vec[0], label=r"$\epsilon_{f,w}=0.6 $"  )
plt.plot(bin_vec[1,0,n_bin_start:n_bin_end]-bin_vec[1,0,n_bin_start],velocity[1,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[1],marker=marker_style[1],ms=2,mew=.1,mfc='none', color=color_vec[1], label=r"$\epsilon_{f,w}=1.0 $" )
plt.plot(bin_vec[2,0,n_bin_start:n_bin_end]-bin_vec[2,0,n_bin_start],velocity[2,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[2],marker=marker_style[2],ms=2,mew=.1,mfc='none', color=color_vec[2], label=r"$\epsilon_{f,w}=1.4 $" )

plt.xlabel(r'$y-y_{w}$', fontsize=12)
plt.ylabel(r'$ \dot{\gamma} (y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='lower center', fontsize=7, edgecolor='k', frameon=False)
plt.locator_params(nbins=6, axis='x')
plt.locator_params(nbins=4, axis='y')
#plt.ylim(-1, 0 )
fig.tight_layout()
fig.savefig('velocityEpsilon.pdf',dpi=fig.dpi)
#plt.show()
##########################################################################################
#shear rate T
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(bin_vec[0,0,n_bin_start:n_bin_end]-bin_vec[0,0,n_bin_start],velocity_grad[0,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[0],marker=marker_style[0],ms=2,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0$" %shear_rate[0] )
plt.plot(bin_vec[0,1,n_bin_start:n_bin_end]-bin_vec[0,1,n_bin_start],velocity_grad[0,1,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[1],marker=marker_style[1],ms=2,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0$" %shear_rate[1])
plt.plot(bin_vec[0,2,n_bin_start:n_bin_end]-bin_vec[0,2,n_bin_start],velocity_grad[0,2,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[2],marker=marker_style[2],ms=2,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0$" %shear_rate[2])
plt.plot(bin_vec[0,3,n_bin_start:n_bin_end]-bin_vec[0,3,n_bin_start],velocity_grad[0,3,n_bin_start:n_bin_end],lw=.5,ls=line_style[3],marker=marker_style[3],ms=2,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0$" %shear_rate[0])

plt.xlabel(r'$y-y_{w}$', fontsize=12)
plt.ylabel(r'$ \dot{\gamma}(y) $', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='lower center', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
plt.locator_params(nbins=6, axis='x') 
plt.locator_params(nbins=4, axis='y') 
plt.ylim(-1,0 )
fig.tight_layout()
fig.savefig('velocity_gradT.pdf',dpi=fig.dpi)
#plt.show()


##########################################################################################
#shear rate Epsilon
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(bin_vec[0,0,n_bin_start:n_bin_end]-bin_vec[0,0,n_bin_start],velocity_grad[0,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[0],marker=marker_style[0],ms=2,mew=.1,mfc='none', color=color_vec[0], label=r"$\epsilon_{f,w}=0.6 $"  )
plt.plot(bin_vec[1,0,n_bin_start:n_bin_end]-bin_vec[1,0,n_bin_start],velocity_grad[1,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[1],marker=marker_style[1],ms=2,mew=.1,mfc='none', color=color_vec[1], label=r"$\epsilon_{f,w}=1.0 $" )
plt.plot(bin_vec[2,0,n_bin_start:n_bin_end]-bin_vec[2,0,n_bin_start],velocity_grad[2,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[2],marker=marker_style[2],ms=2,mew=.1,mfc='none', color=color_vec[2], label=r"$\epsilon_{f,w}=1.4 $" )

plt.xlabel(r'$y-y_{w}$', fontsize=12)
plt.ylabel(r'$ \dot{\gamma} (y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='lower center', fontsize=7, edgecolor='k', frameon=False)
plt.locator_params(nbins=6, axis='x') 
plt.locator_params(nbins=4, axis='y') 
plt.ylim(-1, 0 )
fig.tight_layout()
fig.savefig('velocity_gradEpsilon.pdf',dpi=fig.dpi)
#plt.show()

##########################################################################################
#viscosityT
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(bin_vec[0,0,n_bin_start:n_bin_end]-bin_vec[0,0,n_bin_start],viscosity[0,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[0],marker=marker_style[0],ms=2,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
plt.plot(bin_vec[0,1,n_bin_start:n_bin_end]-bin_vec[0,1,n_bin_start],viscosity[0,1,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[1],marker=marker_style[1],ms=2,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
plt.plot(bin_vec[0,2,n_bin_start:n_bin_end]-bin_vec[0,2,n_bin_start],viscosity[0,2,n_bin_start:n_bin_end],lw=.5,ls=line_style[2],marker=marker_style[2],ms=2,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
plt.plot(bin_vec[0,3,n_bin_start:n_bin_end]-bin_vec[0,3,n_bin_start],viscosity[0,3,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[3],marker=marker_style[3],ms=2,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")

plt.xlabel(r'$y-y_{w}$', fontsize=12)
plt.ylabel(r'$\eta(y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='lower center', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
plt.locator_params(nbins=6, axis='x') 
plt.locator_params(nbins=5, axis='y') 
fig.tight_layout()
fig.savefig('viscosityT.pdf',dpi=fig.dpi)
#plt.show()

##########################################################################################
#viscosity Epsilon
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(bin_vec[0,0,n_bin_start:n_bin_end]-bin_vec[0,0,n_bin_start],viscosity[0,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[0],marker=marker_style[0],ms=2,mew=.1,mfc='none', color=color_vec[0], label=r"$\epsilon_{f,w}=0.6 $")
plt.plot(bin_vec[1,0,n_bin_start:n_bin_end]-bin_vec[1,0,n_bin_start],viscosity[1,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[1],marker=marker_style[1],ms=2,mew=.1,mfc='none', color=color_vec[1], label=r"$\epsilon_{f,w}=1.0 $")
plt.plot(bin_vec[2,0,n_bin_start:n_bin_end]-bin_vec[2,0,n_bin_start],viscosity[2,0,n_bin_start:n_bin_end],lw=.5,ls=line_style[2],marker=marker_style[2],ms=2,mew=.1,mfc='none', color=color_vec[2], label=r"$\epsilon_{f,w}=1.4 $")

plt.xlabel(r'$y-y_{w}$', fontsize=12)
plt.ylabel(r'$\eta(y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='upper center', fontsize=7, edgecolor='k', frameon=False)
plt.locator_params(nbins=6, axis='x') 
plt.locator_params(nbins=5, axis='y') 
fig.tight_layout()
fig.savefig('viscosityEpsilon.pdf',dpi=fig.dpi)
#plt.show()

##########################################################################################
#Temperature
#########################################################################################
fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(bin_vec[0,0,n_bin_start:n_bin_end]-bin_vec[0,0,n_bin_start],temperature[0,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[0],marker=marker_style[0],ms=2,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
plt.plot(bin_vec[0,1,n_bin_start:n_bin_end]-bin_vec[0,1,n_bin_start],temperature[0,1,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[1],marker=marker_style[1],ms=2,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
plt.plot(bin_vec[0,2,n_bin_start:n_bin_end]-bin_vec[0,2,n_bin_start],temperature[0,2,n_bin_start:n_bin_end],lw=.5,ls=line_style[2],marker=marker_style[2],ms=2,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
plt.plot(bin_vec[0,3,n_bin_start:n_bin_end]-bin_vec[0,3,n_bin_start],temperature[0,3,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[3],marker=marker_style[3],ms=2,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")

plt.xlabel(r'$y-y_{w}$', fontsize=12)
plt.ylabel(r'$T(y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='lower center', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
plt.locator_params(nbins=6, axis='x')
plt.locator_params(nbins=5, axis='y')
plt.ylim(0,4.2)
fig.tight_layout()
fig.savefig('temperature0.pdf',dpi=fig.dpi)
#plt.show()

fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(bin_vec[0,0,n_bin_start:n_bin_end]-bin_vec[0,0,n_bin_start],temperature[1,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[0],marker=marker_style[0],ms=2,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
plt.plot(bin_vec[0,1,n_bin_start:n_bin_end]-bin_vec[0,1,n_bin_start],temperature[1,1,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[1],marker=marker_style[1],ms=2,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
plt.plot(bin_vec[0,2,n_bin_start:n_bin_end]-bin_vec[0,2,n_bin_start],temperature[1,2,n_bin_start:n_bin_end],lw=.5,ls=line_style[2],marker=marker_style[2],ms=2,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
plt.plot(bin_vec[0,3,n_bin_start:n_bin_end]-bin_vec[0,3,n_bin_start],temperature[1,3,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[3],marker=marker_style[3],ms=2,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")

plt.xlabel(r'$y-y_{w}$', fontsize=12)
plt.ylabel(r'$T(y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='lower center', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
plt.locator_params(nbins=6, axis='x')
plt.locator_params(nbins=5, axis='y')
plt.ylim(0,4.2)
fig.tight_layout()
fig.savefig('temperature1.pdf',dpi=fig.dpi)
#plt.show()

fig=plt.figure(num=None, figsize=(figWidth,figHeigh), dpi=300, facecolor='w', edgecolor='k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(bin_vec[0,0,n_bin_start:n_bin_end]-bin_vec[0,0,n_bin_start],temperature[2,0,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[0],marker=marker_style[0],ms=2,mew=.1,mfc='none', color=color_vec[0], label=r"$T=1.0 $")
plt.plot(bin_vec[0,1,n_bin_start:n_bin_end]-bin_vec[0,1,n_bin_start],temperature[2,1,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[1],marker=marker_style[1],ms=2,mew=.1,mfc='none', color=color_vec[1], label=r"$T=2.0 $")
plt.plot(bin_vec[0,2,n_bin_start:n_bin_end]-bin_vec[0,2,n_bin_start],temperature[2,2,n_bin_start:n_bin_end],lw=.5,ls=line_style[2],marker=marker_style[2],ms=2,mew=.1,mfc='none', color=color_vec[2], label=r"$T=3.0 $")
plt.plot(bin_vec[0,3,n_bin_start:n_bin_end]-bin_vec[0,3,n_bin_start],temperature[2,3,n_bin_start:n_bin_end] ,lw=.5,ls=line_style[3],marker=marker_style[3],ms=2,mew=.1,mfc='none', color=color_vec[3], label=r"$T=4.0 $")

plt.xlabel(r'$y-y_{w}$', fontsize=12)
plt.ylabel(r'$T(y)$', fontsize=12)
plt.tick_params(axis='both', labelsize=8,color='k' , direction='in')
plt.legend(loc='lower center', fontsize=7, edgecolor='k', frameon=False, ncol=2, columnspacing=1.)
plt.locator_params(nbins=6, axis='x')
plt.locator_params(nbins=5, axis='y')
plt.ylim(0,4.2)
fig.tight_layout()
fig.savefig('temperature2.pdf',dpi=fig.dpi)
#plt.show()

