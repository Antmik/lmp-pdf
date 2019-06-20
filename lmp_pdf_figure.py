#!/home/ar6116/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
import os.path

if os.path.isfile("density.txt"):
	density=np.genfromtxt("density.txt",skip_header=2, dtype=None)
	fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
	#plt.plot(bin_vec,density,lw=1.5,ls='-',label='w=0.001')
	plt.plot(density,lw=1.5,ls='-')
	plt.xlabel(r'y', fontsize=12)
	plt.ylabel(r'$\rho(y)$', fontsize=12)
	plt.tick_params(axis='both', labelsize=10)
	plt.legend()
	fig.tight_layout()
	fig.savefig('rho.pdf',dpi=fig.dpi)
	plt.show()

if os.path.isfile("velocity.txt"):
	velocity=np.genfromtxt("velocity.txt",skip_header=2, dtype=None)
	fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
	plt.plot(velocity,lw=1.5,ls='-')
	plt.xlabel(r'y', fontsize=12)
	plt.ylabel(r'$v_x(y)$', fontsize=12)
	plt.tick_params(axis='both', labelsize=10)
	plt.legend()
	fig.tight_layout()
	fig.savefig('v.pdf',dpi=fig.dpi)
	plt.show()

if os.path.isfile("temperature.txt"):
	temperature=np.genfromtxt("temperature.txt",skip_header=2, dtype=None)
	fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
	plt.plot(temperature,lw=1.5,ls='-')
	plt.xlabel(r'y', fontsize=12)
	plt.ylabel(r'$T(y)$', fontsize=12)
	plt.tick_params(axis='both', labelsize=10)
	plt.legend()
	fig.tight_layout()
	fig.savefig('temp.pdf',dpi=fig.dpi)
	plt.show()

if os.path.isfile("stressK22.txt") and os.path.isfile("stressV22.txt"):
	stressK22=np.genfromtxt("stressK22.txt",skip_header=2, dtype=float)
	stressV22=np.genfromtxt("stressV22.txt",skip_header=2, dtype=float)
	stressT22=stressV22+stressK22
	fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
	plt.plot(stressV22,lw=1.5,ls='-',label=r'$\sigma^v_{yy}(y)$')
	plt.plot(stressK22,lw=1.5,ls='-',label=r'$\sigma^k_{yy}(y)$')
	plt.plot(stressT22,lw=1.5,ls='-',label=r'$\sigma^t_{yy}(y)$')
	plt.xlabel(r'y', fontsize=12)
	plt.ylabel(r'$\sigma(y)$', fontsize=12)
	plt.tick_params(axis='both', labelsize=10)
	plt.legend()
	fig.tight_layout()
	fig.savefig('stress22.pdf',dpi=fig.dpi)
	plt.show()

if os.path.isfile("stressK12.txt") and os.path.isfile("stressV12.txt"):
	stressK12=np.genfromtxt("stressK12.txt",skip_header=2, dtype=float)
	stressV12=np.genfromtxt("stressV12.txt",skip_header=2, dtype=float)
	stressT12=stressV12+stressK12
	fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
	plt.plot(stressV12,lw=1.5,ls='-',label=r'$\sigma^v_{yy}(y)$')
	plt.plot(stressK12,lw=1.5,ls='-',label=r'$\sigma^k_{yy}(y)$')
	plt.plot(stressT12,lw=1.5,ls='-',label=r'$\sigma^t_{yy}(y)$')
	plt.xlabel(r'y', fontsize=12)
	plt.ylabel(r'$\sigma(y)$', fontsize=12)
	plt.tick_params(axis='both', labelsize=10)
	plt.legend()
	fig.tight_layout()
	fig.savefig('stress12.pdf',dpi=fig.dpi)
	plt.show()

if os.path.isfile("viscosity.txt"):
	viscosity=np.genfromtxt("viscosity.txt",skip_header=2, dtype=None)
	fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
	plt.plot(viscosity,lw=1.5,ls='-')
	plt.xlabel(r'y', fontsize=12)
	plt.ylabel(r'$eta_x(y)$', fontsize=12)
	plt.tick_params(axis='both', labelsize=10)
	plt.ylim(0,10)
	plt.legend()
	fig.tight_layout()
	fig.savefig('viscosity.pdf',dpi=fig.dpi)
	plt.show()

	density=np.genfromtxt("density.txt",skip_header=2, dtype=None)
	fig=plt.figure(num=None, figsize=(10*0.39, 8*0.39), dpi=300, facecolor='w', edgecolor='k')
	plt.plot(density,viscosity,lw=1.5,ls='-')
	plt.xlabel(r'y', fontsize=12)
	plt.ylabel(r'$eta_x(y)$', fontsize=12)
	plt.tick_params(axis='both', labelsize=10)
	plt.xlim(0,2)
	plt.ylim(0,15)
	plt.legend()
	fig.tight_layout()
	fig.savefig('viscosity-density.pdf',dpi=fig.dpi)
	plt.show()