import numpy as np
import matplotlib.pyplot as plt


N1_non = np.genfromtxt('non_interacting_5.txt',skip_header=1, skip_footer=151) 	 # read file, only second line
omegar = N1_non[0]						 											 # omega_r from file
rhomin = N1_non[1] 					 												 # rho_min from file
rhomax = N1_non[2]				 													 # rho_max from file
N = N1_non[3]																		 # n from file
Enrgy_non = N1_non[4]																 # eigenvalue from file


wvnr_non = np.genfromtxt('non_interacting_5.txt',skip_header=3) 	 				 # read file, ignore first two lines, get psi

N1 = np.genfromtxt('interacting_5.txt',skip_header=1, skip_footer=151) 	 			 # read file, only second line
Enrgy = N1[4]


wvnr = np.genfromtxt('interacting_5.txt',skip_header=3) 	 						 # read file, ignore first two lines, get psi

rho = np.linspace(rhomin, rhomax, len(wvnr))

plt.plot(rho, wvnr_non)
plt.plot(rho, wvnr)
plt.xlabel(r'$\rho$', fontsize=15)
plt.ylabel(r'$\psi$',fontsize=15)
plt.title(r'Scaled wavefunction for $\omega_r=5.0$',fontsize=15)
plt.legend(['Non-interacting','Interacting'],fontsize=15)
plt.show()


'''
def plot_formatting():
	"""pretty formatting of plots """
	plt.rc('text',usetex=True) # only works with matplotlib version below 1.5.2, comment out elsewise
	axis_font={'family': 'serif','serif':'Computer Modern Roman','size':16}
	plt.rc('font',**axis_font)
	plt.rc('font',weight ='bold')
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)

plot_formatting()
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

plt.subplot(3,1,1)
plt.plot(x1,u1,x1,v1,linewidth=2.0)
plt.legend(['Analytical','Numerical'],loc ='upper right')
plt.setp(plt.gca().get_legend().get_texts(), fontsize=16)
plt.title('n=10',fontsize=16)
plt.xlabel('$x$')
plt.ylabel('$u(x)$')

plt.subplot(3,1,2)
plt.plot(x2,u2,x2,v2,linewidth=2.0)
plt.legend(['Analytical','Numerical'],loc ='upper right')
plt.setp(plt.gca().get_legend().get_texts(), fontsize=16)
plt.title('n=100',fontsize=16)
plt.xlabel('$x$')
plt.ylabel('$u(x)$')

plt.subplot(3,1,3)
plt.plot(x3,u3,x3,v3,linewidth=2.0)
plt.legend(['Analytical','Numerical'],loc ='upper right')
plt.setp(plt.gca().get_legend().get_texts(), fontsize=16)
plt.title('n=1000',fontsize=16)
plt.xlabel('$x$')
plt.ylabel('$u(x)$')


plt.suptitle('Exact and numerical solution of set of linear equations',fontsize=16)
plt.tight_layout()
plt.subplots_adjust(top=0.88)
plt.show()
'''
