from numpy import*
import matplotlib as mpl
mpl.use('TkAgg') # Ensure that the Tkinter backend is used for generating figures
from matplotlib.pyplot import * 

rho_max = array([1.,2.,3.,4.,5.,6.])
known_eig = np.array([3.,7.,11.])
lambda1=array([10.151,3.5296,3.0121,2.9999,2.9997,2.9996])
lambda2=array([39.795,11.168,7.3278,7.0025,6.9987,6.9981])
lambda3=array([89.132,23.524,12.944,11.077,10.997,10.995])



plot(rho_max,lambda1,'b-x',label='Approx. $\lambda=3$')
plot(rho_max,lambda2,'g-x',label='Approx. $\lambda=7$')
plot(rho_max,lambda3,'r-x',label='Approx. $\lambda=11$')


axhline(known_eig[0],color='c',ls='--',label='$\lambda=3$')
axhline(known_eig[1],color='m',ls='--',label='$\lambda=7$')
axhline(known_eig[2],color='k',ls='--',label='$\lambda=11$')
xlabel(r'$\rho_{max}$', fontsize=10)
ylabel('Eigenvalues $\lambda$',fontsize=10)
title(r'Eigenvalues as a function of $\rho_{max}$, $n=170$',fontsize=10)
legend()
show()
