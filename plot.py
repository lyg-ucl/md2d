import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#--------------------------------------------------------------------------------
class plotD:

    def __init__(self, D_t=None, TypeName=None):
       self.D_t = np.array(D_t)
       self.TypeName = TypeName
       self.draw()

    def draw(self):
        layout= (1, 1)
        fig, ax = plt.subplots(*layout,figsize=(8,6))
        plt.subplots_adjust(wspace = 0.15,hspace = 0.3, bottom=0.1, left=0.1,right=0.95, top=0.9)
        plt.style.use('seaborn')

        ax.set_xscale('log', base=2)
        ax.set_xlabel('Time interval (fs)');
        ax.set_ylabel('D (m$^2$/s)')

        for i in range(len(self.D_t[0,:])-1):
            ax.plot(self.D_t[:,0], self.D_t[:,i+1], linewidth=3, label=self.TypeName[i])
    
        ax.legend(loc='best', shadow=True,fontsize='medium',frameon=False)
        plt.savefig('D_tau.pdf',dpi=600)
        plt.show()

class plotsacf:

    def __init__(self, sacf=None ):
       self.sacf = np.array(sacf)
       self.draw()

    def draw(self):
        layout= (1, 1)
        fig, ax = plt.subplots(*layout,figsize=(8,6))
        plt.subplots_adjust(wspace = 0.15,hspace = 0.3, bottom=0.1, left=0.1,right=0.95, top=0.9)
        plt.style.use('seaborn')

        ax.set_xlabel('MD step');
        ax.set_ylabel('SACF (GPa$^2$)')
        ax.set_xscale('log')
        x=list(range(len(self.sacf[:,0])))
        ax.plot(x, self.sacf[:,0], linewidth=3, label="bulk")
        ax.plot(x, self.sacf[:,1], linewidth=3, label="shear")
        ax.legend(loc='best', shadow=True,fontsize='medium',frameon=False)
        plt.savefig('SACF.pdf',dpi=600)
        plt.show()

class plotmsd:

    def __init__(self,msd=None, Niter=None,potim=None,TypeName=None):
        self.msd = np.array(msd)
        self.Niter = Niter
        self.potim = potim
        self.draw()

    def draw(self):
        layout= (1, 1)
        fig, ax = plt.subplots(*layout,figsize=(8,6))
        plt.subplots_adjust(wspace = 0.15,hspace = 0.3, bottom=0.1, left=0.1,right=0.95, top=0.9)
        plt.style.use('seaborn')

        ax.set_xlabel('MD time (fs)');
        ax.set_ylabel('MSD ($\AA^2$)')
        x=np.linspace(0,self.Niter,self.Niter)*self.potim
        ax.plot(x, self.msd, linewidth=3)
        ax.legend(loc='best', shadow=True,fontsize='medium',frameon=False)
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.savefig('msd.pdf',dpi=600)
        plt.show()

#--------------------------------------------------------------------------------
