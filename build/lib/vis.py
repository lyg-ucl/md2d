#--------------------------------------------------------------------------------
# viscosity.py, Dec 2021, Yunguo Li @ USTC
# This script can calculate stress autocorrelation function (SACF) and 
# bulk/shear viscosity
# VASP OUTCAR needed as input file. 
#--------------------------------------------------------------------------------

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftfreq, fftshift, rfft


# Boltzmann Constant in [eV/K]
kb = 8.617332478E-5
# electron volt in [Joule]
ev = 1.60217733E-19
# ev to A^3*GPa
ev2pv = 160.21766208

#--------------------------------------------------------------------------------
class viscosity:

    def __init__(self,File=None):
        if File is None:
            self.outcar = 'OUTCAR'
        else:
            self.outcar = File

        self.potim    = None          # time step of MD
        self.nblock   = None          # trajectory writting interval
        self.Nions    = None          # total no. of ions
        self.volume   = None          # volume of system 
        self.Niter    = 0             # no. of interations
        self.stress   = np.zeros(6)   # stresses each step; 
                                      # direction:XX,YY,ZZ,XY,YZ,ZX
        self.Press    = 0             # average pressure
        self.Temp     = None          # temperature each step
        self.TypeName = ['xx','yy','zz','xy','yz','zx'] 
        # Velocity Autocorrelation Function 
        self.ACF = None
        self.ACF2= None

        self.readoutcar()


    def readoutcar(self):
        print("Reading stresses from OUTCAR ......")
        if os.path.isfile(self.outcar):
            outcar = [line.strip() for line in open(self.outcar)]
            lm=0; lb=0; li=0; lv=0; lk=0; lp=0
            for ll, line in enumerate(outcar):
                if 'POTIM' in line:    # read POTIM 
                    lm=ll
                    self.potim = float(outcar[ll].split()[2]) 
                if 'NBLOCK' in line:    # XDATCAR writting interval
                    lb=ll
                    self.nblock = float(outcar[ll].split()[2].replace(';',''))
                if 'number of ions ' in line:    # read number of ions
                    li=ll
                    self.Nions = float(outcar[ll].split()[11])
                if 'volume of cell' in line:    # read volume
                    lv=ll
                    self.volume = np.array(outcar[ll].split()[4], dtype=float)
                if 'EKIN_LAT=' in line:    # read T
                    self.Temp = np.array(outcar[ll].split()[5], dtype=float)
                    lk=ll
                    self.Niter += 1
                if 'Total+kin.' in line:    # read stresses
                    self.stress = np.array(outcar[ll].split()[1:], dtype=float)
                    lp=ll
                if lm and lb and li and lv and lk and lp:
                    break

            for ll, line in enumerate(outcar):
                if 'EKIN_LAT=' in line:    # read T
                    self.Temp = np.append(self.Temp, 
                                np.array(outcar[ll].split()[5], dtype=float))
                if 'Total+kin.' in line:    # read stresses
                    self.stress = np.append(self.stress, 
                                  np.array(outcar[ll].split()[1:], dtype=float) )
                    self.Niter += 1
        else:
            print("ERROR reading OUTCAR !!!") 

        self.stress = np.reshape(self.stress, (self.Niter, 6)) / 10    # in GPa
        self.Press = np.average(np.average(self.stress[:,0:3]))      # average P

        return


    def getACF(self, method=None):

        #---- Cross-correlation ----
        # in this function, numpy.correlate is used to calculate the ACF
        if method == 'correlate':
            print("Calculating stress auto-correlation function ......")
            self.stress[:,0:3] = self.stress[:,0:3] - self.Press
            # np.correlate(a[:n],a[:m],'full') returns the number of 
            # n+m-1 elements
            self.ACF2 = np.zeros((6, self.Niter*2 - 1))   
            # we have 5 independent shear stresses
            # self.stress in direction:XX,YY,ZZ,XY,YZ,ZX
            self.ACF2[0,:] = np.correlate((self.stress[:,0]+self.stress[:,1]
                             + self.stress[:,2])/3, (self.stress[:,0] 
                             + self.stress[:,1]+self.stress[:,2])/3,'full')
            self.ACF2[1,:] = self.ACF2[0,:]  ; self.ACF2[2,:] = self.ACF2[0,:]
            for i in range(3,6):
                self.ACF2[i,:] = np.correlate(self.stress[:,i],
                                        self.stress[:,i],'full')

            # average: np.correlate(a,a,'full') returns 2*n-1 elements, 
            # which in total performs v(t1+t2)*v(t2) the times of n**2
            self.ACF2[:,:] /=self.Niter
            self.ACF = self.ACF2[:,self.Niter-1:]

        #---- alternatively, in explicite ----
        else:
            print("SACF is calculated straightforwardly ...")
            self.ACF2 = np.zeros((6, self.Niter*2 - 1))
            self.ACF = np.zeros((6, self.Niter))
            for i in range(self.Niter):
                self.ACF[0,i] = self.stress[i,0]*self.stress[0,0]    # XY
                self.ACF[1,i] = self.stress[i,1]*self.stress[0,1]    # YZ
                self.ACF[2,i] = self.stress[i,2]*self.stress[0,2]    # ZX
                self.ACF[3,i] = (np.average(self.stress[i,3:6]) - self.Press) \
                                * (np.average(self.stress[0,3:6])- self.Press)/9
                self.ACF[4,i] = self.ACF[3,i]
                self.ACF[5,i] = self.ACF[3,i]
            # unit GPa^2
            self.ACF2[:,:self.Niter] = self.ACF[:,:]

        return self.ACF
 

    def visco(self, method='fft',acf_method='correlate'):
        print("Calculating viscosity from stress auto-correlation "
              + "function ......")

        # Frequency in THz: omega[0]=0, omega[self.Niter-1]= 
        # fft tranformation has the frequency of Fn=(n-1)*Fs/N for N point 
        # input data, and n<=N/2
        omega = fftfreq( self.Niter, self.potim*self.nblock ) * 1E3
        #print(omega[0],omega[self.Niter-1],omega[self.Niter])

        #---- fft ----
        if method == 'fft':
            if self.ACF is None:
                self.getACF(acf_method)      # (method='correlate') or ()
            vis = np.zeros_like(self.ACF)    # pdos shape
            # viscosity in unit of Pa*S; here self.ACF should be used other than
            #  VSF2, otherwise,
            # vis should be normalised with a factor of 
            # len(self.ACF)/len(self.ACF2)
            vis[:,:] = np.abs(fft(self.ACF[:,:])) \
                       * self.volume/kb/np.average(self.Temp)/ev2pv *1E-6
        #---- integrate ----
        elif method == 'integrate':
            if self.ACF is None:
                self.getACF(acf_method)      # (method='correlate') or ()
            from scipy.integrate import cumtrapz
            vis = cumtrapz(self.ACF, self.potim*self.nblock \
                          * np.arange(1,self.Niter+1)) * \
                           self.volume/kb/np.average(self.Temp)/ev2pv *1E-6 
            # note: integration cumtrapz returns one date less than input.
            return np.arange(2,self.Niter+1), vis

        else:
            print("Chosen spectrum analysis method not found ...")
            sys.exit()
 
        return omega[:self.Niter//2], vis[:,:self.Niter//2]


#-------------------------------------------------------------------------------
