#--------------------------------------------------------------------------------
# diffcoeff.py, Dec 2021, Yunguo Li @ USTC
# This script can calculate diffusion coefficient from Einstein formula 
# VASP XDATCAR and OUTCAR needed as input files. 
#--------------------------------------------------------------------------------

import os, sys
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------------
class diffcoeff:

    def __init__(self, xdatcar="",outcar=""):
        if not os.path.isfile(xdatcar):
            self.xdatcar = 'XDATCAR'
        else:
            self.xdatcar = xdatcar
        if not os.path.isfile(outcar):
            self.outcar = 'OUTCAR'
        else:
            self.outcar = outcar

        # time step of MD
        self.potim = None
        self.nblock = None

        self.readoutcar()

        self.TypeName = None   # list of element name per type
        self.ChemSymb = None   # list of element name per ion
        self.Ntype = None      # no. of types of elements
        self.Nions = None      # total no. of ions  
        self.Nelem = None      # no. of ions for each element
        self.Niter = None      # no. of interations

        # position in Direct Coordinate
        self.position = None
        # position in Cartesian Coordinate
        self.positionC = None
        # Velocity in Angstrom per Femtosecond
        self.readxdat()

        self.mass_and_name_per_ion()
        # Time in femtosecond
        self.Time = self.Niter * self.potim * self.nblock

        self.segments=1
        self.skip=0          # initial MD steps to skip  
        self.skip2=0         # MD steps to skip at the beginning of each segment
        self.msd = None
        self.diffcoeff = None
        self.intercept = None
        self.ddiffcoeff = None

    def mass_and_name_per_ion(self):
        self.ChemSymb = []

        if self.TypeName is None:
            self.TypeName = [chr(i) for i in range(65,91)][:self.Ntype]

        for i in range(self.Ntype):
            self.ChemSymb += [np.tile(self.TypeName[i], self.Nelem[i])]

        self.ChemSymb = np.concatenate(self.ChemSymb)

    def readxdat(self):
        print("Reading XDATCAR ......")

        inp = [line for line in open(self.xdatcar) if line.strip()]
        scale = float(inp[1])
        self.cell = np.array([line.split() for line in inp[2:5]],
                              dtype=float)
        self.cell *= scale
        ta = inp[5].split()
        tb = inp[6].split()
        if ta[0].isalpha():
            self.TypeName = ta
            self.Ntype = len(ta)
            self.Nelem = np.array(tb, dtype=int)
            self.Nions = self.Nelem.sum()
        else:
            print("VASP 4.X Format encountered...")
            self.Nelem = np.array(tb, type=int)
            self.Nions = self.Nelem.sum()
            self.Ntype = len(tb)
            self.TypeName = None

        pos = np.array([line.split() for line in inp[7:]
                        if not line.split()[0].isalpha()],
                        dtype=float)
        self.position = pos.ravel().reshape((-1,self.Nions,3))
        self.Niter = self.position.shape[0]


    def readoutcar(self):
        if os.path.isfile(self.outcar):
            print("Reading OUTCAR ......")
            outcar = [line.strip() for line in open(self.outcar)]
            lp = 0; lb = 0; lm = 0; 
            for ll, line in enumerate(outcar):
                if 'POTIM' in line:
                    lp = ll
                if 'NBLOCK' in line:    # XDATCAR writting interval
                    lb=ll
                    self.nblock = float(outcar[ll].split()[2].replace(';',''))
                if 'Mass of Ions in am' in line:
                    lm = ll + 1
                if lp and lb and lm:
                    break

            # print outcar[lp].split(), lp, lm
            self.potim = float(outcar[lp].split()[2])
            self.mtype = np.array(outcar[lm].split()[2:], dtype=float)
        else:
            print("Reading OUTCAR failed ...")
            self.potim = input("POTIM = ")
            self.nblock = float(input("NBLOCK = "))
            self.mtype = np.array(input("Mass of Ions in am = ").split(), 
                                  dtype=float)

    def getmsd(self):
        # length per segment
        seglength=int(np.floor((self.Niter-self.skip)/self.segments))
        print("Number of segments = ",self.segments)
        print("Segment length = ",seglength)
        # msd for each segment, each type of element, each step, and, x y, z
        self.msd=np.zeros(shape=(self.segments, self.Ntype, seglength, 3))
        for i in range(self.segments):
            # displacement [step,atom,xyz]
            displacement = np.zeros_like(self.position[self.skip 
                           + seglength*i:self.skip+seglength*(i+1),:])
            displacement[0,:,:] = 0
            displacement[1:,:,:] = np.diff(self.position[self.skip 
                           + seglength*i:self.skip+seglength*(i+1),:], axis=0)
            
            # wrap back into box
            displacement[displacement > 0.5] -= 1.0
            displacement[displacement < -0.5] += 1.0

            # vector of accumulated displacement
            displacement = np.cumsum(displacement, axis=0)
            
            displacementC = np.zeros_like(displacement)
            displacementC = np.dot(displacement[:,:,:], self.cell)
            # squared displacement in x, y, z
            displacementC[:,:,:] =  displacementC[:,:,:]**2

            for j in range(self.Ntype):
                if j==0:
                    labelst=0
                else:
                    labelst = np.sum(self.Nelem[:j])
                labeled = np.sum(self.Nelem[:j+1])
                self.msd[i,j,:,:] = np.average(
                                    displacementC[:,labelst:labeled,:], axis=1)


    def getdiff(self):
        self.diffcoeff = np.zeros(shape=(self.segments, self.Ntype, 3))
        self.intercept = np.zeros(shape=(self.segments, self.Ntype, 3))
        self.ddiffcoeff = np.zeros(shape=(self.Ntype, 3))
        for i in range(self.segments):    # loop over segments
            for j in range(self.Ntype):    # loop over type of elements
                for k in range(3):    # loop over x.y,z
                    # time unit in fs, distance unit in A, unit in m^2/s
                    if np.floor((self.Niter-self.skip)/self.segments) < 2*self.skip2 :
                        self.diffcoeff[i,j,k], self.intercept[i,j,k]=np.polyfit(
                        self.potim*self.nblock
                        *np.arange(np.floor((self.Niter-self.skip)/self.segments)),
                        1E-5*self.msd[i,j,:,k], deg=1)
                    else:
                        self.diffcoeff[i,j,k], self.intercept[i,j,k] = np.polyfit(
                        self.potim*self.nblock
                        *np.arange(np.floor((self.Niter-self.skip)/self.segments)
                        - self.skip2), 1E-5*self.msd[i,j,self.skip2:,k], deg=1 )
        # standard deviation of D over segments for x,y,z directions
        self.diffcoeff /= 2 ; 
        self.ddiffcoeff = np.std(self.diffcoeff, axis=0) / 2

    def prt(self):
        for i in range(self.Ntype):
            print('\nD of',self.TypeName[i],'in direction X =',
                  np.average(self.diffcoeff[:,i,0]),'+/-', self.ddiffcoeff[i,0],
                  'm^2/s' )
            print('D of',self.TypeName[i],'in direction Y =',
                  np.average(self.diffcoeff[:,i,1]),'+/-', self.ddiffcoeff[i,1],
                  'm^2/s' )
            print('D of',self.TypeName[i],'in direction Z =', 
                  np.average(self.diffcoeff[:,i,2]),'+/-', self.ddiffcoeff[i,2],
                  'm^2/s' )
            print('Averaged D of',self.TypeName[i],'=', 
                  np.average(self.diffcoeff[:,i,:]), 
                  '+/-', np.average(self.ddiffcoeff[i,:]), 'm^2/s\n' )

#--------------------------------------------------------------------------------
