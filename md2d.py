import os, sys
import numpy as np
import matplotlib.pyplot as plt
import einstein, vis, plot, md2vasp

# Boltzmann Constant in [eV/K]
kb = 8.617332478E-5
# ev to A^3*GPa
ev2pv = 160.21766208

#--------------------------------------------------------------------------------

def md2dexe(xdatcar="",outcar=""):
    xdatcar=input("XDATCAR = ").strip()
    outcar=input("OUTCAR = ").strip()
    diffusion(xdatcar,outcar)

def md2v():
    md2vasp.converter()

class diffusion:

    def __init__(self,xdatcar="",outcar=""):

        self.xdatcar   = xdatcar.strip()
        self.outcar    = outcar.strip()
        self.checkfiles()
        self.loop      = True
        self.choice    = None
        self.tau       = None
        self.D         = None
        self.D3        = None
        self.corr_D    = None
        self.acf       = None
        self.viscosity = None

        while self.loop:
            self.banner()
            self.menu()
            if self.choice==1:
                print("Calculating diffusion coefficient ......")
                try:
                    d.skip = 0; d.skip2=0
                except:
                    d=einstein.diffcoeff(self.xdatcar,self.outcar)
                self.tau       = []
                self.D         = []
                self.D3        = []
                d.skip=int(input("Initial MD steps to drop from the whole"
                                 + " MD run:"))
                d.skip2=int(input("MD steps to skip for each segment:"))

                i=int((d.Niter-d.skip)/2)
                while i >=int((d.Niter-d.skip)/200) :
                    d.segments=i
                    self.tau.append(np.floor((d.Niter-d.skip)/d.segments)
                                    * d.potim*d.nblock)
                    print('Interval (fs) = ', np.floor((d.Niter-d.skip)
                          / d.segments)*d.potim*d.nblock, '\n')
                    d.getmsd(); d.getdiff(); d.prt();
                    self.D.extend(np.average(d.diffcoeff[:,:,:],axis=(0,2)))
                    self.D3.extend(np.average(d.diffcoeff[:,:,:],axis=(0)))
                    i=int(i/2)
                for i in range(int((d.Niter-d.skip)/200),1,-1):
                    d.segments=i; 
                    self.tau.append(np.floor((d.Niter-d.skip)
                                    / d.segments)*d.potim*d.nblock)
                    print('Interval (fs) = ', np.floor((d.Niter-d.skip) 
                          / d.segments)*d.potim*d.nblock,'\n')
                    d.getmsd(); d.getdiff(); d.prt()
                    self.D.extend(np.average(d.diffcoeff[:,:,:],axis=(0,2)))
                    self.D3.extend(np.average(d.diffcoeff[:,:,:],axis=(0)))
                self.D = np.array(self.D); 
                self.tau=np.array(self.tau); 
                self.D3 = np.array(self.D3);
                self.D = self.D.reshape((len(self.tau),d.Ntype))
                self.D3 = self.D3.reshape((len(self.tau),d.Ntype,3))
                np.savetxt("D_tau.dat",np.column_stack((self.tau,self.D)),
                           comments='# MD time (fs)  D (m^2/s)')
                for i in range(d.Ntype):
                    fname = d.TypeName[i]+'_D_tau.dat'
                    np.savetxt(fname,np.c_[self.tau,self.D3[:,i,:]],
                               comments='# MD time (fs)  D (m^2/s)')

            elif self.choice==2:
                print("Calculating viscosities ......")
                v=vis.viscosity(self.outcar)
                acf = v.getACF('correlate')
                for i in range(len(acf)):
                    fname="sacf_of_" + v.TypeName[i] + ".dat"
                    ouf = open(fname, "w")
                    ouf.write("# MD_step  SACF(GPa^2)\n")
                    for j in range(len(acf[0])):
                        ouf.write(str(j) + " " + str(acf[i,j]) + "\n")
                    ouf.close()
                ouf = open("sacf_average.dat", "w")
                ouf.write("# MD_step  bulk_SACF(GPa^2) shear_SACF(GPa^2)\n")
                for j in range(len(acf[0])):
                    ouf.write(str(j) + " " + str(np.average(acf[:3,j])) + " " 
                              + str(np.average(acf[3:6,j])) + "\n")
                ouf.close()
                self.acf = np.column_stack((np.average(acf[:3,:],axis=0),
                                             np.average(acf[3:6,:],axis=0)))

                x,y = v.visco('integrate','correlate'); y = np.transpose(y)
                self.viscosity = np.array(np.column_stack((x,y)))
                for i in range(len(self.acf[:,0])):
                    if self.acf[i,1]<0.0:
                        print("\nCalculated bulk and shear viscosities are: \n"
                              + str(np.average(self.viscosity[i,1:4]))
                              + ", "+str(np.average(self.viscosity[i,4:7]))
                              + "  (Pascal*Second)")
                        dummyB = np.average(self.viscosity[i:,1:4],axis=1); 
                        dummyS = np.average(self.viscosity[i:,4:7],axis=1);
                        dummydB = np.std(dummyB, ddof=1) \
                                / np.sqrt(np.size(dummyB)) * len(self.acf[:,0])/i
                        dummydS = np.std(dummyS, ddof=1) \
                                  / np.sqrt(np.size(dummyS))*len(self.acf[:,0])/i
                        print("\nTheir uncertainties are: \n"+str(dummydB)+", "
                              + str(dummydS)+"  (Pascal*Second)")
                        break
                self.corr_D = np.average(v.Temp)*kb*2.837*ev2pv \
                              / (4*np.pi*dummyS[i]*v.volume**(1/3)) *1E-11 
                print("\nCalculated correction to D is:\n"+str(self.corr_D) 
                      + " +/- "+str(self.corr_D*dummydS/dummyS[i])+"  (m^2/s)")

                ouf = open("viscosity.dat", "w")
                ouf.write("# MD_steps  bulk_ciscosity(Pascal*Second)  "
                          + "shear_viscosity\n")
                for j in range(len(y[0])):
                    ouf.write(str(x[j]) + " " + str(np.average(y[0:3,j])) + " " 
                              + str(np.average(y[3:6,j])) +"\n")
                ouf.close()

            elif self.choice==3:
                try:
                    d.segments=1; d.skip=0 ; d.skip2=0
                except:
                    d=einstein.diffcoeff(self.xdatcar,self.outcar)
                    d.segments=1; d.skip=0 ; d.skip2=0
                d.getmsd()
                np.savetxt("msd.dat",
                           np.squeeze(np.sum(d.msd,axis=3)).transpose() )
                print("Plotting msd ......")
                p2=plot.plotmsd(np.squeeze(np.sum(d.msd,axis=3)).transpose(), 
                                d.Niter,d.potim*d.nblock, d.TypeName)

            elif self.choice==4:
                if os.path.isfile("D_tau.dat") and (self.D is not None):
                    print("Plotting D—τ curve ......")
                    p3=plot.plotD(np.column_stack((self.tau, self.D)), 
                                  d.TypeName)
                else:
                    print(">>> This cannot be done. Step 1 needs to be done "
                          + "firstly. <<<")

            elif self.choice==5:
                if os.path.isfile("sacf_average.dat") and (self.acf is not None):
                    print("Plotting stress auto-correlation function ......")
                    p4=plot.plotsacf(self.acf)
                else:
                    print(">>> This cannot be done. Step 2 needs to be done "
                          + "firstly. <<<")

            else:
                sys.exit("Byebye!!!")

    def checkfiles(self):
        if not os.path.isfile(self.xdatcar):
           sys.exit("!!!XDATCAR does not exist!!!")
        if not os.path.isfile(self.outcar):
            sys.exit("!!!OUTCAR does not exist!!!")

    def menu(self):
        msg = """
-------------------------------------------------------------------------------
-->    1) Diffusion coefficient (D)
-->    2) Viscosities and corretion to D
-->    3) Mean squared displacement plot
-->    4) D-τ plot
-->    5) Stress auto-correlation function plot
-->    6) Quit
-------------------------------------------------------------------------------
    """
        print(msg)
        try:
            inp=int(input("Your option:"))
        except:
            sys.exit("Byebye!!!")
        if isinstance(inp, int) and inp<7 and inp>0:
            self.choice=inp
        else:
            self.choice=inp
            self.loop = False

    def banner(self):
        msg = """
      \\  |  __ \\ ___ \\   __ \\
     |\\/ |  |   |   ) |  |   |
     |   |  |   |  __/   |   |
    _|  _| ____/ _____| ____/    version 1.2.0 @ USTC
        """
        print(msg)
        


#--------------------------------------------------------------------------------
if __name__ == '__main__':
    xdatcar=input("XDATCAR = ").strip()
    outcar=input("OUTCAR = ").strip()
    diffusion(xdatcar,outcar)
