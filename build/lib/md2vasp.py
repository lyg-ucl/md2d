import os
import numpy as np
import MDAnalysis as mda

#--------------------------------------------------------------------------------
class converter:

    def __init__(self):

        self.xdatcar  = 'XDATCAR_converted'
        self.outcar   = 'OUTCAR_converted'
        self.potim    = None      # MD timestep
        self.TypeName = None      # list of element name per type
        self.Ntype    = None      # no. of types of elements
        self.Nions    = None      # total no. of ions
        self.mass     = None      # mass for each type
        self.Nelem    = None      # no. of ions for each element
        self.Niter    = None      # no. of interations
        self.volume   = None
        self.cell     = np.zeros((3,3))
 
        u = self.load()

        self.TypeName, self.Nelem = np.unique(u.atoms.types, return_counts=True)
        self.Nions = u.atoms.n_atoms
        self.potim = u.trajectory.dt
        self.Ntype = len(self.TypeName)

        self.dumpxdatcar(u)
        self.dumpoutcar(u)


    def load(self):
        flag = True
        line = []
        while flag:
            line = input('''
\nPlease provide MD files to be loaded to MDAnalysis. 
Example: ifabp_water.psf ifabp_water_0.pdb rmsfit_ifabp_water_1.dcd
For more info about MDAnalysis, please visit  www.mdanalysis.org
\nInput your file names:
''')
            line = line.split()
            for i in line:
                if not os.path.isfile(i):
                    print('File %s does not exist' % i)
                    flag = False
            if flag == False :
                flag = True
            else:
                flag = False
                print('\nFiles will be converted to %s and %s\n' % (self.xdatcar,self.outcar))
        u=mda.Universe(*line)
        return u

    def dumpxdatcar(self,u):
        ouf = open(self.xdatcar, "w")
        if u.dimensions is None :
            print('>>> No unitcell info provided <<<')
            self.cell[0,0]=max(u.atoms.positions[:,0])-min(u.atoms.positions[:,0])
            self.cell[1,1]=max(u.atoms.positions[:,0])-min(u.atoms.positions[:,1])
            self.cell[2,2]=max(u.atoms.positions[:,0])-min(u.atoms.positions[:,2])
            print('Unitcell vectors will be: \n(in Angstrom)')
            print('a-vector: %f  %f  %f' % tuple(self.cell[0,:]))
            print('b-vector: %f  %f  %f' % tuple(self.cell[1,:]))
            print('c-vector: %f  %f  %f' % tuple(self.cell[2,:]))
            flag=input('Agreed? (y/n)')
            if flag == 'n' or flag == 'no' or flag == 'No' or flag == 'NO':
                line = input('Set new a-vector:'); self.cell[0,:] = np.array(line.split())
                line = input('Set new b-vector:'); self.cell[1,:] = np.array(line.split())
                line = input('Set new c-vector:'); self.cell[2,:] = np.array(line.split())
        else:
            self.cell[0,0] = u.dimensions[0];
            self.cell[1,0] = u.dimensions[1]*np.cos(u.dimensions[5]/180*np.pi)
            self.cell[1,1] = u.dimensions[1]*np.sin(u.dimensions[5]/180*np.pi)
            self.cell[2,0] = u.dimensions[2]*np.cos(u.dimensions[4]/180*np.pi)
            factor = (np.cos(u.dimensions[3]/180*np.pi)-np.cos(u.dimensions[4]/180*np.pi)*np.cos(u.dimensions[5]/180*np.pi))/ \
                     np.sin(u.dimensions[5]/180*np.pi)
            self.cell[2,1] = u.dimensions[2]*factor
            self.cell[2,2] = u.dimensions[2]*np.sqrt(1-np.cos(u.dimensions[4]/180*np.pi)**2-factor**2)
        self.volume = np.dot(self.cell[0,:], np.cross(self.cell[1,:],self.cell[2,:]))

        ouf.write('converted\n')
        ouf.write('1.0\n')
        ouf.write('%16.12f  %16.12f  %16.12f \n' % tuple(self.cell[0,:]))
        ouf.write('%16.12f  %16.12f  %16.12f \n' % tuple(self.cell[1,:]))
        ouf.write('%16.12f  %16.12f  %16.12f \n' % tuple(self.cell[2,:]))
        ouf.write((' '.join(['%s ']*self.TypeName.size)+'\n') % tuple(self.TypeName))
        ouf.write((' '.join(['%s ']*self.Nelem.size)+'\n') % tuple(self.Nelem))

        self.Niter = 0
        self.mass = []
        for ts in u.trajectory:
            self.Niter += 1
            ouf.write('Direct configuration= '+str(self.Niter)+'\n')
            for i in range(self.Ntype) :
                select= 'type ' + self.TypeName[i]
                atm = u.select_atoms(select)
                if self.Niter == 1: self.mass.append(atm.masses[0])
                for j in range(self.Nelem[i]):
                    ouf.write('%12.8f  %12.8f  %12.8f \n' % tuple(np.dot(atm.positions[j,:],np.linalg.inv(self.cell))))
        ouf.close()

    def dumpoutcar(self, u):
        ouf = open(self.outcar, "w")
        print('Timestep will be %f ps' % self.potim)
        flag=str(input('Agreed? (y/n)'))
        if flag == 'n' or flag == 'no' or flag == 'No' or flag == 'NO':
            self.potim = float(input('New timestep in ps :'))
            print('Timestep set as %f ps' % self.potim)
        self.potim = self.potim*1000

        ouf.write('  volume of cell : %f\n' % self.volume)
        ouf.write('   number of dos      NEDOS =    0   number of ions     NIONS =   %d \n' % u.atoms.n_atoms )
        ouf.write('   POTIM  = %f    time-step for ionic-motion\n' % self.potim)
        ouf.write('   NBLOCK =      1;   KBLOCK =    1    inner block; outer block\n')
        self.mass = np.array(self.mass)
        ouf.write('  Mass of Ions in am \n   POMASS = ')
        ouf.write(' '.join(['%6.3f ']*len(self.mass)) % tuple(self.mass))
        ouf.close()

#--------------------------------------------------------------------------------
