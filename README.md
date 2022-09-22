Installing
============

MD2D needs the dependencies:


    numpy
    scipy
    matplotlib
    mdanaysis

It can be installed locally. Under the package directory


    pip install .

It can also be installed from pypi


    pip install md2d

Usage
=====

MD2D can be used as a standalone code


    md2d

or imported as a module


    >>> import md2d
    >>> md2d.diffusion(xdatcar="XDATACR",outcar="OUTCAR")

Files
=====

Input:
------
MD2D reads the output of molecular dynamics simulations in the VASP format. Only two input files will be read: XDATCAR and OUTCAR.

..

    POTIM  = 1.0000    time-step for ionic-motion
    NBLOCK =      1;   KBLOCK =      1    inner block; outer block
    Mass of Ions in am
       POMASS =  55.85 
    volume of cell :     711.71
    kin. lattice  EKIN_LAT=         93.098925  (temperature 7275.16 K)
    Total+kin.     3416.08  3605.40  3585.13   -42.72  -109.39   -38.68
    
..

The keys words “POTIM”, “NBLOCK”, “Mass of Ions in am”, “number of ions“, “volume of cell”, “EKIN_LAT=”, “Total+kin.” are used for searching.

Besides, MD2D uses MDAnalysis to transform the outputs from other molecular dynamics codes to the VASP format. After installation, use the “md2vasp” command to initiate the conversion, or call the “md2v” function.

<img width="340" alt="image" src="https://user-images.githubusercontent.com/37633246/191765094-89fc3353-7b96-4fb7-9c56-2e07f9a3ff6d.png">
<img width="340" alt="image" src="https://user-images.githubusercontent.com/37633246/191765133-9ca6460d-e249-4281-9798-033074058703.png">

Output:
-------
Calculated diffusion coefficients, viscosities, and correction to diffusion coefficient are streamed in standard output. Besides, D-tau, mean squared displacement (MSD) and stress auto-correlation function (SACF) will be written into file as X_D_tau.dat (X is the element name), msd.dat and sacf_average.dat, and their plotsare saved as D_tau.pdf, msd.pdf and SACF.pdf. In each X_D_tau.dat, the first column is the time interval, the left three columns are the Ds for x, y and z directions. The columns in msd.dat corresponds to the MSDs for the elements in the order provided in XDATCAR. The three columns in sacf_average.dat are the MD step, bulk SACF and shear SACF.

How-To Guides
=============

Lauch the code
--------------
Simply import the package within python, and enter the menu by typing:

    import md2d
    md2d.diffusion(XDATCAR", outcar="OUTCAR")
It takes two variables: the names of input files. Default values will be XDATCAR and OUTCAR.

<img width="340" alt="image" src="https://user-images.githubusercontent.com/37633246/191766406-a008d66b-d6ed-4ec0-a294-7cfb03ae51fb.png">

Calculate diffusion coefficient
-------------------------------
Choose option 1 will lead to the calculation of diffusion coefficient. Then, MD2D will exhaust the available time intervals to calculate the diffusion coefficient as a function of time interval. The procedure and results will be streamed on the screen, and the D-tau data will be written into the file X_D_tau.dat. 

Two integer parameters will be asked for, namely, the number of initial MD steps to drop and the number of MD steps to skip for each segment. The first is to drop out the unwanted/unequilibrated initial stage of MD. The second is to exclude the ballistic motion stage from fitting the Einstein relation. The second parameter can be determined from the log-log plot or MSD.

<img width="340" alt="image" src="https://user-images.githubusercontent.com/37633246/191766572-43503b95-89bf-4e4f-9bf4-56a6ca5e47d5.png">

<img width="340" alt="image" src="https://user-images.githubusercontent.com/37633246/191766597-5cc56c2b-6a44-400a-be07-368780b51660.png">

Plot D-tau
----------
If diffusion coefficients have been calculated by doing option 1, then choose the option 4 will plot the D-tau curve.

<img width="340" alt="image" src="https://user-images.githubusercontent.com/37633246/191766678-b9bfc8e0-4e52-4226-af0e-432599cc35d3.png">

If step one hasn’t been done, they a warning message will prompt to do option 1 firstly.

<img width="340" alt="image" src="https://user-images.githubusercontent.com/37633246/191766713-15f1a451-fd7e-4370-ab42-198160b63ae6.png">

Calculate MSD
-------------
The MSD can be calculated and plotted by choosing option 3. The plot usually sees a change of trending, and the length of the first part is the length of ballistic motion that needs to be removed from fitting the Einstein relation.

<img width="340" alt="image" src="https://user-images.githubusercontent.com/37633246/191766798-2a3fc60a-6258-40b1-ba3e-e6933151b4bf.png">

Calculate viscosities and correction to diffusion coefficient
-------------------------------------------------------------
The stress auto-correlation function, viscosities, and correction to diffusion coefficient will be calculated by choosing option 2. The results will be streamed on the screen. 
The correction to diffusion coefficient is best applied to single-component system, e.g., pure Fe liquid or molecular liquid. For complex system like silicate melt, an element- or component-wise correction according to Eq. (4) in Li and Ni, (2022) should be considered.

<img width="340" alt="image" src="https://user-images.githubusercontent.com/37633246/191766846-fd20c376-334c-411e-95d7-c998161580ed.png">

Plot stress auto-correlation function
-------------------------------------
The stress auto-correlation function can be plotted by choosing option 5 after finishing the calculation in option 2. Otherwise, MD2D will ask for doing option 2 firstly.

<img width="340" alt="image" src="https://user-images.githubusercontent.com/37633246/191766944-6937d236-136b-40b9-ac2c-c4a3d5d067c5.png">

<img width="340" alt="image" src="https://user-images.githubusercontent.com/37633246/191766974-9ca7c48f-54c6-45ec-83a5-eda3a9e05224.png">

License & disclaimer
====================
MD2D Copyright (C) University of Science and Technology of China. This file is part of MD2D. 
MD2D is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation version 3 of the License. 
MD2D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 
You should have received a copy of the GNU General Public License along with MD2D. If not, see <http: //www.gnu.org/licenses/>. 
