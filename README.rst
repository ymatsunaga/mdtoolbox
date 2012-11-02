-------------------------------------------------------------------------------------------
MD Toolbox: A MATLAB toolbox for the analysis of molecular dynamics trajectories
-------------------------------------------------------------------------------------------

A collection of MATLAB functions for the analysis of MD trajectories. 

The functions include:

* I/O for topology, coordinates, and trajectory files used for MD simulation
* Least-squares fitting of structures
* Calculation of PMF surface from scattered data
* Anisotropic network model (Elastic network model)
* Other auxiliary functions... such as atom selections

TODOs:

* Principal component analysis
* Hydrogen bond analysis
* WHAM, and MBAR methods
* Some energy calc routines

Installing
-------------------------------------------------------------------------------------------

Download the `zipball <https://github.com/ymatsunaga/mdtoolbox/zipball/master>`_ or 
`tarball <https://github.com/ymatsunaga/mdtoolbox/tarball/master>`_ of the latest codes, 
or just clone this repository.
::

 $ git clone https://github.com/ymatsunaga/mdtoolbox.git

After downloading (and unpacking), invoke ``pathtool`` commmand in MATLAB
and add ``mdtoolbox/function`` to your MATLAB search path. 

Getting started
-------------------------------------------------------------------------------------------

See https://github.com/ymatsunaga/mdtoolbox/wiki

Summary of functions
-------------------------------------------------------------------------------------------

I/O

========================== ==================================================================================================
name                       description
========================== ==================================================================================================
readpdb                    read Protein Data Bank (PDB) file
writepdb                   write Protein Data Bank (PDB) file
readamberparm              read amber parameter/topology file
readambercrd               read amber coordinate/restart file
readamberout               read amber output file
readambertrj               read amber ascii-format trajectory file
readambertrjbox            read amber ascii-format trajectory file including box size
readnetcdf                 read amber netcdf file
writeambercrd              write amber coordinate/restart file
writeambertrj              write amber ascii-format trajectory format file
writenetcdf                write amber netcdf file
readpsf                    read charmm or xplor type Protein Structure File (PSF)
readdcd                    read xplor or charmm (namd) format dcd file
readnamdbin                read namd restart (namdbin) file
readnamdout                read namd output file
writedcd                   write xplor or charmm (namd) format dcd file
writenamdbin               write namd restart (namdbin) file
readfloattrj               read float trajectory file
readgenesisbin             read genesis restart (genesisbin) file
readgenesisout             read genesis output file
readmarblecrd              read marble coordinate/restart file
readmarbletrj              read marble ascii-format trajectory file
writemarbletrj             write marble ascii-format trajectory file
writexplormap              write xplor density format file
========================== ==================================================================================================

calculation

========================== ==================================================================================================
name                       description
========================== ==================================================================================================
calcdihedral               calculate dihedral angles from Cartesian coordinates
calcpmf                    calculate 1D potential of mean force from the scattered 1D-data (using kernel density estimator)
calcpmf2d                  calculate 2D potential of mean force from the scattered 2D-data (using kernel density estimator)
calchistpmf                calculate 1D potential of mean force from the scattered 1D-data (using histogram)
calchistpmf2d              calculate 2D potential of mean force from the scattered 2D-data (using histogram)
calcpairlist               make a pairlist by searching pairs within a cutoff distance
calcpairlist_exhaustive    make a pairlist by searching pairs within a cutoff distance
clusteringbyinformation    clusterize samples according to an information-based criterion
superimpose                least-squares fitting of structures by Kabsch's method
meanstructure              calc average structure by iterative superimposing
anm                        calculate normal modes and anisotropic fluctuations by using Anisotropic Network Model.
transformframe             transform from the Eckart frame to a non-Eckart frame
========================== ==================================================================================================

others

========================== ==================================================================================================
name                       description
========================== ==================================================================================================
formatplot                 fomart the handle properties (fonts, lines, etc.) of the current figure
exportas                   export fig, eps, png, tiff files of the current figure
decenter                   remove the center of mass from coordinates or velocities
orient                     orient molecule using the principal axes of inertia
to3                        convert 1...N atom-index to 1...3N xyz-index
substruct                  create a subset structure from a structure of arrays of same size
searchrange                finds all the atoms within cutoff distance from given atoms
searchrange_exhaustive     finds all the atoms within cutoff distance from given atoms
selectid                   used for atom selection. Finds all the atoms or residues which matches given index
selectname                 used for atom selection. Finds all the atoms or residues which matches given names
selectrange                used for atom selection. Finds all the atoms within cutoff distance from given atoms
kde2d                      fast and accurate state-of-the-art bivariate kernel density estimator
========================== ==================================================================================================

