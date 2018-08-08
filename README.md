# MDToolbox

Note: Now being ported into Julia https://github.com/ymatsunaga/MDToolbox.jl

MDToolbox is a MATLAB/Octave toolbox for statistical analysis of
molecular dynamics (MD) simulation data of biomolecules. It consists
of a collection of functions covering the following types of
scientific computations: 

* I/O for trajectory, coordinate, and topology files used for MD simulation
* Least-squares fitting of structures
* Potential mean force (PMF) or free energy profile from scattered data
* Statistical estimates (WHAM and MBAR methods) from biased data
* Dimensional reductions (Principal Component Analysis, and others)
* Elastic network models (Gaussian and Anisotropic network models)
* Unsupervised learning for refining Markov State Models (Baum-Welch algorithm and others)
* Utility functions, such as atom selections

For more information, see [the documentation](http://mdtoolbox.readthedocs.org/).

## Install

See [the documentation](http://mdtoolbox.readthedocs.io/en/latest/introduction.html#installation-for-matlab).

## Docker

Docker image for MDToolbox is available.
The following single line command starts Octave ready for use with MDToolbox.

```sh
$ docker run -it --rm -v $(pwd):/home/jovyan/work ymatsunaga/octave octave
```

For details, see [the documentation](http://mdtoolbox.readthedocs.io/en/latest/introduction.html#docker-image-for-mdtoolbox).

## Reference

In preparation.

## Files used in our scientific papers

["Minimum Free Energy Path of Ligand-Induced Transition in Adenylate Kinase" PLoS Comput. Biol. (2012)](https://doi.org/10.1371/journal.pcbi.1002555)
* assignvoronoi.m - assigns Voronoi index to trajectory data
* mbar.m - multistate Bennett acceptance ratio method
* calcpca.m  - Principal Component Analysis (PCA)

["Influence of Structural Symmetry on Protein Dynamics" PLoS One (2012)](https://doi.org/10.1371/journal.pone.0050011)
* anmsym.m - anisotropic network model for proteins with circular symmetries
* calcpca.m  - Principal Component Analysis (PCA)

["Sequential data assimilation for single-molecule FRET photon-counting data" J. Chem. Phys. (2015)](https://doi.org/10.1063/1.4921983)
* calcfret.m - generates a FRET photon-count sequence from distance time-series data
* calccontactvector.m - calculates contact map vectors from trajectory data
* calcpca.m  - Principal Component Analysis (PCA)

["Dimensionality of Collective Variables for Describing Conformational Changes of a Multi-Domain Protein" J. Phys. Chem. Lett. (2016)](https://doi.org/10.1021/acs.jpclett.6b00317)
* calcpathcv.m - calculates path collective variable (CV)
* assigntransitionpath.m - detects transition paths from trajectory data

["Energetics and conformational pathways of functional rotation in the multidrug transporter AcrB" eLife (2018)](https://doi.org/10.7554/eLife.31715)
* mbar.m - multistate Bennett acceptance ratio method
* calcpathcv.m - calculates path collective variable (CV)
* superimpose2d.m - superimposes structures by using only translation in XY-space and rotation around Z-axis
* calcmutinf.m - estimates mutual information
* calcgse.m - calculates electrostatic potential from trajectories by using k-space Gaussian split Ewald

"Linking time-series of single-molecule experiments with molecular dynamics simulations by machine learning" submitted
* msmbaumwelchdb.m msmforward.m msmbackward.m msmtransitionmatrix.m - Baum-Welch algorithm with a constraint imposed by the detailed-balance condition
* msmbaumwelchdb_parallel.m msmforward_parallel.m msmbackward_parallel.m - Parallelized version of the above (requires Parallel Computing Toolbox)
* msmgenerate.m - Stochastic simulation according to the constructed Markov State Model (MSM)
* calcqscore.m - calculates Q-score for atomistic model
* calcorientationfactor.m - calculates Orientation factor for FRET dyes
* example/msmbaumwelchdb - test data set for msmbaumwelchdb.m and msmbaumwelchdb_parallel.m explaining the usage of these functions

"Refining Markov State Models for conformational dynamics using ensemble-averaged data and time-series trajectories" submitted
* msmbaumwelchdb.m  msmforward.m msmbackward.m msmtransitionmatrix.m - Baum-Welch algorithm with a constraint imposed by the detailed-balance condition
* msmcountmatrix.m - calculates count matrix from indexed trajectory data
* msmtransitionmatrix.m - reverse Maximum-likelihood estimator for transition matrix from counting matrix
* msmensemble.m - estimates populations from distribution data
* calcdistancevector.m - calculate distance-matrix-based vectors from trajectory data
* calctica.m - time-structure based ICA (tICA)
* example/msmbaumwelchdb - test data set for msmbaumwelchdb.m and msmbaumwelchdb_parallel.m explaining the usage of these functions

## Developer

Yasuhiro Matsunaga ymatsunaga@riken.jp

