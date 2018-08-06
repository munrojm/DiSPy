=====
DiSPy
=====

DiSPy is a utility for applying the distortion symmetry method (DSM) to the calculation of minimum energy pathways using the nudged elastic band (NEB) algorithm. For information on how to run DiSPy, including all relevant necessary input parameters, see the PDF manual in the "docs" folder. Details on the DSM can be found in Ref. [1-2]. Example calculations can be found in the "examples" folder.

Prerequisites and Compilation
=============================

DiSPy is able to run with Python2 or Python3 and requires *numpy*, *spglib* and the *atomic simulation environment (ASE)*. It can be manully installed by typing:

``
python setup.py install
``

Note that a fortran compiler is required. Testing has been done with *gfortran* only. 


References
==========

[1] J.M. Munro et. al. *Discovering minimum energy pathways via distortion symmetry groups*. Phys. Rev. B (2018).

[2] B.K. VanLeeuwen & V. Gopalan. *The antisymmetry of distortions*. Nat. Commun. (2015).





