Implementation of the parabolic Galerkin Boundary Element Method with moving surfaces, in support of our (my and Johannes Tausch's) paper, Quadrature for Parabolic Galerkin Boundary Element Methods with Moving Surfaces.  This code models heat distribution in time and 3D space using C programming, and works for novel 3D shapes such a cubes and polygons.


Necessary libraries: -llapack (quadrature library also necessary, packed as libgaussq.a and doesn't need to be repacked)

Functionality for executable:
make quadjac
quadjac -N=integer valued spacial mesh (5=very coarse, 160=very fine) -M=integer valued temporal mesh (5=very coarse, 160=very fine), -p=integer valued number of quadrature nodes (5=average amount, 10=very fine)
