Necessary libraries: -llapack (quadrature library also necessary, packed as libgaussq.a and doesn't need to be repacked)

Functionality for executable:
make quadjac
quadjac -N=integer valued spacial mesh (5=very coarse, 160=very fine) -M=integer valued temporal mesh (5=very coarse, 160=very fine), -p=integer valued number of quadrature nodes (5=average amount, 10=very fine)
