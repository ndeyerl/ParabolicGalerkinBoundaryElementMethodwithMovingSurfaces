#include <stdio.h>
#include "Jacobi.h"

#ifndef N
#define N 4
#endif

main()
{
    int k;
    double x[N], w[N];
    Jacobi(N, 0.0, 0.0, x, w);
    printf("%d\n", N);
    for (k=0; k<N; k++) printf("%16.9e  %16.9e\n", x[k], w[k]);
}
