/* Exercizer for gaussian integration package.
/**/

#include <math.h>
#include <stdio.h>

#include "Laguerre.h"
#include "Jacobi.h"
#include "Hermite.h"

#define N 8
#define k 5
#define FLOATk 5.0
#define ZERO 0.0
#define SqrtPI 1.77245385090551603

main()
{
    short m;
    double abscis[N];
    double  wt[N];
    double left, right;


    printf("Table 25.8, k= %1d, n=%1d\n\n", k, N);
    Jacobi(N, FLOATk, ZERO, abscis, wt);
    for (m=0; m<N; m++) wt[m] /= (k+1);
    
    for (m=0; m<N; m++) printf("%2d %20.18g %20.18g\n", m, abscis[m], wt[m]);
    getchar();


    printf("\n\n\nTable 25.4, n=9\n\n");
    Radau_Jacobi(4, -0.5, ZERO, abscis, wt, &left);
    left *= 2;
    for (m=0; m<4; m++) abscis[m]= sqrt(abscis[m]);

    printf("%2d %20.18g %20.18g\n", 0, ZERO, left);
    for (m=0; m<4; m++) printf("%2d %20.18g %20.18g\n", m, abscis[m], wt[m]);
    getchar();


    printf("\n\n\nTable 25.9, n= %1d\n\n", N);
    Laguerre(N, ZERO, abscis, wt);

    for (m=0; m<N; m++) printf("%2d %20.18g %20.18g\n", m, abscis[m], wt[m]);
    getchar();


    printf("\n\n\nTable 25.10, n= %1d\n\n", N);
    Hermite(N, ZERO, abscis, wt);

    for (m=0; m<N; m++) wt[m] *= SqrtPI;
    for (m=0; m<N; m++) printf("%2d %20.18g %20.18g\n", m, abscis[m], wt[m]);
    getchar();


    printf("\n\n\nSame, n= %1d\n\n", 2*N);
    Even_Hermite(2*N, ZERO, abscis, wt);
    for (m=0; m<N; m++) wt[m] *= SqrtPI;
    for (m=0; m<N; m++) printf("%2d %20.18g %20.18g\n", m, abscis[m], wt[m]);
    getchar();


    printf("\n\n\nSame, n= %1d\n\n", 9);
    Odd_Hermite(9, ZERO, abscis, wt, &left);
    left *= SqrtPI;
    for (m=0; m<4; m++) wt[m] *= SqrtPI;

    printf("%2d %20.18g %20.18g\n", 0, ZERO, left);
    for (m=0; m<4; m++) printf("%2d %20.18g %20.18g\n", m, abscis[m], wt[m]);
    getchar();


    printf("\n\n\nTable 25.6, n= %1d\n\n", N+1);
    Lobatto_Jacobi(N-1, ZERO, ZERO, abscis, wt, &left, &right);

    for (m=0; m<(N-1); m++) wt[m] *= 2;
    left *= 2;  right *= 2;
    for (m=0; m<(N-1); m++) abscis[m] = (2*abscis[m] - 1);

    printf("   %20.18g %20.18g\n", -1.0, left);
    for (m=0; m<(N-1); m++)
	printf("%2d %20.18g %20.18g\n", m, abscis[m], wt[m]);
    printf("   %20.18g %20.18g\n",  1.0, right);
}
