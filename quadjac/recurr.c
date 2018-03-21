/* How stable is that recurrence relation?
/**/

#include <stdio.h>

main()
{
    static double m[2][2]= { {1.0, 0.0}, {0.0, 1.0} };
    int n;
    double x, alpha, beta;
    double tmp;
    extern double cond();	/* condition number in MAX norm */

    while (!feof(stdin)) {
	printf("x, alpha, beta ? ");
	scanf("%f %f %f%*[^\n]%*1[\n]", &x, &alpha, &beta);

	printf("[ %9f  %9f ]\n[ %9f  %9f ]\n\n",
	m[0][0], m[0][1], m[1][0], m[1][1]);

	for(n=1;; n++) {
	    tmp= n + alpha + beta;
	    m[1][0]/= tmp;  m[1][1]/= tmp;

	    m[1][0]-= m[0][0];  m[1][1]-= m[0][1];

	    m[0][0]+= x*m[1][0];  m[0][1]+= x*m[1][1];

	    tmp= n+alpha;  m[0][0]*= tmp;  m[0][1]*= tmp;
	    tmp= n+beta;   m[1][0]*= tmp;  m[1][1]*= tmp;

	    m[0][0]+= x*m[1][0];  m[0][1]+= x*m[1][1];

	    m[1][0]-= m[0][0];  m[1][1]-= m[0][1];

	    m[0][0]/= n;  m[0][1]/= n;

	    printf("n= %3d  kappa= %f\n", n, cond(m));
	    printf("[ %9f  %9f ]\n[ %9f  %9f ]\n",
	    m[0][0], m[0][1], m[1][0], m[1][1]);
	    if(scanf("%*[^\n]%*1[\n]") == EOF) break;
	}
    }
}

/* A condition number of a 2x2 matrix.
/**/
double cond(m)
double m[2][2];
{
    double n1, n2, delta, tmp;
    extern double fabs();

    delta= m[0][0]*m[1][1] - m[1][0]*m[0][1];

    n1=  fabs(m[0][0]) + fabs(m[0][1]);
    tmp= fabs(m[1][0]) + fabs(m[1][1]);
    if (tmp>n1)  n1= tmp;

    n2=  fabs(m[0][0]) + fabs(m[1][0]);
    tmp= fabs(m[0][1]) + fabs(m[1][1]);
    if (tmp>n2)  n2= tmp;

    if (delta==0)  return 0;
    return n1*n2/delta;
}
