/*
 *  This generates a list with Gauss quadrature points and weights
 *  usage
 *    tPierrePQ -t=TYPE -p=ORDER
 *  where
 *     TYPE  type of Gauss quadrature
 *           0 = Gauss Legendre
 *           1 = Gauss Laguerre
 *           2 = Gauss Hermite
 */
#include <stdio.h>
#include <stdlib.h>
#define SqrtPI 1.77245385090551603



/* Gauss-quadrature points */
void Jacobi(int, double, double, double*, double*); 
void Hermite(int, double, double*, double*);
void Laguerre(int, double, double*, double*);


int main(int nargs, char *argv[]) {
  int i, p=5, type=0;
  double *x, *w; 

  /* parse the command line */
  for ( i=1; i<nargs; i++ )
    if ( argv[i][0] == '-' )
      switch ( argv[i][1] ) {
     case 't': type = atoi( argv[i]+3 );
        break;
      case 'p': p = atoi( argv[i]+3 );
        break;
      }

  x = (double*)malloc( p*sizeof(double) );
  w = (double*)malloc( p*sizeof(double) );

  if ( type==0 ) {
    printf("\nGauss Legendre (for [0,1]), order=%d\n", p);
    Jacobi(p, 0.0, 0.0, x, w); 
  }
  else if ( type==1 ) {
    printf("\nGauss Laguerre (Table 25.9), order=%d\n", p);
    Laguerre(p, 0.0, x, w);
    /* for ( i=0; i<p; i++ ) {
      w[i] *= exp(x[i]);
    } */
  }
  else if ( type==2 ) {
    printf("\nGauss Hermite (Table 25.10), order=%d\n", p);
    Hermite(p, 0.0, x, w);
    for ( i=0; i<p; i++ ) {
      w[i] *= SqrtPI;
    }
    /* for ( i=0; i<p; i++ ) {
      w[i] *= exp(x[i]*x[i]);
      } */
  }
  else {
    printf("\ntype=%d is not supported\n", type);
    exit(1);
  }

  for ( i=0; i<p; i++ ) {
    printf("%18.15lf %18.15lg\n", x[i], w[i] );
  }
  
}
