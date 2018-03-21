#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "panelIA.h"
#include "Laguerre.h"
#include "Jacobi.h"
#include "Hermite.h"


void main(int nargs, char *argv[]) {

  int p=5, d=0;
  double ht = 1.0,  hs = 1.0;
  double *xLeg, *wLeg;
  double **vtx;
  panel *pnls;
  int i, nVtxs, nPnls;
  double solution[8], sumV, sumK;

  /* parse the command line */
  for ( i=1; i<nargs; i++ )
    if ( argv[i][0] == '-' )
      switch ( argv[i][1] ) {
      case 'p': p = atoi( argv[i]+3 );
        break;
      case 'd': d = atoi( argv[i]+3 );
        break;
     case 't': ht = atof( argv[i]+3 );
        break;
     case 'x': hs = atof( argv[i]+3 );
        break;
      }

  xLeg = (double*)malloc( p*sizeof(double) );
  wLeg = (double*)malloc( p*sizeof(double) );
  Jacobi( p, 0.0, 0.0, xLeg, wLeg);

//build vertices
  nVtxs = 3;
  vtx = (double **)calloc(nVtxs, sizeof(double *));
  for(i = 0; i<nVtxs; i++){
    vtx[i] = (double *)calloc(2, sizeof(double));
  }

  vtx[0][0] = 0.0;   vtx[0][1] = 0.0;  
  vtx[1][0] = hs;    vtx[1][1] = 0.5*hs;  
  vtx[2][0] = 0.0;   vtx[2][1] = hs;  

  nPnls = 2;
  pnls = (panel *)malloc(nPnls*sizeof(panel));
  pnls[0].v0 = vtx[0];
  pnls[0].idxv0 = 0; 
  pnls[0].v1 = vtx[1];
  pnls[0].idxv1 = 1;  

  pnls[1].v0 = vtx[1];
  pnls[1].idxv0 = 1; 
  pnls[1].v1 = vtx[2];
  pnls[1].idxv1 = 2;  

  panelIA(d, p, ht,  &pnls[0], &pnls[0], &pnls[0], &pnls[0], xLeg, wLeg, solution);
  sumV = solution[0] + solution[1] + solution[2] + solution[3];
  sumK = solution[4] + solution[5] + solution[6] + solution[7];
  printf("self:  %lf %lf\n", sumV, sumK );

  panelIA(d, p, ht,  &pnls[0], &pnls[0], &pnls[1], &pnls[1], xLeg, wLeg, solution);
  sumV = solution[0] + solution[1] + solution[2] + solution[3];
  sumK = solution[4] + solution[5] + solution[6] + solution[7];
  printf("vtx:   %lf %lf\n", sumV, sumK );
}

