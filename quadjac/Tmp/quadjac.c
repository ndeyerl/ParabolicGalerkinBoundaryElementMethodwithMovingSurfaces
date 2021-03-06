#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "panelIA.h"
#include "Laguerre.h"
#include "Jacobi.h"
#include "Hermite.h"


//implementation: make quadjac
//		./quadjac -p=%whatever desired p value is

void paramcircle(double theta, double *x);
void mparamcircle(double h, int i, double theta, double *x);
void paramsquare(double theta, double *x, int i);
void paramtriangle(double theta, double *x, int i);
void paramellipse(double theta, double *x);
double f1( double *x, double t);
double f2(double *x, double t);
void getfquad(double h, int i, int p, panel *pX0, panel *pX1, double xLeg[], double wLeg[], double solnquad[]);
double l2err(double h, int i, int p, panel *p1, panel *p2, double rhs, double xLeg[], double wLeg[]);
double getfnode(double h, int i, int p, double *x, double *y);
double rsqrt(double *x, double t);
double rprimesqrt(double *x, double t);
double ralph(double *x, double t);
/* LU factorization */
extern void  dgetrf_(int *nRow, int *nCol, double *A, int *lda, int *pvt, int *inf); 
/* LU solve */
extern void  dgetrs_(char *trans, int *n, int *nRhs, double *A, int *lda, int *pvt, double *B, int *ldb, int *inf); 
void MtVadd( int n, int m, double *A, double *x, double *y );
void MtVsub( int n, int m, double *A, double *x, double *y );
double g(double x[], double t, double ny[]);
void getdubdn(double h, int i, panel *p1, panel *p2, double dubdnsoln[]);
double max(double myArray[], int size);

double x0[2];
double alpha; //to change radius r(t), change mparam and getvel 

void main(int nargs, char *argv[]) {

  double h, theta, tmax=1.0, *xLeg, *wLeg, **vtx, xdiff, ydiff, norm, t0, t1, ny[2], fnodal, solnquad[2], *fquad, *fnode, **ftild, *B, *C, **V, **K, solution[8]={0,0,0,0,0,0,0,0}, *innerrhs, *innerdubdn, solnnodal, *z, solndubdn[2], *y, val1, val2, yy, maxerr, l2error, s;

  int p, d, dd, i, j, k, info, nVtx, nPnls, nRhs = 1, test=0, dmax, pnlsX, pnlsY, Xidxv0, Xidxv1, Yidxv0, Yidxv1, idxv0, idxv1, *pvt, n, m;
  char trans='N';
  int verbose=0;

  dmax = 10; //max number of time steps
  nVtx = 5;//number of vertices
  dd = 0;
  p = 5; 
  x0[0] = 0.5;
//  x0[0] = 0.0;
  x0[1] = 0.0;
  alpha = 1.0;

  /* parse the command line */
  for ( i=1; i<nargs; i++ )
    if ( argv[i][0] == '-' )
      switch ( argv[i][1] ) {
     case 't': test = atoi( argv[i]+3 );
        break;
      case 'p': p = atoi( argv[i]+3 ); //number of nodes
        break;
     case 'N': nVtx = atoi( argv[i]+3 ); //max # of vertices
        break;
     case 'M': dmax = atoi( argv[i]+3 ); //max # of time steps
        break;
     case 'T': tmax = atoi( argv[i]+3 ); //time step size numerator
        break;
     case 'd': dd = atoi( argv[i]+3 ); //which matrix to print out
        break;
     case 'v': verbose = atoi( argv[i]+3 ); 
        break;
      }
  nPnls = nVtx; //number of panels, -1 if panels dont wrap around
  h = tmax/dmax; //time step size
  printf("N = %d  M = %d  T = %f  p = %d \n",nVtx, dmax, tmax, p);

//build gauss legendre nodes
  xLeg = (double*)malloc( p*sizeof(double) );
  wLeg = (double*)malloc( p*sizeof(double) );
  Jacobi( p, 0.0, 0.0, xLeg, wLeg);


//allocate memory for vertices
  double ***vtxall;
  vtxall = (double ***)calloc(dmax+1,sizeof(double **));

//allocate space for all panels (all time steps)
    panel **pnls;
    pnls = (panel **)calloc(dmax+2, sizeof(panel*));
    panel *pnlsd;

//allocate memory for velocities
  double **velall;
  velall = (double **)calloc(dmax+2, sizeof(double *));

//build the matrices V, K
  V = (double **)calloc(dmax+1, sizeof(double*));//dmax number of timesteps 
  for(j = 0; j <= dmax; j++){
    V[j] = calloc(nPnls*nPnls, sizeof(double));
  }

  K = (double **)calloc(dmax+1, sizeof(double*));
  for(j = 0; j <= dmax; j++){
    K[j] = calloc(nVtx*nPnls, sizeof(double));
  }

//allocate memory for solution
  double *rhs; //holds dmax+1 of the dmax+1-arrays
  rhs = (double *)calloc(nVtx, sizeof(double));


  double **solnrhs; //holds dmax+1 arrays that are nPnls in length
  solnrhs = (double **)calloc(dmax+1, sizeof(double*));
  for(d = 0; d <= dmax; d++){ //loop over time steps
    innerrhs = (double *)calloc(nVtx, sizeof(double));
    solnrhs[d] = innerrhs;
  }


//allocate space for f, ftilde
  double **ftilde;
  ftilde = (double **)calloc(dmax+1, sizeof(double *));
  double **f;
  f = (double **)calloc(dmax+1, sizeof(double *));

//allocate space for pivots
  pvt  = (int*)calloc(nVtx, sizeof(int));



//allocate space for true solution
  double ** dubdn; //holds dmax+1 arrays that are nPnls in length
  dubdn = (double **)calloc(dmax+1, sizeof(double*));
  for(d = 0; d <= dmax; d++){ //loop over time steps
    innerdubdn = (double *)calloc(nVtx, sizeof(double));
    dubdn[d] = innerdubdn;
  }

//allocate space for errors
  double err[dmax+1];


//build panels
  for(i = 0; i <= dmax; i++){
    //build vertices
    vtx = (double **)calloc(nVtx, sizeof(double *));
    for(k = 0; k<nVtx; k++){
      vtx[k] = (double *)calloc(2, sizeof(double));
    }
    vtxall[i] = vtx;
  
    //geometry
    for(k = 0; k < nVtx; k++){
      theta = 2*M_PI*k/nVtx;
//      paramcircle(theta, vtxall[i][k]); //nonmoving circular geometry
      mparamcircle(h, i, theta, vtxall[i][k]); //moving circular geometry
    }

    //build panels
    pnlsd = (panel *)calloc(nPnls, sizeof(panel));
    for(k=0;k < nPnls; k++){
      pnlsd[k].v1 = vtxall[i][k];
      pnlsd[k].idxv1 = k;  
      pnlsd[k].v0 = vtxall[i][(k+1)%nPnls];
      pnlsd[k].idxv0 = (k+1)%nPnls; 
      xdiff = vtxall[i][k][0]-vtxall[i][(k+1)%nPnls][0];
      ydiff = vtxall[i][k][1]-vtxall[i][(k+1)%nPnls][1];
      norm = sqrt(xdiff*xdiff + ydiff*ydiff);
      pnlsd[k].len = norm;
    }
    pnls[i] = pnlsd;
  }

//solve the system
  for(i = 0; i < dmax; i++){
    if ( verbose ) printf("Time step %d\n",i);
    //build fquad (ftilde) @ step i
    fquad = (double *)calloc(nVtx, sizeof(double));
    for(k = 0; k < nPnls; k++){ //loop over panels (space)
      idxv0 = pnls[i][k].idxv0;  
      idxv1 = pnls[i][k].idxv1;  
      getfquad(h, i, p, &pnls[i][k], &pnls[i+1][k], xLeg, wLeg, solnquad); //moving geometry
      fquad[k] += solnquad[0]+solnquad[1];
    }
    ftilde[i] = fquad;

    //build fnode
    fnode = (double *)calloc(nVtx, sizeof(double)); 
    for(j = 0; j < nVtx; j++){
      solnnodal = getfnode(h, i, p, vtxall[i][j], vtxall[i+1][j]);
      fnode[j] = solnnodal;
    }
    f[i] = fnode;  



    //generate V, K and solve 

    for(m = 0; m <= i; m++){
      d = i-m;
      memset(V[m],0.0,nPnls*nPnls*sizeof(double));
      memset(K[m],0.0,nVtx*nPnls*sizeof(double));
      for(pnlsX = 0; pnlsX < nPnls; pnlsX++){
        for(pnlsY = 0; pnlsY < nPnls; pnlsY++){
          Xidxv0 = pnls[m][pnlsX].idxv0; 
          Xidxv1 = pnls[m][pnlsX].idxv1;
          Yidxv0 = pnls[m][pnlsY].idxv0;
          Yidxv1 = pnls[m][pnlsY].idxv1;
          panelIA(d, p, h, &pnls[i][pnlsX], &pnls[m][pnlsX], &pnls[i][pnlsY], &pnls[m][pnlsY],xLeg, wLeg, solution);
          V[m][pnlsX + nPnls*pnlsY] = solution[0]+solution[1]+solution[2]+solution[3];
          K[m][pnlsX + nPnls*Yidxv0] += solution[4]+solution[6];
          K[m][pnlsX + nPnls*Yidxv1] += solution[5]+solution[7];
        }
      }
      if ( verbose > 2 ) {
        printf("----matrix V[%d]------\n",m);
        int z;
        for(z = 0; z < nPnls; z++){
          for(j = 0; j < nPnls; j++){
            printf("%f ",V[m][z*nPnls+j]);
          } 
          printf("\n");
        }
        
        printf("----matrix K[%d]------\n",m);
        for(z = 0; z < nPnls ; z++){
          for(j = 0; j < nVtx; j++){
            printf("%.11f ",K[m][z*nPnls+j]);
          } 
          printf("\n");
        }
      }
    }//end m-loop

    //history
    for(k = 0; k < nVtx; k++){
      rhs[k] = -.5*ftilde[i][k];  
    } 

    //generate right hand side
    for ( m=0; m < i; m++ ){     
      MtVsub( nPnls, nVtx, V[m], solnrhs[m], rhs);     
    }

    for ( m=0; m <= i; m++ ){  
      MtVadd( nPnls, nVtx, K[m], f[m], rhs);  
    }

    //LU factorize A[0]
    dgetrf_(&nVtx, &nVtx, V[i], &nVtx, pvt, &info); 

    //solver
    dgetrs_( &trans, &nVtx, &nRhs, V[i], &nVtx, pvt, rhs, &nVtx, &info); 

    //save solution
    for(k = 0; k < nVtx; k++){
      solnrhs[i][k] = rhs[k]; 
    }

    //error calculations

    //max error
    z = (double *)calloc(nVtx, sizeof(double));
    for(j = 0; j < nVtx; j++){
      getdubdn(h, i, &pnls[i][j], &pnls[i+1][j], solndubdn); //moving geometry
      //        getdubdn(h, i, &pnls[i][j], &pnls[i][j], solndubdn); //nonmoving geometry
      z[j] += 0.5*(solndubdn[0]+solndubdn[1]);
    }
    dubdn[i] = z;

    if ( verbose>1 ) {
      for(j = 0; j < nPnls; j++){
        printf("solnrhs[%d][%d] = %.16f dubdn = %.16f \n",i,j,solnrhs[i][j], dubdn[i][j]);
      }
    }

    y = (double *)calloc(nVtx, sizeof(double));
    for(j = 0; j < nVtx; j++){
      y[j] = fabs(solnrhs[i][j]-dubdn[i][j]);
    }
    err[i] = max(y,nVtx);


    maxerr = max(err,dmax+1);

  }//end i-loop


  printf("N = %d  M = %d  p = %d\n",nVtx,dmax,p);
  printf("maxerr = %.10f\n",maxerr);

  //l2 error
  l2error = 0.0;
  for(i = 0; i < dmax; i++){
    for(j = 0; j < nPnls; j++){
//      s = l2err(h, i, p, &pnls[i][j], &pnls[i][j], solnrhs[i][j], xLeg, wLeg);//nonmoving geometry
      s = l2err(h, i, p, &pnls[i][j], &pnls[i+1][j], solnrhs[i][j], xLeg, wLeg);//moving geometry
      l2error += s;
    }
  }
  printf("l2 error = %.10f\n",sqrt(l2error));

}


///////////////////////////
//subroutines
///////////////////////////

void paramcircle(double theta, double *x){
  x[0] = cos(theta);
  x[1] = sin(theta);
}

void mparamcircle(double h, int i, double theta, double *x){
  double t, r;
  t = h*i;         //jtdeb: NOT the midpoint 
  r = rsqrt(x, t); 
  x[0] = r*cos(theta); //moving geometry with r(t)
  x[1] = r*sin(theta);
}

void paramsquare(double theta, double *x, int i){
  if(i==1){ //side 1
    x[0] = 1*(1-theta) + 1*theta;
    x[1] = -1*(1-theta) + 1*theta;
  }else if(i==2){ //side 2
    x[0] = 1*(1-theta) + -1*theta;
    x[1] = 1*(1-theta) + 1*theta;
  }else if(i==3){ //side 3
    x[0] = -1*(1-theta) + -1*theta;
    x[1] = 1*(1-theta) + -1*theta;
  }else if(i==4){ //side 4
    x[0] = -1*(1-theta) + 1*theta;
    x[1] = -1*(1-theta) + -1*theta;
  }
}

void paramtriangle(double theta, double *x, int i){
  if(i==1){ //side 1
    x[0] = 1*(1-theta) + 0*theta;
    x[1] = -1*(1-theta) + 1*theta;
  }else if(i==2){ //side 2
    x[0] = 0*(1-theta) + -1*theta;
    x[1] = 1*(1-theta) + -1*theta;
  }else if(i==3){ //side 3
    x[0] = -1*(1-theta) + 1*theta;
    x[1] = -1*(1-theta) + -1*theta;
  }
}

void paramellipse(double theta, double *x){
  x[0] = cos(theta);
  x[1] = 10*sin(theta);
}

// right hand side
double f1( double *x, double t) {
  return t;
} /* fcn */	

double f2(double *x, double t) {
  double r1, r2, r, soln;
  r1 = x[0]-x0[0];
  r2 = x[1]-x0[1];
  r = r1*r1 + r2*r2;
  soln = exp(-r/(4*t))/(4*M_PI*t);
  return soln;
  }

//jtdeb: corrections: Jacobian depends on time
void getfquad(double h, int i, int p, panel *pX0, panel *pX1, double xLeg[], double wLeg[], double solnquad[]){
  int i1, i2;
  double fquad, x[2], tquad, length, *v00, *v10, *v01, *v11, v0t[2], v1t[2];
  v00 = pX0->v0;
  v10 = pX0->v1;
  v01 = pX1->v0;
  v11 = pX1->v1;
  solnquad[0] = solnquad[1] = 0.0;
  for(i1 = 0; i1 < p; i1++){ //time loop
    tquad = h*(i + xLeg[i1]);   
    v0t[0] = v00[0]*(1 - tquad) + v01[0]*tquad; 
    v0t[1] = v00[1]*(1 - tquad) + v01[1]*tquad;
    v1t[0] = v10[0]*(1 - tquad) + v11[0]*tquad;
    v1t[1] = v10[1]*(1 - tquad) + v11[1]*tquad;
    length  = (v0t[0]-v1t[0])*(v0t[0]-v1t[0]);
    length += (v0t[1]-v1t[1])*(v0t[1]-v1t[1]);
    length = sqrt(length);
    for(i2 = 0; i2 < p; i2++){ //space loop
      x[0] =  xLeg[i2]*v0t[0] + (1-xLeg[i2])*v1t[0];
      x[1] =  xLeg[i2]*v0t[1] + (1-xLeg[i2])*v1t[1];  
      //fcn = f(x,t);   
      fquad = f2(x,tquad); 
      fquad *= wLeg[i1]*wLeg[i2]*length;    
      solnquad[0] += fquad*xLeg[i2];
      solnquad[1] += fquad*(1 - xLeg[i2]);
    }
  }
  solnquad[0] *= h;
  solnquad[1] *= h;
}



double l2err(double h, int i, int p, panel *p1, panel *p2, double rhs, double xLeg[], double wLeg[]){
  int i1, i2;
  double t, integ2, integ, xdiff, ydiff, norm;
  double *v00, *v10, *v01, *v11, v0t[2], v1t[2], ny[2], x[2];

  integ = 0.0;
  v00 = p1->v0;
  v10 = p1->v1;
  v01 = p2->v0;
  v11 = p2->v1;

  for(i1 = 0; i1 < p; i1++){//time loop
    t = h*(i + xLeg[i1]);
    v0t[0] = v00[0]*(1 - t) + v01[0]*t; 
    v0t[1] = v00[1]*(1 - t) + v01[1]*t;
    v1t[0] = v10[0]*(1 - t) + v11[0]*t;
    v1t[1] = v10[1]*(1 - t) + v11[1]*t;
    xdiff = v1t[0]-v0t[0];
    ydiff = v1t[1]-v0t[1];
    norm = sqrt(xdiff*xdiff + ydiff*ydiff);
    ny[0] = -ydiff/norm;
    ny[1] = xdiff/norm;
    for(i2 = 0; i2 < p; i2++){//space loop
      x[0] =  xLeg[i2]*v0t[0] + (1-xLeg[i2])*v1t[0];
      x[1] =  xLeg[i2]*v0t[1] + (1-xLeg[i2])*v1t[1];  
      integ2 = g(x, t, ny) - rhs;
      integ += wLeg[i1]*wLeg[i2]*integ2*integ2*norm;    //jtdeb:  norm depends on time
    }
  }
  return h*integ;
}


double getfnode(double h, int i, int p, double *x0, double *x1){
  double tnodal, solnnodal, x[3];                   //jtdeb:  interpolate point
  tnodal = h*(i + 0.5);
  x[0] = 0.5*(x0[0] + x1[0]); 
  x[1] = 0.5*(x0[1] + x1[1]); 
  x[2] = 0.5*(x0[2] + x1[2]); 

  solnnodal = f2(x,tnodal); 
//  printf("solnnodal = %f\n",solnnodal);
  return solnnodal;
}

//r(t) = sqrt(1+t)
//r'(t) = 1/2sqrt(1+t)
double rsqrt(double *x, double t) {
  double soln;
  soln = sqrt(1 + alpha*t);
  return soln;
  }

//r(t) = sqrt(1+t)
//r'(t) = 1/2sqrt(1+t)
double rprimesqrt(double *x, double t) {
  double soln;
  soln = 1/(2*sqrt(1+t));
  return soln;
  }

//r(t) = 1+alpha*t
//r'(t) = alpha
double ralph(double *x, double t) {
  double soln;
  soln = 1.0 + alpha*t;
  return soln;
  }




/* Johannes Tausch
 * Computes the matrix-vector product
 *    y += A*x
 * where A is an n-by-n matrix and adds the result to the vector y
 */
// modified to take in rectangular matrices
void MtVadd( int n, int m, double *A, double *x, double *y){
  int i, j;
//rows n, columns m
  for ( i=0; i<m; i++ ) {
    for ( j=0; j<n; j++ ) {
      y[i] += A[i+n*j]*x[j];
    }
  }
} /* MtVadd */

 /* Johannes Tausch
 * Computes the matrix-vector product
 *    y -= A*x
 * where A is an n-by-n matrix and subtracts the result from the vector y
 */
// modified to take in rectangular matrices
void MtVsub( int n, int m, double *A, double *x, double *y ){
  int i, j;
//rows n, columns m
  for ( i=0; i<m; i++ ) {
    for ( j=0; j<n; j++ ) {
      y[i] -= A[i+n*j]*x[j];
    }
  }
} /* MtVsub */

double g(double x[], double t, double ny[]) {
  double r0, r1, r, rny, soln;
  r0 = x[0]-x0[0];
  r1 = x[1]-x0[1];
  r = r0*r0 + r1*r1;
  rny = r0*ny[0] + r1*ny[1];
  soln = -exp(-r/(4*t))*rny/(4*M_PI*2*t*t);
  return soln;
  }
/*
  fcnVal[0] = exp(-rr/(4*del))/(del*4*M_PI);
  fcnVal[1] = drr*fcnVal[0]/(2*del);
*/


void getdubdn(double h, int i, panel *p1, panel *p2, double dubdnsoln[]){
  double t, xdiff, ydiff, norm;
  double *v00, *v10, *v01, *v11, v0t[2], v1t[2], ny[2];
  v00 = p1->v0;
  v10 = p1->v1;
  v01 = p2->v0;
  v11 = p2->v1;
  t = h*(i + 0.5);
  v0t[0] = v00[0]*(1 - t) + v01[0]*t; 
  v0t[1] = v00[1]*(1 - t) + v01[1]*t;
  v1t[0] = v10[0]*(1 - t) + v11[0]*t;
  v1t[1] = v10[1]*(1 - t) + v11[1]*t;
  xdiff = v1t[0]-v0t[0];
  ydiff = v1t[1]-v0t[1];
  norm = sqrt(xdiff*xdiff + ydiff*ydiff);
  ny[0] = -ydiff/norm;
  ny[1] = xdiff/norm;
  dubdnsoln[0] = 0.0; dubdnsoln[1] = 0.0;
  dubdnsoln[0] = g(v0t,t, ny); 
  dubdnsoln[1] = g(v1t,t, ny);
}


/*http://stackoverflow.com/questions/1690428/finding-max-number-in-an-array-c-programming
*/
double max(double myArray[], int size) {
    double maxValue = myArray[0];
    int i;

    for (i = 1; i < size; ++i) {
        if ( myArray[i] > maxValue ) {
            maxValue = myArray[i];
        }
    }
    return maxValue;
}



