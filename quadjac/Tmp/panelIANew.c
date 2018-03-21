//compile object file:
//gcc -c panelIA.c

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "panelIA.h"

#define TOL 1e-12

void kernel(double r[], double del, double ny[], double fcnVal[]){
  double rr, drr;
  rr = r[0]*r[0] + r[1]*r[1];
  drr = r[0]*ny[0] + r[1]*ny[1];
/* 
  printf("ny[0] = %f ny[1] = %f\n",ny[0],ny[1]);
  printf("drr = %f \n",drr);
  printf("r[0] = %lf r[1] = %lf  s=%lf\n",r[0],r[1],del);
*/

  fcnVal[0] = exp(-rr/(4*del))/(del*4*M_PI);
  fcnVal[1] = drr*fcnVal[0]/(2*del);
}


int nrcommonvtx(panel *p1, panel *p2) {
  if (p1->idxv0==p2->idxv0){

    return 2;//panels are the same
    
  }else if (p1->idxv0==p2->idxv1) {

    return 1;//panels share 1 vertex, pX on left, pY on right
    
  }else if (p1->idxv1==p2->idxv0) {

    return -1;//panels share 1 vertex, pX on right, pY on left

  }else {

    return 0;//panels are not neighbors

  }
}



void intgrnonsing(int d, int p, double h, panel *pX0, panel *pX1, panel *pY0, panel *pY1, double xLeg[], double wLeg[], double solution[]) {
  double *xhat, *yhat, *that, *tauhat, *v00, *v10, *v01, *v11, *w00, *w10, *w01, *w11, r[2], del, jac, fcnVal[2], dfcn, fcn, v0t[2], v1t[2], w0t[2], w1t[2], xdiff, ydiff, ny[2], norma, normb, nVeloc, w0Prm[2], w1Prm[2], shape[4], h2;
  int i, i1, i2, i3, i4;

  h2 = h*h;
  v00 = pX0->v0;
  v10 = pX0->v1;
  v01 = pX1->v0;
  v11 = pX1->v1;
  w00 = pY0->v0;
  w10 = pY0->v1;
  w01 = pY1->v0;
  w11 = pY1->v1;
  w0Prm[0] = w01[0] - w00[0];
  w0Prm[1] = w01[1] - w00[1];
  w1Prm[0] = w11[0] - w10[0];
  w1Prm[1] = w11[1] - w10[1];
  xhat = xLeg;
  yhat = xLeg;
  tauhat = xLeg;
  that = xLeg; 
  for(i3 = 0; i3 < p; i3++) { //that loop
    for(i4 = 0; i4 < p; i4++) { //tauhat loop
      del = h*(d + that[i3] - tauhat[i4]);
      v0t[0] = v00[0]*(1 - that[i3]) + v01[0]*that[i3]; 
      v0t[1] = v00[1]*(1 - that[i3]) + v01[1]*that[i3];
      v1t[0] = v10[0]*(1 - that[i3]) + v11[0]*that[i3];
      v1t[1] = v10[1]*(1 - that[i3]) + v11[1]*that[i3];
      w0t[0] = w00[0]*(1 - tauhat[i4]) + w01[0]*tauhat[i4];
      w0t[1] = w00[1]*(1 - tauhat[i4]) + w01[1]*tauhat[i4];
      w1t[0] = w10[0]*(1 - tauhat[i4]) + w11[0]*tauhat[i4];
      w1t[1] = w10[1]*(1 - tauhat[i4]) + w11[1]*tauhat[i4];
      norma = sqrt((v1t[0]-v0t[0])*(v1t[0]-v0t[0])+(v1t[1]-v0t[1])*(v1t[1]-v0t[1]));
      xdiff = w1t[0]-w0t[0]; 
      ydiff = w1t[1]-w0t[1];
      normb = sqrt(xdiff*xdiff + ydiff*ydiff);
      ny[0] = -ydiff/normb; 
      ny[1] = xdiff/normb; 


      if ( del>TOL ) {
        for(i2 = 0; i2 < p; i2++) { //yhat loop
          nVeloc  = (w0Prm[0]*(1- yhat[i2]) + w0Prm[0]*yhat[i2])*ny[0];
          nVeloc += (w0Prm[1]*(1- yhat[i2]) + w0Prm[1]*yhat[i2])*ny[1];
          nVeloc /= h;
          for(i1 = 0; i1 < p; i1++) { //xhat loop
            r[0] = v1t[0]*(1 - xhat[i1]) + v0t[0]*xhat[i1] - 
              w1t[0]*(1 - yhat[i2]) - w0t[0]*yhat[i2];
            r[1] = v1t[1]*(1 - xhat[i1]) + v0t[1]*xhat[i1] -
              w1t[1]*(1 - yhat[i2]) - w0t[1]*yhat[i2];
            //heat kernel
            kernel(r, del, ny, fcnVal);
            jac = norma*normb;
	    fcn = fcnVal[0]*jac;
            dfcn = (fcnVal[1]-nVeloc)*jac; 
            //weights
            fcn *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
            dfcn *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
            //shape functions
            shape[0] = yhat[i2]*xhat[i1];
            shape[1] = (1 - yhat[i2])*xhat[i1];
            shape[2] = yhat[i2]*(1 - xhat[i1]);
            shape[3] = (1 - yhat[i2])*(1 - xhat[i1]);

            for (i=0; i<4; i++) {
              solution[i] += fcn*shape[i];
              solution[4+i] += dfcn*shape[i];
            }
          }
        }
      }
    }
  }
  for (i=0;i<8; i++ ) solution[i] *= h2;
}


void intgrd1nc1(int d, int p, double h, panel *pX0, panel *pX1, panel *pY0, panel *pY1, double xLeg[], double wLeg[], double solution[], int nc) {
  double *eta, *xi1, *xi2, *xi3, etasq, xi1sq, xi2sq, del1, del2, del3, fcn1, fcn2, fcn3, fcn4, fcn, eta3, r1[2], r2[2], r3[2], rr1, rr2, rr3, r[2], fcnVal[2], dfcn1, dfcn2, dfcn3, dfcn4, rr, xhat1, yhat1, xhat2, yhat2, xhat3, yhat3, jac, *v00, *v10, *v01, *v11, *w00, *w10, *w01, *w11, v0t[2], v1t[2], w0t[2], w1t[2], xdiff, ydiff, ny[2], norma, normb, a[2], b[2];
  double w0Prm[2], w1Prm[2], *dw0, *dw1, nVeloc;
  int i1, i2, i3, i4;

  v00 = pX0->v0;
  v10 = pX0->v1;
  v01 = pX1->v0;
  v11 = pX1->v1;
  w00 = pY0->v0;
  w10 = pY0->v1;
  w01 = pY1->v0;
  w11 = pY1->v1;
  w0Prm[0] = w01[0] - w00[0];
  w0Prm[1] = w01[1] - w00[1];
  w1Prm[0] = w11[0] - w10[0];
  w1Prm[1] = w11[1] - w10[1];
  if(nc==1){
    dw0 = w1Prm;
    dw1 = w0Prm;
  }else if(nc==-1){
    dw0 = w0Prm;
    dw1 = w1Prm;
  }

  eta = xLeg;
  xi1 = xLeg;
  xi2 = xLeg;
  xi3 = xLeg;
  //  rho = norma*norma/(4*h);
  for(i1 = 0; i1 < p; i1++) { //xi1 loop
    for(i2 = 0; i2 < p; i2++) { //xi2 loop
      for(i3 = 0; i3 < p; i3++) { //xi3 loop
	for(i4 = 0; i4 < p; i4++) { //eta loop
	  eta3 = eta[i4]*eta[i4]*eta[i4];
	  del1 = h*eta[i4]*eta[i4]*(xi2[i2]*xi2[i2] + xi3[i3]*xi3[i3]);
	  del2 = h*eta[i4]*eta[i4]*(1 + xi3[i3]*xi3[i3]);
          //t = eta*xi2 tau = eta*xi3
          v0t[0] = v00[0]*(1 - eta[i4]*xi2[i2]) + v01[0]*eta[i4]*xi2[i2]; 
          v0t[1] = v00[1]*(1 - eta[i4]*xi2[i2]) + v01[1]*eta[i4]*xi2[i2];
          v1t[0] = v10[0]*(1 - eta[i4]*xi2[i2]) + v11[0]*eta[i4]*xi2[i2];
          v1t[1] = v10[1]*(1 - eta[i4]*xi2[i2]) + v11[1]*eta[i4]*xi2[i2];
          w0t[0] = w00[0]*(1 - eta[i4]*xi3[i3]) + w01[0]*eta[i4]*xi3[i3];
          w0t[1] = w00[1]*(1 - eta[i4]*xi3[i3]) + w01[1]*eta[i4]*xi3[i3];
          w1t[0] = w10[0]*(1 - eta[i4]*xi3[i3]) + w11[0]*eta[i4]*xi3[i3];
          w1t[1] = w10[1]*(1 - eta[i4]*xi3[i3]) + w11[1]*eta[i4]*xi3[i3];
          if(nc==1){
          a[0] = v1t[0]-v0t[0];
          a[1] = v1t[1]-v0t[1];
          b[0] = w0t[0]-w1t[0];
          b[1] = w0t[1]-w1t[1];
          }else if(nc==-1){
          a[0] = v0t[0] - v1t[0];
          a[1] = v0t[1] - v1t[1];
          b[0] = w1t[0] - w0t[0];
          b[1] = w1t[1] - w0t[1];
          }
	  r1[0] = eta[i4]*a[0] - eta[i4]*xi1[i1]*b[0];
	  r1[1] = eta[i4]*a[1] - eta[i4]*xi1[i1]*b[1];
          xdiff = w1t[0]-w0t[0];
          ydiff = w1t[1]-w0t[1];
          normb = sqrt(xdiff*xdiff + ydiff*ydiff);
          ny[0] = -ydiff/normb;
          ny[1] = xdiff/normb;
          norma = sqrt((v1t[0]-v0t[0])*(v1t[0]-v0t[0])+(v1t[1]-v0t[1])*(v1t[1]-v0t[1]));
          jac = norma*normb*4*h*h*eta[i4]*xi2[i2]*eta[i4]*xi3[i3]*eta3;
          nVeloc  = (dw0[0]*(1-xi1[i1]*eta[i4]) + dw0[0]*xi1[i1]*eta[i4])*ny[0];
          nVeloc += (dw0[1]*(1-xi1[i1]*eta[i4]) + dw0[1]*xi1[i1]*eta[i4])*ny[1];
          nVeloc /= h;
          kernel(r1, del1, ny, fcnVal);
          fcn1 = fcnVal[0]*jac;
          dfcn1 = (fcnVal[1]-nVeloc)*jac;

          //t, tau have same transformation as fcn1
	  r2[0] = eta[i4]*xi1[i1]*a[0] - eta[i4]*b[0];
	  r2[1] = eta[i4]*xi1[i1]*a[1] - eta[i4]*b[1];
          jac = norma*normb*4*h*h*eta[i4]*xi2[i2]*eta[i4]*xi3[i3]*eta3;
          nVeloc  = (dw0[0]*(1-eta[i4]) + dw0[0]*eta[i4])*ny[0];
          nVeloc += (dw0[1]*(1-eta[i4]) + dw0[1]*eta[i4])*ny[1];
          nVeloc /= h;
          kernel(r2, del1, ny, fcnVal);
          fcn2 = fcnVal[0]*jac;
          dfcn2 = (fcnVal[1]-nVeloc)*jac;
 
          //t = eta tau = eta*xi3 (tau same as above->w0t, w1t same as above)
          v0t[0] = v00[0]*(1 - eta[i4]) + v01[0]*eta[i4]; 
          v0t[1] = v00[1]*(1 - eta[i4]) + v01[1]*eta[i4];
          v1t[0] = v10[0]*(1 - eta[i4]) + v11[0]*eta[i4];
          v1t[1] = v10[1]*(1 - eta[i4]) + v11[1]*eta[i4];
          w0t[0] = w00[0]*(1 - eta[i4]*xi3[i3]) + w01[0]*eta[i4]*xi3[i3];
          w0t[1] = w00[1]*(1 - eta[i4]*xi3[i3]) + w01[1]*eta[i4]*xi3[i3];
          w1t[0] = w10[0]*(1 - eta[i4]*xi3[i3]) + w11[0]*eta[i4]*xi3[i3];
          w1t[1] = w10[1]*(1 - eta[i4]*xi3[i3]) + w11[1]*eta[i4]*xi3[i3];
          //b same as above, xdiff, ydiff, normb, ny same as above
          if(nc==1){
          a[0] = v1t[0]-v0t[0];
          a[1] = v1t[1]-v0t[1];
          b[0] = w0t[0]-w1t[0];
          b[1] = w0t[1]-w1t[1];
          }else if(nc==-1){
          a[0] = v0t[0] - v1t[0];
          a[1] = v0t[1] - v1t[1];
          b[0] = w1t[0] - w0t[0];
          b[1] = w1t[1] - w0t[1];
          }
          xdiff = w1t[0]-w0t[0];
          ydiff = w1t[1]-w0t[1];
          normb = sqrt(xdiff*xdiff + ydiff*ydiff);
          ny[0] = -ydiff/normb;
          ny[1] = xdiff/normb;
          norma = sqrt((v1t[0]-v0t[0])*(v1t[0]-v0t[0])+(v1t[1]-v0t[1])*(v1t[1]-v0t[1]));
	  r3[0] = eta[i4]*xi1[i1]*a[0] - eta[i4]*xi2[i2]*b[0];
	  r3[1] = eta[i4]*xi1[i1]*a[1] - eta[i4]*xi2[i2]*b[1];
          jac = norma*normb*4*h*h*eta[i4]*eta[i4]*xi3[i3]*eta3;
          nVeloc  = (dw0[0]*(1-xi2[i2]*eta[i4]) + dw0[0]*xi2[i2]*eta[i4])*ny[0];
          nVeloc += (dw0[1]*(1-xi2[i2]*eta[i4]) + dw0[1]*xi2[i2]*eta[i4])*ny[1];
          nVeloc /= h;
          kernel(r3, del2, ny, fcnVal);
	  fcn3 = fcnVal[0]*jac;
          dfcn3 = (fcnVal[1]-nVeloc)*jac;

          //t = eta*xi3, tau = eta
          v0t[0] = v00[0]*(1 - eta[i4]*xi3[i3]) + v01[0]*eta[i4]*xi3[i3]; 
          v0t[1] = v00[1]*(1 - eta[i4]*xi3[i3]) + v01[1]*eta[i4]*xi3[i3];
          v1t[0] = v10[0]*(1 - eta[i4]*xi3[i3]) + v11[0]*eta[i4]*xi3[i3];
          v1t[1] = v10[1]*(1 - eta[i4]*xi3[i3]) + v11[1]*eta[i4]*xi3[i3];
          w0t[0] = w00[0]*(1 - eta[i4]) + w01[0]*eta[i4];
          w0t[1] = w00[1]*(1 - eta[i4]) + w01[1]*eta[i4];
          w1t[0] = w10[0]*(1 - eta[i4]) + w11[0]*eta[i4];
          w1t[1] = w10[1]*(1 - eta[i4]) + w11[1]*eta[i4];
          if(nc==1){
          a[0] = v1t[0]-v0t[0];
          a[1] = v1t[1]-v0t[1];
          b[0] = w0t[0]-w1t[0];
          b[1] = w0t[1]-w1t[1];
          }else if(nc==-1){
          a[0] = v0t[0] - v1t[0];
          a[1] = v0t[1] - v1t[1];
          b[0] = w1t[0] - w0t[0];
          b[1] = w1t[1] - w0t[1];
          }
          xdiff = w1t[0]-w0t[0];
          ydiff = w1t[1]-w0t[1];
          normb = sqrt(xdiff*xdiff + ydiff*ydiff);
          ny[0] = -ydiff/normb;
          ny[1] = xdiff/normb;
          norma = sqrt((v1t[0]-v0t[0])*(v1t[0]-v0t[0])+(v1t[1]-v0t[1])*(v1t[1]-v0t[1]));
          jac = norma*normb*4*h*h*eta[i4]*xi3[i3]*eta[i4]*eta3;
          kernel(r3, del2, ny, fcnVal);
	  fcn4 = fcnVal[0]*jac;
          dfcn4 = (fcnVal[1]-nVeloc)*jac;

	  //weights, jacobian 
	  fcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  fcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  fcn3 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  fcn4 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn3 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn4 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];

	  //test functions
          xhat1 = eta[i4];
          yhat1 = eta[i4]*xi1[i1];
          xhat2 = eta[i4]*xi1[i1];
          yhat2 = eta[i4];
          xhat3 = eta[i4]*xi1[i1];
          yhat3 = eta[i4]*xi2[i2];

	  solution[0] += fcn1*yhat1*xhat1 + fcn2*yhat2*xhat2 + fcn3*yhat3*xhat3 + fcn4*yhat3*xhat3;
	  solution[1] += fcn1*(1 - yhat1)*xhat1 + fcn2*(1 - yhat2)*xhat2 + fcn3*(1 - yhat3)*xhat3 + fcn4*(1 - yhat3)*xhat3;
	  solution[2] += fcn1*yhat1*(1 - xhat1) + fcn2*yhat2*(1 - xhat2) + fcn3*yhat3*(1 - xhat3) + fcn4*yhat3*(1 - xhat3);
	  solution[3] += fcn1*(1 - yhat1)*(1 - xhat1) + fcn2*(1 - yhat2)*(1 - xhat2) + fcn3*(1 - yhat3)*(1 - xhat3)+ fcn4*(1 - yhat3)*(1 - xhat3); 

	  solution[4] += dfcn1*yhat1*xhat1 + dfcn2*yhat2*xhat2 + dfcn3*yhat3*xhat3 + dfcn4*yhat3*xhat3;
	  solution[5] += dfcn1*(1 - yhat1)*xhat1 + dfcn2*(1 - yhat2)*xhat2 + dfcn3*(1 - yhat3)*xhat3 + dfcn4*(1 - yhat3)*xhat3;
	  solution[6] += dfcn1*yhat1*(1 - xhat1) + dfcn2*yhat2*(1 - xhat2) + dfcn3*yhat3*(1 - xhat3) + dfcn4*yhat3*(1 - xhat3);
	  solution[7] += dfcn1*(1 - yhat1)*(1 - xhat1) + dfcn2*(1 - yhat2)*(1 - xhat2) + dfcn3*(1 - yhat3)*(1 - xhat3)+ dfcn4*(1 - yhat3)*(1 - xhat3); 	

	}
      }
    }
  }
}


/*
 * For the linearly moving panel compute the distance, jacobians and normal
 * Parameters
 *    v00, v10, v01, v11       pX0->v0, pX0->v1, pX1->v0, pX1->v1  vertices of x-panel for t=0 and t=1
 *    w00, w10, w01, w11       pY0->v0, pY0->v1, pY1->v0, pY1->v1  vertices of y-panel for t=0 and t=1
 *    t, tau, x, y             space/time coordinates in [0,1]^4
 *    r                        distance vector  (output)
 *    norma, normb             Jacobians, i.e., lengths (output)
 *    ny                       normal of y-panel (output)
 */
void getVectors(double *v00,double *v10, double *v01,double *v11, double *w00,double *w10, double *w01,double *w11, 
    double t, double tau, double x, double y, double *r, double *norma, double *normb, double *ny){

  double v0t[2], v1t[2], w0t[2], w1t[2], a[2], b[2];

  v0t[0] = v00[0]*(1 - t) + v01[0]*t; 
  v0t[1] = v00[1]*(1 - t) + v01[1]*t; 
  v1t[0] = v10[0]*(1 - t) + v11[0]*t; 
  v1t[1] = v10[1]*(1 - t) + v11[1]*t; 
  w0t[0] = w00[0]*(1 - tau) + w01[0]*tau; 
  w0t[1] = w00[1]*(1 - tau) + w01[1]*tau; 
  w1t[0] = w10[0]*(1 - tau) + w11[0]*tau; 
  w1t[1] = w10[1]*(1 - tau) + w11[1]*tau; 

  r[0] = v0t[0]*(1 - x) + v1t[0]*x;
  r[1] = v0t[1]*(1 - x) + v1t[1]*x;
  r[0] -= w0t[0]*(1 - y) + w1t[0]*y;
  r[1] -= w0t[1]*(1 - y) + w1t[1]*y;

  a[0] = v1t[0] - v0t[0];
  a[1] = v1t[1] - v0t[1];
  b[0] = w1t[0] - w0t[0];
  b[1] = w1t[1] - w0t[1];
  *norma = sqrt(a[0]*a[0] + a[1]*a[1]);
  *normb = sqrt(b[0]*b[0] + b[1]*b[1]);

  ny[0] = -b[1]/(*normb);
  ny[1] = b[0]/(*normb);
} /* getVectors */


void intgrd1nc2(int d, int p, double h, panel *pM, panel *pN, panel *pP, double xLeg[], double wLeg[], double solution[]) {
  double eta1, eta2, eta3, xi, t, tau, s, s2, z, x, y, nVeloc, norma, normb, jacobi, weights; 
  double dw0[2], db[2], r[2], ny[2], v[2], fcnVal[2], shape[4];
  double *v0M, *v1M, *v0N, *v1N, *v0P, *v1P;
  int i1, i2, i3, i;

  v0M = pM->v0;
  v1M = pM->v1;
  v0N = pN->v0;
  v1N = pN->v1;
  v0P = pP->v0;
  v1P = pP->v1;

  dw0[0] =  v0N[0] - v0M[0];
  dw0[1] =  v0N[1] - v0M[1];
  db[0] = v1N[0] - v1M[0] - v0N[0] + v0M[0];
  db[1] = v1N[1] - v1M[1] - v0N[1] + v0M[1];

  for(i = 0; i < p; i++) { 
    xi = xLeg[i]; 
    for(i1 = 0; i1 < p; i1++) { 
      eta1 = xLeg[i1];
      for(i2 = 0; i2 < p; i2++) { 
        eta2 = xLeg[i2];
         for(i3 = 0; i3 < p; i3++) { 
           eta3 = xLeg[i3];
           weights = wLeg[i]*wLeg[i1]*wLeg[i2]*wLeg[i3];
           /* Transformation 1a */
           s = xi*eta1;
           s2 = s*s;
           tau = eta2*s2;
           t = s2 - tau;

           z = xi;
           y = eta3*(1-z);
           x = z + y;
           jacobi = 2*xi*s*s2*(1 - z);
           v[0] = dw0[0] + db[0]*y;
           v[1] = dw0[1] + db[1]*y;
           getVectors(v0N, v1N, v0P, v1P, v0N, v1N, v0M, v1M, t, tau, x, y, r, &norma, &normb, ny);
           jacobi *= norma*normb;
           nVeloc = (v[0]*ny[0] + v[1]*ny[1])/h;
 
           kernel(r, h*s2, ny, fcnVal);
           fcnVal[1] -= nVeloc*fcnVal[0];
           fcnVal[0] *= jacobi*weights;
           fcnVal[1] *= jacobi*weights;

           shape[0] = x*y;
           shape[1] = x*(1 - y);
           shape[2] = y*(1 - x);
           shape[3] = (1 - x)*(1 - y);
       
           solution[0] += fcnVal[0]*shape[0];
           solution[1] += fcnVal[0]*shape[1];
           solution[2] += fcnVal[0]*shape[2];
           solution[3] += fcnVal[0]*shape[3];
           solution[4] += fcnVal[1]*shape[0];
           solution[5] += fcnVal[1]*shape[1];
           solution[6] += fcnVal[1]*shape[2];
           solution[7] += fcnVal[1]*shape[3];

           /* Transformation 1b,   s,tau,t, are the same */
           z = -xi;
           y = -z + eta3*(1+z);
           x = z + y;
           jacobi = 2*xi*s*s2*(1 + z);
           v[0] = dw0[0] + db[0]*y;
           v[1] = dw0[1] + db[1]*y;
           getVectors(v0N, v1N, v0P, v1P, v0N, v1N, v0M, v1M, t, tau, x, y, r, &norma, &normb, ny);
           jacobi *= norma*normb;
           nVeloc = (v[0]*ny[0] + v[1]*ny[1])/h;
 
           kernel(r, h*s2, ny, fcnVal);
           fcnVal[1] -= nVeloc*fcnVal[0];
           fcnVal[0] *= jacobi*weights;
           fcnVal[1] *= jacobi*weights;

           shape[0] = x*y;
           shape[1] = x*(1 - y);
           shape[2] = y*(1 - x);
           shape[3] = (1 - x)*(1 - y);
       
           solution[0] += fcnVal[0]*shape[0];
           solution[1] += fcnVal[0]*shape[1];
           solution[2] += fcnVal[0]*shape[2];
           solution[3] += fcnVal[0]*shape[3];
           solution[4] += fcnVal[1]*shape[0];
           solution[5] += fcnVal[1]*shape[1];
           solution[6] += fcnVal[1]*shape[2];
           solution[7] += fcnVal[1]*shape[3];

           /* Transformation 2a */
           s = xi;
           s2 = s*s;
           tau = eta2*s2;
           t = s2 + tau;

           z = eta1*xi;
           y = eta3*(1-z);
           x = z + y;
           jacobi = 2*xi*s*s2*(1 - z);
           v[0] = dw0[0] + db[0]*y;
           v[1] = dw0[1] + db[1]*y;
           getVectors(v0N, v1N, v0P, v1P, v0N, v1N, v0M, v1M, t, tau, x, y, r, &norma, &normb, ny);
           jacobi *= norma*normb;
           nVeloc = (v[0]*ny[0] + v[1]*ny[1])/h;
 
           kernel(r, h*s2, ny, fcnVal);
           fcnVal[1] -= nVeloc*fcnVal[0];
           fcnVal[0] *= jacobi*weights;
           fcnVal[1] *= jacobi*weights;

           shape[0] = x*y;
           shape[1] = x*(1 - y);
           shape[2] = y*(1-x);
           shape[3] = (1 - x)*(1 - y);
       
           solution[0] += fcnVal[0]*shape[0];
           solution[1] += fcnVal[0]*shape[1];
           solution[2] += fcnVal[0]*shape[2];
           solution[3] += fcnVal[0]*shape[3];
           solution[4] += fcnVal[1]*shape[0];
           solution[5] += fcnVal[1]*shape[1];
           solution[6] += fcnVal[1]*shape[2];
           solution[7] += fcnVal[1]*shape[3];

           /* Transformation 2b,   s,tau,t, are the same */
           z = -eta1*xi;
           y = -z + eta3*(1+z);
           x = z + y;
           jacobi = 2*xi*s*s2*(1 + z);
           v[0] = dw0[0] + db[0]*y;
           v[1] = dw0[1] + db[1]*y;
           getVectors(v0N, v1N, v0P, v1P, v0N, v1N, v0M, v1M, t, tau, x, y, r, &norma, &normb, ny);
           jacobi *= norma*normb;
           nVeloc = (v[0]*ny[0] + v[1]*ny[1])/h;
 
           kernel(r, h*s2, ny, fcnVal);
           fcnVal[1] -= nVeloc*fcnVal[0];
           fcnVal[0] *= jacobi*weights;
           fcnVal[1] *= jacobi*weights;

           shape[0] = x*y;
           shape[1] = x*(1 - y);
           shape[2] = y*(1 - x);
           shape[3] = (1 - x)*(1 - y);
       
           solution[0] += fcnVal[0]*shape[0];
           solution[1] += fcnVal[0]*shape[1];
           solution[2] += fcnVal[0]*shape[2];
           solution[3] += fcnVal[0]*shape[3];
           solution[4] += fcnVal[1]*shape[0];
           solution[5] += fcnVal[1]*shape[1];
           solution[6] += fcnVal[1]*shape[2];
           solution[7] += fcnVal[1]*shape[3];
 

           /* nonsingular part */
           t = 1 - xi*(1 - eta1);
           tau = 1 - xi*eta1;
           y = eta2;
           x = eta3;
           jacobi = xi;
           v[0] = dw0[0] + db[0]*y;
           v[1] = dw0[1] + db[1]*y;
           getVectors(v0N, v1N, v0P, v1P, v0N, v1N, v0M, v1M, t, tau, x, y, r, &norma, &normb, ny);
           jacobi *= norma*normb;
           nVeloc = (v[0]*ny[0] + v[1]*ny[1])/h;
 
           kernel(r, s2, ny, fcnVal);
           fcnVal[1] -= nVeloc*fcnVal[0];
           fcnVal[0] *= jacobi*weights;
           fcnVal[1] *= jacobi*weights;

           shape[0] = x*y;
           shape[1] = x*(1 - y);
           shape[2] = y*(1 - x);
           shape[3] = (1 - x)*(1 - y);
       
           solution[0] += fcnVal[0]*shape[0];
           solution[1] += fcnVal[0]*shape[1];
           solution[2] += fcnVal[0]*shape[2];
           solution[3] += fcnVal[0]*shape[3];
           solution[4] += fcnVal[1]*shape[0];
           solution[5] += fcnVal[1]*shape[1];
           solution[6] += fcnVal[1]*shape[2];
           solution[7] += fcnVal[1]*shape[3];
          }
      }
    }
  }
  for (i = 0; i < 8; i++) solution[i] *= h*h; 


}






void intgrd0nc1(int d, int p, double h, panel *pX0, panel *pX1, panel *pY0, panel *pY1, double xLeg[], double wLeg[], double solution[], int nc) {
  double *omega, *eta, *sigma, *xi1, *xi2, *v00, *v10, *v01, *v11, *w00, *w10, *w01, *w11, etasq, xi2sq, r[2], fcnVal[2], del, jac, fcn1, fcn2, fcn3, xhat1, xhat2, xhat3, yhat1, yhat2, yhat3, dfcn1, dfcn2, dfcn3, v0t[2], v1t[2], w0t[2], w1t[2], xdiff, ydiff, norm, ny[2], norma, normb, a[2], b[2];
  int i1, i2, i3, i4;
  v00 = pX0->v0;
  v10 = pX0->v1;
  v01 = pX1->v0;
  v11 = pX1->v1;
  w00 = pY0->v0;
  w10 = pY0->v1;
  w01 = pY1->v0;
  w11 = pY1->v1;
  omega = xLeg;
  eta = xLeg;
  xi1 = xLeg;
  xi2 = xLeg; 
  sigma = xLeg;
  for(i1 = 0; i1 < p; i1++) { //xi1 loop
    for(i2 = 0; i2 < p; i2++) { //xi2 loop
      for(i3 = 0; i3 < p; i3++) { //eta loop
	for(i4 = 0; i4 < p; i4++) { //sigma loop
          etasq = eta[i3]*eta[i3];
          xi2sq = xi2[i2]*xi2[i2];

          //transf1: t = s + sigma = eta*xi2 + sigma, tau = sigma
          v0t[0] = v00[0]*(1 - etasq*xi2sq - sigma[i4]*(1-etasq*xi2sq)) + v01[0]*(etasq*xi2sq + sigma[i4]*(1-etasq*xi2sq)); 
          v0t[1] = v00[1]*(1 - etasq*xi2sq - sigma[i4]*(1-etasq*xi2sq)) + v01[1]*(etasq*xi2sq + sigma[i4]*(1-etasq*xi2sq));
          v1t[0] = v10[0]*(1 - etasq*xi2sq - sigma[i4]*(1-etasq*xi2sq)) + v11[0]*(etasq*xi2sq + sigma[i4]*(1-etasq*xi2sq));
          v1t[1] = v10[1]*(1 - etasq*xi2sq - sigma[i4]*(1-etasq*xi2sq)) + v11[1]*(etasq*xi2sq + sigma[i4]*(1-etasq*xi2sq));
          w0t[0] = w00[0]*(1 - sigma[i4]*(1-etasq*xi2sq)) + w01[0]*sigma[i4]*(1-etasq*xi2sq);
          w0t[1] = w00[1]*(1 - sigma[i4]*(1-etasq*xi2sq)) + w01[1]*sigma[i4]*(1-etasq*xi2sq);
          w1t[0] = w10[0]*(1 - sigma[i4]*(1-etasq*xi2sq)) + w11[0]*sigma[i4]*(1-etasq*xi2sq);
          w1t[1] = w10[1]*(1 - sigma[i4]*(1-etasq*xi2sq)) + w11[1]*sigma[i4]*(1-etasq*xi2sq);
          if(nc==1){ //later, move this if statement outside of all loops
          a[0] = v1t[0] - v0t[0];
          a[1] = v1t[1] - v0t[1];
          b[0] = w0t[0] - w1t[0];
          b[1] = w0t[1] - w1t[1];
          }else if(nc==-1){
          a[0] = v0t[0] - v1t[0];
          a[1] = v0t[1] - v1t[1];
          b[0] = w1t[0] - w0t[0];
          b[1] = w1t[1] - w0t[1];
          }
          r[0] = eta[i3]*a[0] - eta[i3]*xi1[i1]*b[0];
          r[1] = eta[i3]*a[1] - eta[i3]*xi1[i1]*b[1];
          del = h*etasq*xi2sq;
          xdiff = w1t[0]-w0t[0];
          ydiff = w1t[1]-w0t[1];
          normb = sqrt(xdiff*xdiff + ydiff*ydiff);
          ny[0] = -ydiff/normb;
          ny[1] = xdiff/normb;
          norma = sqrt((v1t[0]-v0t[0])*(v1t[0]-v0t[0])+(v1t[1]-v0t[1])*(v1t[1]-v0t[1]));
          jac = norma*normb*h*h*2*eta[i3]*xi2[i2]*etasq*(1-etasq*xi2sq);
          kernel(r, del, ny, fcnVal);
          fcn1 = fcnVal[0]*jac;
          dfcn1 = fcnVal[1]*jac;

          //transf2: t = s + sigma = eta*xi2 + sigma, tau = sigma (same as transf1)

          r[0] = eta[i3]*xi1[i1]*a[0] - eta[i3]*b[0];
          r[1] = eta[i3]*xi1[i1]*a[1] - eta[i3]*b[1];
          del = h*etasq*xi2sq;
          jac = norma*normb*h*h*2*eta[i3]*xi2[i2]*etasq*(1-etasq*xi2sq);
          kernel(r, del, ny, fcnVal);
          fcn2 = fcnVal[0]*jac;
          dfcn2 = fcnVal[1]*jac;

          //transf3: t = s + sigma = eta + sigma, tau = sigma
          v0t[0] = v00[0]*(1 - etasq - sigma[i4]*(1-etasq)) + v01[0]*(etasq + sigma[i4]*(1-etasq)); 
          v0t[1] = v00[1]*(1 - etasq - sigma[i4]*(1-etasq)) + v01[1]*(etasq + sigma[i4]*(1-etasq));
          v1t[0] = v10[0]*(1 - etasq - sigma[i4]*(1-etasq)) + v11[0]*(etasq + sigma[i4]*(1-etasq));
          v1t[1] = v10[1]*(1 - etasq - sigma[i4]*(1-etasq)) + v11[1]*(etasq + sigma[i4]*(1-etasq));
          w0t[0] = w00[0]*(1 - sigma[i4]*(1-etasq)) + w01[0]*sigma[i4]*(1-etasq);
          w0t[1] = w00[1]*(1 - sigma[i4]*(1-etasq)) + w01[1]*sigma[i4]*(1-etasq);
          w1t[0] = w10[0]*(1 - sigma[i4]*(1-etasq)) + w11[0]*sigma[i4]*(1-etasq);
          w1t[1] = w10[1]*(1 - sigma[i4]*(1-etasq)) + w11[1]*sigma[i4]*(1-etasq);
          if(nc==1){ //later, move this if statement outside of all loops
          a[0] = v1t[0] - v0t[0];
          a[1] = v1t[1] - v0t[1];
          b[0] = w0t[0] - w1t[0];
          b[1] = w0t[1] - w1t[1];
          }else if(nc==-1){
          a[0] = v0t[0] - v1t[0];
          a[1] = v0t[1] - v1t[1];
          b[0] = w1t[0] - w0t[0];
          b[1] = w1t[1] - w0t[1];
          }
          r[0] = eta[i3]*xi1[i1]*a[0] - eta[i3]*xi2[i2]*b[0];
          r[1] = eta[i3]*xi1[i1]*a[1] - eta[i3]*xi2[i2]*b[1];
          del = h*etasq;
          xdiff = w1t[0]-w0t[0];
          ydiff = w1t[1]-w0t[1];
          normb = sqrt(xdiff*xdiff + ydiff*ydiff);
          ny[0] = -ydiff/normb;
          ny[1] = xdiff/normb;
          norma = sqrt((v1t[0]-v0t[0])*(v1t[0]-v0t[0])+(v1t[1]-v0t[1])*(v1t[1]-v0t[1]));
          jac = norma*normb*h*h*2*eta[i3]*etasq*(1-etasq);
          kernel(r, del, ny, fcnVal);
          fcn3 = fcnVal[0]*jac;
          dfcn3 = fcnVal[1]*jac;

          //weights, jacobian 
	  fcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn3 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn3 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];

	  //test functions
          xhat1 = eta[i3];
          yhat1 = eta[i3]*xi1[i1];
          xhat2 = eta[i3]*xi1[i1];
          yhat2 = eta[i3];
          xhat3 = eta[i3]*xi1[i1];
          yhat3 = eta[i3]*xi2[i2];
	  solution[0] += fcn1*yhat1*xhat1 + fcn2*yhat2*xhat2 +fcn3*yhat3*xhat3;
	  solution[1] += fcn1*(1 - yhat1)*xhat1 + fcn2*(1 - yhat2)*xhat2 + fcn3*(1 - yhat3)*xhat3;
	  solution[2] += fcn1*yhat1*(1 - xhat1) + fcn2*yhat2*(1 - xhat2) +fcn3*yhat3*(1 - xhat3);
	  solution[3] += fcn1*(1 - yhat1)*(1 - xhat1) + fcn2*(1 - yhat2)*(1 - xhat2) + fcn3*(1 - yhat3)*(1 - xhat3); 

	  solution[4] += dfcn1*yhat1*xhat1 + dfcn2*yhat2*xhat2 +dfcn3*yhat3*xhat3;
	  solution[5] += dfcn1*(1 - yhat1)*xhat1 + dfcn2*(1 - yhat2)*xhat2 + dfcn3*(1 - yhat3)*xhat3;
	  solution[6] += dfcn1*yhat1*(1 - xhat1) + dfcn2*yhat2*(1 - xhat2) +dfcn3*yhat3*(1 - xhat3);
	  solution[7] += dfcn1*(1 - yhat1)*(1 - xhat1) + dfcn2*(1 - yhat2)*(1 - xhat2) + dfcn3*(1 - yhat3)*(1 - xhat3); 

	}
      }
    }
  }
}



void intgrd0nc2(int p, double h, panel *pX0, panel *pX1, double xLeg[], double wLeg[], double solution[]) {
  double eta1, eta2, eta3, xi, t, tau, s, s2, z, x, y, nVeloc, norma, normb, jacobi, weights; 
  double dw0[2], db[2], r[2], ny[2], v[2], fcnVal[2], shape[4];
  double *v00, *v10, *v01, *v11, *w00, *w10, *w01, *w11;
  int i1, i2, i3, i;
 

  v00 = pX0->v0;
  v10 = pX0->v1;
  v01 = pX1->v0;
  v11 = pX1->v1;

  dw0[0] = v01[0] - v00[0];     //note w00=v00, etc
  dw0[1] = v01[1] - v00[1];
  db[0] = v11[0] - v10[0] - v01[0] + v00[0];
  db[1] = v11[1] - v10[1] - v01[1] + v00[1];

  for(i = 0; i < p; i++) { 
    xi = xLeg[i]; 
    for(i1 = 0; i1 < p; i1++) { 
      eta1 = xLeg[i1];
      for(i2 = 0; i2 < p; i2++) { 
        eta2 = xLeg[i2];
         for(i3 = 0; i3 < p; i3++) { 
           eta3 = xLeg[i3];
           weights = wLeg[i]*wLeg[i1]*wLeg[i2]*wLeg[i3];
           /* Transformation 1a */
           s = xi*eta1;
           s2 = s*s;
           tau = eta2*(1 - s2);
           t = s2 + tau;

           z = xi;
           y = eta3*(1-z);
           x = z + y;
           jacobi = xi*s*(1 - s2)*(1 - z);
           v[0] = dw0[0] + db[0]*y;
           v[1] = dw0[1] + db[1]*y;
           getVectors(v00, v10, v01, v11, v00, v10, v01, v11, t, tau, x, y, r, &norma, &normb, ny);
           jacobi *= norma*normb;
           nVeloc = (v[0]*ny[0] + v[1]*ny[1])/h;
 
           kernel(r, h*s2, ny, fcnVal);
           fcnVal[1] -= nVeloc*fcnVal[0];
           fcnVal[0] *= jacobi*weights;
           fcnVal[1] *= jacobi*weights;

           shape[0] = x*y;
           shape[1] = x*(1 - y);
           shape[2] = y*(1 - x);
           shape[3] = (1 - x)*(1 - y);
       
           solution[0] += fcnVal[0]*shape[0];
           solution[1] += fcnVal[0]*shape[1];
           solution[2] += fcnVal[0]*shape[2];
           solution[3] += fcnVal[0]*shape[3];
           solution[4] += fcnVal[1]*shape[0];
           solution[5] += fcnVal[1]*shape[1];
           solution[6] += fcnVal[1]*shape[2];
           solution[7] += fcnVal[1]*shape[3];

           /* Transformation 1b,   s,tau,t, are the same */
           z = -xi;
           y = -z + eta3*(1+z);
           x = z + y;
           jacobi = xi*s*(1 - s2)*(1 + z);
           v[0] = dw0[0] + db[0]*y;
           v[1] = dw0[1] + db[1]*y;
           getVectors(v00, v10, v01, v11, v00, v10, v01, v11, t, tau, x, y, r, &norma, &normb, ny);
           jacobi *= norma*normb;
           nVeloc = (v[0]*ny[0] + v[1]*ny[1])/h;
 
           kernel(r, h*s2, ny, fcnVal);
           fcnVal[1] -= nVeloc*fcnVal[0];
           fcnVal[0] *= jacobi*weights;
           fcnVal[1] *= jacobi*weights;

           shape[0] = x*y;
           shape[1] = x*(1 - y);
           shape[2] = y*(1 - x);
           shape[3] = (1 - x)*(1 - y);
       
           solution[0] += fcnVal[0]*shape[0];
           solution[1] += fcnVal[0]*shape[1];
           solution[2] += fcnVal[0]*shape[2];
           solution[3] += fcnVal[0]*shape[3];
           solution[4] += fcnVal[1]*shape[0];
           solution[5] += fcnVal[1]*shape[1];
           solution[6] += fcnVal[1]*shape[2];
           solution[7] += fcnVal[1]*shape[3];

           /* Transformation 2a */
           s = xi;
           s2 = s*s;
           tau = eta2*(1 - s2);
           t = s2 + tau;

           z = eta1*xi;
           y = eta3*(1-z);
           x = z + y;
           jacobi = xi*s*(1 - s2)*(1 - z);
           v[0] = dw0[0] + db[0]*y;
           v[1] = dw0[1] + db[1]*y;
           getVectors(v00, v10, v01, v11, v00, v10, v01, v11, t, tau, x, y, r, &norma, &normb, ny);
           jacobi *= norma*normb;
           nVeloc = (v[0]*ny[0] + v[1]*ny[1])/h;
 
           kernel(r, h*s2, ny, fcnVal);
           fcnVal[1] -= nVeloc*fcnVal[0];
           fcnVal[0] *= jacobi*weights;
           fcnVal[1] *= jacobi*weights;

           shape[0] = x*y;
           shape[1] = x*(1 - y);
           shape[2] = y*(1-x);
           shape[3] = (1 - x)*(1 - y);
       
           solution[0] += fcnVal[0]*shape[0];
           solution[1] += fcnVal[0]*shape[1];
           solution[2] += fcnVal[0]*shape[2];
           solution[3] += fcnVal[0]*shape[3];
           solution[4] += fcnVal[1]*shape[0];
           solution[5] += fcnVal[1]*shape[1];
           solution[6] += fcnVal[1]*shape[2];
           solution[7] += fcnVal[1]*shape[3];

           /* Transformation 2b,   s,tau,t, are the same */
           z = -eta1*xi;
           y = -z + eta3*(1+z);
           x = z + y;
           jacobi = xi*s*(1 - s2)*(1 + z);
           v[0] = dw0[0] + db[0]*y;
           v[1] = dw0[1] + db[1]*y;
           getVectors(v00, v10, v01, v11, v00, v10, v01, v11, t, tau, x, y, r, &norma, &normb, ny);
           jacobi *= norma*normb;
           nVeloc = (v[0]*ny[0] + v[1]*ny[1])/h;
 
           kernel(r, h*s2, ny, fcnVal);
           fcnVal[1] -= nVeloc*fcnVal[0];
           fcnVal[0] *= jacobi*weights;
           fcnVal[1] *= jacobi*weights;

           shape[0] = x*y;
           shape[1] = x*(1 - y);
           shape[2] = y*(1 - x);
           shape[3] = (1 - x)*(1 - y);
       
           solution[0] += fcnVal[0]*shape[0];
           solution[1] += fcnVal[0]*shape[1];
           solution[2] += fcnVal[0]*shape[2];
           solution[3] += fcnVal[0]*shape[3];
           solution[4] += fcnVal[1]*shape[0];
           solution[5] += fcnVal[1]*shape[1];
           solution[6] += fcnVal[1]*shape[2];
           solution[7] += fcnVal[1]*shape[3];
         }
      }
    }
  }
  for (i = 0; i < 8; i++) solution[i] *= 2*h*h; 
}





void panelIA(int d, int p, double h, panel *pX0, panel *pX1, panel *pY0, panel *pY1, double xLeg[], double wLeg[], double solution[]) {

  double *v0, *v1,  a[2], *w0, *w1, b[2];
  int nc;
  nc = nrcommonvtx(pX1,pY1);

  memset(solution, 0.0, 8*sizeof(double));

  if(d>1){
    intgrnonsing(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution);
  }else if(d==1){
    if(nc==2){//p and q are the same
      intgrd1nc2(d, p, h, pY0, pX0, pX1, xLeg, wLeg, solution);
    } else if(nc==1) {
      intgrd1nc1(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution, nc);
    } else if(nc==-1) {
      intgrd1nc1(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution, nc);
    } else if(nc==0){   
      intgrnonsing(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution);
    }
  }else if(d==0){
    if(nc==2){
      intgrd0nc2(p, h, pX0, pX1, xLeg, wLeg, solution);   
    } else if(nc==1) {
      intgrd0nc1(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution, nc);
    } else if(nc==-1) {
      intgrd0nc1(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution, nc);
    } else if(nc==0){ 
      intgrnonsing(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution);
    }    
  }
}
