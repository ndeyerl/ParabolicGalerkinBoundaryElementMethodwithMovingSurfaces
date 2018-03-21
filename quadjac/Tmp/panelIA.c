//compile object file:
//gcc -c panelIA.c

#include <math.h>
#include <stdio.h>
#include "panelIA.h"

#define TOL 1e-12

int nrcommonvtx(panel *p1, panel *p2);
void intgrnonsing(int d, int p, double h, panel *pX0, panel *pX1, panel *pY0, panel *pY1, double xLeg[], double wLeg[], double solution[], int nc);
void intgrd1nc1(int d, int p, double h, panel *pX0, panel *pX1, panel *pY0, panel *pY1, double xLeg[], double wLeg[], double solution[], int nc);
void intgrd1nc2(int d, int p, double h, panel *pX0, panel *pX1, panel *pY0, panel *pY1, double xLeg[], double wLeg[], double solution[], int nc);
void intgrd0nc1(int d, int p, double h, panel *pX0, panel *pX1, panel *pY0, panel *pY1, double xLeg[], double wLeg[], double solution[], int nc);
void intgrd0nc2(int d, int p, double h, panel *pX0, panel *pX1, panel *pY0, panel *pY1, double xLeg[], double wLeg[], double solution[], int nc);
void kernel(double r[], double del, double ny[], double fcnVal[]);

void panelIA(int d, int p, double h, panel *pX0, panel *pX1, panel *pY0, panel *pY1, double xLeg[], double wLeg[], double solution[]) {

  double *v0, *v1,  a[2], *w0, *w1, b[2];
  int nc;
  solution[0] = solution[1] = solution[2] = solution[3] = solution[4] = solution[5] = solution[6] = solution[7] = 0.0;
  nc = nrcommonvtx(pX1,pY1);

  if(d>1){
    intgrnonsing(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution, nc);
  }else if(d==1){
    if(nc==2){//p and q are the same
      intgrd1nc2(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution, nc);
    } else if(nc==1) {
      intgrd1nc1(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution, nc);
    } else if(nc==-1) {
      intgrd1nc1(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution, nc);
    } else if(nc==0){   
      intgrnonsing(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution, nc);
    }
  }else if(d==0){
    if(nc==2){
      intgrd0nc2(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution, nc);   
    } else if(nc==1) {
      intgrd0nc1(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution, nc);
    } else if(nc==-1) {
      intgrd0nc1(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution, nc);
    } else if(nc==0){ 
      intgrnonsing(d, p, h, pX0, pX1, pY0, pY1, xLeg, wLeg, solution, nc);
    }    
  }
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


void intgrnonsing(int d, int p, double h, panel *pX0, panel *pX1, panel *pY0, panel *pY1, double xLeg[], double wLeg[], double solution[], int nc) {
  double *xhat, *yhat, *that, *tauhat, *v00, *v10, *v01, *v11, *w00, *w10, *w01, *w11, r[2], del, jac, fcnVal[2], dfcn, fcn, v0t[2], v1t[2], w0t[2], w1t[2], xdiff, ydiff, ny[2], norma, normb;
  int i1, i2, i3, i4;
  v00 = pX0->v0;
  v10 = pX0->v1;
  v01 = pX1->v0;
  v11 = pX1->v1;
  w00 = pY0->v0;
  w10 = pY0->v1;
  w01 = pY1->v0;
  w11 = pY1->v1;
  xhat = xLeg;
  yhat = xLeg;
  tauhat = xLeg;
  that = xLeg; 
  for(i1 = 0; i1 < p; i1++) { //xhat loop
    for(i2 = 0; i2 < p; i2++) { //yhat loop
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

	  r[0] = v1t[0]*(1 - xhat[i1]) + v0t[0]*xhat[i1] - 
	    w1t[0]*(1 - yhat[i2]) - w0t[0]*yhat[i2];
	  r[1] = v1t[1]*(1 - xhat[i1]) + v0t[1]*xhat[i1] -
	    w1t[1]*(1 - yhat[i2]) - w0t[1]*yhat[i2];

	  //heat kernel
          if(del < TOL){
            fcn = 0.0;
            dfcn = 0.0;
          }
          else{
            xdiff = w1t[0]-w0t[0]; 
            ydiff = w1t[1]-w0t[1];
            normb = sqrt(xdiff*xdiff + ydiff*ydiff);
            ny[0] = -ydiff/normb; 
            ny[1] = xdiff/normb; 
            kernel(r, del, ny, fcnVal);
            norma = sqrt((v1t[0]-v0t[0])*(v1t[0]-v0t[0])+(v1t[1]-v0t[1])*(v1t[1]-v0t[1]));
            jac = h*h*norma*normb;
	    fcn = fcnVal[0]*jac;
            dfcn = fcnVal[1]*jac;
          }

	  //weights
	  fcn *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
     
	  //test functions
	  solution[0] += fcn*yhat[i2]*xhat[i1];
	  solution[1] += fcn*(1 - yhat[i2])*xhat[i1];
	  solution[2] += fcn*yhat[i2]*(1 - xhat[i1]);
	  solution[3] += fcn*(1 - yhat[i2])*(1 - xhat[i1]);

	  solution[4] += dfcn*yhat[i2]*xhat[i1];
	  solution[5] += dfcn*(1 - yhat[i2])*xhat[i1];
	  solution[6] += dfcn*yhat[i2]*(1 - xhat[i1]);
	  solution[7] += dfcn*(1 - yhat[i2])*(1 - xhat[i1]);

	}
      }
    }
  }
}

void intgrd1nc1(int d, int p, double h, panel *pX0, panel *pX1, panel *pY0, panel *pY1, double xLeg[], double wLeg[], double solution[], int nc) {
  double *eta, *xi1, *xi2, *xi3, *omega,  etasq, xi1sq, xi2sq, del1, del2, del3, fcn1, fcn2, fcn3, fcn4, fcn, eta3, r1[2], r2[2], r3[2], rr1, rr2, rr3, r[2], fcnVal[2], dfcn1, dfcn2, dfcn3, dfcn4, rr, xhat1, yhat1, xhat2, yhat2, xhat3, yhat3, jac, *v00, *v10, *v01, *v11, *w00, *w10, *w01, *w11, v0t[2], v1t[2], w0t[2], w1t[2], xdiff, ydiff, ny[2], norma, normb, a[2], b[2];
  int i1, i2, i3, i4;
  v00 = pX0->v0;
  v10 = pX0->v1;
  v01 = pX1->v0;
  v11 = pX1->v1;
  w00 = pY0->v0;
  w10 = pY0->v1;
  w01 = pY1->v0;
  w11 = pY1->v1;
  eta = xLeg;
  xi1 = xLeg;
  xi2 = xLeg;
  xi3 = xLeg;
  omega = xLeg; 
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
          kernel(r1, del1, ny, fcnVal);
          fcn1 = fcnVal[0]*jac;
          dfcn1 = fcnVal[1]*jac;

          //t, tau have same transformation as fcn1
	  r2[0] = eta[i4]*xi1[i1]*a[0] - eta[i4]*b[0];
	  r2[1] = eta[i4]*xi1[i1]*a[1] - eta[i4]*b[1];
          jac = norma*normb*4*h*h*eta[i4]*xi2[i2]*eta[i4]*xi3[i3]*eta3;
          kernel(r2, del1, ny, fcnVal);
          fcn2 = fcnVal[0]*jac;
          dfcn2 = fcnVal[1]*jac;

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
          kernel(r3, del2, ny, fcnVal);
	  fcn3 = fcnVal[0]*jac;
          dfcn3 = fcnVal[1]*jac;

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
          dfcn4 = fcnVal[1]*jac;

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

void intgrd1nc2(int d, int p, double h, panel *pX0, panel *pX1, panel *pY0, panel *pY1, double xLeg[], double wLeg[], double solution[], int nc) {
  double *eta, *xi1, *xi2, *omega,  *v00, *v10, *v01, *v11, *w00, *w10, *w01, *w11, etasq, xi1sq, xi2sq, del, fcn1, fcn2, fcn3, xhat1a, xhat1b, xhat2a, xhat2b, xhat3a, xhat3b, xhat9b, yhat1a, yhat1b, yhat2a, yhat2b, yhat3a, yhat3b, fcnVal[2], r[2], jac, dfcn1, dfcn1b, dfcn2, dfcn2b, dfcn3, dfcn3b, fcn1b, fcn2b, fcn3b, v0t[2], v1t[2], w0t[2], w1t[2], xdiff, ydiff, norm, ny[2], norma, normb, a[2];
  int i1, i2, i3, i4;
  v00 = pX0->v0;
  v10 = pX0->v1;
  v01 = pX1->v0;
  v11 = pX1->v1;
  w00 = pY0->v0;
  w10 = pY0->v1;
  w01 = pY1->v0;
  w11 = pY1->v1;
  eta = xLeg;
  xi1 = xLeg;
  xi2 = xLeg;
  omega = xLeg; 
  //rho = norma*norma/(4*h);
  for(i1 = 0; i1 < p; i1++) { //xi1 loop
    for(i2 = 0; i2 < p; i2++) { //xi2 loop
      for(i3 = 0; i3 < p; i3++) { //eta loop
	for(i4 = 0; i4 < p; i4++) { //omega loop
          etasq = eta[i3]*eta[i3];
          xi1sq = xi1[i1]*xi1[i1];
          xi2sq = xi2[i2]*xi2[i2];
//NOTE FOR SELF CASE: does it matter that in the transformation we let a = b and factored it out to transform further?  because now, having a=v0t-v1t and b=w0t-w1t means a=\b due to v0t=\w0t, etc by the transformations on t and tau
//transf1, t = eta, tau = eta*xi1
          v0t[0] = v00[0]*(1 - eta[i3]) + v01[0]*eta[i3]; 
          v0t[1] = v00[1]*(1 - eta[i3]) + v01[1]*eta[i3];
          v1t[0] = v10[0]*(1 - eta[i3]) + v11[0]*eta[i3];
          v1t[1] = v10[1]*(1 - eta[i3]) + v11[1]*eta[i3];
          w0t[0] = w00[0]*(1 - eta[i3]*xi1[i1]) + w01[0]*eta[i3]*xi1[i1];
          w0t[1] = w00[1]*(1 - eta[i3]*xi1[i1]) + w01[1]*eta[i3]*xi1[i1];
          w1t[0] = w10[0]*(1 - eta[i3]*xi1[i1]) + w11[0]*eta[i3]*xi1[i1];
          w1t[1] = w10[1]*(1 - eta[i3]*xi1[i1]) + w11[1]*eta[i3]*xi1[i1];
          a[0] = v0t[0] - v1t[0];
          a[1] = v0t[1] - v1t[1];

          r[0] = a[0]*eta[i3]*xi2[i2];
          r[1] = a[1]*eta[i3]*xi2[i2];
          del = h*etasq*(1 + xi1sq);

          xdiff = v1t[0]-v0t[0];
          ydiff = v1t[1]-v0t[1];
          norma = sqrt(xdiff*xdiff + ydiff*ydiff);
          ny[0] = -ydiff/norma;
          ny[1] = xdiff/norma;

          jac = norma*norma*h*h*4*eta[i3]*eta[i3]*xi1[i1]*etasq*(1-eta[i3]*xi2[i2]);
          kernel(r, del, ny, fcnVal);
          fcn1 = fcnVal[0]*jac;
          dfcn1 = fcnVal[1]*jac;

//transf1b, t = eta, tau = eta*xi1
          r[0] = -a[0]*eta[i3]*xi2[i2];
          r[1] = -a[1]*eta[i3]*xi2[i2];
          kernel(r, del, ny, fcnVal);
          fcn1b = fcnVal[0]*jac;
          dfcn1b = fcnVal[1]*jac;

//transf2a, t = eta*xi1, tau = eta
          v0t[0] = v00[0]*(1 - eta[i3]*xi1[i1]) + v01[0]*eta[i3]*xi1[i1]; 
          v0t[1] = v00[1]*(1 - eta[i3]*xi1[i1]) + v01[1]*eta[i3]*xi1[i1];
          v1t[0] = v10[0]*(1 - eta[i3]*xi1[i1]) + v11[0]*eta[i3]*xi1[i1];
          v1t[1] = v10[1]*(1 - eta[i3]*xi1[i1]) + v11[1]*eta[i3]*xi1[i1];
          w0t[0] = w00[0]*(1 - eta[i3]) + w01[0]*eta[i3];
          w0t[1] = w00[1]*(1 - eta[i3]) + w01[1]*eta[i3];
          w1t[0] = w10[0]*(1 - eta[i3]) + w11[0]*eta[i3];
          w1t[1] = w10[1]*(1 - eta[i3]) + w11[1]*eta[i3];
          a[0] = v0t[0] - v1t[0];
          a[1] = v0t[1] - v1t[1];

          r[0] = a[0]*eta[i3]*xi2[i2];
          r[1] = a[1]*eta[i3]*xi2[i2];
          del = h*etasq*(xi1sq + 1);

          xdiff = v1t[0]-v0t[0];
          ydiff = v1t[1]-v0t[1];
          norma = sqrt(xdiff*xdiff + ydiff*ydiff);
          ny[0] = -ydiff/norma;
          ny[1] = xdiff/norma;

          jac = norma*norma*h*h*4*eta[i3]*xi1[i1]*eta[i3]*etasq*(1-eta[i3]*xi2[i2]);          
          kernel(r, del, ny, fcnVal);
          fcn2 = fcnVal[0]*jac;
          dfcn2 = fcnVal[1]*jac;

//transf2b, t = eta*xi1, tau = eta
          r[0] = -a[0]*eta[i3]*xi2[i2];
          r[1] = -a[1]*eta[i3]*xi2[i2];
          kernel(r, del, ny, fcnVal);
          fcn2b = fcnVal[0]*jac;
          dfcn2b = fcnVal[1]*jac;

//transf3a, t = eta*xi1, tau = eta*xi2
          v0t[0] = v00[0]*(1 - eta[i3]*xi1[i1]) + v01[0]*eta[i3]*xi1[i1]; 
          v0t[1] = v00[1]*(1 - eta[i3]*xi1[i1]) + v01[1]*eta[i3]*xi1[i1];
          v1t[0] = v10[0]*(1 - eta[i3]*xi1[i1]) + v11[0]*eta[i3]*xi1[i1];
          v1t[1] = v10[1]*(1 - eta[i3]*xi1[i1]) + v11[1]*eta[i3]*xi1[i1];
          w0t[0] = w00[0]*(1 - eta[i3]*xi2[i2]) + w01[0]*eta[i3]*xi2[i2];
          w0t[1] = w00[1]*(1 - eta[i3]*xi2[i2]) + w01[1]*eta[i3]*xi2[i2];
          w1t[0] = w10[0]*(1 - eta[i3]*xi2[i2]) + w11[0]*eta[i3]*xi2[i2];
          w1t[1] = w10[1]*(1 - eta[i3]*xi2[i2]) + w11[1]*eta[i3]*xi2[i2];
          a[0] = v0t[0] - v1t[0];
          a[1] = v0t[1] - v1t[1];

          r[0] = a[0]*eta[i3];
          r[1] = a[1]*eta[i3];
          del = h*etasq*(xi1sq + xi2sq);

          xdiff = v1t[0]-v0t[0];
          ydiff = v1t[1]-v0t[1];
          norma = sqrt(xdiff*xdiff + ydiff*ydiff);
          ny[0] = -ydiff/norma;
          ny[1] = xdiff/norma;

          jac = norma*norma*h*h*4*eta[i4]*xi1[i1]*eta[i3]*xi2[i2]*etasq*(1-eta[i3]);
          kernel(r, del, ny, fcnVal);
          fcn3 = fcnVal[0]*jac;
          dfcn3 = fcnVal[1]*jac;

//transf3b, t = eta*xi1, tau = eta*xi2
          r[0] = -a[0]*eta[i3];
          r[1] = -a[1]*eta[i3];
          kernel(r, del, ny, fcnVal);
          fcn3b = fcnVal[0]*jac;
          dfcn3b = fcnVal[1]*jac;

          //weights, jacobian 
	  fcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn3 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  fcn1b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn2b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn3b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn3 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn1b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn2b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn3b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];

	  //test functions
          xhat1a = eta[i3]*xi2[i2] + omega[i4]*(1-eta[i3]*xi2[i2]);
          yhat1a = omega[i4]*(1-eta[i3]*xi2[i2]);
          xhat1b = omega[i4]*(1-eta[i3]*xi2[i2]);
          yhat1b = eta[i3]*xi2[i2] + omega[i4]*(1-eta[i3]*xi2[i2]);
          xhat2a = xhat1a;
          yhat2a = xhat1a;
          xhat2b = xhat1b;
          yhat2b = yhat1b;
          xhat3a = eta[i3] + omega[i4]*(1-eta[i3]);
          yhat3a = omega[i4]*(1-eta[i3]);
          xhat3b = omega[i4]*(1-eta[i3]);
          yhat3b = eta[i3] + omega[i4]*(1-eta[i3]);


	  solution[0] += fcn1*yhat1a*xhat1a + fcn1b*yhat1b*xhat1b + fcn2*yhat2a*xhat2a + fcn2b*yhat2b*xhat2b + fcn3*yhat3a*xhat3a + fcn3b*yhat3b*xhat3b;

	  solution[1] += fcn1*(1 - yhat1a)*xhat1a + fcn1b*(1 - yhat1b)*xhat1b + fcn2*(1 - yhat2a)*xhat2a + fcn2b*(1 - yhat2b)*xhat2b + fcn3*(1 - yhat3a)*xhat3a + fcn3b*(1 - yhat3b)*xhat3b;

	  solution[2] += fcn1*yhat1a*(1 - xhat1a) + fcn1b*yhat1b*(1 - xhat1b) + fcn2*yhat2a*(1 - xhat2a) + fcn2b*yhat2b*(1 - xhat2b) + fcn3*yhat3a*(1 - xhat3a) + fcn3b*yhat3b*(1 - xhat3b);

	  solution[3] += fcn1*(1 - yhat1a)*(1 - xhat1a) + fcn1b*(1 - yhat1b)*(1 - xhat1b) + fcn2*(1 - yhat2a)*(1 - xhat2a) + fcn2b*(1 - yhat2b)*(1 - xhat2b) + fcn3*(1 - yhat3a)*(1 - xhat3a) + fcn3b*(1 - yhat3b)*(1 - xhat3b); 


	  solution[4] += dfcn1*yhat1a*xhat1a + dfcn1b*yhat1b*xhat1b + dfcn2*yhat2a*xhat2a + dfcn2b*yhat2b*xhat2b + dfcn3*yhat3a*xhat3a + dfcn3b*yhat3b*xhat3b;

	  solution[5] += dfcn1*(1 - yhat1a)*xhat1a + dfcn1b*(1 - yhat1b)*xhat1b + dfcn2*(1 - yhat2a)*xhat2a + dfcn2b*(1 - yhat2b)*xhat2b + dfcn3*(1 - yhat3a)*xhat3a + dfcn3b*(1 - yhat3b)*xhat3b;

	  solution[6] += dfcn1*yhat1a*(1 - xhat1a) + dfcn1b*yhat1b*(1 - xhat1b) + dfcn2*yhat2a*(1 - xhat2a) + dfcn2b*yhat2b*(1 - xhat2b) + dfcn3*yhat3a*(1 - xhat3a) + dfcn3b*yhat3b*(1 - xhat3b);

	  solution[7] += dfcn1*(1 - yhat1a)*(1 - xhat1a) + dfcn1b*(1 - yhat1b)*(1 - xhat1b) + dfcn2*(1 - yhat2a)*(1 - xhat2a) + dfcn2b*(1 - yhat2b)*(1 - xhat2b) + dfcn3*(1 - yhat3a)*(1 - xhat3a) + dfcn3b*(1 - yhat3b)*(1 - xhat3b); 

	}
      }
    }
  }
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



void intgrd0nc2(int d, int p, double h, panel *pX0, panel *pX1, panel *pY0, panel *pY1, double xLeg[], double wLeg[], double solution[], int nc) {
  double *omega, *delta, *eta, *xi, *v00, *v10, *v01, *v11, *w00, *w10, *w01, *w11, fcnVal[2], r[2], del, ny[2], jac, fcn1, fcn1b, fcn2, fcn2b, dfcn1, dfcn1b, dfcn2, dfcn2b, xhat1a, xhat1b, xhat2a, xhat2b, yhat1a, yhat1b, yhat2a, yhat2b, norma, normb, xisq, etasq, v0t[2], v1t[2], w0t[2], w1t[2], a[2], b[2], xdiff, ydiff;
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
  delta = xLeg;
  eta = xLeg;
  xi = xLeg; 
  double rr, drr;
  for(i1 = 0; i1 < p; i1++) { //eta loop
    for(i2 = 0; i2 < p; i2++) { //xi loop
      for(i3 = 0; i3 < p; i3++) { //omega loop
	for(i4 = 0; i4 < p; i4++) { //delta loop
          etasq = eta[i1]*eta[i1];
          xisq = xi[i2]*xi[i2];

          //transf 1 t = eta*eta + delta*(1-eta*eta) tau = delta*(1-eta*eta)

          v0t[0] = v00[0]*(1 - etasq - delta[i4]*(1 - etasq)) + v01[0]*(etasq + delta[i4]*(1 - etasq)); 
          v0t[1] = v00[1]*(1 - etasq - delta[i4]*(1 - etasq)) + v01[1]*(etasq + delta[i4]*(1 - etasq));
          v1t[0] = v10[0]*(1 - etasq - delta[i4]*(1 - etasq)) + v11[0]*(etasq + delta[i4]*(1 - etasq));
          v1t[1] = v10[1]*(1 - etasq - delta[i4]*(1 - etasq)) + v11[1]*(etasq + delta[i4]*(1 - etasq));
          w0t[0] = w00[0]*(1 - delta[i4]*(1 - etasq)) + w01[0]*delta[i4]*(1 - etasq);
          w0t[1] = w00[1]*(1 - delta[i4]*(1 - etasq)) + w01[1]*delta[i4]*(1 - etasq);
          w1t[0] = w10[0]*(1 - delta[i4]*(1 - etasq)) + w11[0]*delta[i4]*(1 - etasq);
          w1t[1] = w10[1]*(1 - delta[i4]*(1 - etasq)) + w11[1]*delta[i4]*(1 - etasq);

          a[0] = v0t[0]-v1t[0];
          a[1] = v0t[1]-v1t[1];
          b[0] = w0t[0]-w1t[0];
          b[1] = w0t[1]-w1t[1];


          r[0] = a[0]*eta[i1]*xi[i2];
          r[1] = a[1]*eta[i1]*xi[i2];
          del = h*etasq;
          xdiff = v1t[0]-v0t[0];
          ydiff = v1t[1]-v0t[1];
          norma = sqrt(xdiff*xdiff + ydiff*ydiff);
          ny[0] = -ydiff/norma;
          ny[1] = xdiff/norma;  

          jac = norma*norma*h*h*2*eta[i1]*(1-eta[i1]*xi[i2])*(1-etasq)*eta[i1];
          kernel(r, del, ny, fcnVal);
          fcn1 = fcnVal[0]*jac;
          dfcn1 = fcnVal[1]*jac;

          //transf 1b: t, tau same as transf 1
          r[0] = -a[0]*eta[i1]*xi[i2];
          r[1] = -a[1]*eta[i1]*xi[i2];
          kernel(r, del, ny, fcnVal);
          fcn1b = fcnVal[0]*jac;
          dfcn1b = fcnVal[1]*jac;

          //transf2: t = eta*eta*xi*xi + delta*(1-eta*eta*xi*xi) tau = delta(1-etasq*xisq)
          v0t[0] = v00[0]*(1 - etasq*xisq - delta[i4]*(1 - etasq*xisq)) + v01[0]*(etasq*xisq + delta[i4]*(1 - etasq*xisq)); 
          v0t[1] = v00[1]*(1 - etasq*xisq - delta[i4]*(1 - etasq*xisq)) + v01[1]*(etasq*xisq + delta[i4]*(1 - etasq*xisq));
          v1t[0] = v10[0]*(1 - etasq*xisq - delta[i4]*(1 - etasq*xisq)) + v11[0]*(etasq*xisq + delta[i4]*(1 - etasq*xisq));
          v1t[1] = v10[1]*(1 - etasq*xisq - delta[i4]*(1 - etasq*xisq)) + v11[1]*(etasq*xisq + delta[i4]*(1 - etasq*xisq));
          w0t[0] = w00[0]*(1 - delta[i4]*(1 - etasq*xisq)) + w01[0]*delta[i4]*(1 - etasq*xisq);
          w0t[1] = w00[1]*(1 - delta[i4]*(1 - etasq*xisq)) + w01[1]*delta[i4]*(1 - etasq*xisq);
          w1t[0] = w10[0]*(1 - delta[i4]*(1 - etasq*xisq)) + w11[0]*delta[i4]*(1 - etasq*xisq);
          w1t[1] = w10[1]*(1 - delta[i4]*(1 - etasq*xisq)) + w11[1]*delta[i4]*(1 - etasq*xisq);

          a[0] = v0t[0]-v1t[0];
          a[1] = v0t[1]-v1t[1];
          b[0] = w0t[0]-w1t[0];
          b[1] = w0t[1]-w1t[1];


          r[0] = a[0]*eta[i1];
          r[1] = a[1]*eta[i1];
          del = h*etasq*xisq;

          xdiff = v1t[0]-v0t[0];
          ydiff = v1t[1]-v0t[1];
          norma = sqrt(xdiff*xdiff + ydiff*ydiff);
          ny[0] = -ydiff/norma;
          ny[1] = xdiff/norma;

          jac = norma*norma*h*h*2*eta[i1]*xi[i2]*(1-eta[i1])*(1-etasq*xisq)*eta[i1];
          kernel(r, del, ny, fcnVal);
          fcn2 = fcnVal[0]*jac;
          dfcn2 = fcnVal[1]*jac;

          //transf2b: same as transf 2
          r[0] = -a[0]*eta[i1];
          r[1] = -a[1]*eta[i1];
          kernel(r, del, ny, fcnVal);
          fcn2b = fcnVal[0]*jac;
          dfcn2b = fcnVal[1]*jac;

          //weights, jacobian 
	  fcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn1b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  fcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn2b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn1b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn2b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];

	  //test functions
          xhat1a = eta[i1]*xi[i2] + omega[i3]*(1-eta[i1]*xi[i2]);
          yhat1a = omega[i3]*(1-eta[i1]*xi[i2]);
          xhat1b = omega[i3]*(1-eta[i1]*xi[i2]);
          yhat1b = eta[i1]*xi[i2] + omega[i3]*(1-eta[i1]*xi[i2]);

          xhat2a = eta[i1] + omega[i3]*(1-eta[i1]);
          yhat2a = omega[i3]*(1-eta[i1]);
          xhat2b = omega[i3]*(1-eta[i1]);
          yhat2b = eta[i1] + omega[i3]*(1-eta[i1]);


	  solution[0] += fcn1*yhat1a*xhat1a + fcn1b*yhat1b*xhat1b + fcn2*yhat2a*xhat2a + fcn2b*yhat2b*xhat2b;
	  solution[1] += fcn1*(1 - yhat1a)*xhat1a + fcn1b*(1 - yhat1b)*xhat1b + fcn2*(1 - yhat2a)*xhat2a + fcn2b*(1 - yhat2b)*xhat2b;
	  solution[2] += fcn1*yhat1a*(1 - xhat1a) + fcn1b*yhat1b*(1 - xhat1b) + fcn2*yhat2a*(1 - xhat2a) + fcn2b*yhat2b*(1 - xhat2b);
	  solution[3] += fcn1*(1 - yhat1a)*(1 - xhat1a) + fcn1b*(1 - yhat1b)*(1 - xhat1b) + fcn2*(1 - yhat2a)*(1 - xhat2a) + fcn2b*(1 - yhat2b)*(1 - xhat2b); 

	  solution[4] += dfcn1*yhat1a*xhat1a + dfcn1b*yhat1b*xhat1b + dfcn2*yhat2a*xhat2a + dfcn2b*yhat2b*xhat2b;
	  solution[5] += dfcn1*(1 - yhat1a)*xhat1a + dfcn1b*(1 - yhat1b)*xhat1b + dfcn2*(1 - yhat2a)*xhat2a + dfcn2b*(1 - yhat2b)*xhat2b;
	  solution[6] += dfcn1*yhat1a*(1 - xhat1a) + dfcn1b*yhat1b*(1 - xhat1b) + dfcn2*yhat2a*(1 - xhat2a) + dfcn2b*yhat2b*(1 - xhat2b);
	  solution[7] += dfcn1*(1 - yhat1a)*(1 - xhat1a) + dfcn1b*(1 - yhat1b)*(1 - xhat1b) + dfcn2*(1 - yhat2a)*(1 - xhat2a) + dfcn2b*(1 - yhat2b)*(1 - xhat2b); 

	}
      }
    }
  }
}

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
