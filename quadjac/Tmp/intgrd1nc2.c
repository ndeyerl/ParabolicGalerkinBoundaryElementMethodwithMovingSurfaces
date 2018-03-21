void intgrd1nc2(int d, int p, double h, double *v0, double *v1, double *w0, double *w1, double a[], double b[], double xLeg[], double wLeg[], double solution[]) {
  double *xi1, *xi2, *eta, *omega, norma, normb, rho, eta2, eta3, xi22, xi24, xi12, del1, del2, fcn1, fcn2, fcn3, xhat1a, xhat1b, xhat2a, xhat2b, xhat3a, xhat3b, yhat1a, yhat1b, yhat2a, yhat2b, yhat3a, yhat3b;
  int i1, i2, i3, i4;
  eta = xLeg;
  xi1 = xLeg;
  xi2 = xLeg;
  omega = xLeg; 
  norma = sqrt(a[0]*a[0] + a[1]*a[1]);
  normb = sqrt(b[0]*b[0] + b[1]*b[1]);
  rho = norma*norma/(4*h);
  for(i1 = 0; i1 < p; i1++) { //xi1 loop
    for(i2 = 0; i2 < p; i2++) { //xi2 loop
      for(i3 = 0; i3 < p; i3++) { //eta loop
	for(i4 = 0; i4 < p; i4++) { //omega loop
            eta2 = eta[i3]*eta[i3];
            eta3 = eta2*eta[i3];
            xi22 = xi2[i2]*xi2[i2];
            xi24 = xi22*xi22;
            xi12 = xi1[i1]*xi1[i1];
            del1 = 1+xi12;
            del2 = xi12 + xi22;
            fcn1 = eta3*xi1[i1]*xi2[i2]*exp(-rho*eta2*xi24/del1)*(1-eta2*xi22)/del1;
            fcn2 = fcn1;
            fcn3 = eta3*xi1[i1]*xi2[i2]*exp(-rho*eta2/del2)*(1-eta2)/del2;

fcn2=fcn1=0;
          //weights, jacobian 
	  fcn1 *= 8*h*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn2 *= 8*h*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn3 *= 8*h*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];

	  //test functions
          xhat1a = eta2*xi22 + omega[i4]*(1-eta2*xi22);
          yhat1a = omega[i4]*(1-eta2*xi22);
          xhat1b = 2*eta2*xi22 + omega[i4]*(1-eta2*xi22);
          yhat1b = omega[i4]*(1-eta2*xi22) + eta2*xi22;
          xhat2a = xhat1a;
          yhat2a = xhat1a;
          xhat2b = xhat1b;
          yhat2b = yhat1b;
          xhat3a = eta2 + omega[i4]*(1-eta2);
          yhat3a = omega[i4]*(1-eta2);
          xhat3b = 2*eta2 + omega[i4]*(1-eta2);
          yhat3b = omega[i4]*(1-eta2) + eta2;
 

	  solution[0] += fcn1*yhat1a*xhat1a + fcn1*yhat1b*xhat1b + fcn2*yhat2a*xhat2a + fcn2*yhat2b*xhat2b + fcn3*yhat3a*xhat3a + fcn3*yhat3b*xhat3b;

	  solution[1] += fcn1*(1 - yhat1a)*xhat1a + fcn1*(1 - yhat1b)*xhat1b + fcn2*(1 - yhat2a)*xhat2a + fcn2*(1 - yhat2b)*xhat2b + fcn3*(1 - yhat3a)*xhat3a + fcn3*(1 - yhat3b)*xhat3b;

	  solution[2] += fcn1*yhat1a*(1 - xhat1a) + fcn1*yhat1b*(1 - xhat1b) + fcn2*yhat2a*(1 - xhat2a) + fcn2*yhat2b*(1 - xhat2b) + fcn3*yhat3a*(1 - xhat3a) + fcn3*yhat3b*(1 - xhat3b);

	  solution[3] += fcn1*(1 - yhat1a)*(1 - xhat1a) + fcn1*(1 - yhat1b)*(1 - xhat1b) + fcn2*(1 - yhat2a)*(1 - xhat2a) + fcn2*(1 - yhat2b)*(1 - xhat2b) + fcn3*(1 - yhat3a)*(1 - xhat3a) + fcn3*(1 - yhat3b)*(1 - xhat3b);

	}
      }
    }
  }
}





void intgrd1nc2(int d, int p, double h, double *v0, double *v1, double *w0, double *w1, double a[], double b[], double xLeg[], double wLeg[], double solution[]) {
  double *g1, *g2, *f, *omega, norma, normb, rho, f2, f4, f5, f6, f7, g12, g13, g22, g23, g24, del1, del2, del3, del4, del5, del6, fcn1, fcn2, fcn3, fcn4, fcn5, fcn6, fcn7, fcn8, fcn9, xhat1a, xhat1b, xhat2a, xhat2b, xhat3a, xhat3b, xhat4a, xhat4b, xhat5a, xhat5b, xhat6a, xhat6b, xhat7a, xhat7b, xhat8a, xhat8b, xhat9a, xhat9b, yhat1a, yhat1b, yhat2a, yhat2b, yhat3a, yhat3b, yhat4a, yhat4b, yhat5a, yhat5b, yhat6a, yhat6b, yhat7a, yhat7b, yhat8a, yhat8b, yhat9a, yhat9b;
  int i1, i2, i3, i4;
  f = xLeg;
  g1 = xLeg;
  g2 = xLeg;
  omega = xLeg; 
  norma = sqrt(a[0]*a[0] + a[1]*a[1]);
  normb = sqrt(b[0]*b[0] + b[1]*b[1]);
  rho = norma*norma/(4*h);
  for(i1 = 0; i1 < p; i1++) { //f loop
    for(i2 = 0; i2 < p; i2++) { //g1 loop
      for(i3 = 0; i3 < p; i3++) { //g2 loop
	for(i4 = 0; i4 < p; i4++) { //omega loop
            f2 = f[i1]*f[i1];
            f4 = f2*f2;
            f5 = f4*f[i1];
            f6 = f5*f[i1];
            f7 = f6*f[i1];
            g12 = g1[i2]*g1[i2];
            g13 = g12*g1[i2];
            g22 = g2[i3]*g2[i3];
            g24 = g22*g22;
            del1 = 1 + f2*g12;
            del2 = g12 + g22;
            del3 = 1 + f2;
            del4 = 1 + g22;
            del5 = 1+f2*g22;
            del6 = g22 + 1;
            fcn1 = f7*g1[i2]*g2[i3]*exp(-rho*(f6*g24/del1))*(1-f4*g22)/del1;
            fcn2 = fcn1;
            fcn3 = f5*g1[i2]*g2[i3]*exp(-rho/del2)*(1-f2)/del2;
            fcn4 = f7*g13*g2[i3]*exp(-rho*(f6*g12*g24/del3))*(1-f4*g12*g22)/del3;
            fcn5 = fcn4;
            fcn6 = f5*g13*g2[i3]*exp(-rho*g12/del4)*(1-f2*g12)/del4;
            fcn7 = f7*g13*g2[i3]*exp(-rho*(f6*g12/del5))*(1-f4*g12)/del5;
            fcn8 = fcn7;
            fcn9 = f5*g13*g2[i3]*exp(-rho*g12/del6)*(1-f2*g12)/del6;
          //weights, jacobian 
	  fcn1 *= 8*h*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn2 *= 8*h*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn3 *= 8*h*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  fcn4 *= 8*h*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn5 *= 8*h*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn6 *= 8*h*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  fcn7 *= 8*h*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn8 *= 8*h*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn9 *= 8*h*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  //test functions
          xhat1a = f4*g22 + omega[i4]*(1-f4*g22);
          yhat1a = omega[i4]*(1-f4*g22);
          xhat1b = 2*f4*g22 + omega[i4]*(1-f4*g22);
          yhat1b = omega[i4]*(1-f4*g22) + f4*g22;
          xhat2a = xhat1a;
          yhat2a = xhat1a;
          xhat2b = xhat1b;
          yhat2b = yhat1b;
          xhat3a = f2 + omega[i4]*(1-f2);
          yhat3a = omega[i4]*(1-f2);
          xhat3b = 2*f2 + omega[i4]*(1-f2);
          yhat3b = omega[i4]*(1-f2) + f2;
          xhat4a = f4*g12*g22 + omega[i4]*(1-f4*g12*g22);
          yhat4a = omega[i4]*(1-f4*g12*g22);
          xhat4b = 2*f4*g12*g22 + omega[i4]*(1-f4*g12*g22);
          yhat4b = omega[i4]*(1-f4*g12*g22) + f4*g12*g22;
          xhat5a = xhat4a;
          yhat5a = xhat4a;
          xhat5b = xhat4b;
          yhat5b = yhat4b;
          xhat6a = f2*g12 + omega[i4]*(1-f2*g12);
          yhat6a = omega[i4]*(1-f2*g12);
          xhat6b = 2*f2*g12 + omega[i4]*(1-f2*g12);
          yhat6b = omega[i4]*(1-f2*g12) + f2*g12;
          xhat7a = f4*g12 + omega[i4]*(1-f4*g12);
          yhat7a = omega[i4]*(1-f4*g12);
          xhat7b = 2*f4*g12 + omega[i4]*(1-f4*g12);
          yhat7b = omega[i4]*(1-f4*g12) + f4*g12;
          xhat8a = xhat7a;
          yhat8a = xhat7a;
          xhat8b = xhat7b;
          yhat8b = yhat7b;   
          xhat9a = xhat6a;
          yhat9a = xhat6a;
          xhat9b = xhat6b;
          yhat9b = yhat6b;     

	  solution[0] += fcn1*yhat1a*xhat1a + fcn1*yhat1b*xhat1b + fcn2*yhat2a*xhat2a + fcn2*yhat2b*xhat2b + fcn3*yhat3a*xhat3a + fcn3*yhat3b*xhat3b + fcn4*yhat4a*xhat4a + fcn4*yhat4b*xhat4b + fcn5*yhat5a*xhat5a + fcn5*yhat5b*xhat5b + fcn6*yhat6a*xhat6a + fcn6*yhat6b*xhat6b + fcn7*yhat7a*xhat7a + fcn7*yhat7b*xhat7b + fcn8*yhat8a*xhat8a + fcn8*yhat8b*xhat8b + fcn9*yhat9a*xhat9a + fcn9*yhat9b*xhat9b;

	  solution[1] += fcn1*(1 - yhat1a)*xhat1a + fcn1*(1 - yhat1b)*xhat1b + fcn2*(1 - yhat2a)*xhat2a + fcn2*(1 - yhat2b)*xhat2b + fcn3*(1 - yhat3a)*xhat3a + fcn3*(1 - yhat3b)*xhat3b + fcn4*(1 - yhat4a)*xhat4a + fcn4*(1 - yhat4b)*xhat4b + fcn5*(1 - yhat5a)*xhat5a + fcn5*(1 - yhat5b)*xhat5b + fcn6*(1 - yhat6a)*xhat6a + fcn6*(1 - yhat6b)*xhat6b + fcn7*(1 - yhat7a)*xhat7a + fcn7*(1 - yhat7b)*xhat7b + fcn8*(1 - yhat8a)*xhat8a + fcn8*(1 - yhat8b)*xhat8b + fcn9*(1 - yhat9a)*xhat9a + fcn9*(1 - yhat9b)*xhat9b;

	  solution[2] += fcn1*yhat1a*(1 - xhat1a) + fcn1*yhat1b*(1 - xhat1b) + fcn2*yhat2a*(1 - xhat2a) + fcn2*yhat2b*(1 - xhat2b) + fcn3*yhat3a*(1 - xhat3a) + fcn3*yhat3b*(1 - xhat3b) + fcn4*yhat4a*(1 - xhat4a) + fcn4*yhat4b*(1 - xhat4b) + fcn5*yhat5a*(1 - xhat5a) + fcn5*yhat5b*(1 - xhat5b) + fcn6*yhat6a*(1 - xhat6a) + fcn6*yhat6b*(1 - xhat6b) + fcn7*yhat7a*(1 - xhat7a) + fcn7*yhat7b*(1 - xhat7b) + fcn8*yhat8a*(1 - xhat8a) + fcn8*yhat8b*(1 - xhat8b) + fcn9*yhat9a*(1 - xhat9a) + fcn9*yhat9b*(1 - xhat9b);

	  solution[3] += fcn1*(1 - yhat1a)*(1 - xhat1a) + fcn1*(1 - yhat1b)*(1 - xhat1b) + fcn2*(1 - yhat2a)*(1 - xhat2a) + fcn2*(1 - yhat2b)*(1 - xhat2b) + fcn3*(1 - yhat3a)*(1 - xhat3a) + fcn3*(1 - yhat3b)*(1 - xhat3b) + fcn4*(1 - yhat4a)*(1 - xhat4a) + fcn4*(1 - yhat4b)*(1 - xhat4b) + fcn5*(1 - yhat5a)*(1 - xhat5a) + fcn5*(1 - yhat5b)*(1 - xhat5b) + fcn6*(1 - yhat6a)*(1 - xhat6a) + fcn6*(1 - yhat6b)*(1 - xhat6b) + fcn7*(1 - yhat7a)*(1 - xhat7a) + fcn7*(1 - yhat7b)*(1 - xhat7b) + fcn8*(1 - yhat8a)*(1 - xhat8a) + fcn8*(1 - yhat8b)*(1 - xhat8b) + fcn9*(1 - yhat9a)*(1 - xhat9a) + fcn9*(1 - yhat9b)*(1 - xhat9b); 

	}
      }
    }
  }
}


