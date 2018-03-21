void intgrd0nc2(int d, int p, double h, double *v0, double *v1, double *w0, double *w1, double a[], double b[], double xLeg[], double wLeg[], double solution[]) {
  double *omega, *delta, *eta, *xi, eta2, eta4, xi2, xi4, fcn1, fcn2, xhat1a, xhat1b, xhat2a, xhat2b, yhat1a, yhat1b, yhat2a, yhat2b, norma, rho;
  int i1, i2, i3, i4;
  omega = xLeg;
  delta = xLeg;
  eta = xLeg;
  xi = xLeg; 
  norma = sqrt(a[0]*a[0] + a[1]*a[1]);
  rho = norma*norma/(4*h);
  for(i1 = 0; i1 < p; i1++) { //eta loop
    for(i2 = 0; i2 < p; i2++) { //xi loop
      for(i3 = 0; i3 < p; i3++) { //omega loop
	for(i4 = 0; i4 < p; i4++) { //delta loop
          eta2 = eta[i1]*eta[i1];
          eta4 = eta2*eta2;
          xi2 = xi[i2]*xi[i2];
          xi4 = xi2*xi2;
          fcn1 = eta[i1]*xi[i2]*exp(-rho*xi2*eta4)*(1-eta2)*(1-eta2*xi2);
          fcn2 = eta[i1]*exp(-rho*eta2/xi2)*(1-eta2*xi2)*(1-eta2)/xi[i2];
//fcn1=fcn3=fcn2=0;
//fcn2=fcn4=0;
          //weights, jacobian 
	  fcn1 *= 4*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4]/h;
          fcn2 *= 4*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4]/h;
	  //test functions
          xhat1a = eta2*xi2 + omega[i3]*(1-eta2*xi2);
          yhat1a = omega[i3]*(1-eta2*xi2);
          xhat1b = 2*eta2*xi2 + omega[i3]*(1-eta2*xi2);
          yhat1b = omega[i3]*(1-eta2*xi2) + eta2*xi2;

          xhat2a = eta2 + omega[i3]*(1-eta2);
          yhat2a = omega[i3]*(1-eta2);
          xhat2b = 2*eta2 + omega[i3]*(1-eta2);
          yhat2b = omega[i3]*(1-eta2) + eta2;


	  solution[0] += fcn1*yhat1a*xhat1a + fcn1*yhat1b*xhat1b + fcn2*yhat2a*xhat2a + fcn2*yhat2b*xhat2b;
	  solution[1] += fcn1*(1 - yhat1a)*xhat1a + fcn1*(1 - yhat1b)*xhat1b + fcn2*(1 - yhat2a)*xhat2a + fcn2*(1 - yhat2b)*xhat2b;
	  solution[2] += fcn1*yhat1a*(1 - xhat1a) + fcn1*yhat1b*(1 - xhat1b) + fcn2*yhat2a*(1 - xhat2a) + fcn2*yhat2b*(1 - xhat2b);
	  solution[3] += fcn1*(1 - yhat1a)*(1 - xhat1a) + fcn1*(1 - yhat1b)*(1 - xhat1b) + fcn2*(1 - yhat2a)*(1 - xhat2a) + fcn2*(1 - yhat2b)*(1 - xhat2b); 
	}
      }
    }
  }
}






void intgrd0nc2(int d, int p, double h, double *v0, double *v1, double *w0, double *w1, double a[], double b[], double xLeg[], double wLeg[], double solution[]) {
  double *omega, *delta, *f, *g,  f2, f3, f4, f6, g2, g4, fcn1, fcn2, fcn3, fcn4, xhat1a, xhat1b, xhat2a, xhat2b, xhat3a, xhat3b, xhat4a, xhat4b, yhat1a, yhat1b, yhat2a, yhat2b, yhat3a, yhat3b, yhat4a, yhat4b, norma, rho;
  int i1, i2, i3, i4;
  omega = xLeg;
  delta = xLeg;
  f = xLeg;
  g = xLeg; 
  norma = sqrt(a[0]*a[0] + a[1]*a[1]);
  rho = norma*norma/(4*h);
  for(i1 = 0; i1 < p; i1++) { //f loop
    for(i2 = 0; i2 < p; i2++) { //g loop
      for(i3 = 0; i3 < p; i3++) { //omega loop
	for(i4 = 0; i4 < p; i4++) { //delta loop
          f2 = f[i1]*f[i1];
          f3 = f2*f[i1];
          f4 = f3*f[i1];
          f6 = f3*f3; 
          g2 = g[i2]*g[i2];
          g4 = g2*g2;
          fcn1 = f3*g[i2]*exp(-rho*f6*g4)*(1-f2)*(1-f4*g2);
          fcn2 = f[i1]*exp(-rho*(1/g2))*(1-f4*g2)*(1-f2)/g[i2];
          fcn3 = f3*g[i2]*exp(-rho*f6*g2)*(1-f2*g2)*(1-f4*g2);
          fcn4 = f[i1]*g[i2]*exp(-rho*g2)*(1-f4*g2)*(1-f2*g2);
//fcn1=fcn3=fcn2=0;
//fcn2=fcn4=0;
          //weights, jacobian 
	  fcn1 *= 4*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4]/h;
          fcn2 *= 4*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4]/h;
          fcn3 *= 4*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4]/h;
          fcn4 *= 4*norma*norma*wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4]/h;
	  //test functions
          xhat1a = f4*g2 + omega[i3]*(1-f4*g2);
          yhat1a = omega[i3]*(1-f4*g2);
          xhat1b = 2*f4*g2 + omega[i3]*(1-f4*g2);
          yhat1b = omega[i3]*(1-f4*g2) + f4*g2;

          xhat2a = f2 + omega[i3]*(1-f2);
          yhat2a = omega[i3]*(1-f2);
          xhat2b = 2*f2 + omega[i3]*(1-f2);
          yhat2b = omega[i3]*(1-f2) + f2;

          xhat3a = xhat1a;
          yhat3a = yhat1a;
          xhat3b = xhat1b;
          yhat3b = yhat1b;

          xhat4a = f2*g2 + omega[i3]*(1-f2*g2);
          yhat4a = omega[i3]*(1-f2*g2);
          xhat4b = 2*f2*g2 + omega[i3]*(1-f2*g2);
          yhat4b = omega[i3]*(1-f2*g2) + f2*g2;

	  solution[0] += fcn1*yhat1a*xhat1a + fcn1*yhat1b*xhat1b + fcn2*yhat2a*xhat2a + fcn2*yhat2b*xhat2b + fcn3*yhat3a*xhat3a + fcn3*yhat3b*xhat3b + fcn4*yhat4a*xhat4a + fcn4*yhat4b*xhat4b;
	  solution[1] += fcn1*(1 - yhat1a)*xhat1a + fcn1*(1 - yhat1b)*xhat1b + fcn2*(1 - yhat2a)*xhat2a + fcn2*(1 - yhat2b)*xhat2b + fcn3*(1 - yhat3a)*xhat3a + fcn3*(1 - yhat3b)*xhat3b + fcn4*(1 - yhat4a)*xhat4a + fcn4*(1 - yhat4b)*xhat4b;
	  solution[2] += fcn1*yhat1a*(1 - xhat1a) + fcn1*yhat1b*(1 - xhat1b) + fcn2*yhat2a*(1 - xhat2a) + fcn2*yhat2b*(1 - xhat2b) + fcn3*yhat3a*(1 - xhat3a) + fcn3*yhat3b*(1 - xhat3b) + fcn4*yhat4a*(1 - xhat4a) + fcn4*yhat4b*(1 - xhat4b);
	  solution[3] += fcn1*(1 - yhat1a)*(1 - xhat1a) + fcn1*(1 - yhat1b)*(1 - xhat1b) + fcn2*(1 - yhat2a)*(1 - xhat2a) + fcn2*(1 - yhat2b)*(1 - xhat2b) + fcn3*(1 - yhat3a)*(1 - xhat3a) + fcn3*(1 - yhat3b)*(1 - xhat3b) + fcn4*(1 - yhat4a)*(1 - xhat4a) + fcn4*(1 - yhat4b)*(1 - xhat4b); 
	}
      }
    }
  }
}
