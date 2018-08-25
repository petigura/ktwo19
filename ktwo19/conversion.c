// compute x, y, z, vx, vy, vz array from a, e, i, omega, Omega, f, mass array
double *keptostate(double a, double e, double i, double lo, double bo, double f, double m) { 

  const int sofd = SOFD;

  double mass = m; 
  if (isnan(m)) {
    printf("mass error\n");
    exit(0);
  }
  double r = a*(1.0-pow(e,2))/(1.0+e*cos(f));
  double x0 = r*cos(f);
  double y0 = r*sin(f);
  double vx0 = -sqrt(mass/(a*(1.0-pow(e,2)))) * sin(f);
  double vy0 = sqrt(mass/(a*(1.0-pow(e,2)))) * (e+cos(f));

  double *statearr = malloc(6*sofd);  
  double *state1arr, *state2arr;
  state1arr = rotatekep(x0, y0, i, lo, bo); 
  state2arr = rotatekep(vx0, vy0, i, lo, bo); 
  memcpy(statearr,state1arr,3*sofd);
  memcpy(statearr+3,state2arr,3*sofd);
  free(state1arr);free(state2arr);

  return statearr;

}



// p is a vector that has for each planet: (p, T0, eparameter1, eparameter2, inclination, bigomega, mass, radius)
// state is a matrix of the (x, y, z, vx, vy, vz) values in a stellar-centric coordinate system of each planet. 

double ***dsetup2 (double *p, const int npl){
  const int sofd = SOFD;
  const int sofds = SOFDS;

  int pperplan = PPERPLAN;
  int pstar = PSTAR;

  double epoch = EPOCH;

  double brightstar=0;
  double bsum=0;
  int i;

  double ms = p[npl*pperplan+0];
  double rstar = p[npl*pperplan+1];
  double c1 = p[npl*pperplan+2];
  double c2 = p[npl*pperplan+3];
  double dilute = p[npl*pperplan+4];

  double bigg = 1.0e0; //Newton's constant
  double ghere = G; //2.9591220363e-4; 
  double jos = 1.0/MSOMJ;  //9.545e-4; //M_jup/M_sol

  double *mp = malloc(npl*sofd);
  double *mpjup = malloc(npl*sofd);
  double *msys = malloc((npl+1)*sofd);
  msys[0] = ms;

  double *a = malloc(npl*sofd);
  double *e = malloc(npl*sofd);
  double *inc = malloc(npl*sofd);
  double *bo = malloc(npl*sofd);
  double *lo = malloc(npl*sofd);
  double *lambda = malloc(npl*sofd);
  double *f = malloc(npl*sofd);

  for (i=0;i<npl; i++) {
    if (SQRTE) {
      e[i] = pow( sqrt(pow(p[i*pperplan+2],2)+pow(p[i*pperplan+3],2)), 2);
    } else {
      e[i] = sqrt(pow(p[i*pperplan+2],2)+pow(p[i*pperplan+3],2));
    }
    inc[i] = p[i*pperplan+4]*M_PI/180;
    bo[i] = p[i*pperplan+5]*M_PI/180;
    lo[i] = atan2(p[i*pperplan+3],p[i*pperplan+2]);
    mp[i]= p[i*pperplan+6];

    mpjup[i] = mp[i]*jos;       //          ; M_Jup
    msys[i+1] = msys[i]+mpjup[i];
    a[i] = cbrt(ghere*(msys[i+1])) * pow(cbrt(p[i*pperplan+0]),2) * pow(cbrt(2*M_PI),-2);

    double pomega = bo[i]+lo[i];
    double lambda0 = getlambda( (M_PI/2-lo[i]), e[i], pomega);
    double m0 = lambda0-pomega;
    double me = m0 + 2*M_PI*(epoch - p[i*pperplan+1])/p[i*pperplan+0];
    double mepomega = me+pomega;
    double *lambdaepoint = pushpi(&mepomega,1);
    double lambdae = lambdaepoint[0];
    f[i] = getf(lambdae, e[i], pomega);

  }

  double **state = malloc(npl*sofds);
  for (i=0; i<npl;i++) {
    state[i] = keptostate(a[i],e[i],inc[i],lo[i],bo[i],f[i],(ghere*msys[i+1]));
    int j;
    for (j=0; j<6; j++) {
      state[i][j] = -state[i][j];
    }
  }


  // removed extraneous code

}




