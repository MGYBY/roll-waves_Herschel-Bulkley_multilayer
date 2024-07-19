/* A collection of Riemann solvers for the Saint-Venant system 
* power-law model
 */

#define epsilon 1e-30

static double mmo (double a, double b)
{
  if (a>0 && b>0) {
    if (a>b) return b;
    if (b>a) return a;
  }
  else if (a<0 && b<0) {
    if (a>b) return a;
    if (b>a) return b;
  }
  else {
    return 0.0;
  }
}

void hllc (double hm, double hp, double um, double up, double Delta,
	   double * fh, double * fq, double * dtmax)
{
  // eigenvalues based on power-law system:
  // lambda1 = u beta + sqrt(alpha h + u^2 (beta^2-beta))
  // lambda2 = u beta
  // lambda3 = u beta - sqrt(alpha h + u^2 (beta^2-beta))
  // double cm = sqrt (G*hm), cp = sqrt (G*hp);
  double cm = sqrt (alphaCoeff*hm+um*um*(betaCoeff*betaCoeff-betaCoeff));
  double cp = sqrt (alphaCoeff*hp+up*up*(betaCoeff*betaCoeff-betaCoeff));
  double ustar = betaCoeff*(um + up)/2. + cm - cp;
  double cstar = (cm + cp)/2. + betaCoeff*(um - up)/4.;
  double SL = hm == 0. ? betaCoeff*up - 2.*cp : min (betaCoeff*um - cm, ustar - cstar);
  double SR = hp == 0. ? betaCoeff*um + 2.*cm : max (betaCoeff*up + cp, ustar + cstar);

  if (0. <= SL) {
    // *fh = um*hm;
    // *fq = hm*(um*um + G*hm/2.);
    *fh = um*hm;
    *fq = hm*(betaCoeff*um*um + alphaCoeff*hm/2.);
  }
  else if (0. >= SR) {
    // *fh = up*hp;
    // *fq = hp*(up*up + G*hp/2.);
    *fh = up*hp;
    *fq = hp*(betaCoeff*up*up + alphaCoeff*hp/2.);
  }
  else {
    // double fhm = um*hm;
    // double fum = hm*(um*um + G*hm/2.);
    double fhm = um*hm;
    double fum = hm*(betaCoeff*um*um + alphaCoeff*hm/2.);
    // double fhp = up*hp;
    // double fup = hp*(up*up + G*hp/2.);
    double fhp = up*hp;
    double fup = hp*(betaCoeff*up*up + alphaCoeff*hp/2.);
    *fh = (SR*fhm - SL*fhp + SL*SR*(hp - hm))/(SR - SL);
    *fq = (SR*fum - SL*fup + SL*SR*(hp*up - hm*um))/(SR - SL);
  }

  double a = max(fabs(SL), fabs(SR));
  if (a > epsilon) {
    double dt = CFL*Delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
}

void kurganov (double hm, double hp, double um, double up, double Delta,
	       double * fh, double * fq, double * dtmax)
{
  double cp = sqrt(G*hp+sq(up)*(betaCoeff*betaCoeff-betaCoeff)), cm = sqrt(G*hm+sq(um)*(betaCoeff*betaCoeff-betaCoeff));
  double ap = max(betaCoeff*up + cp, betaCoeff*um + cm); ap = max(ap, 0.);
  double am = min(betaCoeff*up - cp, betaCoeff*um - cm); am = min(am, 0.);
  double qm = hm*um, qp = hp*up;
  double a = max(ap, -am);
  if (a > epsilon) {
    *fh = (ap*qm - am*qp + ap*am*(hp - hm))/(ap - am); // (4.5) of [1]
    *fq = (ap*(betaCoeff*qm*um + G*sq(hm)/2.) - am*(betaCoeff*qp*up + G*sq(hp)/2.) +
	    ap*am*(qp - qm))/(ap - am);
    double dt = CFL*Delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
  else
    *fh = *fq = 0.;
}

void kurganovSharp (double hm, double hp, double um, double up, double Delta,
	       double * fh, double * fq, double * dtmax)
{
  double cp = sqrt(G*hp+sq(up)*(betaCoeff*betaCoeff-betaCoeff)), cm = sqrt(G*hm+sq(um)*(betaCoeff*betaCoeff-betaCoeff));
  double ap = max(betaCoeff*up + cp, betaCoeff*um + cm); ap = max(ap, 0.);
  double am = min(betaCoeff*up - cp, betaCoeff*um - cm); am = min(am, 0.);
  double qm = hm*um, qp = hp*up;
  double a = max(ap, -am);
  if (a > epsilon) {
    double wint = (ap*hp-am*hm-(qp-qm))/(ap-am);
    double qCorr = mmo((hp-wint)/(ap-am), (wint-hm)/(ap-am));
    *fh = (ap*qm - am*qp)/(ap - am) + (ap*am*((hp - hm)/(ap - am)-qCorr)); // (4.5) of [1]
    wint = (ap*qp-am*qm-((betaCoeff*qp*up + G*pow(hp ,2.0)/2.)-(betaCoeff*qm*um + G*pow(hm ,2.0)/2.)))/(ap-am);
    qCorr = mmo((qp-wint)/(ap-am), (wint-qm)/(ap-am));
    *fq = (ap*(betaCoeff*qm*um + G*pow(hp ,2.0)/2.) - am*(betaCoeff*qp*up + G*pow(hp ,2.0)/2.))/(ap - am) + (ap*am*((qp - qm)/(ap - am)-qCorr));

    //*fh = (ap*qm - am*qp + ap*am*(hp - hm))/(ap - am); // (4.5) of [1]
    //*fq = (ap*(betaCoeff*qm*um + G*sq(hm)/2.) - am*(betaCoeff*qp*up + G*sq(hp)/2.) +
	    //ap*am*(qp - qm))/(ap - am);
    double dt = CFL*Delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
  else
    *fh = *fq = 0.;
}

void hlle (double hm, double hp, double um, double up, double Delta,
	   double * fh, double * fq, double * dtmax)
{
  // Roe average
  double uhat = (sqrt(hm)*um+sqrt(hp)*up)/(sqrt(hm)+sqrt(hp));
  double cm = sqrt (G*hm+sq(um)*(betaCoeff*betaCoeff-betaCoeff)), cp = sqrt (G*hp+sq(up)*(betaCoeff*betaCoeff-betaCoeff));
  double chat = sqrt(G*(hp+hm)/2.0+sq(uhat)*(betaCoeff*betaCoeff-betaCoeff));
  double SL = min (betaCoeff*um - cm, betaCoeff*uhat - chat); SL = min(SL, betaCoeff*up - cp);
  double SR = max (betaCoeff*up + cp, betaCoeff*uhat + chat); SR = max(SR, betaCoeff*um + cm);

  if (0. <= SL) {
    *fh = um*hm;
    *fq = hm*(betaCoeff*um*um + G*hm/2.);
  }
  else if (0. >= SR) {
    *fh = up*hp;
    *fq = hp*(betaCoeff*up*up + G*hp/2.);
  }
  else {
    double fhm = um*hm;
    double fum = hm*(betaCoeff*um*um + G*hm/2.);
    double fhp = up*hp;
    double fup = hp*(betaCoeff*up*up + G*hp/2.);
    *fh = (SR*fhm - SL*fhp + SL*SR*(hp - hm))/(SR - SL);
    *fq = (SR*fum - SL*fup + SL*SR*(hp*up - hm*um))/(SR - SL);
  }

  double a = max(fabs(SL), fabs(SR));
  if (a > epsilon) {
    double dt = CFL*Delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
}