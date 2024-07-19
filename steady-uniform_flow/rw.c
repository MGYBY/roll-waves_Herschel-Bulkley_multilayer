/**
# Steady-uniform flow (normal flow) simulation
Based on dimensional variables.

## Problem

# Code
*/
// #include "grid/cartesian1D.h"
#include "grid/multigrid1D.h"
#include "./hydroNN.h"
#include "layered/remap.h"

// Rheological properties
#define tauY 3.4765171
#define muN 0.10
#define nHB 1.00

#define rhoFluid 1000.0
#define gravity 9.81
#define sinTheta 0.06
#define cosTheta (pow(1.0-sinTheta*sinTheta,0.50))
#define So (sinTheta/cosTheta)
#define gPrime (gravity*cosTheta)

#define froude 0.97601
#define normalDepth 0.019688057
#define normalVel (froude*pow(gPrime*normalDepth,0.50))

#define DOMAINLENGTH (20.0*normalDepth)

#define simTime 150.0

/**
 
### Non Newtonian viscosity
 The definition of viscosity as a function of shear:
 */

// double nu_eq(double shear,double pipi){
double nu_eq(double shear){
  double nu_eq;
  nu_eq = muN/rhoFluid*pow(sqrt(sq(shear) + sq(1.e-10)), (nHB-1.0))+tauY/rhoFluid/(sqrt(sq(shear) + sq(1.e-10)));
  return nu_eq;
}
/**



*/

char s[80];

int main() {
  L0 = DOMAINLENGTH;
  periodic (right);
  // G  = gravity;
  // change to g' for rotated axis frame
  G = gPrime;
  N  = 16;
  nl = 32;
  nu = 1.; // dummy

  // CFL_H = 0.40;
  CFL = 0.40;

  run();
}

/**
## Initialization  */
event init (i = 0) {
  /**
  We initialize *h*. */
  foreach(){
    // periodic BC cannot use realistic topo
    // zb[] = -x*sinTheta;
    zb[] = 0.0;
    foreach_layer() {
      u.x[] = 0.0;
      h[] = normalDepth/nl;
    }
    }
}

  /**
## Output
  We print the elevation  */
event acceleration (i++) {
  foreach(){
    foreach_layer() {
      u.x[] += dt*gravity*sinTheta;
    }
    }
    boundary ((scalar *){u});
}

/**
save the hight the flux and the yield surface as a function of time
*/ 
event output  (t = 0; t <= simTime; t+=5.0){
  sprintf (s, "slice-%g.txt", t);
  FILE * fp = fopen (s, "w"); 
  foreach (serial) {
    double H = 0.;
    foreach_layer() {
      H += h[];
      fprintf (fp, "%g %g %g \n", x, H, u.x[]);
    }
  }
  fclose(fp);
}
