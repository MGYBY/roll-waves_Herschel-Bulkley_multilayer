/**
# Steady-uniform flow (normal flow) simulation
Based on dimensional variables.

## Problem

# Code
*/
// #include "grid/cartesian1D.h"
// #include "grid/bitree.h"
#include "grid/quadtree.h"
#include "./saint-venantNN.h"
// #include "./output_vts.h"
#include "vtk.h"

// Rheological properties
#define bParam 0.30
#define nHB 0.33333
#define froude 0.50
#define grav (1.0/pow(froude,2.0))
#define rParam (pow(((nHB/(nHB+1.0))*pow(grav,(1.0/nHB))*pow((1.0-bParam),((nHB+1.0)/nHB))*(1.0-(nHB/(2.0*nHB+1.0))*(1.0-bParam))), nHB))

#define normalDepth 1.00
// #define normalDepth 0.095263786
#define normalVel 1.00
// #define alphaCoeff (tauY/(rhoFluid*gravity*normalDepth*sinTheta))
#define normalMaxVel (normalVel/(1.0-nHB/(2.0*nHB+1.0)*(1.0-bParam)))
#define colorbarMax (normalMaxVel*3.50)

#define piVal 3.14159265
#define disAmp 0.20
// #define hBC(t) (t<=distPeriod/2.0 ? (normalDepth*(1.0+disAmp*sin(2. * pi * t / distPeriod))) : (normalDepth))

#define DOMAINLENGTH (200.0)
#define PLOTRANGEMAX (5.00*normalDepth)

#define distWL 2.0
#define refLen (distWL*1.0)
#define Ly (distWL*16.0)

#define MAXLEVEL 14
#define MINLEVEL 2
#define INITLEVEL 8

#define simTime 250.0
#define outputInterval 0.25

#define visRegThre 1.20e-18
#define yieldSurfThre 2.e-11

#define hError (normalDepth/300.0)
#define hmaxError (normalDepth/500.0)
#define uError (normalVel/300.0)
#define hGradError (0.01)
#define uxGradError (0.01)

/**
 
### Non Newtonian viscosity
 The definition of viscosity as a function of shear:
 */

// double nu_eq(double shear,double pipi){
double nu_eq(double shear){
  double nu_eq;
  nu_eq = (rParam*pow(sqrt(sq(shear)), (nHB)) + (bParam/pow(froude,2.0)))/(sqrt(sq(shear)+visRegThre));
//   nu_eq = rParam*pow(sqrt(sq(shear) + sq(visRegThre)), (nHB-1.0)) + (bParam/pow(froude,2.0))*1.0/(sqrt(sq(shear) + sq(visRegThre)));
  // nu_eq = muN/rhoFluid*pow(sqrt(sq(shear) + sq(1.e-10)), (nHB-1.0))+tauY/rhoFluid/(sqrt(sq(shear) + sq(1.e-10)))*(1.-exp((-1.)*papaCoeff*sqrt(sq(shear) + sq(1.e-10))));
  return nu_eq;
}

double hbProfile (double z, double ht, double umaxSurf)
{
  // general normal-flow profile of Herschel-Bulkley fluids
  double zo;
//   zo = ht*(1.-bParam);
  zo = normalDepth*(1.-bParam);
  // zo = normalDepth*alphaCoeff;
  // if (z<=(normalDepth*alphaCoeff)) {
  if (z<=zo) {
    // shear zone
    return ((1.0-pow((1.0-z/zo),((nHB+1.0)/nHB)))*umaxSurf);
  }
  else {
    // plug zone
    return (umaxSurf);
  }
}

char s[40], vtkName[40], gerrisName[40];
scalar depthGrad[], uAve[], yieldSurf[], hmax[];

int main() {
  L0 = DOMAINLENGTH;
  // G  = gravity;
  // change to g' for rotated axis frame
  G = grav;
  init_grid(1 << (INITLEVEL));
  nl = 64;
  nu = 1.; // dummy

  CFL = 0.399;

  // slightly higher resolution
  theta = 1.50;

  periodic (right);

  run();
}

/**
## Initialization  */
event init (i = 0) {
  refine(y>=(Ly-1.1*DOMAINLENGTH/pow(2, MAXLEVEL)) && y<=(Ly+1.1*DOMAINLENGTH/pow(2, MAXLEVEL)) && level<MAXLEVEL);
  refine(fabs(y)<=(distWL*2.0) && x<=(distWL*2.05) && level<MAXLEVEL);

  // for 2-D case 40x3m
  mask(y > Ly ? top : none);

  /**
  We initialize *h*. */
  foreach(){
//     double totalDepth = pow((x-(distWL/2.0+distWL/4.0)), 2.0)+pow(y, 2.0)<=pow((distWL/4.0), 2.0) ? normalDepth*(1.0+disAmp*pow((pow((distWL/4.0), 2.0)-(pow((x-(distWL/2.0+distWL/4.0)), 2.0)+pow(y, 2.0))), 0.50)/(distWL/4.0)) : normalDepth;
    double totalDepth = pow((x-(distWL/2.0+distWL/4.0)), 2.0)+pow((y), 2.0)<=pow((distWL/4.0), 2.0) ? normalDepth*(1.0+disAmp*pow((pow((distWL/4.0), 2.0)-(pow((x-(distWL/2.0+distWL/4.0)), 2.0)+pow((y), 2.0))), 0.50)/(distWL/4.0)) : normalDepth;
    double maxTotalVel = pow((x-(distWL/2.0+distWL/4.0)), 2.0)+pow((y), 2.0)<=pow((distWL/4.0), 2.0) ? pow(totalDepth, 0.50)*(normalMaxVel/normalVel) : pow(normalDepth, 0.50)*(normalMaxVel/normalVel);
    double z = zb[];
    // periodic BC cannot use realistic topo
    // zb[] = -x*sinTheta;
    zb[] = 0.0;
    h[] = totalDepth;

    hGradMag[] = 0.0;
    uGradMag[] = 0.0;

    hmax[] = h[];

    uAve[] = 0.0;

//     depthGrad[] = fabs(h[-1]-h[1])/(2.*Delta);

      for (int l = 0; l < nl; l++) {
        z += layer[l]*h[]/2.;
        u = ul[l];
        u.x[] = hbProfile(z, totalDepth, maxTotalVel);
        // u.x[] = normalVel;
        u.y[] = 0.00;
        // TODO: close follow Gaussian test case
//         ul.n[left] = t<=distPeriod/2.0 ? dirichlet(hbProfileBoundary(point, _s, z, distDepth)) : dirichlet(hbProfileBoundary(point, _s, z, normalDepth));
//         u.n[right] = neumann(0);
        z += layer[l]*h[]/2.;

        uAve[] += pow((u.x[]*u.x[]+u.y[]*u.y[]),0.50)*layer[l];
      }
    }
}

  /**
## Output
  We print the elevation  */
event acceleration (i++) {
  foreach(){
//     depthGrad[] = fabs(h[-1]-h[1])/(2.*Delta);
    uAve[] = 0.0;
    double zCoord = zb[];
    double uVertGrad = 0.0;
    yieldSurf[] = 0.0;

      for (int l = 0; l < nl; l++) {
        zCoord += layer[l]*h[]*0.50;
        u = ul[l];
        u.x[] += dt*grav;

        uAve[] += pow((u.x[]*u.x[]+u.y[]*u.y[]),0.50)*layer[l];

//         if (l>0 && l<(nl-1)) {
//         // middle layers
//         vector um = ul[l-1] ;
//         vector up = ul[l+1] ;
//         uVertGrad = (up.x[]-um.x[])/(layer[l-1]*h[]*0.50+layer[l+1]*h[]*0.50+layer[l]*h[]);
//       }
//       else if (l==(nl-1))
//       {
//         // top layer
//         vector um = ul[l-1] ;
//         vector up = ul[l] ;
//         uVertGrad = (up.x[]-um.x[])/(layer[l-1]*h[]*0.50+layer[l]*h[]*1.50);
//       }
//       else
//       {
//         // bottom layer
//         vector um = ul[l] ;
//         vector up = ul[l+1] ;
//         uVertGrad = (up.x[]-um.x[])/(layer[l+1]*h[]*0.50+layer[l]*h[]*1.50);
//       }
//
// //       if (uVertGrad<yieldSurfThre && yieldSurf[]<=0 && bParam>1.e-5)
//       if (uVertGrad<yieldSurfThre && yieldSurf[]<=0 )
//       {
//         yieldSurf[] = zCoord;
//       }

      zCoord += layer[l]*h[]*0.50;
      }

      if (h[]>hmax[]) {
        hmax[]=h[];
      }

  }
    boundary ((scalar *){u});
}

/*
 * AMR here
 *
 */
event adapt1 (i++) {
  // adapt_wavelet({h, depthGrad, uAve}, (double[]){normalDepth/200.0, 0.01, normalVel/200.0}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // a more reasonable AMR criteria
//   adapt_wavelet({h, depthGrad, uAve, yieldSurf}, (double[]){normalDepth/300.0, 0.007, normalVel/300.0, normalDepth/150.0}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);

  adapt_wavelet({h, hmax, uAve, hGradMag, uGradMag}, (double[]){hError, hmaxError, uError, hGradError, uxGradError}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);

//      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);

  refine(x>=(distWL/2.0/2.50) && x<=distWL+refLen && fabs(y)<=(refLen*0.50) && t<=(0.50*(distWL/normalVel)) && level<MAXLEVEL);
//   refine(x<=(9.9*DOMAINLENGTH/pow(2,MAXLEVEL)) && level<MAXLEVEL);
//   refine(x>=(DOMAINLENGTH-9.9*DOMAINLENGTH/pow(2,MAXLEVEL)) && level<MAXLEVEL);
}

event hmax(i+=20)
{
  //    double maxDepth = normalDepth;
     double maxDepthLocX = 0.0;
  //    double maxDepthVel = 0.0;
  //    double minDepth = normalDepth;
  //    double minDepthLocX = 0.0;
	 // double max2LocX = 0.0;

     stats s1 = statsf (h);
     stats s2 = statsf (uAve);

     FILE *fp2 = fopen("maxMinDepth", "a+");

     // FR: maxDepth, minDepth, H, waveLength
     // foreach ()
     // {
     //           if (y<=2.01*(DOMAINLENGTH/pow(2, MAXLEVEL)) && h[]>h[1] && h[]>h[-1] && x>maxDepthLocX && h[]>(normalDepth*(1.0+0.250*disAmp))){
     //                maxDepth = h[];
     //                maxDepthLocX = x;
     //                maxDepthVel = uAve[];
     //           }
     // }

     foreach ()
     {
       if (h[]==s1.max)
       {
         maxDepthLocX = x;
      }
    }

// 	 foreach ()
// 	 {
// // 		if (h[]<h[1] && h[]<h[-1] && h[-1]<h[-2] && h[1]<h[2] && x>minDepthLocX && x<maxDepthLocX ) {
//        if (h[]<h[1] && h[]<h[-1] && h[-1]<h[-2] && h[1]<h[2] && x>minDepthLocX && h[]<(normalDepth*(1.0+disAmp/10.0)) ) {
//                     minDepthLocX = x;
//                     minDepth = h[];
//        }
//
// // 		else if (h[]>h[1] && h[]>h[-1] && h[1]>h[2] && h[-1]>h[-2] && h[]>1.1 && x<maxDepthLocX && x>max2LocX) {
// //                     max2LocX = x;
// //                }
// 	 }
//      fprintf(fp2, "%g %13.9f %13.9f %13.9f %13.9f %13.9f \n", t, maxDepthLocX, maxDepth, minDepth, (maxDepth-minDepth), (maxDepthLocX-max2LocX));
//     fprintf(fp2, "%g %13.9f %13.9f %13.9f %13.9f %13.9f \n", t, maxDepthLocX, maxDepth, maxDepthVel, minDepth, (maxDepth-minDepth));
    fprintf(fp2, "%g %13.9f %13.9f %13.9f \n", t, maxDepthLocX, s1.max, s2.max );
    fclose(fp2);
}

/**
save the hight the flux and the yield surface as a function of time
*/ 
event output  (t = 0; t <= simTime; t+=outputInterval){
  // sprintf (s, "slice-%g.txt", t);
  // sprintf(vtkName, "outVTK-%g.vtk", t);
  sprintf(gerrisName, "snapshot-%g.gfs", t);
  // FILE * fp1 = fopen (s, "w");
  // FILE * fp3 = fopen(vtkName, "w");
  FILE * fp4 = fopen(gerrisName, "w");

//   foreach(serial){
//     double zCoord = zb[];
//     double uVertGrad = 0.0;
//     for (int l = 0; l < nl; l++) {
//       zCoord += layer[l]*h[]*0.50;
//       u = ul[l];
//
// //       if (l>0 && l<(nl-1)) {
// //         // middle layers
// //         vector um = ul[l-1] ;
// //         vector up = ul[l+1] ;
// //         uVertGrad = (up.x[]-um.x[])/(layer[l-1]*h[]*0.50+layer[l+1]*h[]*0.50+layer[l]*h[]);
// //       }
// //       else if (l==(nl-1))
// //       {
// //         // top layer
// //         vector um = ul[l-1] ;
// //         vector up = ul[l] ;
// //         uVertGrad = (up.x[]-um.x[])/(layer[l-1]*h[]*0.50+layer[l]*h[]*1.50);
// //       }
// //       else
// //       {
// //         // bottom layer
// //         vector um = ul[l] ;
// //         vector up = ul[l+1] ;
// //         uVertGrad = (up.x[]-um.x[])/(layer[l+1]*h[]*0.50+layer[l]*h[]*1.50);
// //       }
//
//       // output yield surface only for HB or BP
//       // if(bParam>1.e-5)
// //         fprintf (fp1, "%g %g %g %g \n", x, zCoord, u.x[], uVertGrad);
//       fprintf (fp1, "%g %g %g %g %g \n", x, y, zCoord, u.x[], u.y[]);
//
//       zCoord += layer[l]*h[]*0.50;
//     }
//   }
//   fclose(fp1);

  sprintf (s, "depth-%g.txt", t);
  FILE * fp2 = fopen (s, "w");
  foreach (serial) {
//     fprintf (fp2, "%g %g %g %g \n", x, h[], uAve[], yieldSurf[]);
    fprintf (fp2, "%g %g %g %g \n", x, y, h[], uAve[]);
  }
  fclose(fp2);

  // output_vtk((scalar *) {h, uAve}, 1 << MAXLEVEL, (FILE *) fp3, false);
  // fclose(fp3);
  output_gfs(fp4, list = {h, uAve}, translate = true);
  fclose(fp4);
}

// event moviemaker (t += 12) {
//   // output_ppm (h, map = jet, linear = true, n = 400, file = "depth.mp4");
//   // map = radial, map = cool_warm
//   // scalar l[];
//   // foreach()
//   //   l[] = level;
//   // output_ppm (l, map = cool_warm, min = 4, max = LEVEL, n = 400,
// 	 //      file = "level.mp4");
//
//   view (map = jet, fov = 45, width = 600, height = 600, samples = 1);
//   clear();
//   squares ("h", linear = true);
//   save ("h_animation.mp4");
// }
