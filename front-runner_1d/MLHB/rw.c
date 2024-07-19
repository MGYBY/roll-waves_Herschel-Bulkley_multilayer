/**
# Steady-uniform flow (normal flow) simulation
Based on dimensional variables.

## Problem

# Code
*/
// #include "grid/cartesian1D.h"
#include "grid/bitree.h"
#include "./saint-venantNN.h"

// Rheological properties
#define bParam 0.40
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
#define refLen (distWL*4.0)

#define MAXLEVEL 14
#define MINLEVEL 3
#define INITLEVEL 14

#define simTime 80.0
#define outputInterval 0.25

#define visRegThre 1.0e-12
#define yieldSurfThre 2.e-11

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
  zo = ht*(1.-bParam);
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

char s[80];
scalar depthGrad[], uAve[], yieldSurf[];

int main() {
  L0 = DOMAINLENGTH;
  // G  = gravity;
  // change to g' for rotated axis frame
  G = grav;
  init_grid(1 << (INITLEVEL));
  nl = 64;
  nu = 1.; // dummy

  CFL = 0.399;

  periodic (right);

  run();
}

// double uleft (Point point, scalar s, double ht)
// {
//   double zo;
//   double zc = zb[];
//   zo = ht*(1.-bParam);
//   // zo = normalDepth*alphaCoeff;
//   // if (z<=(normalDepth*alphaCoeff)) {
//   for (int l = 0; l < s.l; l++)
//     zc += layer[l]*h[];
//   zc += layer[s.l]*h[]/2.;
//   if (zc<=zo) {
//     // shear zone
//     return ((1.0-pow((1.0-z/zo),((nHB+1.0)/nHB)))*normalMaxVel);
//   }
//   else {
//     // plug zone
//     return (normalMaxVel);
//   }
// }

/**
## Initialization  */
event init (i = 0) {
  /**
  We initialize *h*. */
  foreach(){
    double totalDepth = (x<=distWL && x>=(distWL/2.0)) ? normalDepth*(1.0+disAmp*sin(2. * pi * (x-(distWL/2.0)) / distWL)) : normalDepth;
    double maxTotalVel = (x<=distWL && x>=(distWL/2.0)) ? pow(normalDepth*(1.0+disAmp*sin(2. * pi * (x-(distWL/2.0)) / distWL)), 0.50)*(normalMaxVel/normalVel) : pow(normalDepth, 0.50)*(normalMaxVel/normalVel);
    double z = zb[];
    // periodic BC cannot use realistic topo
    // zb[] = -x*sinTheta;
    zb[] = 0.0;
    h[] = totalDepth;

    uAve[] = 0.0;

    depthGrad[] = fabs(h[-1]-h[1])/(2.*Delta);

      for (int l = 0; l < nl; l++) {
        z += layer[l]*h[]/2.;
        u = ul[l];
        u.x[] = hbProfile(z, totalDepth, maxTotalVel);
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
    depthGrad[] = fabs(h[-1]-h[1])/(2.*Delta);
    uAve[] = 0.0;
    double zCoord = zb[];
    double uVertGrad = 0.0;
    yieldSurf[] = 0.0;

      for (int l = 0; l < nl; l++) {
        zCoord += layer[l]*h[]*0.50;
        u = ul[l];
        u.x[] += dt*grav;

        uAve[] += pow((u.x[]*u.x[]+u.y[]*u.y[]),0.50)*layer[l];

        if (l>0 && l<(nl-1)) {
        // middle layers
        vector um = ul[l-1] ;
        vector up = ul[l+1] ;
        uVertGrad = (up.x[]-um.x[])/(layer[l-1]*h[]*0.50+layer[l+1]*h[]*0.50+layer[l]*h[]);
      }
      else if (l==(nl-1))
      {
        // top layer
        vector um = ul[l-1] ;
        vector up = ul[l] ;
        uVertGrad = (up.x[]-um.x[])/(layer[l-1]*h[]*0.50+layer[l]*h[]*1.50);
      }
      else
      {
        // bottom layer
        vector um = ul[l] ;
        vector up = ul[l+1] ;
        uVertGrad = (up.x[]-um.x[])/(layer[l+1]*h[]*0.50+layer[l]*h[]*1.50);
      }

//       if (uVertGrad<yieldSurfThre && yieldSurf[]<=0 && bParam>1.e-5)
      if (uVertGrad<yieldSurfThre && yieldSurf[]<=0 )
      {
        yieldSurf[] = zCoord;
      }

      zCoord += layer[l]*h[]*0.50;
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
  adapt_wavelet({h, depthGrad, uAve, yieldSurf}, (double[]){normalDepth/300.0, 0.007, normalVel/300.0, normalDepth/150.0}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);

  refine(x<=refLen && t<=(4.0*(distWL/normalVel)) && level<MAXLEVEL);
  refine(x<=(9.9*DOMAINLENGTH/pow(2,MAXLEVEL)) && level<MAXLEVEL);
  refine(x>=(DOMAINLENGTH-9.9*DOMAINLENGTH/pow(2,MAXLEVEL)) && level<MAXLEVEL);
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
  sprintf (s, "slice-%g.txt", t);
  FILE * fp1 = fopen (s, "w");
  foreach(serial){
    double zCoord = zb[];
    double uVertGrad = 0.0;
    for (int l = 0; l < nl; l++) {
      zCoord += layer[l]*h[]*0.50;
      u = ul[l];

      if (l>0 && l<(nl-1)) {
        // middle layers
        vector um = ul[l-1] ;
        vector up = ul[l+1] ;
        uVertGrad = (up.x[]-um.x[])/(layer[l-1]*h[]*0.50+layer[l+1]*h[]*0.50+layer[l]*h[]);
      }
      else if (l==(nl-1))
      {
        // top layer
        vector um = ul[l-1] ;
        vector up = ul[l] ;
        uVertGrad = (up.x[]-um.x[])/(layer[l-1]*h[]*0.50+layer[l]*h[]*1.50);
      }
      else
      {
        // bottom layer
        vector um = ul[l] ;
        vector up = ul[l+1] ;
        uVertGrad = (up.x[]-um.x[])/(layer[l+1]*h[]*0.50+layer[l]*h[]*1.50);
      }

      // output yield surface only for HB or BP
      // if(bParam>1.e-5)
        fprintf (fp1, "%g %g %g %g \n", x, zCoord, u.x[], uVertGrad);

      zCoord += layer[l]*h[]*0.50;
    }
  }
  fclose(fp1);

  sprintf (s, "depth-%g.txt", t);
  FILE * fp2 = fopen (s, "w");
  foreach (serial) {
    fprintf (fp2, "%g %g %g %g \n", x, h[], uAve[], yieldSurf[]);
  }
  fclose(fp2);
}

/**
Post-processing visualization module
*/ 
// void setup (FILE * fp)
// {
//   // FIXME: other customized color schemes
//   fprintf (fp,
// 	   "set pm3d map interpolate 2,2\n"
// // 	   "# jet colormap\n"
// // 	   "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
// // 	   " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
// // 	   " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
// //        "# rainbow colormap\n"
//        "set palette rgb 33,13,10\n"
//       //  "load 'turbo.pal'\n"
// //        "# green-red-violet colormap\n"
// //        "set palette rgb 3,11,6\n"
// 	   "unset key\n"
// 	   "set cbrange [0:%g]\n"
// 	   "set xlabel 'x'\n"
// 	   "set ylabel 'height'\n"
// 	   "set xrange [%g:%g]\n"
// 	   "set yrange [0.0:%g]\n"
//      "set lmargin at screen 0.05\n"
//      "set rmargin at screen 0.95\n", colorbarMax, 0.0, DOMAINLENGTH, PLOTRANGEMAX
// 	   );
// }

void setup (FILE * fp)
{
  // FIXME: other customized color schemes
  fprintf (fp,
	   "set pm3d map interpolate 2,2\n"
// 	   "# jet colormap\n"
// 	   "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
// 	   " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
// 	   " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
//        "# rainbow colormap\n"
       "set palette rgb 33,13,10\n"
      //  "load 'turbo.pal'\n"
//        "# green-red-violet colormap\n"
//        "set palette rgb 3,11,6\n"
	   "unset key\n"
	   "set cbrange [0:]\n"
	   "set xlabel 'x'\n"
	   "set ylabel 'height'\n"
	   "set xrange [%g:%g]\n"
	   "set yrange [0.0:]\n"
     "set lmargin at screen 0.05\n"
     "set rmargin at screen 0.95\n" , 0.0,DOMAINLENGTH
	   );
}

void plot (FILE * fp)
{
  fprintf (fp,
	   "set title 't = %g'\n"
	   "sp '-' u 1:2:3\n", t);
  foreach (serial) {
    double z = zb[];
    // fprintf (fp, "%g %g %g\n", x, max(0.5, z), u.x[]);
    for (int l = 0; l < nl; l++) {
      z += layer[l]*h[]*0.50;
      u = ul[l];
      // z += h[]/2.0;
      fprintf (fp, "%g %g %g \n", x, z, u.x[]);
      // z += h[]/2.0;
      z += layer[l]*h[]*0.50;
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  //  fprintf (fp, "pause 1\n");
  fflush (fp);  
}

event gnuplot (t += outputInterval; t <= simTime)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  fprintf (fp,
	   "set term pngcairo font \",10\" size 1000,313\n"
	   "set output 'plot-%g.png'\n", t);
  if (i == 0)
    setup (fp);
  plot (fp);
}

// event moviemaker (t = end) {
//     system ("for f in plot-*.png; do convert $f ppm:- && rm -f $f; done | "
// 	    "ppm2mp4 movie.mp4");
// }
