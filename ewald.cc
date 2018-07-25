/*  Simple implementation of Ewald Summation
    Shane Fogerty, July 2017
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "ewald.hpp"

#define PI   3.14159265358979323846264338328

using namespace std;


//TODO: load params from argument of script, or input deck?
//Define params
double alpha = 8.0;
double L = 1.0; //System length, assume cubic
double rcutoff = 0.49;//direct space cutoff less than 0.5*L
#define kmax 40

//Define params related to lattice
#define np 8 //Number of particles (cube)
double eps = 1.0; //Dielectric constant

//---------Set up arrays--------
double x[np], y[np], z[np]; //Particle coords
double q[np];               //Particle charges
static double expsum_re[2*kmax+1][2*kmax+1][2*kmax+1];
static double expsum_im[2*kmax+1][2*kmax+1][2*kmax+1];
static double G[kmax+1][kmax+1][kmax+1];

//-----------------------------


int main(void) {

   //-----Set up NaCl lattice-----
   //Should load this from input args or input deck
   int i,j,k;
   double sp = 0.5*L;//spacing
   double E = 0.0;//total energy
   double Er,Ek,Es;//energy components
   int dimx,dimy,dimz;//lattice dimensions
   dimx=dimy=dimz=2;//for a 2x2x2 cube of 8 particles

   //NaCl lattice alternating +1 and -1 charges in box
   for (i=0; i<dimx; i++) {
       for (j=0; j<dimy; j++) {
           for (k=0; k<dimz; k++) {
               x[i+j*dimy+k*dimy*dimz] = i*sp;
               y[i+j*dimy+k*dimy*dimz] = j*sp;
               z[i+j*dimy+k*dimy*dimz] = k*sp;
               q[i+j*dimy+k*dimy*dimz] = pow(-1.0,i+j*(dimy+1.0)+k*(dimy+1.0)*(dimz+1.0));
           }
       }
   }
   //-----------------------------

   Er = realComponent(q,x,y,z);   
   cout << Er;
   cout << " real component \n";
   Ek =  recipComponent(q,x,y,z);
   cout << Ek;
   cout << " k-space component \n";
   Es =  selfComponent(q);
   cout << -Es;
   cout << " self component \n";
   E += Er + Ek - Es;
   cout << E;
   cout << " total energy \n";
   cout << L*E/np;
   cout << " Madelung \n";
return 0;

}


//TODO:Use this
void Setup_Ewald(double len, int numparticles, double realcutoff, 
                 int kspacecutoff, double valuealpha)
{


}

double realComponent(double *q, double *x, double *y, double *z)
{
   double xdist,ydist,zdist,r,sum;
   int i,j;

   sum=0.0;
   //TODO: could loop through this faster I think
   for (i=0; i<np-1; i++) {
       for (j=i+1; j<np; j++) {
          xdist = x[i] - x[j];
          ydist = y[i] - y[j];
          zdist = z[i] - z[j];
          r = sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
          if (r <= rcutoff) {
             sum += q[i]*q[j]*erfc(alpha*r)/r;
          } 
       }
   }
   return 0.5*sum;
} 

double recipComponent(double *q, double *x, double *y, double *z)
{
   double inf, kr, ksum = 0.0;
   double realpart, imagpart, inside, cosine, sine;
   int i,kx,ky,kz;

   //Zero out expsum
   for (kx=-kmax; kx<=kmax; kx++){
       for (ky=-kmax; ky<=kmax; ky++){
           for (kz=-kmax; kz<=kmax; kz++){
              expsum_re[kx+kmax][ky+kmax][kz+kmax] = 0.0;
              expsum_im[kx+kmax][ky+kmax][kz+kmax] = 0.0;
           }
       }
   }
   
   influenceFunction();//compute influence funtion
   //expsum is cos(2*PI*k.r/L) + i*sin(2*PI*k.r/L)
   for (kx=-kmax; kx<=kmax; kx++) {
       for (ky=-kmax; ky<=kmax; ky++) {
           for (kz=-kmax; kz<=kmax; kz++) {
               kr = (kx*kx) + (ky*ky) + (kz*kz); 
               for (i=0; i<np; i++) {
                   if (kr < kmax*kmax) {
                      inside = 2.0 * PI * (kx*x[i] + ky*y[i] + kz*z[i]) / L;
                      expsum_re[kx+kmax][ky+kmax][kz+kmax] += q[i]*cos(inside);
                      expsum_im[kx+kmax][ky+kmax][kz+kmax] -= q[i]*sin(inside);
                   }
               }
                                
               for (i=0; i<np; i++) {
                   if (kr < kmax*kmax) {
                      inf = G[abs(kx)][abs(ky)][abs(kz)];//influence function
                      inside = 2.0*PI*(kx*x[i]+ky*y[i]+kz*z[i])/L;
                      cosine = cos(inside);
                      sine = sin(inside);
                      realpart = expsum_re[kx+kmax][ky+kmax][kz+kmax] * cosine;
                      imagpart = expsum_im[kx+kmax][ky+kmax][kz+kmax] * sine;
                      ksum += q[i]*inf*(realpart-imagpart);
                   }
               }

           }
       }
   }
   
   return (ksum*L/(4.0*PI));
}


double selfComponent(double *q)
{
   int i;
   double sum = 0.0;
   for (i=0; i<np; i++) {
       sum += q[i]*q[i];
   }
   return alpha*sum/sqrt(PI);

}

void influenceFunction(void) {
   int kx,ky,kz;
   double kr, f1, f2;
   f1 = 2.0/(L*L);
   f2 = PI/(alpha*L);
   for (kx=0; kx<=kmax; kx++){
       for (ky=0; ky<=kmax; ky++){
           for (kz=0; kz<=kmax; kz++){
               kr = (kx*kx) + (ky*ky) + (kz*kz);
               G[kx][ky][kz] = f1*(exp(-(f2*f2)*kr))/kr;
           }
       }
   }
   G[0][0][0] = 0.0;
}
