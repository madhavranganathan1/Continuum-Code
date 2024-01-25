#include "initial.h"  /*This file contains some included header files, the running parameters, the global variables and declarations of subroutines */
int on,off,aa,bb;
#include "runge_kutta.c"   /*This is the routine that carries out the time-step integration */

#include "EValuesRead.c"   /*Eigenvalues and Eigenvectors are evaluated outside using a Mathematica program and a read using these routines */
#include "EvectorRead.c"

#include "linear_OneTCaladd.c"
void main(void){
int i,j,k;
double r,x[nx],y[ny];

double rnum[nx][ny],Rsum,Ravg,rn;

srand(time(NULL));
begin =clock();
//double test_height[nx][ny]; 
//
for(i=0;i<nx;i++){ x[i]=i*dx; kx[i]=i*((2*M_PI)/LX); }
for(i=0;i<ny;i++){ y[i]=i*dy; ky[i]=i*((2*M_PI)/LY); }
//(below is right one)
for(i=0;i<(nx/2);i++){wx1[i]=2*i*M_PI/LX;}
wx1[nx/2]=0;
for(i=(nx/2+1);i<nx;i++){wx1[i]=((i-nx)*2*M_PI)/(double)LX;}

for(i=0;i<(ny/2);i++){wy1[i]=(i*2*M_PI)/(double)LX;}
wy1[ny/2]=0;
for(i=(ny/2+1);i<ny;i++){wy1[i]=((i-ny)*2*M_PI)/(double)LY;}

for(i=0;i<(nx/2);i++){wx2[i]=2*i*M_PI/(double)LX;}
for(i=(nx/2);i<nx;i++){wx2[i]=((i-nx)*2*M_PI)/(double)LX;}

for(i=0;i<(ny/2);i++){wy2[i]=(i*2*M_PI)/LY;}
for(i=(ny/2);i<ny;i++){wy2[i]=((i-ny)*2*M_PI)/(double)LY;}

/*Fourier Transforms of the height, the substrate profile, derivatives of height, displacement fields at different orders are declared here */

dft  = fftw_plan_dft_2d(nx,ny,&in[0][0],&out[0][0],-1,FFTW_MEASURE);
idft = fftw_plan_dft_2d(nx,ny,&wh[0][0],&in2[0][0],+1,FFTW_MEASURE);

hxidft = fftw_plan_dft_2d(nx,ny,&whx[0][0],&hxin2[0][0],+1,FFTW_MEASURE);
hyidft = fftw_plan_dft_2d(nx,ny,&why[0][0],&hyin2[0][0],+1,FFTW_MEASURE);
hxxidft = fftw_plan_dft_2d(nx,ny,&whxx[0][0],&hxxin2[0][0],+1,FFTW_MEASURE);
hxyidft = fftw_plan_dft_2d(nx,ny,&whxy[0][0],&hxyin2[0][0],+1,FFTW_MEASURE);
hyyidft = fftw_plan_dft_2d(nx,ny,&whyy[0][0],&hyyin2[0][0],+1,FFTW_MEASURE);

Hxxidft = fftw_plan_dft_2d(nx,ny,&WHxx[0][0],&Hxxin2[0][0],+1,FFTW_MEASURE);
Hxyidft = fftw_plan_dft_2d(nx,ny,&WHxy[0][0],&Hxyin2[0][0],+1,FFTW_MEASURE);
Hyyidft = fftw_plan_dft_2d(nx,ny,&WHyy[0][0],&Hyyin2[0][0],+1,FFTW_MEASURE);

u1x_1idft= fftw_plan_dft_2d(nx,ny,&wu1x_1[0][0],&u1x_1in2[0][0],+1,FFTW_MEASURE);
u2x_1idft= fftw_plan_dft_2d(nx,ny,&wu2x_1[0][0],&u2x_1in2[0][0],+1,FFTW_MEASURE);
u3x_1idft= fftw_plan_dft_2d(nx,ny,&wu3x_1[0][0],&u3x_1in2[0][0],+1,FFTW_MEASURE);

u1y_1idft= fftw_plan_dft_2d(nx,ny,&wu1y_1[0][0],&u1y_1in2[0][0],+1,FFTW_MEASURE);
u2y_1idft= fftw_plan_dft_2d(nx,ny,&wu2y_1[0][0],&u2y_1in2[0][0],+1,FFTW_MEASURE);
u3y_1idft= fftw_plan_dft_2d(nx,ny,&wu3y_1[0][0],&u3y_1in2[0][0],+1,FFTW_MEASURE);

u1z_1idft= fftw_plan_dft_2d(nx,ny,&wu1z_1[0][0],&u1z_1in2[0][0],+1,FFTW_MEASURE);
u2z_1idft= fftw_plan_dft_2d(nx,ny,&wu2z_1[0][0],&u2z_1in2[0][0],+1,FFTW_MEASURE);
u3z_1idft= fftw_plan_dft_2d(nx,ny,&wu3z_1[0][0],&u3z_1in2[0][0],+1,FFTW_MEASURE);

u1_1corrdft= fftw_plan_dft_2d(nx,ny,&u1_1corr[0][0],&wu1_1corr[0][0],-1,FFTW_MEASURE);
u2_1corrdft= fftw_plan_dft_2d(nx,ny,&u2_1corr[0][0],&wu2_1corr[0][0],-1,FFTW_MEASURE);
u3_1corrdft= fftw_plan_dft_2d(nx,ny,&u3_1corr[0][0],&wu3_1corr[0][0],-1,FFTW_MEASURE);

u1z_1corrdft = fftw_plan_dft_2d(nx,ny,&u1z_1corr[0][0],&u1z_1corrout[0][0],-1,FFTW_MEASURE);
u2z_1corrdft = fftw_plan_dft_2d(nx,ny,&u2z_1corr[0][0],&u2z_1corrout[0][0],-1,FFTW_MEASURE);
u3z_1corrdft = fftw_plan_dft_2d(nx,ny,&u3z_1corr[0][0],&u3z_1corrout[0][0],-1,FFTW_MEASURE);

ss31_1corridft = fftw_plan_dft_2d(nx,ny,&wss31_1corr[0][0],&ss31_1corrin2[0][0],+1,FFTW_MEASURE);
ss32_1corridft = fftw_plan_dft_2d(nx,ny,&wss32_1corr[0][0],&ss32_1corrin2[0][0],+1,FFTW_MEASURE);
ss33_1corridft = fftw_plan_dft_2d(nx,ny,&wss33_1corr[0][0],&ss33_1corrin2[0][0],+1,FFTW_MEASURE);

u2_2idft= fftw_plan_dft_2d(nx,ny,&wu2_2[0][0],&u2_2in2[0][0],+1,FFTW_MEASURE);
u1_2idft= fftw_plan_dft_2d(nx,ny,&wu1_2[0][0],&u1_2in2[0][0],+1,FFTW_MEASURE);
u3_2idft= fftw_plan_dft_2d(nx,ny,&wu3_2[0][0],&u3_2in2[0][0],+1,FFTW_MEASURE);

ss31_2dft= fftw_plan_dft_2d(nx,ny,&ss31_2[0][0],&wss31_2[0][0],-1,FFTW_MEASURE);
ss32_2dft= fftw_plan_dft_2d(nx,ny,&ss32_2[0][0],&wss32_2[0][0],-1,FFTW_MEASURE);
ss33_2dft= fftw_plan_dft_2d(nx,ny,&ss33_2[0][0],&wss33_2[0][0],-1,FFTW_MEASURE);

els_temp_2dft= fftw_plan_dft_2d(nx,ny,&els_temp_2[0][0],&wels_temp_2[0][0],-1,FFTW_MEASURE);


ells_1idft= fftw_plan_dft_2d(nx,ny,&wels_1[0][0],&ells_1[0][0],+1,FFTW_MEASURE);
ells_1corridft= fftw_plan_dft_2d(nx,ny,&wels_1corr[0][0],&ells_1corr[0][0],+1,FFTW_MEASURE);
ells_2idft= fftw_plan_dft_2d(nx,ny,&wels_2[0][0],&ells_2[0][0],+1,FFTW_MEASURE);

ells_idft= fftw_plan_dft_2d(nx,ny,&wels[0][0],&ells[0][0],+1,FFTW_MEASURE);
gdft = fftw_plan_dft_2d(nx,ny,&gin[0][0],&gout[0][0],-1,FFTW_MEASURE);
gidft = fftw_plan_dft_2d(nx,ny,&gout[0][0],&gin2[0][0],+1,FFTW_MEASURE);

//================================================
//Rtemp6idft= fftw_plan_dft_3d(nx,ny,3,&wRtemp6[0][0][0],&Rtemp6_in2[0][0][0],+1,FFTW_MEASURE);
//Ltemp6idft= fftw_plan_dft_3d(nx,ny,3,&wLtemp6[0][0][0],&Ltemp6_in2[0][0][0],+1,FFTW_MEASURE);

Rtemp6_1idft=fftw_plan_dft_2d(nx,ny,&wRtemp6_1[0][0],&Rtemp6_1in2[0][0],+1,FFTW_MEASURE);
Rtemp6_2idft=fftw_plan_dft_2d(nx,ny,&wRtemp6_2[0][0],&Rtemp6_2in2[0][0],+1,FFTW_MEASURE);
Rtemp6_3idft=fftw_plan_dft_2d(nx,ny,&wRtemp6_3[0][0],&Rtemp6_3in2[0][0],+1,FFTW_MEASURE);

Rtemp61_1dft=fftw_plan_dft_2d(nx,ny,&Rtemp61_1[0][0],&wRtemp61_1[0][0],-1,FFTW_MEASURE);
Rtemp61_2dft=fftw_plan_dft_2d(nx,ny,&Rtemp61_2[0][0],&wRtemp61_2[0][0],-1,FFTW_MEASURE);
Rtemp61_3dft=fftw_plan_dft_2d(nx,ny,&Rtemp61_3[0][0],&wRtemp61_3[0][0],-1,FFTW_MEASURE);

Ltemp6_1idft=fftw_plan_dft_2d(nx,ny,&wLtemp6_1[0][0],&Ltemp6_1in2[0][0],+1,FFTW_MEASURE);
Ltemp6_2idft=fftw_plan_dft_2d(nx,ny,&wLtemp6_2[0][0],&Ltemp6_2in2[0][0],+1,FFTW_MEASURE);
Ltemp6_3idft=fftw_plan_dft_2d(nx,ny,&wLtemp6_3[0][0],&Ltemp6_3in2[0][0],+1,FFTW_MEASURE);

Ltemp61_1dft=fftw_plan_dft_2d(nx,ny,&Ltemp61_1[0][0],&wLtemp61_1[0][0],-1,FFTW_MEASURE);
Ltemp61_2dft=fftw_plan_dft_2d(nx,ny,&Ltemp61_2[0][0],&wLtemp61_2[0][0],-1,FFTW_MEASURE);
Ltemp61_3dft=fftw_plan_dft_2d(nx,ny,&Ltemp61_3[0][0],&wLtemp61_3[0][0],-1,FFTW_MEASURE);

//================================================

testdft = fftw_plan_dft_2d(nx,ny,&testin[0][0],&testout[0][0],-1,FFTW_MEASURE);
testidft = fftw_plan_dft_2d(nx,ny,&testout[0][0],&testin2[0][0],+1,FFTW_MEASURE);


zeta_dft =  fftw_plan_dft_2d(nx,ny,&zeta_t[0][0],&zeta_out[0][0],-1,FFTW_MEASURE);
EValuesRead();
EvectorRead();
int l,m;

//=============================================================
//Input is read from the following files below. The initial morphology is taken as pit patterned.
//

FILE *input, *input1,*in, *pit, *wpit;
in =fopen("input_height.dat","r");
input = fopen("initial_str.dat","w");
input1 = fopen("initial_str_array.dat","w");
pit= fopen("initial_pit_morpd.dat","w");
wpit=fopen("ft_initial_pit.dat","w");
double number;

for(i=0;i<nx;i++){ for(j=0;j<ny;j++){
      zeta_t[i][j][0]=  2.0*zeta_A-zeta_A*exp(-(cos(2.0*M_PI*(i)/zeta_l)+cos(2.0*M_PI*(j)/zeta_l)+2.0)/zeta_w);         zeta_t[i][j][1]=0;

  }}                       //Interface variations defined

fftw_execute(zeta_dft);

for(i=0;i<nx;i++){ for(j=0;j<ny;j++){for(k=0;k<2;k++){if(fabs(zeta_out[i][j][k])<1.00*pow(10,-10)){ zeta_out[i][j][k]=0;}}}}


for(i=0;i<nx;i++){ for(j=0;j<ny;j++){ wzeta_t[i][j][0]=zeta_out[i][j][0]; wzeta_t[i][j][1]=zeta_out[i][j][1];}}

for(i=0;i<nx;i++){ for(j=0;j<ny;j++){ wzeta[i][j]= wzeta_t[i][j][0]+I*wzeta_t[i][j][1];}}

//h0_avg= .2;

for(i=0;i<nx;i++){
  for(j=0;j<ny;j++){
     r=(rand()/(double)RAND_MAX)*.01+h0_avg+zeta_t[i][j][0];

     hight_t[i][j][0]=r;
     hight_t[i][j][1]=0.00;
     hight[i][j]=hight_t[i][j][0];
     test_height[i][j]=r;

}}
   
for(i=0;i<nx;i++){ for(j=0;j<ny;j++){fprintf(wpit, "%e\t%e\t\t%e\t%e\n",wzeta_t[i][j][0],wzeta_t[i][j][1],creal(wzeta[i][j]),cimag(wzeta[i][j]));}}

for(i=0;i<nx;i++){
  for(j=0;j<ny;j++){

     fprintf(input,"%d\t%d\t%12f\t%12f\t%12f\n",i,j,hight[i][j],hight_t[i][j][0],hight_t[i][j][1]);
     fprintf(pit,"%d\t%d\t%12f\n",i,j,zeta_t[i][j][0]);
  }
  fprintf(input,"\n");
     fprintf(input1,"\n");
      fprintf(wpit,"\n");
}
fclose(input);
fclose(input1);

linear_OneTCaladd(); /*This routine calculates some quantities which can be calculated using only the initial conditions. These are used each time step to calculate the additional terms in the height field */  

runge_kutta(); /*The main calculation of the program is the time-integration using the Runge-Kutta scheme in pseudospectral space. This program will call several other routines.  */

end = clock();
time_spent =(double)(end - begin)/CLOCKS_PER_SEC;
printf("execution_time=%lf\n",time_spent/(3600.00));

fftw_destroy_plan(dft);
fftw_destroy_plan(idft);

fftw_destroy_plan(hxdft);
fftw_destroy_plan(hxidft);

fftw_destroy_plan(hydft);
fftw_destroy_plan(hyidft);

fftw_destroy_plan(hxxdft);
fftw_destroy_plan(hxxidft);

fftw_destroy_plan(hxydft);
fftw_destroy_plan(hxyidft);

fftw_destroy_plan(hyydft);
fftw_destroy_plan(hyyidft);

////fftw_destroy_plan(Hxxdft);
//fftw_destroy_plan(Hxxidft);

////fftw_destroy_plan(Hxydft);
//fftw_destroy_plan(Hxyidft);

////fftw_destroy_plan(Hyydft);
//fftw_destroy_plan(Hyyidft);
/*
fftw_destroy_plan(fn1dft);
fftw_destroy_plan(fn2dft);
fftw_destroy_plan(fn3dft);

fftw_destroy_plan(Hfn1idft);
fftw_destroy_plan(Hfn2idft);
fftw_destroy_plan(Hfn3idft);
*/
fftw_destroy_plan(gdft);
fftw_destroy_plan(gidft);

fftw_destroy_plan(u1_1idft);
fftw_destroy_plan(u2_1idft);
fftw_destroy_plan(u3_1idft);

 
fftw_destroy_plan(u1x_1idft);
fftw_destroy_plan(u1y_1idft);
fftw_destroy_plan(u1z_1idft);

fftw_destroy_plan(u2x_1idft);
fftw_destroy_plan(u2y_1idft);
fftw_destroy_plan(u2z_1idft);

fftw_destroy_plan(u3x_1idft);
fftw_destroy_plan(u3y_1idft);
fftw_destroy_plan(u3z_1idft);

fftw_destroy_plan(u1_2idft);
fftw_destroy_plan(u2_2idft);
fftw_destroy_plan(u3_2idft);

fftw_destroy_plan(ss31_2dft);
fftw_destroy_plan(ss32_2dft);
fftw_destroy_plan(ss33_2dft);

fftw_destroy_plan(els_temp_2dft);
void fftw_cleanup(void);
}
