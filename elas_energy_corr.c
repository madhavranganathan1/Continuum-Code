#include "evolution_eqn_elas.c"
void elas_energy_corr(void){
int i,j,k,l,m,c;
double ui[nx][ny];
for(i=0;i<nx;i++){for(j=0;j<ny;j++){
    wss31_1[i][j][0]=ss_0*whx[i][j][0];               wss31_1[i][j][1]=ss_0*whx[i][j][1];   // First order stress 
    wss32_1[i][j][0]=ss_0*why[i][j][0];               wss32_1[i][j][1]=ss_0*why[i][j][1];
    wss33_1[i][j][0]=0;                              wss33_1[i][j][1]=0;                       }}

for(i=0;i<nx;i++){for(j=0;j<ny;j++){
/*
wss31_1corr[i][j][0]=(lkfg[i][j][0][0][0]*wss31_1[i][j][0]-lkfg[i][j][0][0][1]*wss31_1[i][j][1])+(lkfg[i][j][0][1][0]*wss32_1[i][j][0]-lkfg[i][j][0][1][1]*wss32_1[i][j][1])+(lkfg[i][j][0][2][0]*wss33_1[i][j][0]-lkfg[i][j][0][2][1]*wss33_1[i][j][1]);

wss31_1corr[i][j][1]=(lkfg[i][j][0][0][0]*wss31_1[i][j][1]+lkfg[i][j][0][0][1]*wss31_1[i][j][0])+(lkfg[i][j][0][1][0]*wss32_1[i][j][1]+lkfg[i][j][0][1][1]*wss32_1[i][j][0])+(lkfg[i][j][0][2][0]*wss33_1[i][j][1]+lkfg[i][j][0][2][1]*wss33_1[i][j][0]);

wss32_1corr[i][j][0]=(lkfg[i][j][1][0][0]*wss31_1[i][j][0]-lkfg[i][j][1][0][1]*wss31_1[i][j][1])+(lkfg[i][j][1][1][0]*wss32_1[i][j][0]-lkfg[i][j][1][1][1]*wss32_1[i][j][1])+(lkfg[i][j][1][2][0]*wss33_1[i][j][0]-lkfg[i][j][1][2][1]*wss33_1[i][j][1]);

wss32_1corr[i][j][1]=(lkfg[i][j][1][0][0]*wss31_1[i][j][1]+lkfg[i][j][1][0][1]*wss31_1[i][j][0])+(lkfg[i][j][1][1][0]*wss32_1[i][j][1]+lkfg[i][j][1][1][1]*wss32_1[i][j][0])+(lkfg[i][j][1][2][0]*wss33_1[i][j][1]+lkfg[i][j][1][2][1]*wss33_1[i][j][0]);

wss33_1corr[i][j][0]=(lkfg[i][j][2][0][0]*wss31_1[i][j][0]-lkfg[i][j][2][0][1]*wss31_1[i][j][1])+(lkfg[i][j][2][1][0]*wss32_1[i][j][0]-lkfg[i][j][2][1][1]*wss32_1[i][j][1])+(lkfg[i][j][2][2][0]*wss33_1[i][j][0]-lkfg[i][j][2][2][1]*wss33_1[i][j][1]);

wss33_1corr[i][j][1]=(lkfg[i][j][2][0][0]*wss31_1[i][j][1]+lkfg[i][j][2][0][1]*wss31_1[i][j][0])+(lkfg[i][j][2][1][0]*wss32_1[i][j][1]+lkfg[i][j][2][1][1]*wss32_1[i][j][0])+(lkfg[i][j][2][2][0]*wss33_1[i][j][1]+lkfg[i][j][2][2][1]*wss33_1[i][j][0]);

*/
wss31_1corr[i][j][0]=(lkfg[i][j][0][0][0]*wss31_1[i][j][0]-lkfg[i][j][0][0][1]*wss31_1[i][j][1])+(lkfg[i][j][0][1][0]*wss32_1[i][j][0]-lkfg[i][j][0][1][1]*wss32_1[i][j][1])+(lkfg[i][j][0][2][0]*wss33_1[i][j][0]-lkfg[i][j][0][2][1]*wss33_1[i][j][1])+wss1z_corr[i][j][0];

wss31_1corr[i][j][1]=(lkfg[i][j][0][0][0]*wss31_1[i][j][1]+lkfg[i][j][0][0][1]*wss31_1[i][j][0])+(lkfg[i][j][0][1][0]*wss32_1[i][j][1]+lkfg[i][j][0][1][1]*wss32_1[i][j][0])+(lkfg[i][j][0][2][0]*wss33_1[i][j][1]+lkfg[i][j][0][2][1]*wss33_1[i][j][0])+wss1z_corr[i][j][1];

wss32_1corr[i][j][0]=(lkfg[i][j][1][0][0]*wss31_1[i][j][0]-lkfg[i][j][1][0][1]*wss31_1[i][j][1])+(lkfg[i][j][1][1][0]*wss32_1[i][j][0]-lkfg[i][j][1][1][1]*wss32_1[i][j][1])+(lkfg[i][j][1][2][0]*wss33_1[i][j][0]-lkfg[i][j][1][2][1]*wss33_1[i][j][1])+wss2z_corr[i][j][0];

wss32_1corr[i][j][1]=(lkfg[i][j][1][0][0]*wss31_1[i][j][1]+lkfg[i][j][1][0][1]*wss31_1[i][j][0])+(lkfg[i][j][1][1][0]*wss32_1[i][j][1]+lkfg[i][j][1][1][1]*wss32_1[i][j][0])+(lkfg[i][j][1][2][0]*wss33_1[i][j][1]+lkfg[i][j][1][2][1]*wss33_1[i][j][0])+wss2z_corr[i][j][1];

wss33_1corr[i][j][0]=(lkfg[i][j][2][0][0]*wss31_1[i][j][0]-lkfg[i][j][2][0][1]*wss31_1[i][j][1])+(lkfg[i][j][2][1][0]*wss32_1[i][j][0]-lkfg[i][j][2][1][1]*wss32_1[i][j][1])+(lkfg[i][j][2][2][0]*wss33_1[i][j][0]-lkfg[i][j][2][2][1]*wss33_1[i][j][1])+wss3z_corr[i][j][0];

wss33_1corr[i][j][1]=(lkfg[i][j][2][0][0]*wss31_1[i][j][1]+lkfg[i][j][2][0][1]*wss31_1[i][j][0])+(lkfg[i][j][2][1][0]*wss32_1[i][j][1]+lkfg[i][j][2][1][1]*wss32_1[i][j][0])+(lkfg[i][j][2][2][0]*wss33_1[i][j][1]+lkfg[i][j][2][2][1]*wss33_1[i][j][0])+wss3z_corr[i][j][1];

}}

fftw_execute(ss31_1corridft);
fftw_execute(ss32_1corridft);
fftw_execute(ss33_1corridft);


for(i=0;i<nx;i++){for(j=0;j<ny;j++){
ss31_1corr[i][j][0] = (ss31_1corrin2[i][j][0]/(nx*ny))*(hight[i][j]-hbar);             ss31_1corr[i][j][1] = 0;
ss32_1corr[i][j][0] = (ss32_1corrin2[i][j][0]/(nx*ny))*(hight[i][j]-hbar);             ss32_1corr[i][j][1] = 0;
ss33_1corr[i][j][0] = (ss33_1corrin2[i][j][0]/(nx*ny))*(hight[i][j]-hbar);             ss33_1corr[i][j][1] = 0;
}}

for(i=0;i<nx;i++){for(j=0;j<ny;j++){
/*
wu1_1[i][j][0]=(fg[i][j][0][0][0]*wss31_1[i][j][0]-fg[i][j][0][0][1]*wss31_1[i][j][1])+(fg[i][j][0][1][0]*wss32_1[i][j][0]-fg[i][j][0][1][1]*wss32_1[i][j][1])+(fg[i][j][0][2][0]*wss33_1[i][j][0]-fg[i][j][0][2][1]*wss33_1[i][j][1]);

wu1_1[i][j][1]=(fg[i][j][0][0][0]*wss31_1[i][j][1]+fg[i][j][0][0][1]*wss31_1[i][j][0])+(fg[i][j][0][1][0]*wss32_1[i][j][1]+fg[i][j][0][1][1]*wss32_1[i][j][0])+(fg[i][j][0][2][0]*wss33_1[i][j][1]+fg[i][j][0][2][1]*wss33_1[i][j][0]);

wu2_1[i][j][0]=(fg[i][j][1][0][0]*wss31_1[i][j][0]-fg[i][j][1][0][1]*wss31_1[i][j][1])+(fg[i][j][1][1][0]*wss32_1[i][j][0]-fg[i][j][1][1][1]*wss32_1[i][j][1])+(fg[i][j][1][2][0]*wss33_1[i][j][0]-fg[i][j][1][2][1]*wss33_1[i][j][1]);

wu2_1[i][j][1]=(fg[i][j][1][0][0]*wss31_1[i][j][1]+fg[i][j][1][0][1]*wss31_1[i][j][0])+(fg[i][j][1][1][0]*wss32_1[i][j][1]+fg[i][j][1][1][1]*wss32_1[i][j][0])+(fg[i][j][1][2][0]*wss33_1[i][j][1]+fg[i][j][1][2][1]*wss33_1[i][j][0]);

wu3_1[i][j][0]=(fg[i][j][2][0][0]*wss31_1[i][j][0]-fg[i][j][2][0][1]*wss31_1[i][j][1])+(fg[i][j][2][1][0]*wss32_1[i][j][0]-fg[i][j][2][1][1]*wss32_1[i][j][1])+(fg[i][j][2][2][0]*wss33_1[i][j][0]-fg[i][j][2][2][1]*wss33_1[i][j][1]);

wu3_1[i][j][1]=(fg[i][j][2][0][0]*wss31_1[i][j][1]+fg[i][j][2][0][1]*wss31_1[i][j][0])+(fg[i][j][2][1][0]*wss32_1[i][j][1]+fg[i][j][2][1][1]*wss32_1[i][j][0])+(fg[i][j][2][2][0]*wss33_1[i][j][1]+fg[i][j][2][2][1]*wss33_1[i][j][0]);

*/
wu1_1[i][j][0]=((fg[i][j][0][0][0]*wss31_1[i][j][0]-fg[i][j][0][0][1]*wss31_1[i][j][1])+(fg[i][j][0][1][0]*wss32_1[i][j][0]-fg[i][j][0][1][1]*wss32_1[i][j][1])+(fg[i][j][0][2][0]*wss33_1[i][j][0]-fg[i][j][0][2][1]*wss33_1[i][j][1]))+wu_1zeta[i][j][0][0];

wu1_1[i][j][1]=((fg[i][j][0][0][0]*wss31_1[i][j][1]+fg[i][j][0][0][1]*wss31_1[i][j][0])+(fg[i][j][0][1][0]*wss32_1[i][j][1]+fg[i][j][0][1][1]*wss32_1[i][j][0])+(fg[i][j][0][2][0]*wss33_1[i][j][1]+fg[i][j][0][2][1]*wss33_1[i][j][0]))+wu_1zeta[i][j][0][1];

wu2_1[i][j][0]=((fg[i][j][1][0][0]*wss31_1[i][j][0]-fg[i][j][1][0][1]*wss31_1[i][j][1])+(fg[i][j][1][1][0]*wss32_1[i][j][0]-fg[i][j][1][1][1]*wss32_1[i][j][1])+(fg[i][j][1][2][0]*wss33_1[i][j][0]-fg[i][j][1][2][1]*wss33_1[i][j][1]))+wu_1zeta[i][j][1][0];

wu2_1[i][j][1]=((fg[i][j][1][0][0]*wss31_1[i][j][1]+fg[i][j][1][0][1]*wss31_1[i][j][0])+(fg[i][j][1][1][0]*wss32_1[i][j][1]+fg[i][j][1][1][1]*wss32_1[i][j][0])+(fg[i][j][1][2][0]*wss33_1[i][j][1]+fg[i][j][1][2][1]*wss33_1[i][j][0]))+wu_1zeta[i][j][1][1];

wu3_1[i][j][0]=((fg[i][j][2][0][0]*wss31_1[i][j][0]-fg[i][j][2][0][1]*wss31_1[i][j][1])+(fg[i][j][2][1][0]*wss32_1[i][j][0]-fg[i][j][2][1][1]*wss32_1[i][j][1])+(fg[i][j][2][2][0]*wss33_1[i][j][0]-fg[i][j][2][2][1]*wss33_1[i][j][1]))+wu_1zeta[i][j][2][0];

wu3_1[i][j][1]=((fg[i][j][2][0][0]*wss31_1[i][j][1]+fg[i][j][2][0][1]*wss31_1[i][j][0])+(fg[i][j][2][1][0]*wss32_1[i][j][1]+fg[i][j][2][1][1]*wss32_1[i][j][0])+(fg[i][j][2][2][0]*wss33_1[i][j][1]+fg[i][j][2][2][1]*wss33_1[i][j][0]))+wu_1zeta[i][j][2][1];

}}

for(i=0;i<nx;i++){for(j=0;j<ny;j++){
         wu1x_1[i][j][0]=-wx1[i]*wu1_1[i][j][1];        wu1x_1[i][j][1]=wx1[i]*wu1_1[i][j][0];

         wu1y_1[i][j][0]=-wy1[j]*wu1_1[i][j][1];        wu1y_1[i][j][1]=wy1[j]*wu1_1[i][j][0];
/*
wu1z_1[i][j][0]=(kfg[i][j][0][0][0]*wss31_1[i][j][0]-kfg[i][j][0][0][1]*wss31_1[i][j][1])+(kfg[i][j][0][1][0]*wss32_1[i][j][0]-kfg[i][j][0][1][1]*wss32_1[i][j][1])+(kfg[i][j][0][2][0]*wss33_1[i][j][0]-kfg[i][j][0][2][1]*wss33_1[i][j][1]);      

wu1z_1[i][j][1]=(kfg[i][j][0][0][0]*wss31_1[i][j][1]+kfg[i][j][0][0][1]*wss31_1[i][j][0])+(kfg[i][j][0][1][0]*wss32_1[i][j][1]+kfg[i][j][0][1][1]*wss32_1[i][j][0])+(kfg[i][j][0][2][0]*wss33_1[i][j][1]+kfg[i][j][0][2][1]*wss33_1[i][j][0]);
*/
wu1z_1[i][j][0]=(kfg[i][j][0][0][0]*wss31_1[i][j][0]-kfg[i][j][0][0][1]*wss31_1[i][j][1])+(kfg[i][j][0][1][0]*wss32_1[i][j][0]-kfg[i][j][0][1][1]*wss32_1[i][j][1])+(kfg[i][j][0][2][0]*wss33_1[i][j][0]-kfg[i][j][0][2][1]*wss33_1[i][j][1])+wu1zeta_corr[i][j][0][0];

wu1z_1[i][j][1]=(kfg[i][j][0][0][0]*wss31_1[i][j][1]+kfg[i][j][0][0][1]*wss31_1[i][j][0])+(kfg[i][j][0][1][0]*wss32_1[i][j][1]+kfg[i][j][0][1][1]*wss32_1[i][j][0])+(kfg[i][j][0][2][0]*wss33_1[i][j][1]+kfg[i][j][0][2][1]*wss33_1[i][j][0])+wu1zeta_corr[i][j][0][1];


         wu2x_1[i][j][0]=-wx1[i]*wu2_1[i][j][1];        wu2x_1[i][j][1]=wx1[i]*wu2_1[i][j][0];

         wu2y_1[i][j][0]=-wy1[j]*wu2_1[i][j][1];        wu2y_1[i][j][1]=wy1[j]*wu2_1[i][j][0];
/*
wu2z_1[i][j][0]=(kfg[i][j][1][0][0]*wss31_1[i][j][0]-kfg[i][j][1][0][1]*wss31_1[i][j][1])+(kfg[i][j][1][1][0]*wss32_1[i][j][0]-kfg[i][j][1][1][1]*wss32_1[i][j][1])+(kfg[i][j][1][2][0]*wss33_1[i][j][0]-kfg[i][j][1][2][1]*wss33_1[i][j][1]);

wu2z_1[i][j][1]=(kfg[i][j][1][0][0]*wss31_1[i][j][1]+kfg[i][j][1][0][1]*wss31_1[i][j][0])+(kfg[i][j][1][1][0]*wss32_1[i][j][1]+kfg[i][j][1][1][1]*wss32_1[i][j][0])+(kfg[i][j][1][2][0]*wss33_1[i][j][1]+kfg[i][j][1][2][1]*wss33_1[i][j][0]);
*/
wu2z_1[i][j][0]=(kfg[i][j][1][0][0]*wss31_1[i][j][0]-kfg[i][j][1][0][1]*wss31_1[i][j][1])+(kfg[i][j][1][1][0]*wss32_1[i][j][0]-kfg[i][j][1][1][1]*wss32_1[i][j][1])+(kfg[i][j][1][2][0]*wss33_1[i][j][0]-kfg[i][j][1][2][1]*wss33_1[i][j][1])+wu1zeta_corr[i][j][1][0];

wu2z_1[i][j][1]=(kfg[i][j][1][0][0]*wss31_1[i][j][1]+kfg[i][j][1][0][1]*wss31_1[i][j][0])+(kfg[i][j][1][1][0]*wss32_1[i][j][1]+kfg[i][j][1][1][1]*wss32_1[i][j][0])+(kfg[i][j][1][2][0]*wss33_1[i][j][1]+kfg[i][j][1][2][1]*wss33_1[i][j][0])+wu1zeta_corr[i][j][1][1];


         wu3x_1[i][j][0]=-wx1[i]*wu3_1[i][j][1];        wu3x_1[i][j][1]=wx1[i]*wu3_1[i][j][0];

         wu3y_1[i][j][0]=-wy1[j]*wu3_1[i][j][1];        wu3y_1[i][j][1]=wy1[j]*wu3_1[i][j][0];
/*
wu3z_1[i][j][0]=(kfg[i][j][2][0][0]*wss31_1[i][j][0]-kfg[i][j][2][0][1]*wss31_1[i][j][1])+(kfg[i][j][2][1][0]*wss32_1[i][j][0]-kfg[i][j][2][1][1]*wss32_1[i][j][1])+(kfg[i][j][2][2][0]*wss33_1[i][j][0]-kfg[i][j][2][2][1]*wss33_1[i][j][1]);

wu3z_1[i][j][1]=(kfg[i][j][2][0][0]*wss31_1[i][j][1]+kfg[i][j][2][0][1]*wss31_1[i][j][0])+(kfg[i][j][2][1][0]*wss32_1[i][j][1]+kfg[i][j][2][1][1]*wss32_1[i][j][0])+(kfg[i][j][2][2][0]*wss33_1[i][j][1]+kfg[i][j][2][2][1]*wss33_1[i][j][0]);
*/
wu3z_1[i][j][0]=(kfg[i][j][2][0][0]*wss31_1[i][j][0]-kfg[i][j][2][0][1]*wss31_1[i][j][1])+(kfg[i][j][2][1][0]*wss32_1[i][j][0]-kfg[i][j][2][1][1]*wss32_1[i][j][1])+(kfg[i][j][2][2][0]*wss33_1[i][j][0]-kfg[i][j][2][2][1]*wss33_1[i][j][1])+wu1zeta_corr[i][j][2][0];

wu3z_1[i][j][1]=(kfg[i][j][2][0][0]*wss31_1[i][j][1]+kfg[i][j][2][0][1]*wss31_1[i][j][0])+(kfg[i][j][2][1][0]*wss32_1[i][j][1]+kfg[i][j][2][1][1]*wss32_1[i][j][0])+(kfg[i][j][2][2][0]*wss33_1[i][j][1]+kfg[i][j][2][2][1]*wss33_1[i][j][0])+wu1zeta_corr[i][j][2][1];

}}

         fftw_execute(u1x_1idft);
         fftw_execute(u1y_1idft);
         fftw_execute(u2x_1idft);
         fftw_execute(u2y_1idft);
         fftw_execute(u3x_1idft);
         fftw_execute(u3y_1idft); 
         fftw_execute(u1z_1idft);
         fftw_execute(u2z_1idft);
         fftw_execute(u3z_1idft);

for(i=0;i<nx;i++){for(j=0;j<ny;j++){
         u1x_1[i][j][0]=u1x_1in2[i][j][0]/(nx*ny);              u1x_1[i][j][1]=u1x_1in2[i][j][1]/(nx*ny);
         u1y_1[i][j][0]=u1y_1in2[i][j][0]/(nx*ny);              u1y_1[i][j][1]=u1y_1in2[i][j][1]/(nx*ny);
         u1z_1[i][j][0]=u1z_1in2[i][j][0]/(nx*ny);              u1z_1[i][j][1]=u1z_1in2[i][j][1]/(nx*ny);

         u2x_1[i][j][0]=u2x_1in2[i][j][0]/(nx*ny);              u2x_1[i][j][1]=u2x_1in2[i][j][1]/(nx*ny);
         u2y_1[i][j][0]=u2y_1in2[i][j][0]/(nx*ny);              u2y_1[i][j][1]=u2y_1in2[i][j][1]/(nx*ny);
         u2z_1[i][j][0]=u2z_1in2[i][j][0]/(nx*ny);              u2z_1[i][j][1]=u2z_1in2[i][j][1]/(nx*ny);

         u3x_1[i][j][0]=u3x_1in2[i][j][0]/(nx*ny);              u3x_1[i][j][1]=u3x_1in2[i][j][1]/(nx*ny);
         u3y_1[i][j][0]=u3y_1in2[i][j][0]/(nx*ny);              u3y_1[i][j][1]=u3y_1in2[i][j][1]/(nx*ny);
         u3z_1[i][j][0]=u3z_1in2[i][j][0]/(nx*ny);              u3z_1[i][j][1]=u3z_1in2[i][j][1]/(nx*ny);
                                                                                      }}


for(i=0;i<nx;i++){for(j=0;j<ny;j++){
u1_1corr[i][j][0] = (u1z_1[i][j][0]*(hight_t[i][j][0]-hbar))-(u1z_1[i][j][1]*hight_t[i][j][1]);
u1_1corr[i][j][1] = (u1z_1[i][j][0]*hight_t[i][j][1])+(u1z_1[i][j][1]*(hight_t[i][j][0]-hbar));

u2_1corr[i][j][0] = (u2z_1[i][j][0]*(hight_t[i][j][0]-hbar))-(u2z_1[i][j][1]*hight_t[i][j][1]);
u2_1corr[i][j][1] = (u2z_1[i][j][0]*hight_t[i][j][1])+(u2z_1[i][j][1]*(hight_t[i][j][0]-hbar));

}}


fftw_execute(u1_1corrdft);
fftw_execute(u2_1corrdft);

for(i=0;i<nx;i++){for(j=0;j<ny;j++){
wu1x_1corr[i][j][0] = -wx1[i]*wu1_1corr[i][j][1];        wu1x_1corr[i][j][1] = wx1[i]*wu1_1corr[i][j][0];
wu2y_1corr[i][j][0] = -wy1[j]*wu2_1corr[i][j][1];        wu2y_1corr[i][j][1] = wy1[j]*wu2_1corr[i][j][0];

}}



for(i=0;i<nx;i++){for(j=0;j<ny;j++){        
         ss11_1[i][j][0]=(c11*u1x_1[i][j][0]+c12*u2y_1[i][j][0]+c12*u3z_1[i][j][0]);
         ss11_1[i][j][1]=(c11*u1x_1[i][j][1]+c12*u2y_1[i][j][1]+c12*u3z_1[i][j][1]);

         ss12_1[i][j][0]=c44*(u1y_1[i][j][0]+u2x_1[i][j][0]);
        ss12_1[i][j][1]=c44*(u1y_1[i][j][1]+u2x_1[i][j][1]);
 
         ss22_1[i][j][0]=c12*u1x_1[i][j][0]+c11*u2y_1[i][j][0]+c12*u3z_1[i][j][0];
         ss22_1[i][j][1]=c12*u1x_1[i][j][1]+c11*u2y_1[i][j][1]+c12*u3z_1[i][j][1];

         ss31_2[i][j][0]=(c11*u1x_1[i][j][0]+c12*u2y_1[i][j][0]+c12*u3z_1[i][j][0])*hx[i][j]+c44*(u1y_1[i][j][0]+u2x_1[i][j][0])*hy[i][j]-ss31_1corr[i][j][0];
         ss31_2[i][j][1]=(c11*u1x_1[i][j][1]+c12*u2y_1[i][j][1]+c12*u3z_1[i][j][1])*hx[i][j]+c44*(u1y_1[i][j][1]+u2x_1[i][j][1])*hy[i][j];

         ss32_2[i][j][0]=c44*(u1y_1[i][j][0]+u2x_1[i][j][0])*hx[i][j]+(c12*u1x_1[i][j][0]+c11*u2y_1[i][j][0]+c12*u3z_1[i][j][0])*hy[i][j]-ss32_1corr[i][j][0];
         ss32_2[i][j][1]=c44*(u1y_1[i][j][1]+u2x_1[i][j][1])*hx[i][j]+(c12*u1x_1[i][j][1]+c11*u2y_1[i][j][1]+c12*u3z_1[i][j][1])*hy[i][j];
 
         ss33_2[i][j][0]=ss_0*pow(hx[i][j],2)+ss_0*pow(hy[i][j],2)-ss33_1corr[i][j][0];
         ss33_2[i][j][1]=0;    
                                                                                       }}

         fftw_execute(ss31_2dft);
         fftw_execute(ss32_2dft);
         fftw_execute(ss33_2dft);


for(i=0;i<nx;i++){for(j=0;j<ny;j++){
/*
wu1_2[i][j][0]=(fg[i][j][0][0][0]*wss31_2[i][j][0]-fg[i][j][0][0][1]*wss31_2[i][j][1])+(fg[i][j][0][1][0]*wss32_2[i][j][0]-fg[i][j][0][1][1]*wss32_2[i][j][1])+(fg[i][j][0][2][0]*wss33_2[i][j][0]-fg[i][j][0][2][1]*wss33_2[i][j][1]);

wu1_2[i][j][1]=(fg[i][j][0][0][0]*wss31_2[i][j][1]+fg[i][j][0][0][1]*wss31_2[i][j][0])+(fg[i][j][0][1][0]*wss32_2[i][j][1]+fg[i][j][0][1][1]*wss32_2[i][j][0])+(fg[i][j][0][2][0]*wss33_2[i][j][1]+fg[i][j][0][2][1]*wss33_2[i][j][0]);

wu2_2[i][j][0]=(fg[i][j][1][0][0]*wss31_2[i][j][0]-fg[i][j][1][0][1]*wss31_2[i][j][1])+(fg[i][j][1][1][0]*wss32_2[i][j][0]-fg[i][j][1][1][1]*wss32_2[i][j][1])+(fg[i][j][1][2][0]*wss33_2[i][j][0]-fg[i][j][1][2][1]*wss33_2[i][j][1]);

wu2_2[i][j][1]=(fg[i][j][1][0][0]*wss31_2[i][j][1]+fg[i][j][1][0][1]*wss31_2[i][j][0])+(fg[i][j][1][1][0]*wss32_2[i][j][1]+fg[i][j][1][1][1]*wss32_2[i][j][0])+(fg[i][j][1][2][0]*wss33_2[i][j][1]+fg[i][j][1][2][1]*wss33_2[i][j][0]);
*/
wu3_2[i][j][0]=(fg[i][j][2][0][0]*wss31_2[i][j][0]-fg[i][j][2][0][1]*wss31_2[i][j][1])+(fg[i][j][2][1][0]*wss32_2[i][j][0]-fg[i][j][2][1][1]*wss32_2[i][j][1])+(fg[i][j][2][2][0]*wss33_2[i][j][0]-fg[i][j][2][2][1]*wss33_2[i][j][1]);

wu3_2[i][j][1]=(fg[i][j][2][0][0]*wss31_2[i][j][1]+fg[i][j][2][0][1]*wss31_2[i][j][0])+(fg[i][j][2][1][0]*wss32_2[i][j][1]+fg[i][j][2][1][1]*wss32_2[i][j][0])+(fg[i][j][2][2][0]*wss33_2[i][j][1]+fg[i][j][2][2][1]*wss33_2[i][j][0]);


wu1_2[i][j][0]=(fg[i][j][0][0][0]*wss31_2[i][j][0]-fg[i][j][0][0][1]*wss31_2[i][j][1])+(fg[i][j][0][1][0]*wss32_2[i][j][0]-fg[i][j][0][1][1]*wss32_2[i][j][1])+(fg[i][j][0][2][0]*wss33_2[i][j][0]-fg[i][j][0][2][1]*wss33_2[i][j][1])+wu_2zeta[i][j][0][0];

wu1_2[i][j][1]=(fg[i][j][0][0][0]*wss31_2[i][j][1]+fg[i][j][0][0][1]*wss31_2[i][j][0])+(fg[i][j][0][1][0]*wss32_2[i][j][1]+fg[i][j][0][1][1]*wss32_2[i][j][0])+(fg[i][j][0][2][0]*wss33_2[i][j][1]+fg[i][j][0][2][1]*wss33_2[i][j][0])+wu_2zeta[i][j][0][1];

wu2_2[i][j][0]=(fg[i][j][1][0][0]*wss31_2[i][j][0]-fg[i][j][1][0][1]*wss31_2[i][j][1])+(fg[i][j][1][1][0]*wss32_2[i][j][0]-fg[i][j][1][1][1]*wss32_2[i][j][1])+(fg[i][j][1][2][0]*wss33_2[i][j][0]-fg[i][j][1][2][1]*wss33_2[i][j][1])+wu_2zeta[i][j][1][0];

wu2_2[i][j][1]=(fg[i][j][1][0][0]*wss31_2[i][j][1]+fg[i][j][1][0][1]*wss31_2[i][j][0])+(fg[i][j][1][1][0]*wss32_2[i][j][1]+fg[i][j][1][1][1]*wss32_2[i][j][0])+(fg[i][j][1][2][0]*wss33_2[i][j][1]+fg[i][j][1][2][1]*wss33_2[i][j][0])+wu_2zeta[i][j][1][1];

                                                                                                              }}
     

for(i=0;i<nx;i++){for(j=0;j<ny;j++){
         wu1x_2[i][j][0]=-wx1[i]*wu1_2[i][j][1];        wu1x_2[i][j][1]=wx1[i]*wu1_2[i][j][0];
         wu2y_2[i][j][0]=-wy1[j]*wu2_2[i][j][1];        wu2y_2[i][j][1]=wy1[j]*wu2_2[i][j][0];

                                                                                                                }}

//els_0 need not to be added to the the complete expression of elasticity cz it would become zero after differentiation
//each wels_1 & wels_2 are devided by els_0 to make them in terms of reduced unit  
//els_0 is calculated in linear_OneTCaladd.c
//         els_0=ss_0*sn_0;


for(i=0;i<nx;i++){for(j=0;j<ny;j++){
         wels_1[i][j][0]=2*.5*ss_0*(wu1x_1[i][j][0]+wu2y_1[i][j][0]);
         wels_1[i][j][1]=2*.5*ss_0*(wu1x_1[i][j][1]+wu2y_1[i][j][1]);
 
         wels_1corr[i][j][0]=2*.5*ss_0*(wu1x_1corr[i][j][0]+wu2y_1corr[i][j][0]);
         wels_1corr[i][j][1]=2*.5*ss_0*(wu1x_1corr[i][j][1]+wu2y_1corr[i][j][1]);


//(below is previous one)
//els_temp_2[i][j][0]=els_0+((.5*(ss11_1[i][j][0]*u1x_1[i][j][0]-ss11_1[i][j][1]*u1x_1[i][j][1]))+(.5*(ss12_1[i][j][0]*u1y_1[i][j][0]-ss12_1[i][j][1]*u1y_1[i][j][1]))+(.5*(ss12_1[i][j][0]*u2x_1[i][j][0]-ss12_1[i][j][1]*u2x_1[i][j][1]))+(.5*(ss22_1[i][j][0]*u2y_1[i][j][0]-ss22_1[i][j][1]*u2y_1[i][j][1]))+(ss_0*hx[i][j]*(u3x_1[i][j][0]+u1z_1[i][j][0]))+(ss_0*hy[i][j]*(u3y_1[i][j][0]+u2z_1[i][j][0])))/(2*(1+vv)*els_0);

els_temp_2[i][j][0]=/*els_0+*/((.5*(ss11_1[i][j][0]*u1x_1[i][j][0]-ss11_1[i][j][1]*u1x_1[i][j][1]))+(.5*(ss12_1[i][j][0]*u1y_1[i][j][0]-ss12_1[i][j][1]*u1y_1[i][j][1]))+(.5*(ss12_1[i][j][0]*u2x_1[i][j][0]-ss12_1[i][j][1]*u2x_1[i][j][1]))+(.5*(ss22_1[i][j][0]*u2y_1[i][j][0]-ss22_1[i][j][1]*u2y_1[i][j][1]))+.5*(ss_0*hx[i][j]*(u3x_1[i][j][0]+u1z_1[i][j][0]))+.5*(ss_0*hy[i][j]*(u3y_1[i][j][0]+u2z_1[i][j][0])));

//els_temp_2[i][j][1]=0;

els_temp_2[i][j][1]=/*els_0+*/((.5*(ss11_1[i][j][0]*u1x_1[i][j][1]+ss11_1[i][j][1]*u1x_1[i][j][0]))+(.5*(ss12_1[i][j][0]*u1y_1[i][j][1]+ss12_1[i][j][1]*u1y_1[i][j][0]))+(.5*(ss12_1[i][j][0]*u2x_1[i][j][1]+ss12_1[i][j][1]*u2x_1[i][j][0]))+(.5*(ss22_1[i][j][0]*u2y_1[i][j][1]+ss22_1[i][j][1]*u2y_1[i][j][0]))+.5*(ss_0*hx[i][j]*(u3x_1[i][j][1]+u1z_1[i][j][1]))+.5*(ss_0*hy[i][j]*(u3y_1[i][j][1]+u2z_1[i][j][1])));

                                                                                                                }}

         fftw_execute(els_temp_2dft);

for(i=0;i<nx;i++){for(j=0;j<ny;j++){
 wels_2[i][j][0]=2*.5*((ss_0*wu1x_2[i][j][0])+(ss_0*wu2y_2[i][j][0]))+(wels_temp_2[i][j][0]);
 wels_2[i][j][1]=2*.5*((ss_0*wu1x_2[i][j][1])+(ss_0*wu2y_2[i][j][1]))+(wels_temp_2[i][j][1]);                          

 wels[i][j][0]=(wels_1[i][j][0]+(wels_1corr[i][j][0])+wels_2[i][j][0])/(2*(1+vv)*els_0);
 wels[i][j][1]=(wels_1[i][j][1]+(wels_1corr[i][j][1])+wels_2[i][j][1])/(2*(1+vv)*els_0);

                                                                                                                }}

//fftw_execute(ells_1idft);
//fftw_execute(ells_1corridft);
//fftw_execute(ells_2idft);

fftw_execute(ells_idft);
evolution_eqn_elas();
}
