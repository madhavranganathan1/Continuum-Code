#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <complex.h>
#include <time.h>
#define LX 128 //dim of system
#define LY 128 //dim of system
#define dx 1 //grid sepration
#define dy 1 //grid sepration
#define nx 128 //no of modes
#define ny 128 //no of modes
#define nxd1 nx/3
#define nxd2 2*nx/3+2
#define nyd1 ny/3
#define nyd2 2*ny/3+2
#define M  16   //no of points for complex mean
#define tot_step 500 //Reduced total timesteps for quick running of code 
#define Nprint 50 //Print output every Nprint steps 
//#define tot_step 50000 
#define dt .01  //time step
#define yf 1
#define cw 0.09//.058//.07//.05
#define dw 0.03129//.0083//.01//.01//.005
#define a100 .48//.583463
#define a105 .02//.0307489
#define n100 10
#define n105 /*20*/15
#define e100 /*.01*/.1 
#define e105 /*.0006*/.001
#define flux 0.00
#define vv 0.2664//0.28
#define c11 1.32//1.42//1.4265//1.511 
#define c12 .4815//.51575//.53675
#define c44 .70/*0.904250 .226062 0.271275 0.44102 0.452577  1.8085 1.35637.7234  0.293881 0.135637 1.94850.487174  0.491996 .746006*/
#define check 200000000
#define ms .02 //misit for 50% Ge in film
#define hbar .032//.105
#define zeta_A 0.0
#define zeta_l 32.0
#define zeta_w 0.2
#define h0_avg 0.6
#define N_G 0.004
// zenar anisotropy ratio is 1.001
//=============================================

double dyn105_1temp[nx][ny],dyn105_1[nx][ny],dyn105_2temp[nx][ny],dyn105_2[nx][ny],dyn105_3temp[nx][ny],dyn105_3[nx][ny],dyn105_4temp[nx][ny],dyn105_4[nx][ny],dyn100_temp[nx][ny],dyn100[nx][ny];
double mu100_temp[nx][ny],mu105_1temp[nx][ny],mu105_2temp[nx][ny],mu105_3temp[nx][ny],mu105_4temp[nx][ny],muyn_temp[nx][ny],muyn[nx][ny];
double mu100[nx][ny],mu105_1[nx][ny],mu105_2[nx][ny],mu105_3[nx][ny],mu105_4[nx][ny];
double mu105_1t1[nx][ny],mu105_1t2[nx][ny];
double jj[nx][ny],kk1[nx][ny],kk2[nx][ny],ff[nx][ny];
double ss00_[nx][ny];
fftw_complex cals_temp[nx][ny];
//=============================================

double fg[nx][ny][3][3][2],kfg[nx][ny][3][3][2],kkfg[nx][ny][3][3][2],lkfg[nx][ny][3][3][2];
double ss_0,sn_0,els_0;
double begin,end,time_spent;
double hight[nx][ny],wx1[nx],wy1[ny],wx2[nx],wy2[ny],kx[nx],ky[ny];
fftw_complex hight_t[nx][ny],Nhight_t[nx][ny];
fftw_complex whx[nx][ny],hx_t[nx][ny],why[nx][ny],hy_t[nx][ny],whxx[nx][ny],hxx_t[nx][ny],whxy[nx][ny],hxy_t[nx][ny],whyy[nx][ny],hyy_t[nx][ny],WHxx[nx][ny],WHxy[nx][ny],WHyy[nx][ny],Hxx[nx][ny],Hxy[nx][ny],Hyy[nx][ny],fn1[nx][ny],fn2[nx][ny],fn3[nx][ny],wfn1[nx][ny],wfn2[nx][ny],wfn3[nx][ny],WHfn1[nx][ny],WHfn2[nx][ny],WHfn3[nx][ny];
double hx[nx][ny],hy[nx][ny],hxx[nx][ny],hyy[nx][ny],hxy[nx][ny],p[nx][ny],k[nx][ny],W[nx][ny],Hfn1[nx][ny],Hfn2[nx][ny],Hfn3[nx][ny];
fftw_complex ss31_1[nx][ny],ss32_1[nx][ny],ss33_1[nx][ny];

fftw_complex in[nx][ny],out[nx][ny],in2[nx][ny],hin2[nx][ny],hxin[nx][ny],hxout[nx][ny],hxin2[nx][ny],hyin[nx][ny],hyout[nx][ny],hyin2[nx][ny],hxxin[nx][ny],hxxout[nx][ny],hxxin2[nx][ny],hxyin[nx][ny],hxyout[nx][ny],hxyin2[nx][ny],hyyin[nx][ny],hyyout[nx][ny],hyyin2[nx][ny],Hxxout[nx][ny],Hxxin2[nx][ny],Hxyout[nx][ny],Hxyin2[nx][ny],Hyyout[nx][ny],Hyyin2[nx][ny],gin[nx][ny],gout[nx][ny],gin2[nx][ny],wh[nx][ny],wcals[nx][ny],Hfn1out[nx][ny],Hfn2out[nx][ny],Hfn3out[nx][ny],Hfn1in2[nx][ny],Hfn2in2[nx][ny],Hfn3in2[nx][ny];

fftw_plan dft,idft,hxdft,hxidft,hydft,hyidft,hxxdft,hxxidft,hxydft,hxyidft,hyydft,hyyidft,gdft,gidft,Hxxidft,Hxyidft,Hyyidft,Hfn1idft,Hfn2idft,Hfn3idft,fn1dft,fn2dft,fn3dft;

fftw_plan u1x_1idft,u2x_1idft,u3x_1idft,u1y_1idft,u2y_1idft,u3y_1idft,u1z_1idft,u2z_1idft,ss31_2dft,ss32_2dft,ss33_2dft,els_temp_2dft;
fftw_complex L[nx][ny],root[M],k1[nx][ny],k2[nx][ny],k3[nx][ny],Nk1[nx][ny],Nk2[nx][ny],Nk3[nx][ny],wh_1[nx][ny],Nv[nx][ny];
double E[nx][ny],E2[nx][ny],Q[nx][ny],f1[nx][ny],f2[nx][ny],wg[nx][ny],f3[nx][ny];
double Q_1[nx][ny][M],f1_1[nx][ny][M],f2_1[nx][ny][M],f3_1[nx][ny][M];
fftw_complex LR_1[nx][ny][M],LR_2[nx][ny][M],LR[nx][ny][M];
//dft  = fftw_plan_dft_2d(nx,ny,&in[0][0],&out[0][0],-1,FFTW_MEASURE);
//idft = fftw_plan_dft_2d(nx,ny,&out[0][0],&in2[0][0],+1,FFTW_MEASURE);

//double yh[nx][ny],yn[nx][ny];

double yhn_hx[nx][ny],yhn_hy[nx][ny],yhn_hxhx[nx][ny],yhn_hxhy[nx][ny],yhn_hyhy[nx][ny],yhn_h[nx][ny];
double yhn100_hx[nx][ny],yhn105_1_hx[nx][ny],yhn105_2_hx[nx][ny],yhn105_3_hx[nx][ny],yhn105_4_hx[nx][ny];
double yhn100_hy[nx][ny],yhn105_1_hy[nx][ny],yhn105_2_hy[nx][ny],yhn105_3_hy[nx][ny],yhn105_4_hy[nx][ny];
double yhn100_hxhx[nx][ny],yhn105_1_hxhx[nx][ny],yhn105_2_hxhx[nx][ny],yhn105_3_hxhx[nx][ny],yhn105_4_hxhx[nx][ny];
double yhn100_hxhy[nx][ny],yhn105_1_hxhy[nx][ny],yhn105_2_hxhy[nx][ny],yhn105_3_hxhy[nx][ny],yhn105_4_hxhy[nx][ny],yhn105_11_hxhy[nx][ny],yhn105_21_hxhy[nx][ny],yhn105_31_hxhy[nx][ny],yhn105_41_hxhy[nx][ny];
double yhn100_hyhy[nx][ny],yhn105_1_hyhy[nx][ny],yhn105_2_hyhy[nx][ny],yhn105_3_hyhy[nx][ny],yhn105_4_hyhy[nx][ny];
double rr105_1_[nx][ny],rr105_2_[nx][ny],rr105_3_[nx][ny],rr105_4_[nx][ny];
double ss105_1_[nx][ny],ss105_2_[nx][ny],ss105_3_[nx][ny],ss105_4_[nx][ny];

double rr100_[nx][ny],ss100_[nx][ny];
double t1[nx][ny],t2[nx][ny],t3[nx][ny],t4[nx][ny],t5[nx][ny],t6[nx][ny],t7[nx][ny],t8[nx][ny];
//***************************************************************************************************
fftw_complex wss31_1[nx][ny],wss32_1[nx][ny],wss33_1[nx][ny];
fftw_complex wu1_1[nx][ny],wu2_1[nx][ny],wu3_1[nx][ny],u1_1[nx][ny],u2_1[nx][ny],u3_1[nx][ny];
fftw_complex wu1x_1[nx][ny],wu1y_1[nx][ny],wu1z_1[nx][ny],u1x_1out[nx][ny],u1y_1out[nx][ny],u1z_1out[nx][ny],u1x_1in2[nx][ny],u1y_1in2[nx][ny],u1z_1in2[nx][ny],u1x_1[nx][ny],u1y_1[nx][ny],u1z_1[nx][ny];
fftw_complex wu2x_1[nx][ny],wu2y_1[nx][ny],wu2z_1[nx][ny],u2x_1out[nx][ny],u2y_1out[nx][ny],u2z_1out[nx][ny],u2x_1in2[nx][ny],u2y_1in2[nx][ny],u2z_1in2[nx][ny],u2x_1[nx][ny],u2y_1[nx][ny],u2z_1[nx][ny];
fftw_complex wu3x_1[nx][ny],wu3y_1[nx][ny],wu3z_1[nx][ny],u3x_1out[nx][ny],u3y_1out[nx][ny],u3z_1out[nx][ny],u3x_1in2[nx][ny],u3y_1in2[nx][ny],u3z_1in2[nx][ny],u3x_1[nx][ny],u3y_1[nx][ny],u3z_1[nx][ny];
fftw_complex ss11_1[nx][ny],ss12_1[nx][ny],ss22_1[nx][ny];
fftw_complex ss31_2[nx][ny],ss32_2[nx][ny],ss33_2[nx][ny],ss31_2out[nx][ny],ss32_2out[nx][ny],ss33_2out[nx][ny],wss31_2[nx][ny],wss32_2[nx][ny],wss33_2[nx][ny];
fftw_complex wu1_2[nx][ny],wu2_2[nx][ny],wu3_2[nx][ny];
fftw_complex wu1x_2[nx][ny],wu2y_2[nx][ny];
fftw_complex els_temp_2[nx][ny],wels_temp_2out[nx][ny],wels_temp_2[nx][ny], wels_2[nx][ny],wels_1[nx][ny],wels[nx][ny];
fftw_plan u1x_1idft,u1y_1idft,u1z_1idft,u2x_1idft,u2y_1idft,u2z_1idft,u3x_1idft,u3y_1idft,u3z_1idft;
fftw_plan ss31_2dft,ss32_2dft,ss33_2dft;
fftw_plan els_temp_2dft;

fftw_complex u1_1out[nx][ny],u2_1out[nx][ny],u3_1out[nx][ny],u1_1in2[nx][ny],u2_1in2[nx][ny],u3_1in2[nx][ny],u1_1[nx][ny],u2_1[nx][ny],u3_1[nx][ny],u1_2out[nx][ny],u2_2out[nx][ny],u3_2out[nx][ny],u1_2in2[nx][ny],u2_2in2[nx][ny],u3_2in2[nx][ny],u1_2[nx][ny],u2_2[nx][ny],u3_2[nx][ny];

fftw_plan u1_1idft,u2_1idft,u3_1idft,u1_2idft,u2_2idft,u3_2idft;
/*
u1_1idft= fftw_plan_dft_2d(nx,ny,&u1_1out[0][0],&u1_1in2[0][0],+1,FFTW_MEASURE);
u2_1idft= fftw_plan_dft_2d(nx,ny,&u2_1out[0][0],&u2_1in2[0][0],+1,FFTW_MEASURE);
u3_1idft= fftw_plan_dft_2d(nx,ny,&u3_1out[0][0],&u3_1in2[0][0],+1,FFTW_MEASURE);
*/
fftw_complex els_temp[nx][ny],els[nx][ny];
fftw_plan els_idft;
double els_tot,els_avg[tot_step],surf_tot,surf_avg[tot_step];
//**************************************************************************************************
fftw_complex wss31_1corr[nx][ny],wss32_1corr[nx][ny],wss33_1corr[nx][ny],ss31_1corr[nx][ny],ss32_1corr[nx][ny],ss33_1corr[nx][ny],ss31_1corrin2[nx][ny],ss32_1corrin2[nx][ny],ss33_1corrin2[nx][ny];
fftw_plan ss31_1corridft,ss32_1corridft,ss33_1corridft;

fftw_complex u1_1corr[nx][ny],u2_1corr[nx][ny],u3_1corr[nx][ny],u1_1corrout[nx][ny],u2_1corrout[nx][ny],u3_1corrout[nx][ny],wu1_1corr[nx][ny],wu2_1corr[nx][ny],wu3_1corr[nx][ny];
fftw_complex wu1zz_1[nx][ny],wu2zz_1[nx][ny],wu3zz_1[nx][ny],u1zz_1in2[nx][ny],u2zz_1in2[nx][ny],u3zz_1in2[nx][ny],u1zz_1[nx][ny],u2zz_1[nx][ny],u3zz_1[nx][ny],u1z_1corrout[nx][ny],u2z_1corrout[nx][ny],u3z_1corrout[nx][ny],u1z_1corr[nx][ny],u2z_1corr[nx][ny],u3z_1corr[nx][ny]; 

fftw_plan u1_1corrdft,u2_1corrdft,u3_1corrdft,u1zz_1idft,u2zz_1idft,u3zz_1idft,u1z_1corrdft,u2z_1corrdft,u3z_1corrdft;

fftw_complex wu1x_1corr[nx][ny],wu1y_1corr[nx][ny],wu1z_1corr[nx][ny],wu2x_1corr[nx][ny],wu2y_1corr[nx][ny],wu2z_1corr[nx][ny],wu3x_1corr[nx][ny],wu3y_1corr[nx][ny],wu3z_1corr[nx][ny];

fftw_complex wels_1corr[nx][ny];

double sn33_0;
fftw_complex wels_1cont[nx][ny],wels_1corr_cont[nx][ny],wels_2cont[nx][ny],wu3z_2[nx][ny];
fftw_plan ells_1idft,ells_1corridft,ells_2idft,ells_idft;
fftw_complex wels_1[nx][ny],wels_1corr[nx][ny],wels_2[nx][ny],ells_1[nx][ny],ells_1corr[nx][ny],ells_2[nx][ny],wels[nx][ny],ells[nx][ny];
//************************************************************************************************
fftw_plan ss31_1corrdft,ss32_1corrdft,ss33_1corrdft;
fftw_complex ss31_1corrin21[nx][ny],ss32_1corrin21[nx][ny],ss33_1corrin21[nx][ny],ss31_11corr[nx][ny],ss32_11corr[nx][ny],ss33_11corr[nx][ny];
double ss32_2temp[nx][ny];
double test_height[nx][ny];
//fftw_complex whcon[nx][ny],hcon_in2[nx][ny],hcon[nx][ny],whxcon[nx][ny],hxcon_in2[nx][ny],hxcon[nx][ny],whycon[nx][ny],hycon_in2[nx][ny],hycon[nx][ny];
//fftw_complex wu1x_1con[nx][ny],wu2x_1con[nx][ny],wu3x_1con[nx][ny],u1x_1conin2[nx][ny],u2x_1conin2[nx][ny],u3x_1conin2[nx][ny],u1x_1con[nx][ny],u2x_1con[nx][ny],u3x_1con[nx][ny];
//fftw_complex wu1y_1con[nx][ny],wu2y_1con[nx][ny],wu3y_1con[nx][ny],u1y_1conin2[nx][ny],u2y_1conin2[nx][ny],u3y_1conin2[nx][ny],u1y_1con[nx][ny],u2y_1con[nx][ny],u3y_1con[nx][ny];
//fftw_complex wu1z_1con[nx][ny],wu2z_1con[nx][ny],wu3z_1con[nx][ny],u1z_1conin2[nx][ny],u2z_1conin2[nx][ny],u3z_1conin2[nx][ny],u1z_1con[nx][ny],u2z_1con[nx][ny],u3z_1con[nx][ny];
//fftw_complex ss11_1con[nx][ny],ss12_1con[nx][ny],ss22_1con[nx][ny];

//fftw_plan hconidft,hxconidft,hyconidft;
//fftw_plan u1z_1conidft,u2z_1conidft,u3z_1conidft;
//fftw_plan u1y_1conidft,u2y_1conidft,u3y_1conidft;
//fftw_plan u1x_1conidft,u2x_1conidft,u3x_1conidft;
//fftw_plan els_temp_2dftcon;
//fftw_complex els_temp_2con[nx][ny],wels_temp_2outcon[nx][ny],wels_temp_2con[nx][ny];

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Interface Variations
double r1,u_0;
double complex CMatSum0,CMatSum1,CMatSum2,CMatSum3,CMatSum4,CDet1,CDet2;
double complex CValpha[nx][ny][3][6],CVdd[nx][ny][3][6];

fftw_complex zeta_t[nx][ny],zeta_out[nx][ny],wzeta_t[nx][ny];
fftw_plan zeta_dft;
double complex wzeta[nx][ny];
double v[nx][ny][6][2],alpha[nx][ny][3][6][2],Ialpha[nx][ny][3][3][2],dd[nx][ny][3][6][2],Idd[nx][ny][3][3][2],praa[nx][ny][3][6][2],Iaatemp[nx][ny][3][3][2],prdd[nx][ny][3][6][2],praaMprdd[nx][ny][3][6][2],aMd[nx][ny][3][3][2];
double complex  CDet_alpha[nx][ny],CIalpha_test[nx][ny][3][3],Cdd_test[nx][ny][3][3],CIdd_test[nx][ny][3][3];
double UU[nx][ny][3][2][2][2],ddUU[nx][ny][3][2][2][2],UU1[nx][ny][3][3][2],ddUU1[nx][ny][3][3][2],ddUU2[nx][ny][3][3][2],RR[nx][ny][3][2][2][2],RR1[nx][ny][3][3][2],RR2[nx][ny][3][3][2];
double Co1[nx][ny][3][3][3][2];

double complex Cv[nx][ny][6],Calpha[nx][ny][3][6],CIalpha[nx][ny][3][3],Cdd[nx][ny][3][6],CIdd[nx][ny][3][3],CpraaIdd[nx][ny][3][3],Cpraa[nx][ny][3][6],CIaatemp[nx][ny][3][3],Cprdd[nx][ny][3][6],CpraaMprdd[nx][ny][3][6],CaMd[nx][ny][3][3];

double complex CvpraaIdd[nx][ny][3][3],CvprddIdd[nx][ny][3][3];
double complex CUU[nx][ny][3][2][2],CddUU[nx][ny][3][2][2],CUU1[nx][ny][3][3],CddUU1[nx][ny][3][3],CddUU2[nx][ny][3][3],CRR[nx][ny][3][2][2],CRR1[nx][ny][3][3],CRR2[nx][ny][3][3];
double complex CCo1[nx][ny][3][3][3],CCo2[nx][ny][3][3][3],CCo21[nx][ny][3][3],CCo22[nx][ny][3][3]/*,CaCCo21[nx][ny][3][3]*/,CDetdd[nx][ny],CDetaMd[nx][ny];
double complex CCo21_temp[nx][ny][3][3][3],CCo22_temp[nx][ny][3][3],CCo[nx][ny][6][3];
//double complex CCotest[nx][ny][3][3],CaCCotest[nx][ny][3][3],CaCCotestTr[nx][ny][3][3],CCo22Trtest[nx][ny][3][3];
double complex Ct1[nx][ny][3][3];

double complex wp1[nx][ny][3],Cwu_1zeta[nx][ny][3];
fftw_complex wu_1zeta[nx][ny][3];
fftw_complex testin[nx][ny],testout[nx][ny],testin2[nx][ny];
fftw_plan testdft,testidft;

//Second Order
double complex CXX21[nx][ny][3][3],CXX22[nx][ny][3][3],CXX[nx][ny][6][3];
double complex CRtemp1[nx][ny][3][6],CRtemp2[nx][ny][3][6],CRtemp3[nx][ny][3][3],CRtemp4[nx][ny][3],CRtemp5[nx][ny][3],CRtemp6[nx][ny][3];
double complex CLtemp1[nx][ny][3][6],CLtemp2[nx][ny][3][6],CLtemp3[nx][ny][3][3],CLtemp4[nx][ny][3],CLtemp5[nx][ny][3],CLtemp6[nx][ny][3];
//fftw_complex wRtemp6[nx][ny][3], Rtemp6_in2[nx][ny][3],Rtemp6[nx][ny][3],wLtemp6[nx][ny][3], Ltemp6_in2[nx][ny][3],Ltemp6[nx][ny][3];
//fftw_plan Rtemp6idft,Ltemp6idft;

fftw_complex wRtemp6_1[nx][ny],wRtemp6_2[nx][ny],wRtemp6_3[nx][ny],Rtemp6_1in2[nx][ny],Rtemp6_2in2[nx][ny],Rtemp6_3in2[nx][ny];
fftw_complex Rtemp61_1[nx][ny],Rtemp61_2[nx][ny],Rtemp61_3[nx][ny],wRtemp61_1[nx][ny],wRtemp61_2[nx][ny],wRtemp61_3[nx][ny];
fftw_plan Rtemp6_1idft,Rtemp6_2idft,Rtemp6_3idft,Rtemp61_1dft,Rtemp61_2dft,Rtemp61_3dft;
double complex CRtemp7[nx][ny][3],CRtemp8[nx][ny][3]; 

fftw_complex wLtemp6_1[nx][ny],wLtemp6_2[nx][ny],wLtemp6_3[nx][ny],Ltemp6_1in2[nx][ny],Ltemp6_2in2[nx][ny],Ltemp6_3in2[nx][ny];
fftw_complex Ltemp61_1[nx][ny],Ltemp61_2[nx][ny],Ltemp61_3[nx][ny],wLtemp61_1[nx][ny],wLtemp61_2[nx][ny],wLtemp61_3[nx][ny];
fftw_plan Ltemp6_1idft,Ltemp6_2idft,Ltemp6_3idft,Ltemp61_1dft,Ltemp61_2dft,Ltemp61_3dft;
double complex CLtemp7[nx][ny][3],CLtemp8[nx][ny][3];
double complex CLR[nx][ny][3];

double complex Ct1corr[nx][ny][3][3],Cwu1zeta_corr[nx][ny][3];
double complex Ct2[nx][ny][3][3],Cwssz_corr[nx][ny][3];
double complex Cwu_2zeta[nx][ny][3];

fftw_complex wss1z_corr[nx][ny],wss2z_corr[nx][ny],wss3z_corr[nx][ny];
double wu1zeta_corr[nx][ny][3][2],wu_2zeta[nx][ny][3][2];

//fftw_complex ss1z_corrIn2[nx][ny],ss2z_corrIn2[nx][ny],ss3z_corrIn2[nx][ny];
//fftw_plan ss1z_corridft,ss2z_corridft,ss3z_corridft;

//fftw_plan ss1zeta_corrdft,ss2zeta_corrdft,ss3zeta_corrdft;
//fftw_complex ss1zeta_corr[nx][ny],ss2zeta_corr[nx][ny],ss3zeta_corr[nx][ny],wss1zeta_corr[nx][ny],wss2zeta_corr[nx][ny],wss3zeta_corr[nx][ny];

//double complex Cwu_2zeta[nx][ny][3];

//fftw_plan wu11zeta_corridft,wu12zeta_corridft,wu13zeta_corridft;
//fftw_complex  wu11zeta_corr[nx][ny],wu12zeta_corr[nx][ny],wu13zeta_corr[nx][ny];
//fftw_complex  u11zeta_corrIn2[nx][ny],u12zeta_corrIn2[nx][ny],u13zeta_corrIn2[nx][ny];
//fftw_complex u1_1corrzeta[nx][ny],u2_1corrzeta[nx][ny],u3_1corrzeta[nx][ny];

//double complex CAA[nx][ny][6],CVaA[nx][ny][3],CVddA[nx][ny][3];
//double complex CBBtemp1[nx][ny][3],CBB[nx][ny][3],CVaB[nx][ny][3],CVddB[nx][ny][3];
//double complex CVaA_M_aB[nx][ny][3],CVdA_M_dB[nx][ny][3];
//double complex CDis2t1[nx][ny][3],CStr2[nx][ny][3];
//========
//double wVaA_M_aB[nx][ny][3][2],wVdA_M_dB[nx][ny][3][2],VaA_M_aB[nx][ny][3][2],VdA_M_dB[nx][ny][3][2];
//=======

//fftw_plan VaA_M_aBidft,VdA_M_dBidft;
//CwRtemp71[nx][ny],CwRtemp72[nx][ny],CwRtemp73[nx][ny],CRtemp71[nx][ny],CRtemp72[nx][ny],CRtemp73[nx][ny];

//ftw_complex wRtemp71[nx][ny],wRtemp72[nx][ny],wRtemp73[nx][ny],Rtemp71[nx][ny],Rtemp72[nx][ny],Rtemp73[nx][ny],Rtemp81[nx][ny],Rtemp82[nx][ny];
//fftw_plan Rtemp71idft,Rtemp72idft,Rtemp73idft;
//fftw_complex wzx[nx][ny],zxin2[nx][ny],wzy[nx][ny],zyin2[nx][ny];
//double zx[nx][ny][2],zy[nx][ny][2];
//fftw_plan zxidft,zyidft;



//"""""""""""""""""""""""""""""""
double mu[nx][ny],cals_mu[nx][ny];
fftw_complex zeta_corr[nx][ny];
