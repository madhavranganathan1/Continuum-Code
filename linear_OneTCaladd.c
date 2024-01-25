void linear_OneTCaladd(void){
// Several quantities used in the correction to the height field are calculated in this routine. 

u_0  = -((c11+2*c12)/c11)*ms;
ss_0 = ((c11-c12)*(c11+2*c12)/c11)*ms;
sn_0 = ms;
sn33_0 = -2*(c12/c11)*ms;
els_0=ss_0*sn_0;

int i,j,k,l,m,c;
for(i=0;i<nx;i++){   for(j=0;j<ny;j++){  for(k=0;k<6;k++){ Cv[i][j][k] = v[i][j][k][0] + v[i][j][k][1]*I; }}}        //[EigenValues]


for(i=0;i<nx;i++){   for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<6;l++){
       Calpha[i][j][k][l] = alpha[i][j][k][l][0] +  alpha[i][j][k][l][1] * I ;                            }}}}      //[EigenVectors]       


for(i=0;i<nx;i++){   for(j=0;j<ny;j++){  for(k=0;k<6;k++){
Cdd[i][j][0][k] = c44*(Cv[i][j][k]*Calpha[i][j][0][k]+I*wx2[i]*Calpha[i][j][2][k]);
Cdd[i][j][1][k] = c44*(Cv[i][j][k]*Calpha[i][j][1][k]+I*wy2[j]*Calpha[i][j][2][k]);
Cdd[i][j][2][k] = c12*(I*wx2[i]*Calpha[i][j][0][k]+I*wy2[j]*Calpha[i][j][1][k])+c11*Cv[i][j][k]*Calpha[i][j][2][k];
}}}                                                                                                                  //[ddValues]

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){     CDet1 = 0 + 0 * I;        CDet2 = 0 + 0 * I;                                           
           for(c=0;c<3;c++){
                     CDet1  = CDet1 + Calpha[i][j][0][c]*(Calpha[i][j][1][(c+1)%3]*Calpha[i][j][2][(c+2)%3]-Calpha[i][j][1][(c+2)%3]*Calpha[i][j][2][(c+1)%3]);
                     CDet2  = CDet2 + Cdd[i][j][0][c]*(Cdd[i][j][1][(c+1)%3]*Cdd[i][j][2][(c+2)%3]-Cdd[i][j][1][(c+2)%3]*Cdd[i][j][2][(c+1)%3]);
                      }
               CDet_alpha[i][j]=CDet1;       CDetdd[i][j] = CDet2;                                                      
                  }}                                                                                              // Determinant of [EigenVectors] and [ddValues]  


for(i=0;i<nx;i++){  for(j=0;j<ny;j++){for(k=0;k<3;k++){  for(l=0;l<3;l++){
CIalpha[i][j][l][k]=((Calpha[i][j][(k+1)%3][(l+1)%3] * Calpha[i][j][(k+2)%3][(l+2)%3]) - (Calpha[i][j][(k+1)%3][(l+2)%3]*Calpha[i][j][(k+2)%3][(l+1)%3]))/CDet_alpha[i][j];
CIdd[i][j][l][k]   =((   Cdd[i][j][(k+1)%3][(l+1)%3] *    Cdd[i][j][(k+2)%3][(l+2)%3]) -    (Cdd[i][j][(k+1)%3][(l+2)%3]*   Cdd[i][j][(k+2)%3][(l+1)%3]))/CDetdd[i][j];
}}}}                                                                                                              // Inverse of [EigenVectors] and [ddValues]

for(i=0;i<3;i++){  for(j=0;j<3;j++){
CIalpha[0][0][i][j] = 0 + I * 0;   CIdd[0][0][i][j] = 0 + I * 0;
}}

//=========================
for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<6;l++){
                              if(fabs(creal(CIalpha[i][j][k][l]))<1.00*pow(10,-12)){ CIalpha[i][j][k][l]= 0 + cimag(CIalpha[i][j][k][l]) * I; }
                              if(fabs(cimag(CIalpha[i][j][k][l]))<1.00*pow(10,-12)){ CIalpha[i][j][k][l]= creal(CIalpha[i][j][k][l])+ 0 * I; }                  }}}}

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<6;l++){
                              if(fabs(creal(CIdd[i][j][k][l]))<1.00*pow(10,-12)){ CIdd[i][j][k][l]= 0 + cimag(CIdd[i][j][k][l]) * I; }
                              if(fabs(cimag(CIdd[i][j][k][l]))<1.00*pow(10,-12)){ CIdd[i][j][k][l]= creal(CIdd[i][j][k][l])+ 0 * I; }                  }}}}



for(i=0;i<nx;i++){   for(j=0;j<ny;j++){  for(k=0;k<3;k++){ for(l=0;l<6;l++){ CValpha[i][j][k][l] = Cv[i][j][l] * Calpha[i][j][k][l];}}}}   //[EigenValues*EigenVector]
for(i=0;i<nx;i++){   for(j=0;j<ny;j++){  for(k=0;k<3;k++){ for(l=0;l<6;l++){    CVdd[i][j][k][l] = Cv[i][j][l] * Cdd[i][j][k][l];}}}}      //[EigenVatues*ddValues]


//==========================================================
for(i=0;i<nx;i++){   for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){ CMatSum0 = 0 + 0 * I;
                     for(c=0;c<3;c++){CMatSum0=CMatSum0+Calpha[i][j][k][c]*CIdd[i][j][c][l];}          
                     CpraaIdd[i][j][k][l]=CMatSum0;                                                      }}}}//[EigenVectors].[Inverse ddValues]


for(i=0;i<nx;i++){   for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){ CMatSum0 = 0 + 0 * I;
                     for(c=0;c<3;c++){CMatSum0=CMatSum0+CValpha[i][j][k][c]*CIdd[i][j][c][l];}
                     CvpraaIdd[i][j][k][l]=CMatSum0;                                                      }}}}//[EigenValues * EigenVectors].[Inverse ddValues]


for(i=0;i<nx;i++){   for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){ CMatSum0 = 0 + 0 * I;
                     for(c=0;c<3;c++){CMatSum0=CMatSum0+CVdd[i][j][k][c]*CIdd[i][j][c][l];}
                     CvprddIdd[i][j][k][l]=CMatSum0;                                                      }}}}//[(EigenValues * ddValues)].[Inverse ddValues]

for(i=0;i<nx;i++){   for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){
  fg[i][j][k][l][0]   = creal(CpraaIdd[i][j][k][l]);           fg[i][j][k][l][1]= cimag(CpraaIdd[i][j][k][l]);
  kfg[i][j][k][l][0]  = creal(CvpraaIdd[i][j][k][l]);         kfg[i][j][k][l][1]= cimag(CvpraaIdd[i][j][k][l]);
  lkfg[i][j][k][l][0] = creal(CvprddIdd[i][j][k][l]);        lkfg[i][j][k][l][1]= cimag(CvprddIdd[i][j][k][l]);
}}}}

//==============================================================
for(k=0;k<nx;k++){for(i=0;i<3;i++){ for(j=0;j<3;j++){
fg[nx/2][k][i][j][0]=fg[0][k][i][j][0];     fg[nx/2][k][i][j][1]=fg[0][k][i][j][1];
fg[k][nx/2][i][j][0]=fg[k][0][i][j][0];     fg[k][nx/2][i][j][1]=fg[k][0][i][j][1];

kfg[nx/2][k][i][j][0]=kfg[0][k][i][j][0];   kfg[nx/2][k][i][j][1]=kfg[0][k][i][j][1];
kfg[k][nx/2][i][j][0]=kfg[k][0][i][j][0];   kfg[k][nx/2][i][j][1]=kfg[k][0][i][j][1];

lkfg[nx/2][k][i][j][0]=lkfg[0][k][i][j][0];   lkfg[nx/2][k][i][j][1]=lkfg[0][k][i][j][1];
lkfg[k][nx/2][i][j][0]=kfg[k][0][i][j][0];   lkfg[k][nx/2][i][j][1]=lkfg[k][0][i][j][1];

}}}

for(i=0;i<nx;i++){for(j=0;j<ny;j++){for(k=0;k<3;k++){for(l=0;l<3;l++){ for(m=0;m<2;m++){if(fabs(fg[i][j][k][l][m])<1.00*pow(10,-12)){ fg[i][j][k][l][m]=0;}}}}}}

for(i=0;i<nx;i++){for(j=0;j<ny;j++){for(k=0;k<3;k++){for(l=0;l<3;l++){ for(m=0;m<2;m++){if(fabs(kfg[i][j][k][l][m])<1.00*pow(10,-12)){kfg[i][j][k][l][m]=0;}}}}}}

for(i=0;i<nx;i++){for(j=0;j<ny;j++){for(k=0;k<3;k++){for(l=0;l<3;l++){ for(m=0;m<2;m++){if(fabs(lkfg[i][j][k][l][m])<1.00*pow(10,-12)){lkfg[i][j][k][l][m]=0;}}}}}}


//==================================================================================================
//This section is used for calculating inverse.

for(i=0;i<nx;i++){   for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){
              CMatSum0 = 0 + 0 * I;     CMatSum1 = 0 + 0 * I;    CMatSum2 = 0 + 0 * I;     CMatSum3 = 0 + 0 * I;
//                 printf("\n");
              for(c=0;c<3;c++){
//                 printf("%d  %d %d %d %d \t %e  %e \n",i,j,k,l,c,CMatSum0);
                 CMatSum0=CMatSum0+CIalpha[i][j][k][c]*Calpha[i][j][c][l];                                 //[Inverse EigenVector](3,3).[EigenVector](3,6)
//                 printf("%d  %d %d %d %d \t %e  %e \t %e %e\t %e %e \n",i,j,k,l,c,CMatSum0,CIalpha[i][j][k][c],Calpha[i][j][c][l]);
                 CMatSum1=CMatSum1+CIdd[i][j][k][c]*Cdd[i][j][c][l];                                       //[Inverse ddValues](3,3).[ddValues](3,6)
                 CMatSum2=CMatSum2+CIalpha[i][j][k][c]*Calpha[i][j][c][l+3];
                 CMatSum3=CMatSum3+CIdd[i][j][k][c]*Cdd[i][j][c][l+3];
               }
         Cpraa[i][j][k][l] = CMatSum0;         Cpraa[i][j][k][l+3]=CMatSum2;
         Cprdd[i][j][k][l] = CMatSum1;         Cprdd[i][j][k][l+3]=CMatSum3;                                 }}}}

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<6;l++){
                              if(fabs(creal(Cpraa[i][j][k][l]))<1.00*pow(10,-12)){ Cpraa[i][j][k][l]= 0 + cimag(Cpraa[i][j][k][l]) * I; }
                              if(fabs(cimag(Cpraa[i][j][k][l]))<1.00*pow(10,-12)){ Cpraa[i][j][k][l]= creal(Cpraa[i][j][k][l])+ 0 * I; }                  }}}}


for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<6;l++){
                              if(fabs(creal(Cprdd[i][j][k][l]))<1.00*pow(10,-12)){ Cprdd[i][j][k][l]= 0 + cimag(Cprdd[i][j][k][l]) * I; }
                              if(fabs(cimag(Cprdd[i][j][k][l]))<1.00*pow(10,-12)){ Cprdd[i][j][k][l]= creal(Cprdd[i][j][k][l])+ 0 * I; }                  }}}}


for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){   for(l=0;l<6;l++){    
                                   CpraaMprdd[i][j][k][l]=Cpraa[i][j][k][l]-Cprdd[i][j][k][l];     }}}} 


for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){   for(l=0;l<3;l++){   
                                   CaMd[i][j][k][l]=Cpraa[i][j][k][l+3]-Cprdd[i][j][k][l+3];     }}}} 


for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){   for(l=0;l<2;l++){   for(m=0;m<2;m++){

       CUU[i][j][k][l][m]= Cdd[i][j][(k+1)%3][l]*Cdd[i][j][(k+2)%3][m+1]-Cdd[i][j][(k+2)%3][l]*Cdd[i][j][(k+1)%3][m+1];
       CRR[i][j][k][l][m]= CaMd[i][j][(k+1)%3][l]*CaMd[i][j][(k+2)%3][m+1]-CaMd[i][j][(k+2)%3][l]*CaMd[i][j][(k+1)%3][m+1];
}}}}}

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){   
                           
                             CUU1[i][j][k][0]=CUU[i][j][k][0][0];                   CRR1[i][j][k][0]=CRR[i][j][k][0][0];
                             CUU1[i][j][k][1]=CUU[i][j][k][0][1];                   CRR1[i][j][k][1]=CRR[i][j][k][0][1];
                             CUU1[i][j][k][2]=CUU[i][j][k][1][1];                   CRR1[i][j][k][2]=CRR[i][j][k][1][1];    }}}


for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){   for(l=0;l<3;l++){   CMatSum0= 0+0*I;
           for(c=0;c<3;c++){ CMatSum0 =  CMatSum0+Cdd[i][j][c][k+3]*CUU1[i][j][c][l];  }
                             CddUU1[i][j][k][l]=CMatSum0;                                                            }}}}


for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  
                      
                            CRR2[i][j][k][0]= CRR1[i][j][k][2];                     CddUU2[i][j][k][0]= CddUU1[i][j][k][2];
                            CRR2[i][j][k][1]=-CRR1[i][j][k][1];                     CddUU2[i][j][k][1]=-CddUU1[i][j][k][1];
                            CRR2[i][j][k][2]= CRR1[i][j][k][0];                     CddUU2[i][j][k][2]= CddUU1[i][j][k][0];   }}}


for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){  for(m=0;m<3;m++){
                            CCo1[i][j][k][l][m]=CRR2[i][j][k][m]*CddUU2[i][j][m][l];                                }}}}}


for(i=0;i<nx;i++){  for(j=0;j<ny;j++){     CDet1 = 0 + 0 * I;        CDet2 = 0 + 0 * I;  
           for(c=0;c<3;c++){
                     CDet2  = CDet2 + CaMd[i][j][0][c]*(CaMd[i][j][1][(c+1)%3]*CaMd[i][j][2][(c+2)%3]-CaMd[i][j][1][(c+2)%3]*CaMd[i][j][2][(c+1)%3]); 
                      }
                    CDetaMd[i][j] = CDet2;                                                       }}                                          //Determinant Calculation

//===================================

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){  for(m=0;m<3;m++){  CCo21_temp[i][j][k][l][m]=-CCo1[i][j][l][k][m]/(CDetdd[i][j]*CDetaMd[i][j]); }}}}}
for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){  CCo22_temp[i][j][k][l]=CRR2[i][j][l][k]/CDetaMd[i][j];   }}}}//Inverse Cal

for(k=0;k<3;k++){  for(l=0;l<3;l++){  for(m=0;m<3;m++){ CCo21_temp[0][0][k][l][m]=0+I*0; }}}
for(k=0;k<3;k++){  for(l=0;l<3;l++){ CCo22_temp[0][0][k][l]=0+I*0;   }}

//==================================
for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){  CMatSum0 = 0 + 0 * I;
                    for(m=0;m<3;m++){ CMatSum0 = CMatSum0 + CCo21_temp[i][j][k][l][m]*cexp(-Cv[i][j][m]*h0_avg);}
                    CCo21[i][j][k][l] = CMatSum0;                                                                                }}}}       //Inverse cal

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){ CCo22[i][j][k][l]= CCo22_temp[i][j][k][l]*cexp(-Cv[i][j][k]*h0_avg);}}}} // Inverse cal

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){  CCo[i][j][k][l]=CCo21[i][j][k][l];
                                                                             CCo[i][j][k+3][l]=CCo22[i][j][k][l];                }}}}      // combining two 3*3 inverse matrices to get a 6*3 matrix

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){  CMatSum0 = 0 + 0 * I;
                    for(c=0;c<6;c++){CMatSum0 = CMatSum0+Calpha[i][j][k][c]*CCo[i][j][c][l];}
                    Ct1[i][j][k][l]=CMatSum0;                                                                               }}}}           //multiplication to Eigen Vector Matrix   

for(k=0;k<3;k++){  for(l=0;l<3;l++){ Ct1[0][0][k][l]=0+I*0;   }}

for(i=0;i<nx;i++){   for(j=0;j<ny;j++){
wp1[i][j][0]=-(CIdd[i][j][0][0]*wx1[i]+CIdd[i][j][0][1]*wy1[j])*I*wzeta[i][j]*ss_0- u_0*CIalpha[i][j][0][2]*wzeta[i][j];
wp1[i][j][1]=-(CIdd[i][j][1][0]*wx1[i]+CIdd[i][j][1][1]*wy1[j])*I*wzeta[i][j]*ss_0- u_0*CIalpha[i][j][1][2]*wzeta[i][j];
wp1[i][j][2]=-(CIdd[i][j][2][0]*wx1[i]+CIdd[i][j][2][1]*wy1[j])*I*wzeta[i][j]*ss_0- u_0*CIalpha[i][j][2][2]*wzeta[i][j]; }}


//==================================================
//This section will also be copied to diff subroutine.

for(i=0;i<nx;i++){   for(j=0;j<ny;j++){  for(k=0;k<3;k++){ CMatSum0 = 0 + 0 * I;
                     for(c=0;c<3;c++){CMatSum0=CMatSum0+Ct1[i][j][k][c]*wp1[i][j][c];}
                     Cwu_1zeta[i][j][k]=CMatSum0;                                                      }}}

for(k=0;k<nx;k++){for(i=0;i<3;i++){Cwu_1zeta[nx/2][k][i]=Cwu_1zeta[0][k][i];     Cwu_1zeta[k][nx/2][i]=Cwu_1zeta[k][0][i];     }}

for(i=0;i<nx;i++){for(j=0;j<ny;j++){for(k=0;k<3;k++){ wu_1zeta[i][j][k][0]= creal(Cwu_1zeta[i][j][k]);   wu_1zeta[i][j][k][1]= cimag(Cwu_1zeta[i][j][k]);    }}}


// added to "1st order displacement correction term"

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){  CMatSum0 = 0 + 0 * I;
                    for(c=0;c<6;c++){CMatSum0 = CMatSum0+CValpha[i][j][k][c]*CCo[i][j][c][l];}
                    Ct1corr[i][j][k][l]=CMatSum0;                                                                               }}}}           //multiplication to Eigen Vector Matrix   
for(k=0;k<3;k++){  for(l=0;l<3;l++){ Ct1corr[0][0][k][l]=0+I*0;   }}

for(i=0;i<nx;i++){   for(j=0;j<ny;j++){  for(k=0;k<3;k++){ CMatSum0 = 0 + 0 * I;
                     for(c=0;c<3;c++){CMatSum0=CMatSum0+Ct1corr[i][j][k][c]*wp1[i][j][c];}
                     Cwu1zeta_corr[i][j][k]=CMatSum0;                                                      }}}

for(k=0;k<nx;k++){for(i=0;i<3;i++){Cwu1zeta_corr[nx/2][k][i]=Cwu1zeta_corr[0][k][i];     Cwu1zeta_corr[k][nx/2][i]=Cwu1zeta_corr[k][0][i];     }}

for(i=0;i<nx;i++){for(j=0;j<ny;j++){for(k=0;k<3;k++){ wu1zeta_corr[i][j][k][0]= creal(Cwu1zeta_corr[i][j][k]);   wu1zeta_corr[i][j][k][1]= cimag(Cwu1zeta_corr[i][j][k]);    }}}


for(i=0;i<nx/2;i++){for(j=0;j<ny/2;j++){for(k=0;k<3;k++){
                         if(fabs(wu1zeta_corr[i][j][k][0]-wu1zeta_corr[nx-i][ny-j][k][0])<=1.00*pow(10,-12)){ wu1zeta_corr[nx-i][ny-j][k][0] = wu1zeta_corr[i][j][k][0];}
                         if(fabs(wu1zeta_corr[i][j][k][1]+wu1zeta_corr[nx-i][ny-j][k][1])<=1.00*pow(10,-12)){ wu1zeta_corr[nx-i][ny-j][k][1] =-wu1zeta_corr[i][j][k][1];}
}}}

// resulting terms are used in file elas_energy_corr.c .

//=======================================
//added to "1st order stress correction term":
for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){  CMatSum0 = 0 + 0 * I;
                    for(c=0;c<6;c++){CMatSum0 = CMatSum0+CVdd[i][j][k][c]*CCo[i][j][c][l];}
                    Ct2[i][j][k][l]=CMatSum0;                                                                               }}}}           //multiplication to Eigen Vector Matrix   
for(k=0;k<3;k++){  for(l=0;l<3;l++){ Ct2[0][0][k][l]=0+I*0;   }}

for(i=0;i<nx;i++){   for(j=0;j<ny;j++){  for(k=0;k<3;k++){ CMatSum0 = 0 + 0 * I;
                     for(c=0;c<3;c++){CMatSum0=CMatSum0+Ct2[i][j][k][c]*wp1[i][j][c];}
                     Cwssz_corr[i][j][k]=CMatSum0;                                                      }}}

for(i=0;i<nx;i++){   for(j=0;j<ny;j++){
wss1z_corr[i][j][0]=creal(Cwssz_corr[i][j][0]);    wss1z_corr[i][j][1]=cimag(Cwssz_corr[i][j][0]);
wss2z_corr[i][j][0]=creal(Cwssz_corr[i][j][1]);    wss2z_corr[i][j][1]=cimag(Cwssz_corr[i][j][1]);
wss3z_corr[i][j][0]=creal(Cwssz_corr[i][j][2]);    wss3z_corr[i][j][1]=cimag(Cwssz_corr[i][j][2]);}}

// resulting terms are used in file elas_energy_corr.c .
//===================================================

//                            ++++++++++++++++++Second Order variables:++++++++++++++++++


//double complex CXX21[nx][ny][3][3],CXX22[nx][ny][3][3],CXX[nx][ny][6][3];


for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){  CMatSum0 = 0 + 0 * I;
                    for(m=0;m<3;m++){ CMatSum0 = CMatSum0 + CCo21_temp[i][j][k][l][m]*cexp(-(Cv[i][j][m]+Cv[i][j][k])*h0_avg);}
                    CXX21[i][j][k][l] = CMatSum0;                                                                                }}}}       //Inverse cal

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){ CXX22[i][j][k][l]= CCo22_temp[i][j][k][l];}}}} // Inverse cal

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<3;l++){ CXX[i][j][k][l]=CXX21[i][j][k][l];

                                                                           CXX[i][j][k+3][l]=CXX22[i][j][k][l];                }}}}      // combining two 3*3 inverse matrices to get a 6*3 matrix

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){ for(l=0;l<6;l++){   CMatSum0 = 0 + 0 * I;
                    for(c=0;c<3;c++){CMatSum0 = CMatSum0+CVdd[i][j][k][c]*Cpraa[i][j][c][l];}
                    CRtemp1[i][j][k][l]=CMatSum0;                                                                               }}}}  

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){ for(l=0;l<6;l++){CRtemp2[i][j][k][l]=CVdd[i][j][k][l]-CRtemp1[i][j][k][l];  }}}}


for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<6;l++){
                              if(fabs(creal(CRtemp2[i][j][k][l]))<1.00*pow(10,-12)){ CRtemp2[i][j][k][l]= 0 + cimag(CRtemp2[i][j][k][l]) * I; }
                              if(fabs(cimag(CRtemp2[i][j][k][l]))<1.00*pow(10,-12)){ CRtemp2[i][j][k][l]= creal(CRtemp2[i][j][k][l])+ 0 * I; }                  }}}}

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){ for(l=0;l<3;l++){   CMatSum0 = 0 + 0 * I;
                    for(c=0;c<3;c++){CMatSum0 = CMatSum0+CRtemp2[i][j][k][c+3]*CXX22[i][j][c][l];}
                    CRtemp3[i][j][k][l]=CMatSum0;                                                                               }}}}     

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){    CMatSum0 = 0 + 0 * I;
                    for(c=0;c<3;c++){CMatSum0 = CMatSum0+CRtemp3[i][j][k][c]*wp1[i][j][c];}
                    CRtemp4[i][j][k]=CMatSum0;                                                                               }}}          

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){    CMatSum0 = 0 + 0 * I;
                    for(c=0;c<3;c++){CMatSum0 = CMatSum0+CVdd[i][j][k][c]*CIalpha[i][j][c][2];}
                    CRtemp5[i][j][k]=CMatSum0*u_0*wzeta[i][j];                                                               }}}         

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){ CRtemp6[i][j][k]=CRtemp4[i][j][k]-CRtemp5[i][j][k]; }}}

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){
wRtemp6_1[i][j][0]=creal(CRtemp6[i][j][0]);                          wRtemp6_1[i][j][1]=cimag(CRtemp6[i][j][0]);   
wRtemp6_2[i][j][0]=creal(CRtemp6[i][j][1]);                          wRtemp6_2[i][j][1]=cimag(CRtemp6[i][j][1]);   
wRtemp6_3[i][j][0]=creal(CRtemp6[i][j][2]);                          wRtemp6_3[i][j][1]=cimag(CRtemp6[i][j][2]);   }}

         fftw_execute(Rtemp6_1idft);
         fftw_execute(Rtemp6_2idft);
         fftw_execute(Rtemp6_3idft);

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){
Rtemp61_1[i][j][0]=Rtemp6_1in2[i][j][0]/(nx*ny)*zeta_t[i][j][0];                  Rtemp61_1[i][j][1]= Rtemp6_1in2[i][j][1]/(nx*ny)*zeta_t[i][j][0];
Rtemp61_2[i][j][0]=Rtemp6_2in2[i][j][0]/(nx*ny)*zeta_t[i][j][0];                  Rtemp61_2[i][j][1]= Rtemp6_2in2[i][j][1]/(nx*ny)*zeta_t[i][j][0];
Rtemp61_3[i][j][0]=Rtemp6_3in2[i][j][0]/(nx*ny)*zeta_t[i][j][0];                  Rtemp61_3[i][j][1]= Rtemp6_3in2[i][j][1]/(nx*ny)*zeta_t[i][j][0];}}


         fftw_execute(Rtemp61_1dft);
         fftw_execute(Rtemp61_2dft);
         fftw_execute(Rtemp61_3dft);

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){
CRtemp7[i][j][0]= wRtemp61_1[i][j][0]+I*wRtemp61_1[i][j][1];
CRtemp7[i][j][1]= wRtemp61_2[i][j][0]+I*wRtemp61_2[i][j][1];
CRtemp7[i][j][2]= wRtemp61_3[i][j][0]+I*wRtemp61_3[i][j][1]+ (pow(wx2[i],2)+pow(wy2[j],2))*wzeta[i][j]*ss_0 ;  }}

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){    CMatSum0 = 0 + 0 * I;
                    for(c=0;c<3;c++){CMatSum0 = CMatSum0+CIdd[i][j][k][c]*CRtemp7[i][j][c];}
                    CRtemp8[i][j][k]=CMatSum0;                                                                               }}}       


for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){ for(l=0;l<6;l++){   CMatSum0 = 0 + 0 * I;
                    for(c=0;c<3;c++){CMatSum0 = CMatSum0+CValpha[i][j][k][c]*Cpraa[i][j][c][l];}
                    CLtemp1[i][j][k][l]=CMatSum0;                                                                               }}}}

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){ for(l=0;l<6;l++){CLtemp2[i][j][k][l]=CValpha[i][j][k][l]-CLtemp1[i][j][k][l];  }}}}


for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){  for(l=0;l<6;l++){
                              if(fabs(creal(CLtemp2[i][j][k][l]))<1.00*pow(10,-12)){ CLtemp2[i][j][k][l]= 0 + cimag(CLtemp2[i][j][k][l]) * I; }
                              if(fabs(cimag(CLtemp2[i][j][k][l]))<1.00*pow(10,-12)){ CLtemp2[i][j][k][l]= creal(CLtemp2[i][j][k][l])+ 0 * I; }                  }}}}

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){ for(l=0;l<3;l++){   CMatSum0 = 0 + 0 * I;
                    for(c=0;c<3;c++){CMatSum0 = CMatSum0+CLtemp2[i][j][k][c+3]*CXX22[i][j][c][l];}
                    CLtemp3[i][j][k][l]=CMatSum0;                                                                               }}}}

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){    CMatSum0 = 0 + 0 * I;
                    for(c=0;c<3;c++){CMatSum0 = CMatSum0+CLtemp3[i][j][k][c]*wp1[i][j][c];}
                    CLtemp4[i][j][k]=CMatSum0;                                                                               }}}  

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){    CMatSum0 = 0 + 0 * I;
                    for(c=0;c<3;c++){CMatSum0 = CMatSum0+CValpha[i][j][k][c]*CIalpha[i][j][c][2];}
                    CLtemp5[i][j][k]=CMatSum0*u_0*wzeta[i][j];                                                               }}}   

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){ CLtemp6[i][j][k]=CLtemp4[i][j][k]-CLtemp5[i][j][k]; }}}
for(i=0;i<nx;i++){  for(j=0;j<ny;j++){

wLtemp6_1[i][j][0]=creal(CLtemp6[i][j][0]);                          wLtemp6_1[i][j][1]=cimag(CLtemp6[i][j][0]);
wLtemp6_2[i][j][0]=creal(CLtemp6[i][j][1]);                          wLtemp6_2[i][j][1]=cimag(CLtemp6[i][j][1]);
wLtemp6_3[i][j][0]=creal(CLtemp6[i][j][2]);                          wLtemp6_3[i][j][1]=cimag(CLtemp6[i][j][2]);   }}

         fftw_execute(Ltemp6_1idft);
         fftw_execute(Ltemp6_2idft);
         fftw_execute(Ltemp6_3idft);

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){
Ltemp61_1[i][j][0]=Ltemp6_1in2[i][j][0]/(nx*ny)*zeta_t[i][j][0];                  Ltemp61_1[i][j][1]= Ltemp6_1in2[i][j][1]/(nx*ny)*zeta_t[i][j][0];
Ltemp61_2[i][j][0]=Ltemp6_2in2[i][j][0]/(nx*ny)*zeta_t[i][j][0];                  Ltemp61_2[i][j][1]= Ltemp6_2in2[i][j][1]/(nx*ny)*zeta_t[i][j][0];
Ltemp61_3[i][j][0]=Ltemp6_3in2[i][j][0]/(nx*ny)*zeta_t[i][j][0];                  Ltemp61_3[i][j][1]= Ltemp6_3in2[i][j][1]/(nx*ny)*zeta_t[i][j][0];}}

//for(i=0;i<nx;i++){   for(j=0;j<ny;j++){printf("%2d  %2d  %2d\t\t %e  %e  %e  %e\n",i,j,k,Rtemp61_3[i][j][0],Rtemp61_3[i][j][1],Ltemp61_3[i][j][0],Ltemp61_3[i][j][1]);}
//printf("\n");}
  
         fftw_execute(Ltemp61_1dft);
         fftw_execute(Ltemp61_2dft);
         fftw_execute(Ltemp61_3dft);

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){
CLtemp7[i][j][0]= wLtemp61_1[i][j][0]+I*wLtemp61_1[i][j][1];
CLtemp7[i][j][1]= wLtemp61_2[i][j][0]+I*wLtemp61_2[i][j][1];
CLtemp7[i][j][2]= wLtemp61_3[i][j][0]+I*wLtemp61_3[i][j][1];  }}

for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){    CMatSum0 = 0 + 0 * I;
                    for(c=0;c<3;c++){CMatSum0 = CMatSum0+CIalpha[i][j][k][c]*CLtemp7[i][j][c];}
                    CLtemp8[i][j][k]=CMatSum0;                                                                               }}}              



for(i=0;i<nx;i++){  for(j=0;j<ny;j++){  for(k=0;k<3;k++){CLR[i][j][k]=-CLtemp8[i][j][k]+CRtemp8[i][j][k];  }}}

for(i=0;i<nx;i++){   for(j=0;j<ny;j++){  for(k=0;k<3;k++){ CMatSum0 = 0 + 0 * I;
                     for(c=0;c<3;c++){CMatSum0=CMatSum0+Ct1[i][j][k][c]*CLR[i][j][c];}
                     Cwu_2zeta[i][j][k]=CMatSum0;                                                      }}}

for(k=0;k<nx;k++){for(i=0;i<3;i++){Cwu_2zeta[nx/2][k][i]=Cwu_2zeta[0][k][i];     Cwu_2zeta[k][nx/2][i]=Cwu_2zeta[k][0][i];     }}

for(i=0;i<nx;i++){for(j=0;j<ny;j++){for(k=0;k<3;k++){ wu_2zeta[i][j][k][0]= creal(Cwu_2zeta[i][j][k]);   wu_2zeta[i][j][k][1]= cimag(Cwu_2zeta[i][j][k]);    }}}


for(i=0;i<nx/2;i++){for(j=0;j<ny/2;j++){for(k=0;k<3;k++){
                                                          if(fabs(wu_2zeta[i][j][k][0]-wu_2zeta[nx-i][ny-j][k][0])<=1.00*pow(10,-12)){ wu_2zeta[nx-i][ny-j][k][0] = wu_2zeta[i][j][k][0];}
                                                          if(fabs(wu_2zeta[i][j][k][1]+wu_2zeta[nx-i][ny-j][k][1])<=1.00*pow(10,-12)){ wu_2zeta[nx-i][ny-j][k][1] =-wu_2zeta[i][j][k][1];}
}}}


}
