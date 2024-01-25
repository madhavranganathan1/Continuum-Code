/* This is a 4th order pseudo-spectral Runge-Kutta time-step integrator */
int ttime;
#include "partial_deri.c"
void runge_kutta(void){
int i,j,k;

int count=-1;
FILE *output[55],*foutput[55];
FILE *output_mu[25],*output_surf[25],*output_elas_1[25],*output_elas_2[25],*Ginput2;
char file_name[256],file_name_mu[256],file_name_surf[256],file_name_elas_1[256],file_name_elas_2[256],file_name_f[256];

double Gsh[nx][ny],Gavgh,Gavgh2,Gavgsh,Grou[tot_step],Gsum,Gsum2;


for(i=0;i<nx;i++){ for(j=0;j<ny;j++){
L[i][j][0]=-(pow(wx2[i],4)+(2*pow(wx2[i],2)*pow(wy2[j],2))+pow(wy2[j],4));        L[i][j][1]=0.00;
                                           E[i][j]=exp(dt*L[i][j][0]);        E2[i][j]=exp(dt*L[i][j][0]/2.00);
}}
for(i=0;i<M;i++){root[i][0]=cos((M_PI*((double)(i+1)-.5))/(double)M);         root[i][1]=sin((M_PI*((double)(i+1)-.5))/(double)M);}

for(i=0;i<nx;i++){for(j=0;j<ny;j++){for(k=0;k<M;k++){
                                           LR_1[i][j][k][0]=dt*L[i][j][0];       LR_1[i][j][k][1]=L[i][j][1];
                                           LR_2[i][j][k][0]=root[k][0];       LR_2[i][j][k][1]=root[k][1];
                      LR[i][j][k][0]=LR_1[i][j][k][0]+LR_2[i][j][k][0];       LR[i][j][k][1]=LR_1[i][j][k][1]+LR_2[i][j][k][1];
                                                                       }}}      

for(i=0;i<nx;i++){ for(j=0;j<ny;j++){ for(k=0;k<M;k++){

          Q_1[i][j][k]=(LR[i][j][k][0]*(exp(LR[i][j][k][0]/2.0)*cos(LR[i][j][k][1]/2.0)-1)+LR[i][j][k][1]*(exp(LR[i][j][k][0]/2.0)*sin(LR[i][j][k][1]/2.0)))/(pow(LR[i][j][k][0],2)+pow(LR[i][j][k][1],2));

          f1_1[i][j][k]=(((-4-LR[i][j][k][0]+exp(LR[i][j][k][0])*cos(LR[i][j][k][1])*(4-3*LR[i][j][k][0]+pow(LR[i][j][k][0],2)-pow(LR[i][j][k][1],2))+(exp(LR[i][j][k][0])*sin(LR[i][j][k][1])*(3*LR[i][j][k][1]-2*LR[i][j][k][0]*LR[i][j][k][1])))*(pow(LR[i][j][k][0],3)-3*LR[i][j][k][0]*pow(LR[i][j][k][1],2)))+((-LR[i][j][k][1]+exp(LR[i][j][k][0])*sin(LR[i][j][k][1])*(4-3*LR[i][j][k][0]+pow(LR[i][j][k][0],2)-pow(LR[i][j][k][1],2))+exp(LR[i][j][k][0])*cos(LR[i][j][k][1])*(-3*LR[i][j][k][1]+2*LR[i][j][k][0]*LR[i][j][k][1]))*(-pow(LR[i][j][k][1],3)+3*pow(LR[i][j][k][0],2)*LR[i][j][k][1])))/(pow(LR[i][j][k][0],6)+3*pow(LR[i][j][k][0],2)*pow(LR[i][j][k][1],4)+3*pow(LR[i][j][k][1],2)*pow(LR[i][j][k][0],4)+pow(LR[i][j][k][1],6));

            f2_1[i][j][k]=(((2+LR[i][j][k][0]+exp(LR[i][j][k][0])*cos(LR[i][j][k][1])*(LR[i][j][k][0]-2)-LR[i][j][k][1]*exp(LR[i][j][k][0])*sin(LR[i][j][k][1]))*(pow(LR[i][j][k][0],3)-3*LR[i][j][k][0]*pow(LR[i][j][k][1],2)))+(LR[i][j][k][1]+(LR[i][j][k][0]-2)*exp(LR[i][j][k][0])*sin(LR[i][j][k][1])+LR[i][j][k][1]*exp(LR[i][j][k][0])*cos(LR[i][j][k][1]))*(-pow(LR[i][j][k][1],3)+3*pow(LR[i][j][k][0],2)*LR[i][j][k][1]))/(pow(LR[i][j][k][0],6)+3*pow(LR[i][j][k][0],2)*pow(LR[i][j][k][1],4)+3*pow(LR[i][j][k][1],2)*pow(LR[i][j][k][0],4)+pow(LR[i][j][k][1],6));


            f3_1[i][j][k]=(((-4-3*LR[i][j][k][0]-pow(LR[i][j][k][0],2)+pow(LR[i][j][k][1],2)+exp(LR[i][j][k][0])*cos(LR[i][j][k][1])*(4-LR[i][j][k][0])+LR[i][j][k][1]*exp(LR[i][j][k][0])*sin(LR[i][j][k][1]))*(pow(LR[i][j][k][0],3)-3*LR[i][j][k][0]*pow(LR[i][j][k][1],2)))+((-3*LR[i][j][k][1]-2*LR[i][j][k][0]*LR[i][j][k][1]-LR[i][j][k][1]*exp(LR[i][j][k][0])*cos(LR[i][j][k][1])+(4-LR[i][j][k][0])*exp(LR[i][j][k][0])*sin(LR[i][j][k][1]))*(3*LR[i][j][k][1]*pow(LR[i][j][k][0],2)-pow(LR[i][j][k][1],3))))/(pow(LR[i][j][k][0],6)+3*pow(LR[i][j][k][0],2)*pow(LR[i][j][k][1],4)+3*pow(LR[i][j][k][1],2)*pow(LR[i][j][k][0],4)+pow(LR[i][j][k][1],6));

                    }}}

for(i=0;i<nx;i++){for(j=0;j<ny;j++){ Q[i][j]=0.0; f1[i][j]=0.0;  f2[i][j]=0.0; f3[i][j]=0.0; }}

for(i=0;i<nx;i++){for(j=0;j<ny;j++){for(k=0;k<M;k++){
        Q[i][j]=Q[i][j]+Q_1[i][j][k];

        f1[i][j]=f1[i][j]+f1_1[i][j][k];

        f2[i][j]=f2[i][j]+f2_1[i][j][k];

        f3[i][j]=f3[i][j]+f3_1[i][j][k];
        }}}
for(i=0;i<nx;i++){for(j=0;j<ny;j++){Q[i][j]=dt*(Q[i][j]/(double)M);
                                    f1[i][j]=dt*(f1[i][j]/(double)M);
                                    f2[i][j]=dt*(f2[i][j]/(double)M);
                                    f3[i][j]=dt*(f3[i][j]/(double)M);}}
for(i=0;i<nx;i++){for(j=0;j<ny;j++){ wg[i][j]=-(pow(wx2[i],2)+pow(wy2[j],2));   }}

for(ttime=0;ttime<=tot_step;ttime++){

       	/////////////////////////////UNIFORM NOISE/////////////////////////////////////////////////////////////////
 /*
   double rn,noise[nx][ny];
   srand(time(NULL));
   for(i=0;i<nx;i++){
     for(j=0;j<ny;j++){
        rn=((double)rand()/(double)RAND_MAX)*0.008 + hight_t[i][j][0];
        //noise[i][j]= rn;
        //Nhight_t[i][j][0]=noise[i][j] + hight_t[i][j][0];
          Nhight_t[i][j][0]=rn;
          Nhight_t[i][j][1]=hight_t[i][j][1];

   }}
*/
///////////////////////Gaussian Noise//////////////

  static double V1, V2, S;
  static int phase = 0;
  double X, Y[nx][ny],U1,U2,noise[nx][ny];


  for (i=0; i<nx; i++){for(j=0;j<ny;j++){
        if(phase == 0) {
       // for (i=0, i<nx, i++){
                do {
                        U1 = (rand() /(double) RAND_MAX);
                        U2 = (rand() /(double) RAND_MAX);

                        V1 = 2 * U1 - 1;
                        V2 = 2 * U2 - 1;
                        S = V1 * V1 + V2 * V2;
                       }
               while(S >= 1.0 || S == 0.0);
                X = V1 * sqrt(-2 * log(S) / S);
        } else
                X = V2 * sqrt(-2 * log(S) / S);
              phase = 1 - phase;
         noise[i][j]=X*N_G;
         Nhight_t[i][j][0]= noise[i][j] + flux + hight_t[i][j][0];
         Nhight_t[i][j][1]=hight_t[i][j][1];
  }}


//////////////////////////////////////////////////////////////////////////////////////////////
  for(i=0;i<nx;i++){ for(j=0;j<ny;j++){  in[i][j][0] = Nhight_t[i][j][0] ;        in[i][j][1] = Nhight_t[i][j][1]; }}

  fftw_execute(dft);

  for(i=0;i<nx;i++){ for(j=0;j<ny;j++){  wh[i][j][0] = out[i][j][0];            wh[i][j][1] = out[i][j][1]; }}
  
  for(i=0;i<nx;i++){  for(j=0;j<ny;j++){
         wh_1[i][j][0]=wh[i][j][0];                                       wh_1[i][j][1]=wh[i][j][1];}}
//           printf("wh[0][0][0]=%lf\twh[0][0][1]=%lf\n\n",i,j,wh[0][0][0],i,j,wh[0][0][1]);
//          printf("out[0][0][0]=%lf\tout[0][0][1]=%lf\n\n",i,j,out[0][0][0],i,j,out[0][0][1]);
  partial_deri();

  for(i=0;i<nx;i++){  for(j=0;j<ny;j++){
//          printf("wcals[%d][%d][0]=%lf\n",i,j,wcals[i][j][0]);
         Nv[i][j][0]=wg[i][j]*wcals[i][j][0];                             Nv[i][j][1]=wg[i][j]*wcals[i][j][1];
//          printf("Nv[%d][%d][0]=%lf\tNv[%d][%d][1]=%lf@@@@@\n",i,j,Nv[i][j][0],Nv[i][j][1]);

          k1[i][j][0]=E2[i][j]*wh_1[i][j][0]+Q[i][j]*Nv[i][j][0];          k1[i][j][1]=E2[i][j]*wh_1[i][j][1]+Q[i][j]*Nv[i][j][1];
                                         wh[i][j][0]=k1[i][j][0];          wh[i][j][1]=k1[i][j][1];
//          printf("wh[%d][%d][0]=%lf\twh[%d][%d][1]=%lf@@@@@\n",i,j,wh[i][j][0],i,j,wh[i][j][1]);
   }}
   partial_deri();

   for(i=0;i<nx;i++){  for(j=0;j<ny;j++){
//          printf("wcals[%d][%d][0]=%lf\n",i,j,wcals[i][j][0]);
          Nk1[i][j][0]=wg[i][j]*wcals[i][j][0];                            Nk1[i][j][1]=wg[i][j]*wcals[i][j][1];

          k2[i][j][0]=E2[i][j]*wh_1[i][j][0]+Q[i][j]*Nv[i][j][0];          k2[i][j][1]=E2[i][j]*wh_1[i][j][1]+Q[i][j]*Nv[i][j][1];
          wh[i][j][0]=k2[i][j][0];                                         wh[i][j][1]=k2[i][j][1];
//          printf("wh[%d][%d][0]=%lf\twh[%d][%d][1]=%lf@@@@@\n",i,j,wh[i][j][0],wh[i][j][1]);
    }}

    partial_deri();

   for(i=0;i<nx;i++){  for(j=0;j<ny;j++){      
//          printf("wcals[%d][%d][0]=%lf@@@@@@\n",i,j,wcals[i][j][0]);
          Nk2[i][j][0]=wg[i][j]*wcals[i][j][0];                            Nk2[i][j][1]=wg[i][j]*wcals[i][j][1];

          k3[i][j][0]=E2[i][j]*k1[i][j][0]+Q[i][j]*(2*Nk2[i][j][0]-Nv[i][j][0]);
          k3[i][j][1]=E2[i][j]*k1[i][j][1]+Q[i][j]*(2*Nk2[i][j][1]-Nv[i][j][1]);
          wh[i][j][0]=k3[i][j][0];                                         wh[i][j][1]=k3[i][j][1];
//          printf("wh[%d][%d][0]=%lf\twh[%d][%d][1]=%lf@@@@@\n",i,j,wh[i][j][0],i,j,wh[i][j][1]);
     }}

     partial_deri(); /* This is the subroutine which calls the new surface function which calls the elastic energy */

     for(i=0;i<nx;i++){  for(j=0;j<ny;j++){
//          printf("%d\t%d\t%lf\n",i,j,wcals[i][j][0]);  
//          printf("wcals[%d][%d][0]=%lf\n@@@@@",i,j,wcals[i][j][0]);
          Nk3[i][j][0]=wg[i][j]*wcals[i][j][0];                            Nk3[i][j][1]=wg[i][j]*wcals[i][j][1];

          wh[i][j][0]=(E[i][j]*wh_1[i][j][0])+(Nv[i][j][0]*f1[i][j])+(2*(Nk1[i][j][0]+Nk2[i][j][0])*f2[i][j])+(Nk3[i][j][0]*f3[i][j]);
//          printf("%d\t%d\t%lf\n",i,j,wh[i][j][0]);
          wh[i][j][1]=(E[i][j]*wh_1[i][j][1])+(Nv[i][j][1]*f1[i][j])+(2*(Nk1[i][j][1]+Nk2[i][j][1])*f2[i][j])+(Nk3[i][j][1]*f3[i][j]);
          
          //printf("\n***********************************************\n");
          //printf("\n***********************************************\n");
          //          printf("wh[%d][%d][0]=%lf\twh[%d][%d][1]=%lf@@@@@@\n",i,j,wh[i][j][0],i,j,wh[i][j][1]);
         
          //printf("\n");
     }}         
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     for(i=0;i<nx;i++){ for(j=0;j<ny;j++){
              out[i][j][0] = wh[i][j][0];     out[i][j][1] = wh[i][j][1];
//             printf("out[%d][%d][0]=%lf\tout[%d][%d][1]=%lf\n",i,j,out[i][j][0],out[i][j][1]);
                                                                                                  }}
     fftw_execute(idft);

     for(i=0;i<nx;i++){ for(j=0;j<ny;j++){
         hight_t[i][j][0]=in2[i][j][0]/(nx*ny);       // hight_t[i][j][1]=in2[i][j][1]/(nx*ny);
         hight_t[i][j][1]=0;    
         hight[i][j]=hight_t[i][j][0];                                                           }}


   for(i=0;i<nx;i++){for(j=0;j<ny;j++){out[i][j][0] = wh[i][j][0];     out[i][j][1] = wh[i][j][1];}}
   fftw_execute(idft);
          
   if((ttime%Nprint)==0){
     // for(i=0;i<nx;i++){for(j=0;j<ny;j++){out[i][j][0] = wh[i][j][0];     out[i][j][1] = wh[i][j][1];}}
     //     fftw_execute(idft);

       printf("ttime=%d\n",ttime);
       count=count+1;
       sprintf(file_name,"data_%d.dat",count);
//   sprintf(file_name_mu,"data_%d_mu.dat",count);
//    sprintf(file_name_f,"fdata_%d.dat",count); 
//   sprintf(file_name_surf,"data_%d_surf.dat",count);
//   sprintf(file_name_elas_1,"data_%d_elas_1.dat",count);
//   sprintf(file_name_elas_2,"data_%d_elas_2.dat",count);
   
        output[count]=fopen(file_name,"w");
//   output_mu[count]=fopen(file_name_mu,"w");
//   foutput[count]=fopen(file_name_f,"w");
//   output_surf[count]=fopen(file_name_surf,"w");
//   output_elas_1[count]=fopen(file_name_elas_1,"w");
//   output_elas_2[count]=fopen(file_name_elas_2,"w");

        for(i=0;i<nx;i++){
             for(j=0;j<ny;j++){
                 fprintf(output[count],"%d\t%d\t%lf\n",i,j,in2[i][j][0]/(nx*ny));
             }
             fprintf(output[count],"\n");
        }
        fclose(output[count]);
  
   } /*End of if loop for printing output files */
 
   for(i=0; i<nx; i++){for(j=0; j<ny; j++){ Gsh[i][j]= pow((in2[i][j][0]/(nx*ny)), 2);}}

        Gsum = 0.0;
        for (i=0; i<nx; i++){for(j=0; j<ny; j++){Gsum= Gsum + in2[i][j][0]/(nx*ny);}}

        Gsum2=0.0;
        for (i=0; i<nx; i++){for(j=0; j<ny; j++){Gsum2= Gsum2 + Gsh[i][j];}}

        Gavgh = Gsum/(nx*ny);
        Gavgh2 = pow(Gavgh,2);
        Gavgsh = Gsum2/(nx*ny);

        Grou[ttime]= sqrt(Gavgsh-Gavgh2);

} /* End of loop over ttime */

Ginput2 = fopen("roug.dat","w");

for(i=0; i<=tot_step; i++){
           fprintf(Ginput2,"%d\t%lf\n", i, Grou[i]);
                }
fclose (Ginput2);

}
