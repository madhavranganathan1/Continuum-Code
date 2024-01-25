//#include "surface_energy.c"
//#include "new_surface_energy.c"
#include "new_surface.c"
void partial_deri(void){
int i,j;
double hi[nx][ny];

         fftw_execute(idft);

        for(i=0;i<nx;i++){ for(j=0;j<ny;j++){
         hight_t[i][j][0]=in2[i][j][0]/(nx*ny);        hight_t[i][j][1]=in2[i][j][1]/(nx*ny);
         hight[i][j]=hight_t[i][j][0];}}


      for(i=0;i<nx;i++){ for(j=0;j<ny;j++){
         whx[i][j][0]=-wx1[i]*wh[i][j][1];      whx[i][j][1]=wx1[i]*wh[i][j][0];
         why[i][j][0]=-wy1[j]*wh[i][j][1];      why[i][j][1]=wy1[j]*wh[i][j][0];
         whxx[i][j][0]=-pow(wx2[i],2)*wh[i][j][0];      whxx[i][j][1]=-pow(wx2[i],2)*wh[i][j][1];
         whxy[i][j][0]=-wx2[i]*wy2[j]*wh[i][j][0];      whxy[i][j][1]=-wx2[i]*wy2[j]*wh[i][j][1];
         whyy[i][j][0]=-pow(wy2[j],2)*wh[i][j][0];      whyy[i][j][1]=-pow(wy2[j],2)*wh[i][j][1];
         }}

         fftw_execute(hxidft);
         fftw_execute(hyidft);
         fftw_execute(hxxidft);
         fftw_execute(hxyidft);
         fftw_execute(hyyidft);        

        for(i=0;i<nx;i++){ for(j=0;j<ny;j++){
         hx[i][j] = hxin2[i][j][0]/(nx*ny);
         hy[i][j] = hyin2[i][j][0]/(nx*ny);
        hxx[i][j] = hxxin2[i][j][0]/(nx*ny);
        hxy[i][j] = hxyin2[i][j][0]/(nx*ny);
        hyy[i][j] = hyyin2[i][j][0]/(nx*ny);        

         p[i][j]=1+pow(hx[i][j],2)+pow(hy[i][j],2);
         k[i][j]=-(((1+pow(hy[i][j],2))*hxx[i][j])-(2*hx[i][j]*hy[i][j]*hxy[i][j])+((1+pow(hx[i][j],2))*hyy[i][j]))/pow(p[i][j],1.5);                                                                                                          }}


    for(i=0;i<nx;i++){ for(j=0;j<ny;j++){
        WHxx[i][j][0]=((wx2[i]*wx2[i])/W[i][j])*wh[i][j][0];  WHxx[i][j][1]=((wx2[i]*wx2[i])/W[i][j])*wh[i][j][1];
        WHxy[i][j][0]=((wx2[i]*wy2[j])/W[i][j])*wh[i][j][0];  WHxy[i][j][1]=((wx2[i]*wy2[j])/W[i][j])*wh[i][j][1];
        WHyy[i][j][0]=((wy2[j]*wy2[j])/W[i][j])*wh[i][j][0];   WHyy[i][j][1]=((wy2[j]*wy2[j])/W[i][j])*wh[i][j][1];
}}

        fftw_execute(Hxxidft);
        fftw_execute(Hxyidft);
        fftw_execute(Hyyidft);

     for(i=0;i<nx;i++){ for(j=0;j<ny;j++){
        Hxx[i][j][0]=Hxxin2[i][j][0]/(nx*ny);                   Hxx[i][j][1]=Hxxin2[i][j][1]/(nx*ny);
        Hxy[i][j][0]=Hxyin2[i][j][0]/(nx*ny);                   Hxy[i][j][1]=Hxyin2[i][j][1]/(nx*ny);
        Hyy[i][j][0]=Hyyin2[i][j][0]/(nx*ny);                   Hyy[i][j][1]=Hyyin2[i][j][1]/(nx*ny);       
}}
new_surface();
}
