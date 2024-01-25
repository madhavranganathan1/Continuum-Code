fftw_complex cals[nx][ny];
void evolution_eqn_elas(void){
//double cals1[nx][ny],cals2[nx][ny],cals3[nx][ny];
int i,j;
//flux=0;
        for(i=0;i<nx;i++){for(j=0;j<ny;j++){

//     cals[i][j][0]=-(yhn[i][j]*(hxx[i][j]+hyy[i][j]))+((1-.5*(pow(hx[i][j],2)+pow(hy[i][j],2)))*yhn_h[i][j]);
        cals[i][j][0]=-yhn[i][j]*ff[i][j]+((1-.5*ss00_[i][j])*yhn_h[i][j])+muyn[i][j];
        cals[i][j][1]=0.00;
//=======================================================================
//following variable can be used to compare surface energies defined in files "new_surface.c" and "new_surface_energy.c"
//        cals_temp[i][j][0]=-(yhn[i][j]*(hxx[i][j]+hyy[i][j]))+((1-.5*(pow(hx[i][j],2)+pow(hy[i][j],2)))*yhn_h[i][j])-(hxx[i][j]*yhn_hxhx[i][j])-(2*hxy[i][j]*yhn_hxhy[i][j])-(hyy[i][j]*yhn_hyhy[i][j])-(2*(hx[i][j]*hxx[i][j]+hy[i][j]*hxy[i][j])*yhn_hx[i][j])-(2*(hx[i][j]*hxy[i][j]+hy[i][j]*hyy[i][j])*yhn_hy[i][j]);
//        cals_temp[i][j][1]=0.00;
//       printf("%d  %d\t\t%e\t\t%e\t%e\n",i,j,cals[i][j][0],cals_temp[i][j][0],cals[i][j][0]-cals_temp[i][j][0]); 
//=======================================================================


//        cals[i][j][0]=0.00;
//""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
//This segment is only for chemical potential calculation:
//        cals_mu[i][j]=-((1.00+yhn[i][j])*(hxx[i][j]+hyy[i][j]))+((1-.5*(pow(hx[i][j],2)+pow(hy[i][j],2)))*yhn_h[i][j]);
//        mu[i][j]=cals_mu[i][j]+ells[i][j][0]/(nx*ny);

//""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        gin[i][j][0]=cals[i][j][0];                gin[i][j][1]=cals[i][j][1];    }}

        fftw_execute(gdft);
     
       for(i=0;i<nx;i++){for(j=0;j<ny;j++){
       wcals[i][j][0]=gout[i][j][0]+wels[i][j][0];              wcals[i][j][1]=gout[i][j][1]+wels[i][j][1];   
//       printf("%d  %d\t\t%e\t\t%e\n",i,j,wels[i][j][0],wels[i][j][1]);

//      printf("cals[%d][%d][0]=%lf\tcals[%d][%d][1]=%lf\twcals[%d][%d][0]=%lf\twcals[%d][%d][1]=%lf\n",i,j,cals[i][j][0],i,j,cals[i][j][1],i,j,wcals[i][j][0],i,j,wcals[i][j][1]);                                                         
//       printf("out[%d][%d][0]=%lf\tout[%d][%d][1]=%lf*********\n",i,j,out[i][j][0],i,j,out[i][j][1]);
       }}
      }
