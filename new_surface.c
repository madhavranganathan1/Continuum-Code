double yhn[nx][ny];
#include "elas_energy_corr.c"

void new_surface(void){
int i,j;
double yh[nx][ny],yn[nx][ny];
    for(i=0;i<nx;i++){ for(j=0;j<ny;j++){
//============================================================================
//Wetting Term 
//          yh[i][j] = cw*exp(-hight[i][j]/dw);
          yh[i][j]=cw*exp(-(hight[i][j]-zeta_t[i][j][0])/dw);
//          yh[i][j]=cw*exp(-(hight[i][j]-zeta_t[i][j][0])*(1+sqrt(zeta_corr[i][j][0]))/dw);

       yhn_h[i][j] = -yh[i][j]/dw; //first order derivative wrt h

//       yhn_h[i][j] = -(cw*exp(-hight[i][j]/dw))/dw;
//          yh[i][j]=cw*exp(-(hight[i][j]-zeta_t[i][j][0])/dw);

       ff[i][j] = hxx[i][j]+hyy[i][j];
       ss00_[i][j]  = pow(hx[i][j],2)+pow(hy[i][j],2);

//==================================================================================================================
/*
//Open this in case of isotropic surface energy.

      yn[i][j] = 0;
      yhn[i][j]= yh[i][j]+ yn[i][j];
      muyn[i][j] = 0;
*/
//===================================================================================================================
//Open this in case of ANisotropic surface energy.

//Surface Anisotropy Variables
        ss100_[i][j] = ss00_[i][j]+e100;
      ss105_1_[i][j] = ss00_[i][j]+.04+e105-.4*hy[i][j];
      ss105_2_[i][j] = ss00_[i][j]+.04+e105+.4*hy[i][j];
      ss105_3_[i][j] = ss00_[i][j]+.04+e105-.4*hx[i][j];
      ss105_4_[i][j] = ss00_[i][j]+.04+e105+.4*hx[i][j];

////        ss100_[i][j] = pow(hx[i][j],2)+pow(hy[i][j],2)+e100;
        rr100_[i][j] = exp(-n100*sqrt(ss100_[i][j]));
//        rr100_[i][j] = (1-n100*sqrt(ss100_[i][j]));

////      ss105_1_[i][j] = pow(hx[i][j],2)+pow((hy[i][j]-.2),2)+e105;
      rr105_1_[i][j] = exp(-n105*sqrt(ss105_1_[i][j]));
//      rr105_1_[i][j] = (1-n105*sqrt(ss105_1_[i][j]));  


////      ss105_2_[i][j] = pow(hx[i][j],2)+pow((hy[i][j]+.2),2)+e105;
      rr105_2_[i][j] = exp(-n105*sqrt(ss105_2_[i][j]));
//      rr105_2_[i][j] = (1-n105*sqrt(ss105_2_[i][j]));

////      ss105_3_[i][j] = pow((hx[i][j]-.2),2)+pow(hy[i][j],2)+e105;      
      rr105_3_[i][j] = exp(-n105*sqrt(ss105_3_[i][j]));
//      rr105_3_[i][j] = (1-n105*sqrt(ss105_3_[i][j]));

////      ss105_4_[i][j] = pow((hx[i][j]+.2),2)+pow(hy[i][j],2)+e105;
      rr105_4_[i][j] = exp(-n105*sqrt(ss105_4_[i][j]));
//      rr105_4_[i][j] = (1-n105*sqrt(ss105_4_[i][j]));

////        rr100_[i][j] = exp(-n100*sqrt(pow(hx[i][j],2)+pow(hy[i][j],2)+e100));
////      rr105_1_[i][j] = exp(-n105*sqrt(pow(hx[i][j],2)+pow((hy[i][j]-.2),2)+e105));
////      rr105_2_[i][j] = exp(-n105*sqrt(pow(hx[i][j],2)+pow((hy[i][j]+.2),2)+e105));
////      rr105_3_[i][j] = exp(-n105*sqrt(pow((hx[i][j]-.2),2)+pow(hy[i][j],2)+e105));
////      rr105_4_[i][j] = exp(-n105*sqrt(pow((hx[i][j]+.2),2)+pow(hy[i][j],2)+e105));

// ==================================================
//Surface Anisotropy Defined
      yn[i][j] = -a100*rr100_[i][j]-a105*rr105_1_[i][j]-a105*rr105_2_[i][j]-a105*rr105_3_[i][j]-a105*rr105_4_[i][j]+a100*exp(-n100*sqrt(e100))+4.0*a105*exp(-n105*sqrt(e105));

      yhn[i][j]= yh[i][j]+ yn[i][j];
//===================================================
//Surface Anisotropy Derivative

       jj[i][j] = hxx[i][j]*pow(hx[i][j],2)+2*hxy[i][j]*hx[i][j]*hy[i][j]+hyy[i][j]*pow(hy[i][j],2);
      kk1[i][j] = .4*(hx[i][j]*hxy[i][j]+hy[i][j]*hyy[i][j]);
      kk2[i][j] = .4*(hx[i][j]*hxx[i][j]+hy[i][j]*hxy[i][j]);


      mu100[i][j] = a100*n100*rr100_[i][j]*(jj[i][j]*(1+n100*sqrt(ss100_[i][j]))-(2*jj[i][j]+ff[i][j])*ss100_[i][j])/pow(ss100_[i][j],1.5);

      mu105_1[i][j] = a105*n105*rr105_1_[i][j]*((jj[i][j]-kk1[i][j]+.04*hyy[i][j])*(1+n105*sqrt(ss105_1_[i][j]))-(2*jj[i][j]-kk1[i][j]+ff[i][j])*ss105_1_[i][j])/pow(ss105_1_[i][j],1.5);
      mu105_2[i][j] = a105*n105*rr105_2_[i][j]*((jj[i][j]+kk1[i][j]+.04*hyy[i][j])*(1+n105*sqrt(ss105_2_[i][j]))-(2*jj[i][j]+kk1[i][j]+ff[i][j])*ss105_2_[i][j])/pow(ss105_2_[i][j],1.5);

      mu105_3[i][j] = a105*n105*rr105_3_[i][j]*((jj[i][j]-kk2[i][j]+.04*hxx[i][j])*(1+n105*sqrt(ss105_3_[i][j]))-(2*jj[i][j]-kk2[i][j]+ff[i][j])*ss105_3_[i][j])/pow(ss105_3_[i][j],1.5);
      mu105_4[i][j] = a105*n105*rr105_4_[i][j]*((jj[i][j]+kk2[i][j]+.04*hxx[i][j])*(1+n105*sqrt(ss105_4_[i][j]))-(2*jj[i][j]+kk2[i][j]+ff[i][j])*ss105_4_[i][j])/pow(ss105_4_[i][j],1.5);

      muyn[i][j] = mu100[i][j]+mu105_1[i][j]+mu105_2[i][j]+mu105_3[i][j]+mu105_4[i][j];


//muyn_temp[i][j] = mu100[i][j]+mu105_1[i][j]+mu105_2[i][j]+mu105_3[i][j]+mu105_4[i][j];
//muyn_temp[i][j] = -(hxx[i][j]*yhn_hxhx[i][j])-(2*hxy[i][j]*yhn_hxhy[i][j])-(hyy[i][j]*yhn_hyhy[i][j])-(2*(hx[i][j]*hxx[i][j]+hy[i][j]*hxy[i][j])*yhn_hx[i][j])-(2*(hx[i][j]*hxy[i][j]+hy[i][j]*hyy[i][j])*yhn_hy[i][j]);

//printf("%2d  %2d\t%lf\t%lf\t\t%e\n",i,j,muyn[i][j],muyn_temp[i][j],muyn[i][j]-muyn_temp[i][j]);

        }
       }
elas_energy_corr();
  }

