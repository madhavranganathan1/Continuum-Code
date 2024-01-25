double anum11raa,anum12raa,anum13raa,anum21raa,anum22raa,anum23raa,anum31raa,anum32raa,anum33raa;
double anum11rab,anum12rab,anum13rab,anum21rab,anum22rab,anum23rab,anum31rab,anum32rab,anum33rab;
double anum11rba,anum12rba,anum13rba,anum21rba,anum22rba,anum23rba,anum31rba,anum32rba,anum33rba;
double anum11rbb,anum12rbb,anum13rbb,anum21rbb,anum22rbb,anum23rbb,anum31rbb,anum32rbb,anum33rbb;
double anum11iaa,anum12iaa,anum13iaa,anum21iaa,anum22iaa,anum23iaa,anum31iaa,anum32iaa,anum33iaa;
double anum11iab,anum12iab,anum13iab,anum21iab,anum22iab,anum23iab,anum31iab,anum32iab,anum33iab;
double anum11iba,anum12iba,anum13iba,anum21iba,anum22iba,anum23iba,anum31iba,anum32iba,anum33iba;
double anum11ibb,anum12ibb,anum13ibb,anum21ibb,anum22ibb,anum23ibb,anum31ibb,anum32ibb,anum33ibb;

void EvectorRead(void){
int i,j;
 
FILE *ina11raa,*ina12raa,*ina13raa,*ina21raa,*ina22raa,*ina23raa,*ina31raa,*ina32raa,*ina33raa;
FILE *ina11rab,*ina12rab,*ina13rab,*ina21rab,*ina22rab,*ina23rab,*ina31rab,*ina32rab,*ina33rab;
FILE *ina11rba,*ina12rba,*ina13rba,*ina21rba,*ina22rba,*ina23rba,*ina31rba,*ina32rba,*ina33rba;
FILE *ina11rbb,*ina12rbb,*ina13rbb,*ina21rbb,*ina22rbb,*ina23rbb,*ina31rbb,*ina32rbb,*ina33rbb;

FILE *ina11iaa,*ina12iaa,*ina13iaa,*ina21iaa,*ina22iaa,*ina23iaa,*ina31iaa,*ina32iaa,*ina33iaa;
FILE *ina11iab,*ina12iab,*ina13iab,*ina21iab,*ina22iab,*ina23iab,*ina31iab,*ina32iab,*ina33iab;
FILE *ina11iba,*ina12iba,*ina13iba,*ina21iba,*ina22iba,*ina23iba,*ina31iba,*ina32iba,*ina33iba;
FILE *ina11ibb,*ina12ibb,*ina13ibb,*ina21ibb,*ina22ibb,*ina23ibb,*ina31ibb,*ina32ibb,*ina33ibb;

ina11raa=fopen("Evector/a11raa.dat","r");
ina12raa=fopen("Evector/a12raa.dat","r");
ina13raa=fopen("Evector/a13raa.dat","r");
ina21raa=fopen("Evector/a21raa.dat","r");
ina22raa=fopen("Evector/a22raa.dat","r");
ina23raa=fopen("Evector/a23raa.dat","r");
ina31raa=fopen("Evector/a31raa.dat","r");
ina32raa=fopen("Evector/a32raa.dat","r");
ina33raa=fopen("Evector/a33raa.dat","r");

ina11rab=fopen("Evector/a11rab.dat","r");
ina12rab=fopen("Evector/a12rab.dat","r");
ina13rab=fopen("Evector/a13rab.dat","r");
ina21rab=fopen("Evector/a21rab.dat","r");
ina22rab=fopen("Evector/a22rab.dat","r");
ina23rab=fopen("Evector/a23rab.dat","r");
ina31rab=fopen("Evector/a31rab.dat","r");
ina32rab=fopen("Evector/a32rab.dat","r");
ina33rab=fopen("Evector/a33rab.dat","r");

ina11rba=fopen("Evector/a11rba.dat","r");
ina12rba=fopen("Evector/a12rba.dat","r");
ina13rba=fopen("Evector/a13rba.dat","r");
ina21rba=fopen("Evector/a21rba.dat","r");
ina22rba=fopen("Evector/a22rba.dat","r");
ina23rba=fopen("Evector/a23rba.dat","r");
ina31rba=fopen("Evector/a31rba.dat","r");
ina32rba=fopen("Evector/a32rba.dat","r");
ina33rba=fopen("Evector/a33rba.dat","r");

ina11rbb=fopen("Evector/a11rbb.dat","r");
ina12rbb=fopen("Evector/a12rbb.dat","r");
ina13rbb=fopen("Evector/a13rbb.dat","r");
ina21rbb=fopen("Evector/a21rbb.dat","r");
ina22rbb=fopen("Evector/a22rbb.dat","r");
ina23rbb=fopen("Evector/a23rbb.dat","r");
ina31rbb=fopen("Evector/a31rbb.dat","r");
ina32rbb=fopen("Evector/a32rbb.dat","r");
ina33rbb=fopen("Evector/a33rbb.dat","r");

ina11iaa=fopen("Evector/a11iaa.dat","r");
ina12iaa=fopen("Evector/a12iaa.dat","r");
ina13iaa=fopen("Evector/a13iaa.dat","r");
ina21iaa=fopen("Evector/a21iaa.dat","r");
ina22iaa=fopen("Evector/a22iaa.dat","r");
ina23iaa=fopen("Evector/a23iaa.dat","r");
ina31iaa=fopen("Evector/a31iaa.dat","r");
ina32iaa=fopen("Evector/a32iaa.dat","r");
ina33iaa=fopen("Evector/a33iaa.dat","r");

ina11iab=fopen("Evector/a11iab.dat","r");
ina12iab=fopen("Evector/a12iab.dat","r");
ina13iab=fopen("Evector/a13iab.dat","r");
ina21iab=fopen("Evector/a21iab.dat","r");
ina22iab=fopen("Evector/a22iab.dat","r");
ina23iab=fopen("Evector/a23iab.dat","r");
ina31iab=fopen("Evector/a31iab.dat","r");
ina32iab=fopen("Evector/a32iab.dat","r");
ina33iab=fopen("Evector/a33iab.dat","r");

ina11iba=fopen("Evector/a11iba.dat","r");
ina12iba=fopen("Evector/a12iba.dat","r");
ina13iba=fopen("Evector/a13iba.dat","r");
ina21iba=fopen("Evector/a21iba.dat","r");
ina22iba=fopen("Evector/a22iba.dat","r");
ina23iba=fopen("Evector/a23iba.dat","r");
ina31iba=fopen("Evector/a31iba.dat","r");
ina32iba=fopen("Evector/a32iba.dat","r");
ina33iba=fopen("Evector/a33iba.dat","r");

ina11ibb=fopen("Evector/a11ibb.dat","r");
ina12ibb=fopen("Evector/a12ibb.dat","r");
ina13ibb=fopen("Evector/a13ibb.dat","r");
ina21ibb=fopen("Evector/a21ibb.dat","r");
ina22ibb=fopen("Evector/a22ibb.dat","r");
ina23ibb=fopen("Evector/a23ibb.dat","r");
ina31ibb=fopen("Evector/a31ibb.dat","r");
ina32ibb=fopen("Evector/a32ibb.dat","r");
ina33ibb=fopen("Evector/a33ibb.dat","r");

for(i=0;i<nx/2;i++){
   for(j=0;j<ny/2;j++){

      fscanf(ina11raa,"%lf",&anum11raa);
      alpha[i][j][0][0][0]=anum11raa;
//      printf("%e\n",num11r1);
      

      fscanf(ina11iaa,"%lf",&anum11iaa);
     alpha[i][j][0][0][1]=anum11iaa;
//      printf("%e\n",num2);

      fscanf(ina12raa,"%lf",&anum12raa);
      alpha[i][j][0][1][0]=anum12raa;
 //      printf("%e\t",num3);

      fscanf(ina12iaa,"%lf",&anum12iaa);
      alpha[i][j][0][1][1]=anum12iaa;

      fscanf(ina13raa,"%lf",&anum13raa);
      alpha[i][j][0][2][0]=anum13raa;

      fscanf(ina13iaa,"%lf",&anum13iaa);
      alpha[i][j][0][2][1]=anum13iaa;

      fscanf(ina21raa,"%lf",&anum21raa);
      alpha[i][j][1][0][0]=anum21raa;

      fscanf(ina21iaa,"%lf",&anum21iaa);
      alpha[i][j][1][0][1]=anum21iaa;

      fscanf(ina22raa,"%lf",&anum22raa);
      alpha[i][j][1][1][0]=anum22raa;

      fscanf(ina22iaa,"%lf",&anum22iaa);
      alpha[i][j][1][1][1]=anum22iaa;

      fscanf(ina23raa,"%lf",&anum23raa);
      alpha[i][j][1][2][0]=anum23raa;

      fscanf(ina23iaa,"%lf",&anum23iaa);
      alpha[i][j][1][2][1]=anum23iaa;

      fscanf(ina31raa,"%lf",&anum31raa);
      alpha[i][j][2][0][0]=anum31raa;

      fscanf(ina31iaa,"%lf",&anum31iaa);
      alpha[i][j][2][0][1]=anum31iaa;


      fscanf(ina32raa,"%lf",&anum32raa);
      alpha[i][j][2][1][0]=anum32raa;

      fscanf(ina32iaa,"%lf",&anum32iaa);
      alpha[i][j][2][1][1]=anum32iaa;

      fscanf(ina33raa,"%lf",&anum33raa);
      alpha[i][j][2][2][0]=anum33raa;

      fscanf(ina33iaa,"%lf",&anum33iaa);
      alpha[i][j][2][2][1]=anum33iaa;


}}

for(i=nx/2;i<nx;i++){
  for(j=ny/2;j<ny;j++){
//for(i=0;i<nx/2;i++){
//   for(j=0;j<ny/2;j++){

      fscanf(ina11rbb,"%lf",&anum11rbb);
      alpha[i][j][0][0][0]=anum11rbb;
   //   printf("%e\n",num1);
      

      fscanf(ina11ibb,"%lf",&anum11ibb);
      alpha[i][j][0][0][1]=anum11ibb;
//      printf("%e\n",num2);

      fscanf(ina12rbb,"%lf",&anum12rbb);
      alpha[i][j][0][1][0]=anum12rbb;
 //      printf("%e\t",num3);

      fscanf(ina12ibb,"%lf",&anum12ibb);
      alpha[i][j][0][1][1]=anum12ibb;

      fscanf(ina13rbb,"%lf",&anum13rbb);
      alpha[i][j][0][2][0]=anum13rbb;

      fscanf(ina13ibb,"%lf",&anum13ibb);
      alpha[i][j][0][2][1]=anum13ibb;

      fscanf(ina21rbb,"%lf",&anum21rbb);
      alpha[i][j][1][0][0]=anum21rbb;

      fscanf(ina21ibb,"%lf",&anum21ibb);
      alpha[i][j][1][0][1]=anum21ibb;

      fscanf(ina22rbb,"%lf",&anum22rbb);
      alpha[i][j][1][1][0]=anum22rbb;

      fscanf(ina22ibb,"%lf",&anum22ibb);
      alpha[i][j][1][1][1]=anum22ibb;

      fscanf(ina23rbb,"%lf",&anum23rbb);
      alpha[i][j][1][2][0]=anum23rbb;

      fscanf(ina23ibb,"%lf",&anum23ibb);
      alpha[i][j][1][2][1]=anum23ibb;

      fscanf(ina31rbb,"%lf",&anum31rbb);
      alpha[i][j][2][0][0]=anum31rbb;

      fscanf(ina31ibb,"%lf",&anum31ibb);
      alpha[i][j][2][0][1]=anum31ibb;


      fscanf(ina32rbb,"%lf",&anum32rbb);
      alpha[i][j][2][1][0]=anum32rbb;

      fscanf(ina32ibb,"%lf",&anum32ibb);
      alpha[i][j][2][1][1]=anum32ibb;

      fscanf(ina33rbb,"%lf",&anum33rbb);
      alpha[i][j][2][2][0]=anum33rbb;

      fscanf(ina33ibb,"%lf",&anum33ibb);
      alpha[i][j][2][2][1]=anum33ibb;


}}
for(i=nx/2;i<nx;i++){
   for(j=0;j<ny/2;j++){


      fscanf(ina11rba,"%lf",&anum11rba);
      alpha[i][j][0][0][0]=anum11rba;
//      printf("%e\n",num11r1);
      

      fscanf(ina11iba,"%lf",&anum11iba);
      alpha[i][j][0][0][1]=anum11iba;
//      printf("%e\n",num2);

      fscanf(ina12rba,"%lf",&anum12rba);
      alpha[i][j][0][1][0]=anum12rba;
 //      printf("%e\t",num3);

      fscanf(ina12iba,"%lf",&anum12iba);
      alpha[i][j][0][1][1]=anum12iba;

      fscanf(ina13rba,"%lf",&anum13rba);
      alpha[i][j][0][2][0]=anum13rba;

      fscanf(ina13iba,"%lf",&anum13iba);
      alpha[i][j][0][2][1]=anum13iba;

      fscanf(ina21rba,"%lf",&anum21rba);
      alpha[i][j][1][0][0]=anum21rba;

      fscanf(ina21iba,"%lf",&anum21iba);
      alpha[i][j][1][0][1]=anum21iba;

      fscanf(ina22rba,"%lf",&anum22rba);
      alpha[i][j][1][1][0]=anum22rba;

      fscanf(ina22iba,"%lf",&anum22iba);
      alpha[i][j][1][1][1]=anum22iba;

      fscanf(ina23rba,"%lf",&anum23rba);
      alpha[i][j][1][2][0]=anum23rba;

      fscanf(ina23iba,"%lf",&anum23iba);
      alpha[i][j][1][2][1]=anum23iba;

      fscanf(ina31rba,"%lf",&anum31rba);
      alpha[i][j][2][0][0]=anum31rba;

      fscanf(ina31iba,"%lf",&anum31iba);
      alpha[i][j][2][0][1]=anum31iba;


      fscanf(ina32rba,"%lf",&anum32rba);
      alpha[i][j][2][1][0]=anum32rba;

      fscanf(ina32iba,"%lf",&anum32iba);
      alpha[i][j][2][1][1]=anum32iba;

      fscanf(ina33rba,"%lf",&anum33rba);
      alpha[i][j][2][2][0]=anum33rba;

      fscanf(ina33iba,"%lf",&anum33iba);
      alpha[i][j][2][2][1]=anum33iba;


}}

for(i=0;i<nx/2;i++){
   for(j=ny/2;j<ny;j++){

      fscanf(ina11rab,"%lf",&anum11rab);
      alpha[i][j][0][0][0]=anum11rab;
//      printf("%e \t alpha[%d][%d][0][0][0] \t %e\n",anum11rab,i,j,alpha[i][j][0][0][0]);
            

      fscanf(ina11iab,"%lf",&anum11iab);
      alpha[i][j][0][0][1]=anum11iab;
//      printf("%e\n",num2);


      fscanf(ina12rab,"%lf",&anum12rab);
      alpha[i][j][0][1][0]=anum12rab;
 //      printf("%e\t",num3);
//      printf("%e \t alpha[%d][%d][0][1][0] \t %e\n",anum12rab,i,j,alpha[i][j][0][1][0]);


      fscanf(ina12iab,"%lf",&anum12iab);
      alpha[i][j][0][1][1]=anum12iab;

      fscanf(ina13rab,"%lf",&anum13rab);
      alpha[i][j][0][2][0]=anum13rab;

      fscanf(ina13iab,"%lf",&anum13iab);
      alpha[i][j][0][2][1]=anum13iab;

      fscanf(ina21rab,"%lf",&anum21rab);
      alpha[i][j][1][0][0]=anum21rab;

      fscanf(ina21iab,"%lf",&anum21iab);
      alpha[i][j][1][0][1]=anum21iab;
   

      fscanf(ina22rab,"%lf",&anum22rab);
      alpha[i][j][1][1][0]=anum22rab;

      fscanf(ina22iab,"%lf",&anum22iab);
      alpha[i][j][1][1][1]=anum22iab;

      fscanf(ina23rab,"%lf",&anum23rab);
      alpha[i][j][1][2][0]=anum23rab;

      fscanf(ina23iab,"%lf",&anum23iab);
      alpha[i][j][1][2][1]=anum23iab;

      fscanf(ina31rab,"%lf",&anum31rab);
      alpha[i][j][2][0][0]=anum31rab;

      fscanf(ina31iab,"%lf",&anum31iab);
      alpha[i][j][2][0][1]=anum31iab;


      fscanf(ina32rab,"%lf",&anum32rab);
      alpha[i][j][2][1][0]=anum32rab;

      fscanf(ina32iab,"%lf",&anum32iab);
      alpha[i][j][2][1][1]=anum32iab;

      fscanf(ina33rab,"%lf",&anum33rab);
      alpha[i][j][2][2][0]=anum33rab;

      fscanf(ina33iab,"%lf",&anum33iab);
      alpha[i][j][2][2][1]=anum33iab;


}}
int m;


//=====This part of the code is commented out. This will be used when the system anisotropy is smaller than 1.=======================================
/*
printf("========\nWARNING : Eigen vectors are defined for Anisotropy < 1.00\n=======\n");
printf("========\nc44 = %lf\n=======\n",c44);

for(i=0;i<nx;i++){
   for(j=0;j<ny;j++){
       for(m=0;m<2;m++){
      alpha[i][j][0][3][m]=alpha[i][j][0][0][m];
      alpha[i][j][0][4][m]=alpha[i][j][0][1][m];
      alpha[i][j][0][5][m]=alpha[i][j][0][2][m];

      alpha[i][j][1][3][m]=alpha[i][j][1][0][m];
      alpha[i][j][1][4][m]=alpha[i][j][1][1][m];
      alpha[i][j][1][5][m]=alpha[i][j][1][2][m];

      alpha[i][j][2][3][m]=-alpha[i][j][2][0][m];
      alpha[i][j][2][4][m]=-alpha[i][j][2][1][m];
      alpha[i][j][2][5][m]=-alpha[i][j][2][2][m];

}}}

for(i=0;i<nx;i++){
alpha[i][0][1][4][0]=alpha[i][0][1][4][0];
alpha[i][0][1][4][1]=-alpha[i][0][1][4][1];
}
*/
//===============================================

//===============================================

printf("========\nINFO : Eigen vectors are defined for Anisotropy > 1.00\n=======\n");
printf("========\nc44 = %lf\n=======\n",c44);
for(i=0;i<nx;i++){
   for(j=0;j<ny;j++){
       for(m=0;m<2;m++){
      alpha[i][j][0][3][m]=-alpha[i][j][0][0][m];
      alpha[i][j][0][4][m]=-alpha[i][j][0][1][m];
      alpha[i][j][0][5][m]= alpha[i][j][0][2][m];

      alpha[i][j][1][3][m]=-alpha[i][j][1][0][m];
      alpha[i][j][1][4][m]=-alpha[i][j][1][1][m];
      alpha[i][j][1][5][m]= alpha[i][j][1][2][m];

      alpha[i][j][2][3][m]= alpha[i][j][2][0][m];
      alpha[i][j][2][4][m]= alpha[i][j][2][1][m];
      alpha[i][j][2][5][m]=-alpha[i][j][2][2][m];

}}}

for(i=0;i<nx;i++){
     for(m=0;m<2;m++){
      alpha[i][0][0][3][m]=-alpha[i][0][0][3][m];
      alpha[i][0][0][4][m]=-alpha[i][0][0][4][m];
      alpha[i][0][0][5][m]=-alpha[i][0][0][5][m];

      alpha[i][0][1][3][m]=-alpha[i][0][1][3][m];
      alpha[i][0][1][4][m]=-alpha[i][0][1][4][m];
      alpha[i][0][1][5][m]=-alpha[i][0][1][5][m];

      alpha[i][0][2][3][m]=-alpha[i][0][2][3][m];
      alpha[i][0][2][4][m]=-alpha[i][0][2][4][m];
      alpha[i][0][2][5][m]=-alpha[i][0][2][5][m];
}}

fclose(ina11raa);
fclose(ina12raa);
fclose(ina13raa);
fclose(ina21raa);
fclose(ina22raa);
fclose(ina23raa);
fclose(ina31raa);
fclose(ina32raa);
fclose(ina33raa);

fclose(ina11iaa);
fclose(ina12iaa);
fclose(ina13iaa);
fclose(ina21iaa);
fclose(ina22iaa);
fclose(ina23iaa);
fclose(ina31iaa);
fclose(ina32iaa);
fclose(ina33iaa);

fclose(ina11rab);
fclose(ina12rab);
fclose(ina13rab);
fclose(ina21rab);
fclose(ina22rab);
fclose(ina23rab);
fclose(ina31rab);
fclose(ina32rab);
fclose(ina33rab);

fclose(ina11iab);
fclose(ina12iab);
fclose(ina13iab);
fclose(ina21iab);
fclose(ina22iab);
fclose(ina23iab);
fclose(ina31iab);
fclose(ina32iab);
fclose(ina33iab);

fclose(ina11rba);
fclose(ina12rba);
fclose(ina13rba);
fclose(ina21rba);
fclose(ina22rba);
fclose(ina23rba);
fclose(ina31rba);
fclose(ina32rba);
fclose(ina33rba);

fclose(ina11iba);
fclose(ina12iba);
fclose(ina13iba);
fclose(ina21iba);
fclose(ina22iba);
fclose(ina23iba);
fclose(ina31iba);
fclose(ina32iba);
fclose(ina33iba);

fclose(ina11rbb);
fclose(ina12rbb);
fclose(ina13rbb);
fclose(ina21rbb);
fclose(ina22rbb);
fclose(ina23rbb);
fclose(ina31rbb);
fclose(ina32rbb);
fclose(ina33rbb);

fclose(ina11ibb);
fclose(ina12ibb);
fclose(ina13ibb);
fclose(ina21ibb);
fclose(ina22ibb);
fclose(ina23ibb);
fclose(ina31ibb);
fclose(ina32ibb);
fclose(ina33ibb);

}
