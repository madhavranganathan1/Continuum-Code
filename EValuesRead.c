/* Real and imaginary part of eigenvalues of a 6X6 matrix for each value of kx and ky are are read and stored in arrays. */

double vnum1raa,vnum2raa,vnum3raa;
double vnum1rab,vnum2rab,vnum3rab;
double vnum1rba,vnum2rba,vnum3rba;
double vnum1rbb,vnum2rbb,vnum3rbb;
double vnum1iaa,vnum2iaa,vnum3iaa;
double vnum1iab,vnum2iab,vnum3iab;
double vnum1iba,vnum2iba,vnum3iba;
double vnum1ibb,vnum2ibb,vnum3ibb;

void EValuesRead(void){
int i,j;

FILE *ind1raa,*ind2raa,*ind3raa;
FILE *ind1rab,*ind2rab,*ind3rab;
FILE *ind1rba,*ind2rba,*ind3rba;
FILE *ind1rbb,*ind2rbb,*ind3rbb;

FILE *ind1iaa,*ind2iaa,*ind3iaa;
FILE *ind1iab,*ind2iab,*ind3iab;
FILE *ind1iba,*ind2iba,*ind3iba;
FILE *ind1ibb,*ind2ibb,*ind3ibb;

/*The data files containing the tables of eigenvalues is generated by a Mathematica program */

ind1raa=fopen("Evalues/v1raa.dat","r");
ind2raa=fopen("Evalues/v2raa.dat","r");
ind3raa=fopen("Evalues/v3raa.dat","r");

ind1rab=fopen("Evalues/v1rab.dat","r");
ind2rab=fopen("Evalues/v2rab.dat","r");
ind3rab=fopen("Evalues/v3rab.dat","r");

ind1rba=fopen("Evalues/v1rba.dat","r");
ind2rba=fopen("Evalues/v2rba.dat","r");
ind3rba=fopen("Evalues/v3rba.dat","r");

ind1rbb=fopen("Evalues/v1rbb.dat","r");
ind2rbb=fopen("Evalues/v2rbb.dat","r");
ind3rbb=fopen("Evalues/v3rbb.dat","r");

ind1iaa=fopen("Evalues/v1iaa.dat","r");
ind2iaa=fopen("Evalues/v2iaa.dat","r");
ind3iaa=fopen("Evalues/v3iaa.dat","r");

ind1iab=fopen("Evalues/v1iab.dat","r");
ind2iab=fopen("Evalues/v2iab.dat","r");
ind3iab=fopen("Evalues/v3iab.dat","r");

ind1iba=fopen("Evalues/v1iba.dat","r");
ind2iba=fopen("Evalues/v2iba.dat","r");
ind3iba=fopen("Evalues/v3iba.dat","r");

ind1ibb=fopen("Evalues/v1ibb.dat","r");
ind2ibb=fopen("Evalues/v2ibb.dat","r");
ind3ibb=fopen("Evalues/v3ibb.dat","r");

for(i=0;i<nx/2;i++){
   for(j=0;j<ny/2;j++){

      fscanf(ind1raa,"%lf",&vnum1raa);
        v[i][j][0][0]=vnum1raa;
//      printf("%e\n",num1r1d;
      

      fscanf(ind1iaa,"%lf",&vnum1iaa);
        v[i][j][0][1]=vnum1iaa;
//      printf("%e\n",num2);

      fscanf(ind2raa,"%lf",&vnum2raa);
        v[i][j][1][0]=vnum2raa;
 //      printf("%e\t",num3);

      fscanf(ind2iaa,"%lf",&vnum2iaa);
        v[i][j][1][1]=vnum2iaa;

      fscanf(ind3raa,"%lf",&vnum3raa);
        v[i][j][2][0]=vnum3raa;

      fscanf(ind3iaa,"%lf",&vnum3iaa);
        v[i][j][2][1]=vnum3iaa;


}} 

for(i=nx/2;i<nx;i++){
  for(j=ny/2;j<ny;j++){
//for(i=0;i<nx/2;i++){
//   for(j=0;j<ny/2;j++){

      fscanf(ind1rbb,"%lf",&vnum1rbb);
        v[i][j][0][0]=vnum1rbb;
   //   printf("%e\n",num1);
      

      fscanf(ind1ibb,"%lf",&vnum1ibb);
        v[i][j][0][1]=vnum1ibb;
//      printf("%e\n",num2);

      fscanf(ind2rbb,"%lf",&vnum2rbb);
        v[i][j][1][0]=vnum2rbb;
 //      printf("%e\t",num3);

      fscanf(ind2ibb,"%lf",&vnum2ibb);
        v[i][j][1][1]=vnum2ibb;

      fscanf(ind3rbb,"%lf",&vnum3rbb);
        v[i][j][2][0]=vnum3rbb;

      fscanf(ind3ibb,"%lf",&vnum3ibb);
        v[i][j][2][1]=vnum3ibb;


}}
for(i=nx/2;i<nx;i++){
   for(j=0;j<ny/2;j++){


      fscanf(ind1rba,"%lf",&vnum1rba);
        v[i][j][0][0]=vnum1rba;
//      printf("%e\n",num1r1d;
      

      fscanf(ind1iba,"%lf",&vnum1iba);
        v[i][j][0][1]=vnum1iba;
//      printf("%e\n",num2);

      fscanf(ind2rba,"%lf",&vnum2rba);
        v[i][j][1][0]=vnum2rba;
 //      printf("%e\t",num3);

      fscanf(ind2iba,"%lf",&vnum2iba);
        v[i][j][1][1]=vnum2iba;

      fscanf(ind3rba,"%lf",&vnum3rba);
        v[i][j][2][0]=vnum3rba;

      fscanf(ind3iba,"%lf",&vnum3iba);
        v[i][j][2][1]=vnum3iba;


}}


for(i=0;i<nx/2;i++){
   for(j=ny/2;j<ny;j++){


      fscanf(ind1rab,"%lf",&vnum1rab);
        v[i][j][0][0]=vnum1rab;
//      printf("%e\n",num1r1d;
      

      fscanf(ind1iab,"%lf",&vnum1iab);
        v[i][j][0][1]=vnum1iab;
//      printf("%e\n",num2);

      fscanf(ind2rab,"%lf",&vnum2rab);
       v[i][j][1][0]=vnum2rab;
 //      printf("%e\t",num3);

      fscanf(ind2iab,"%lf",&vnum2iab);
        v[i][j][1][1]=vnum2iab;

      fscanf(ind3rab,"%lf",&vnum3rab);
        v[i][j][2][0]=vnum3rab;

      fscanf(ind3iab,"%lf",&vnum3iab);
        v[i][j][2][1]=vnum3iab;


}}

/* The Matrix symmetry ensures that if lambda is an eigenvalue, then -lambda is also an eigenvalue */

int m;
for(i=0;i<nx;i++){
   for(j=0;j<ny;j++){
      for(m=0;m<2;m++){
          v[i][j][3][m]= -v[i][j][0][m];
          v[i][j][4][m]= -v[i][j][1][m];
          v[i][j][5][m]= -v[i][j][2][m];

}}}



fclose(ind1raa);
fclose(ind2raa);
fclose(ind3raa);

fclose(ind1iaa);
fclose(ind2iaa);
fclose(ind3iaa);

fclose(ind1rab);
fclose(ind2rab);
fclose(ind3rab);

fclose(ind1iab);
fclose(ind2iab);
fclose(ind3iab);

fclose(ind1rba);
fclose(ind2rba);
fclose(ind3rba);

fclose(ind1iba);
fclose(ind2iba);
fclose(ind3iba);

fclose(ind1rbb);
fclose(ind2rbb);
fclose(ind3rbb);

fclose(ind1ibb);
fclose(ind2ibb);
fclose(ind3ibb);

}
