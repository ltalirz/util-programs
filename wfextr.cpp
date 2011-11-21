#include<stdio>
#include<stdlib>
#include<string>
#include<fftw3>
#include<math>
#include<complex>

#include"mystruct.h"
#include"myio.h"
#include"myreplica.h"
#include"mydata_analysis.h"
//#include"mystm.h"
#include"mygrid_tools.h"


int main( int argc, char **argv) {
	FILE *fline;
	char eigfile[FILENAME_MAX],lista_wf[FILENAME_MAX],hfile[FILENAME_MAX],nomewf[1000],fileline[1000];
int a,b,c,i,j,k,nyh,rt=0,nstates=0,ncf_up=0,ncf_down,r=0,*listaci_up,*listacj_up,*listaci_down,*listacj_down,extnz_top,extnz_bot,u=0,extrn,h,t;
struct eig eig;
struct cube hartree,wf,extcube;
double *rho_up,*rho_down,*hav,***rrval;
double heigth,sl,hvac,Vb,exp_u,eignow,extrl,m1,m2;
fftw_complex *outrho_up,*outrho_down;
fftw_plan plan_forward;
double complex C,Cexp;


if(argc != 6 || strcmp(argv[1],"--help") == 0)
{
printf("The program needs the following arguments:\n");
printf(" 1) the output file of the CP2K calculations containing the eigenvalues\n 2) the file with the list of wavefunction cubes you wish to extrapolate\n 3) the hartree cube\n");
printf(" 4) the distance in a.u. between the extrapolation plane and the outermost atom\n 5) how many a.u. you want to extrapolate\n");
exit(1);
}

strcpy(eigfile,argv[1]);
strcpy(lista_wf,argv[2]);
strcpy(hfile,argv[3]);
heigth   = atof(argv[4]);
extrl    = atof(argv[5]);


//reading eigenvalues
eig    = read_eig(eigfile);

//reading hartree
hartree = read_cube(hfile);
hav=malloc(hartree.nz*sizeof(double));

//opening hartree and averaging on a plane
FILE *pip;
pip=fopen("hartreeout","w");
for(i=0;i<hartree.nz;i++)
 {
 hav[i]=avplane(hartree,i);
 fprintf(pip,"%lf\n",hav[i]);
 }
fclose(pip);

//find top surface
sl=getmax(hartree.z,hartree.nat);
k=(int)((sl+heigth)/hartree.dx[2][2]); //setting the surface at height a.u. on top of last atom
hvac=hav[k];
extrn=(int)(extrl/hartree.dx[2][2]);
extnz_top=extrn;
//find bottom
sl=getmin(hartree.z,hartree.nat);
t=(int)((sl-heigth)/hartree.dx[2][2]); //setting the surface at height a.u. on top of last atom
if(t<0)
{
printf("not enough space left under bottom surface for the extrapolation.\n Use a program to put the slab at the center of the cube file. Quitting\n");
exit(1);
}

extrn=(int)(extrl/hartree.dx[2][2]);
extnz_bot=extrn-t;




//put Vinf to 0
for(i=0;i<eig.norb;i++)
eig.val[i]-=hvac;
eig.Fermi-=hvac;
hvac=0;


//reading states  cube files
FILE *fw;
fw=fopen(lista_wf,"r");

//
rho_up=malloc(hartree.nx*hartree.ny*sizeof(double));
rho_down=malloc(hartree.nx*hartree.ny*sizeof(double));
outrho_up   = fftw_malloc ( sizeof ( fftw_complex ) * hartree.nx*hartree.ny*1);
outrho_down = fftw_malloc ( sizeof ( fftw_complex ) * hartree.nx*hartree.ny*1);
int count = 0;
int dum=1;
int ki=0,kj=0;
double vec1,vec2,scal=1/(double)(hartree.nx*hartree.ny);
char nameext[FILENAME_MAX];

extnz_bot=t-extrn;
extnz_top=k+extrn;

k=k-(t-extrn); // this shifts the upper plane accordingly
t=extrn;

while(fscanf(fw,"%s",&nomewf)!=EOF)
{

wf=read_cube(nomewf);
wf=shift_cube(wf,extnz_bot,extnz_top,"for_extrapolation");


// fourier transforming at plane k
//printf("...computing forward fft ---> fft(rho-up)..... \n");
 for(i=0;i<wf.nx;i++)
 for(j=0;j<wf.ny;j++)
 rho_up[0+(j+wf.ny*i)]=wf.val[i][j][k];



plan_forward = fftw_plan_dft_r2c_3d ( wf.nx, wf.ny, dum,  rho_up, outrho_up, FFTW_ESTIMATE );
fftw_execute ( plan_forward );
fftw_destroy_plan ( plan_forward );

//printf("...computing forward fft ---> fft(rho-down).....\n");
 
 for(i=0;i<wf.nx;i++)
 for(j=0;j<wf.ny;j++)
 rho_down[0+(j+wf.ny*i)]=wf.val[i][j][t];

plan_forward = fftw_plan_dft_r2c_3d ( wf.nx, wf.ny, dum,  rho_down, outrho_down, FFTW_ESTIMATE );
fftw_execute ( plan_forward );
fftw_destroy_plan ( plan_forward );


if(count==0)
 {
  rrval=malloc(wf.nx*sizeof(double**));

   for(i=0;i<wf.nx;i++)
   {
   rrval[i]=malloc(wf.ny*sizeof(double*));
    for(j=0;j<wf.ny;j++)
    rrval[i][j]=malloc(wf.nz*sizeof(double));
   }
 }

//choosing the state to be extrapolated
eignow=eig.val[eig.norb-abs(wf.state)-1];

ncf_up=0;
r=0;
//reducing the number of Fourier components
 for(i=0;i<wf.nx;i=i+1)
 for(j=0;j<wf.ny;j=j+1)
 {
 if(cabs(outrho_up[j+wf.ny*i][0]+_Complex_I*outrho_up[j+wf.ny*i][1])/cabs(outrho_up[0][0]+_Complex_I*outrho_up[0][1])>5E-2)
  ncf_up++;
 
 }

//creating lists for the coefficients to be retained
listaci_up=malloc(wf.nx*wf.ny*sizeof(int));
listacj_up=malloc(wf.nx*wf.ny*sizeof(int));

//filling up lists
for(i=0;i<wf.nx;i=i+1)
for(j=0;j<wf.ny;j=j+1)
{
 if(cabs(outrho_up[j+wf.ny*i][0]+_Complex_I*outrho_up[j+wf.ny*i][1])/cabs(outrho_up[0][0]+_Complex_I*outrho_up[0][1])>5E-2)
 {
  //listacf[r]=j+wf.ny*i;
  listaci_up[r]=i;
  listacj_up[r]=j;
  r++;
 }
}
printf("retaining %d Fourier coefficients (top). Total number is %d\n",ncf_up,wf.nx*wf.ny);


///-------------------------DOWN------------------------------


r=0;
ncf_down=0;
//reducing the number of Fourier components
 for(i=0;i<wf.nx;i=i+1)
 for(j=0;j<wf.ny;j=j+1)
 {
 if(cabs(outrho_down[j+wf.ny*i][0]+_Complex_I*outrho_down[j+wf.ny*i][1])/cabs(outrho_down[0][0]+_Complex_I*outrho_down[0][1])>5E-2)
  ncf_down++;
 }

//creating lists for the coefficients to be retained
listaci_down=malloc(wf.nx*wf.ny*sizeof(int));
listacj_down=malloc(wf.nx*wf.ny*sizeof(int));

//filling up lists
for(i=0;i<wf.nx;i=i+1)
for(j=0;j<wf.ny;j=j+1)
{
 if(cabs(outrho_down[j+wf.ny*i][0]+_Complex_I*outrho_down[j+wf.ny*i][1])/cabs(outrho_down[0][0]+_Complex_I*outrho_down[0][1])>5E-2)
 {
  //listacf[r]=j+wf.ny*i;
  listaci_down[r]=i;
  listacj_down[r]=j;
  r++;
 }
}
printf("retaining %d Fourier coefficients (bottom). Total number is %d\n\n",ncf_down,wf.nx*wf.ny);




sprintf(nameext,"ext_%d.cube",wf.state);

// doing extrapolation
 for(u=0;u<wf.nz;u++)
 {
  //printf("we are at %d and k=%d\n",u,k);
  if(u<=k && u>=t) //copy the prexisting value
   {
     for(i=0;i<wf.nx;i=i+1)
     for(j=0;j<wf.ny;j=j+1)
     rrval[i][j][u]=wf.val[i][j][u];
   }
  else if(u>k)  //extrapolate
   {
    //printf("entering ext %d k is %d and top is %d\n",u,k,wf.nz);
    for(i=0;i<wf.nx;i=i+1)
     {
      
    for(j=0;j<wf.ny;j=j+1)
       {
       rrval[i][j][u]=0;
       C=0;
       
       for(ki=0;ki<ncf_up;ki++)
        {
        
        vec1=((2*M_PI)/wf.cell[0][0])*listaci_up[ki];
        if(listaci_up[ki]>wf.nx/2)
        vec1=((2*M_PI)/wf.cell[0][0])*(listaci_up[ki]-wf.nx);

        vec2=((2*M_PI)/wf.cell[1][1])*listacj_up[ki];
        if(listacj_up[ki]>wf.ny/2)
        vec2=((2*M_PI)/wf.cell[1][1])*(listacj_up[ki]-wf.ny);

        
        Cexp            =          _Complex_I*(vec1*i*wf.dx[0][0]+vec2*j*wf.dx[1][1]);    
        exp_u           =           exp(-sqrt(2*(hvac-eignow)+vec1*vec1+vec2*vec2)*(u-k)*wf.dx[2][2]);   
        rrval[i][j][u] +=           creal(scal*(outrho_up[(listacj_up[ki]+wf.ny*listaci_up[ki])][0] + _Complex_I*outrho_up[(listacj_up[ki]+wf.ny*listaci_up[ki])][1])*cexp(Cexp)*exp_u);
        //printf("value is %6.6le %d %d %d %lf %lf %d %d %d\n",rrval[i][j][u],i,j,u,sqrt(2*(hvac-eignow)+vec1*vec1+vec2*vec2),vec1,ki,listaci_up[ki],listacj_up[ki]);
        } // end for ki
      
       }//end for j
     }//end for i

    }//end else up
   else if(u<t)  //extrapolate
   {
    //printf("entering ext %d k is %d and top is %d\n",u,k,wf.nz);
    for(i=0;i<wf.nx;i=i+1)
     {
      
    for(j=0;j<wf.ny;j=j+1)
       {
       rrval[i][j][u]=0;
       C=0;
       
       for(ki=0;ki<ncf_down;ki++)
        {
        
        vec1=((2*M_PI)/wf.cell[0][0])*listaci_down[ki];
        if(listaci_down[ki]>wf.nx/2)
        vec1=((2*M_PI)/wf.cell[0][0])*(listaci_down[ki]-wf.nx);

        vec2=((2*M_PI)/wf.cell[1][1])*listacj_down[ki];
        if(listacj_down[ki]>wf.ny/2)
        vec2=((2*M_PI)/wf.cell[1][1])*(listacj_down[ki]-wf.ny);

        
        Cexp            =          _Complex_I*(vec1*i*wf.dx[0][0]+vec2*j*wf.dx[1][1]);    
        exp_u           =           exp(-sqrt(2*(hvac-eignow)+vec1*vec1+vec2*vec2)*(t-u)*wf.dx[2][2]);   
        rrval[i][j][u] +=           creal(scal*(outrho_down[(listacj_down[ki]+wf.ny*listaci_down[ki])][0] + _Complex_I*outrho_down[(listacj_down[ki]+wf.ny*listaci_down[ki])][1])*cexp(Cexp)*exp_u);
        //printf("value is %6.6le %d %d %d %lf %lf %d %d %d\n",rrval[i][j][u],i,j,u,sqrt(2*(hvac-eignow)+vec1*vec1+vec2*vec2),vec1,ki,listaci_up[ki],listacj_up[ki]);
        } // end for ki
      
       }//end for j
     }//end for i

    }//end else down

//exit(1);
   }//end for u

sprintf(fileline,"line_%d.dat",wf.state);
fline=fopen(fileline,"w");
for(i=0;i<wf.nz;i++)
{
 for(j=0;j<wf.nx;j++)
 for(h=0;h<wf.ny;h++)
 {
 m1+=rrval[j][h][i]*scal;
 m2+=wf.val[j][h][i]*scal;
 }
 fprintf(fline,"%6.6lf  %6.6le  %6.6le\n",wf.nz*i,m1,m2);
m1=0;
m2=0;
}
fclose(fline);

extcube = create_cube(wf.nat,wf.nx,wf.ny,wf.nz,wf.dx,rrval,wf.x,wf.y,wf.z,wf.state);
print_cube( extcube, nameext, "wf");

count++;

 }// end while

return(0);
}
