#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DIM 3

typedef double matrix[3][3];
typedef float  fmatrix[3][3];

float norm(float a[])
{
  return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

int main(int argc, char **argv)
{

  int n,d,i,j,k,flnr;
  int bdouble;
  int nx,ny,nz;
  float com[3],scale[3];
  float dx[3],length,dbin;
  int bin[3],ndx;
  float *HISTO,*RPRES,*A1PRES,*A2PRES,*C01PRES,*C10PRES;

  double  box[3][3];
  float   fbox[3][3];
  matrix *localpressure;
  matrix *slice;
  fmatrix tmpmat;
  double max[3];
  double min[3];
  double visualfilter=10000;  /*You might want to filter signal for rasmol visualation*/

  float r_comp[3][3],R[3][3],xylength;

  FILE *fp,*out,*out3D,*outRASMOL;

  /* READ THE BINARY FILE */
  
  for(flnr=1;flnr<argc;flnr++){
    fp=fopen(argv[flnr],"r");
    
    fread(&bdouble,sizeof(int),1,fp);
    
    if(bdouble)
      {
	fread(box,sizeof(double),9,fp);    
      }
    else
      {
	fread(fbox,sizeof(float),9,fp);
	for(i=0;i<3;i++)
	  for(j=0;j<3;j++)
	    box[i][j]=fbox[i][j];
      }
    
    fread(&nx,sizeof(int),1,fp);
    fread(&ny,sizeof(int),1,fp);
    fread(&nz,sizeof(int),1,fp);
    
    if(flnr==1){
      localpressure=malloc(sizeof(matrix)*nx*ny*nz);
      slice = malloc(sizeof(matrix)*nz);
    }
    
    if(bdouble)
      {
	fread(localpressure,sizeof(matrix),nx*ny*nz,fp);
      }
    else
      {
	for(k=0;k<nx*ny*nz;k++)
	  {
	    fread(tmpmat,sizeof(tmpmat),1,fp);
	    for(i=0;i<3;i++)
	      for(j=0;j<3;j++)
		localpressure[k][i][j]+=tmpmat[i][j];
	  }
      }
    
    fclose(fp);
    
  }

  /* AVERAGE OVER FILES */
  /*printf("%i",argc-1);*/
  for(k=0;k<nx*ny*nz;k++)
    {
      for(i=0;i<3;i++)
	for(j=0;j<3;j++)
	  localpressure[k][i][j]=localpressure[k][i][j]/(argc-1);
    }
  
  /* PRINT 3D PRESSURE FIELD */
  
  out3D=fopen("3Dpp.dat","w");
  
  for(i=0;i<nx;i++)
    for(j=0;j<ny;j++)   
      for(k=0;k<nz;k++)
	{
	  fprintf(out3D,"%2d %d %2d  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f\n",
		  i,j,k,
		  localpressure[i*ny*nz+j*nz+k][0][0],
		  localpressure[i*ny*nz+j*nz+k][0][1],
		  localpressure[i*ny*nz+j*nz+k][0][2],
		  localpressure[i*ny*nz+j*nz+k][1][0],
		  localpressure[i*ny*nz+j*nz+k][1][1],
		  localpressure[i*ny*nz+j*nz+k][1][2],
		  localpressure[i*ny*nz+j*nz+k][2][0],
		  localpressure[i*ny*nz+j*nz+k][2][1],
		  localpressure[i*ny*nz+j*nz+k][2][2]);
	}
  
     
  /* PRINT THE RASMOL FILE FOR VISUALATION */
    
  outRASMOL=fopen("3Dppras.mol","w");
  
  for(k=0;k<nx*ny*nz;k++)
    {
      for(i=0;i<3;i++){	
	if(localpressure[k][i][i] > max[i])
	  max[i] = localpressure[k][i][i];
	if(localpressure[k][i][i] < min[i])
	  min[i] = localpressure[k][i][i];
      }
    }
  
  for(i=0;i<nx;i++)
    for(j=0;j<ny;j++)   
      for(k=0;k<nz;k++)
	{
	  
	  if(localpressure[i*ny*nz+j*nz+k][0][0]>(-1)*visualfilter && localpressure[i*ny*nz+j*nz+k][0][0]< visualfilter)
	    fprintf(outRASMOL,"%-6s%5u  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f \n","ATOM ",(i+j+k+1)%10000,"CA","ALA",' ',(i+j+k+1)%10000,i*1.0,j*1.0,k*1.0,1.0,(localpressure[i*ny*nz+j*nz+k][0][0]+max[0])/(max[0]-min[0]));
	  
	}
  fclose(outRASMOL);
  fclose(out3D);
  
  /*AVERAGE OVER SPHERICAL ANGLE*/

  bin[0]=nx;   /* x bins */
  bin[1]=ny;   /* y bins */
  bin[2]=nz;   /* z bins */ 
  
  com[0]=bin[0]/2;
  com[1]=bin[1]/2;
  com[2]=bin[2]/2;

  printf("Your origin of the coordinate system is set into voxel %i %i %i (should be COM of your spherical system) \n If you want to change this go to the line 149-151 in the source code \n ",(int)com[0],(int)com[1],(int)com[2]);

  dbin=norm(com)/(nx);
  
  HISTO = (float*)malloc(nx*sizeof(float));
  RPRES = (float*)malloc(nx*sizeof(float));
  A1PRES = (float*)malloc(nx*sizeof(float));
  A2PRES = (float*)malloc(nx*sizeof(float));
  /*  C01PRES = (float*)malloc(nx*sizeof(float));
      C10PRES = (float*)malloc(nx*sizeof(float));*/

  for(i=0;i<nx;i++)
    {
      RPRES[i]  = 0;
      A1PRES[i] = 0;
      A2PRES[i] = 0;      
      HISTO[i]  = 0;
    }
 

  for(i=0;i<nx;i++)
    for(j=0;j<ny;j++)   
      for(k=0;k<nz;k++)
	{
	  dx[0]= i-com[0];
	  dx[1]= j-com[1];
	  dx[2]= k-com[2];
	  
	  length = norm(dx);
	  
	  for(d=0;d<DIM;d++)
	    {   
	      R[0][d]=dx[d]/length;
	      for(n=0;n<DIM;n++)
		r_comp[d][n]=0;
	    }
	  
	  if(sqrt(dx[0]*dx[0]+dx[1]*dx[1])>0.000000001)
	    xylength=1/sqrt(dx[0]*dx[0]+dx[1]*dx[1]);
	  R[1][0]=(-1)*dx[1]*xylength;
	  R[1][1]=dx[0]*xylength;
	  R[1][2]=0;
  
	  R[2][0]=dx[0]*dx[2]*xylength/length;
	  R[2][1]=dx[1]*dx[2]*xylength/length;
	  R[2][2]=(-1)/(xylength*length);
  
  
	  for(n=0;n<DIM;n++){
	    /*diagonal*/
	    for(d=0;d<DIM;d++){
	      r_comp[n][n] += R[n][d]*R[n][d]*localpressure[i*ny*nz+j*nz+k][d][d];
	    }
	    
	    /*off diagonal*/
	    r_comp[n][n] += R[n][0]*R[n][1]*(localpressure[i*ny*nz+j*nz+k][0][1] + localpressure[i*ny*nz+j*nz+k][1][0]);
	    r_comp[n][n] += R[n][1]*R[n][2]*(localpressure[i*ny*nz+j*nz+k][1][2] + localpressure[i*ny*nz+j*nz+k][2][1]);
	    r_comp[n][n] += R[n][0]*R[n][2]*(localpressure[i*ny*nz+j*nz+k][0][2] + localpressure[i*ny*nz+j*nz+k][2][0]); 
	  }


	  /* ONLY DIAGONAL COMPONENTS IN SPHERICAL COORDINATE SYSTEM ARE CALCULATED, 
	     BECAUSE NON-DIAGONALS SHOULD BE ZERO IN SPHERICAL SYMMETRIC SYSTEM. 
	     IF YOU WANT TO GET NON-DIAGONAL COMPONENTS YOU CAN USE EQUATIONS BELOW FOR [0,1] and [1,0] COMPONENTS BY REMOVING
	     COMMENTS FROM HERE AND OTHER PARTS OF THE CODE WHERE r_comp[0][1] and r_comp[1][0] COMPONENTS OCCUR. */ 

	  /*  r_comp[0][1]=localpressure[i*ny*nz+j*nz+k][0][0]*R[0][0]*R[1][0] + localpressure[i*ny*nz+j*nz+k][1][0]*R[0][1]*R[1][0] + localpressure[i*ny*nz+j*nz+k][2][0]*R[0][2]*R[1][0] + localpressure[i*ny*nz+j*nz+k][0][1]*R[0][0]*R[1][1] + localpressure[i*ny*nz+j*nz+k][1][1]*R[0][1]*R[1][1] +  localpressure[i*ny*nz+j*nz+k][2][1]*R[0][2]*R[1][1] + localpressure[i*ny*nz+j*nz+k][0][2]*R[0][0]* R[1][2] + localpressure[i*ny*nz+j*nz+k][1][2]*R[0][1]*R[1][2] + localpressure[i*ny*nz+j*nz+k][2][2]*R[0][2]*R[1][2];
	  
	      r_comp[1][0]= localpressure[i*ny*nz+j*nz+k][0][0]*R[0][0]*R[1][0] + localpressure[i*ny*nz+j*nz+k][0][1]*R[0][1]*R[1][0] + localpressure[i*ny*nz+j*nz+k][0][2]*R[0][2]*R[1][0] + localpressure[i*ny*nz+j*nz+k][1][0]*R[0][0]*R[1][1] + localpressure[i*ny*nz+j*nz+k][1][1]*R[0][1]*R[1][1] +  localpressure[i*ny*nz+j*nz+k][1][2]*R[0][2]*R[1][1] + localpressure[i*ny*nz+j*nz+k][2][0]*R[0][0]*R[1][2] + localpressure[i*ny*nz+j*nz+k][2][1]*R[0][1]*R[1][2] + localpressure[i*ny*nz+j*nz+k][2][2]*R[0][2]*R[1][2];*/
	  
	  ndx=(int)(length/dbin);
	  
	  if(ndx < nx)
	    {
	      RPRES[ndx] += r_comp[0][0];
	      A1PRES[ndx] += r_comp[1][1];
	      A2PRES[ndx] += r_comp[2][2];
	      /*  C01PRES[ndx] += r_comp[0][1];
		  C10PRES[ndx] += r_comp[1][0];*/
	      
	      HISTO[ndx]++;
	    }
	  
	}
  


  out=fopen("rad_pres.xvg","w");

  fprintf(out,"# Angular and radial components of the pressure field averaged over angle \n");
  fprintf(out,"# r      P_theta   P_phi   P_r \n ");
  
  for(n=1;n<nx;n++)
    if(HISTO[n]!=0)
      fprintf(out,"%f %f %f %f \n",(float)(n*dbin),A1PRES[n]/HISTO[n],A2PRES[n]/HISTO[n],RPRES[n]/HISTO[n]);
  
  
  
  fclose(out);
 
  free(HISTO);
  free(RPRES);
  free(A1PRES);
  free(A2PRES);
  
}
