#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DIM 3

typedef double matrix[3][3];
typedef float  fmatrix[3][3];

float norm(float a[])
{
  return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

int main(int argc, char **argv)
{
  int i,j,k,flnr;
  int bdouble;
  int nx,ny,nz;
  int d,dd;
  float   fbox[3][3];
  double  box[3][3];
  double max[3];
  double min[3];
  double dbin;
  
  matrix *localpressure;
  matrix *slice;
  fmatrix tmpmat;
  matrix average;
  FILE *fp,*outRASMOL,*out3D,*out2D,*outCYL,*outRTENS,*outKAPPAC0;

  int bin[3],ndx,bSKIP;
  float dx[3],length, P;
  float r_comp[3],com[3];
  double ***cil_p=NULL;

  double visualfilter=10000;     /*YOU MAY WANT TO FILTER VALUES FOR RASMOL VISUALATION*/
  double scale = 0.3;           /*THIS IS THE SIZE OF YOUR VOXEL */

  double tottension,cylradius=0,sum1;
  double zcenter=-1,bilayerTOP=2,bilayerBOT=-4,zcenterTOP=0,zcenterBOT=-2;  /*LOCATION OF UPPER AND LOWER MONOLAYER DEFINED*/
  
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
  
   
    /* AVERAGE OVER PLANAR, "THE LATERAL PRESSURE PROFILE" */

    out2D=fopen("2Dlpp.dat","w");

    for(d=0;d<3;d++)
      for(dd=0;dd<3;dd++)
	average[d][dd]=0;
    
    for(k=0;k<nz;k++)
      {
        for(d=0;d<3;d++)
	  for(dd=0;dd<3;dd++)
	    slice[k][d][dd]=0;
        
        for(i=0;i<nx;i++)
	  for(j=0;j<ny;j++)
	    for(d=0;d<3;d++)
	      for(dd=0;dd<3;dd++)
		slice[k][d][dd] += localpressure[i*ny*nz+j*nz+k][d][dd]/(nx*ny);
	
        for(d=0;d<3;d++)
	  for(dd=0;dd<3;dd++)
	    average[d][dd] += slice[k][d][dd]/nz;
	
        fprintf(out2D,"%2d %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
		k,
               slice[k][0][0],slice[k][0][1],slice[k][0][2],
               slice[k][1][0],slice[k][1][1],slice[k][1][2],               
               slice[k][2][0],slice[k][2][1],slice[k][2][2]);
      }
    printf("TOTAL PRESSURE TENSOR COMPONENTS: \n");
    printf("%12.6s  %12.6s  %12.6s  %12.6s  %12.6s  %12.6s  %12.6s  %12.6s  %12.6s\n",
	   "Pxx","Pxy","Pxz","Pyx","Pyy","Pyz","Pzx","Pzy","Pzz");
    printf("%12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
	   average[0][0],average[0][1],average[0][2],
	   average[1][0],average[1][1],average[1][2],
	   average[2][0],average[2][1],average[2][2]);

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
    fclose(out2D);


    /* AVERAGE OVER CYLINDRICAL ANGLE */


    bin[0]=nx;   /* x bins */
    bin[1]=ny;   /* y bins */
    bin[2]=nz;   /* z bins */ 
    
    com[0]=bin[0]/2;
    com[1]=bin[1]/2;
    com[2]=bin[2]/2;

    cil_p = (double***)malloc(bin[2]*sizeof(double));  /* matrix for clyndrical symmetry  */
    for(i=0;i<bin[2];i++)
      {
	cil_p[i] = (double**)malloc(bin[0]*sizeof(double));
	
	for(j=0;j<bin[0];j++)
	  cil_p[i][j] = (double*)malloc(DIM*sizeof(double));
	
      }
    
    for(i=0;i<nx;i++)
      for(k=0;k<nz;k++)
	for(d=0;d<DIM;d++)
	  cil_p[k][i][d]=0;
	  
    

    for(i=0;i<nx;i++)
      for(j=0;j<ny;j++)   
	for(k=0;k<nz;k++)
	  {
	    
	    dx[0]= i-com[0];
	    dx[1]= j-com[1];
	    dx[2]= 0;

	    length = norm(dx);

	    for(d=0;d<DIM;d++)
	      {   
		if(length !=0)
		  dx[d]=dx[d]/length;
		r_comp[d]=0;
	      }


	     r_comp[0] = dx[0]*dx[0]*localpressure[i*ny*nz+j*nz+k][0][0]+dx[1]*dx[1]*localpressure[i*ny*nz+j*nz+k][1][1];

	     /*off diagonal*/
	     r_comp[0] += dx[0]*dx[1]*(localpressure[i*ny*nz+j*nz+k][0][1] + localpressure[i*ny*nz+j*nz+k][1][0]) ;
	     
	     r_comp[2] = localpressure[i*ny*nz+j*nz+k][2][2];
	     
	     
	     /*angular*/
	     
	     r_comp[1]=(localpressure[i*ny*nz+j*nz+k][0][0] + localpressure[i*ny*nz+j*nz+k][1][1] +  localpressure[i*ny*nz+j*nz+k][2][2]  -  r_comp[0] - r_comp[2]);

	     bSKIP = 1;
	     /* don't go beyond periodic boundaries for radius*/
	     for(d=0;d<DIM-1;d++)
	       if(length > (bin[d] - com[d]))
		 bSKIP = 0;
	    
	     /*IN THIS VERSION AVERAGE OVER ANGLE IS TAKEN OVER SHELLS WHICH THICKNESS
	       IS HALF OF THE SIZE OF THE VOXEL. THAT IS WHY THERE IS 0.5 BELOW. 
	       NOT SURE IF THIS IS GOOD OR NOT? IF YOU WANT THE SHELLS AS THICK AS
	       VOXELS CHANGE 0.5 -> 1.0  BELOW*/

	     dbin=0.5/bin[0];	     
	     ndx=(int) floor(length/(nx*dbin));
	     
	    
	     
	     if(bSKIP == 1)
	       if(ndx >= 0 && ndx < bin[0])
		 {
		   cil_p[k][ndx][0] += r_comp[0] + r_comp[1];
		   cil_p[k][ndx][1] += r_comp[2]; 
		   cil_p[k][ndx][2]++; 
		 }
	     
	  }

    outCYL=fopen("pres_prof_cil.xvg","w");
    outRTENS=fopen("tension_radius.xvg","w");
    outKAPPAC0=fopen("kappac0_radius.xvg","w");

    double tension[nx]; 
    double tensionDOWN[nx]; 
    double tensionUP[nx];
    double kappac0[nx];
    double kappac0DOWN[nx];
    double kappac0UP[nx];
    
    for(j=0;j<bin[0];j++){
      tension[j]=0;
      tensionDOWN[j]=0;
      tensionUP[j]=0;
      kappac0[j]=0;
      kappac0DOWN[j]=0;
      kappac0UP[j]=0;
    }

    printf("\n You have now set the following bilayer location:\n center: %f \n upper boundary: %f \n lower boundary: %f \n center of upper leaflet: %f \n center of lower leaflet: %f \n", zcenter,bilayerTOP,bilayerBOT,zcenterTOP,zcenterBOT);
    printf("To change these settings look line 42 in the source code \n");

    zcenter=zcenter/scale+com[2];
    bilayerTOP=bilayerTOP/scale+com[2];
    bilayerBOT=bilayerBOT/scale+com[2];
    zcenterTOP=zcenterTOP/scale+com[2];
    zcenterBOT=zcenterBOT/scale+com[2];

    /* printf("%f",zcenter);*/

    for(i=0;i<bin[2];i++)
    {
      for(j=0;j<0.5*bin[0]/(bin[0]*dbin);j++)
	{
	  P=0;
	  if(cil_p[i][j][2]!=0)
	    P=(0.5*cil_p[i][j][0]-cil_p[i][j][1])/cil_p[i][j][2];
	  /*   P=cil_p[i][j][0]/cil_p[i][j][2];*/

	  /*IF YOU WANT TO FILTER FLUCTUATIONS FOR VISUALATION YOU CAN REMOVE THE COMMENTS BELOW*/
	  /*if(P>300)
	    P=299;
	  if(P<-700)
	  P=-699;	 */
	
	  /* SURFACE TENSION AS A FUNTION OF RADIUS */
  
	  tension[j]+=P*scale; 
	  if(i>zcenter && i<bilayerTOP)
	    tensionUP[j]+=P*scale;
	  if(i<zcenter && i>bilayerBOT)
	    tensionDOWN[j]+=P*scale;
	  
	  /* \KAPPA c_0 AS A FUNTION OF RADIUS */
	  
	  if(i<bilayerTOP && i>bilayerBOT)
	    kappac0[j]+=(i-zcenter)*P*scale*scale;
	  if(i<zcenter && i>bilayerBOT)
	    kappac0DOWN[j]+=(i-zcenterBOT)*P*scale*scale;
	  if(i<bilayerTOP && i>zcenter)
	    kappac0UP[j]+=(i-zcenterTOP)*P*scale*scale;
	  

	  fprintf(outCYL,"%f %f %f\n",(float)(scale)*(i - com[2]),(float)(scale*j*nx*dbin),(float)P);
	  
	  
	}
      
      fprintf(outCYL,"\n");
          
    }

    fclose(outCYL);

    fprintf(outRTENS,"#%s %s %s %s \n","r       ","bilayer  ","LOWERleaflet","UPPERleaflet"); 
    fprintf(outKAPPAC0,"#%s %s %s %s \n","r       ","bilayer  ","LOWERleaflet","UPPERleaflet");

    for(i=0;i<0.5*bin[0]/(bin[0]*dbin);i++){
      fprintf(outRTENS,"%f %f %f %f \n",(float)(scale*i*nx*dbin),(-1)*tension[i]*0.1,(-1)*tensionDOWN[i]*0.1,(-1)*tensionUP[i]*0.1); 
      fprintf(outKAPPAC0,"%f %f %f %f \n",(float)(scale*i*nx*dbin),kappac0[i],(-1)*kappac0DOWN[i],kappac0UP[i]); 
      
      if((float)(scale*i*nx*dbin)>cylradius){
	tottension+=tension[i]*0.1;
	sum1++;
      }
    }
    
    /* printf("# total tension %f \n",(-1)*tottension/sum1);*/

    fclose(outRTENS);
    fclose(outKAPPAC0);
    
    return 0;
    
}
