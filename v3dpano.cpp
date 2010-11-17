#include "global.h"
#include "v3dpano.h"
#include "geopano.h"
#include "bestfit.h"

/*
v3Dpano::v3Dpano(gpano_t *self)
{
  m_bandIm=self->imPano[0]->nChannels;
  m_wIm= self->imPano[0]->width;
  m_hIm=self->imPano[0]->height;
  m_stepIm=self->imPano[0]->widthStep;

  m_bandRa=self->raPano[0]->nChannels;
  m_wRa=self->raPano[0]->width;
  m_hRa=self->raPano[0]->height;
  m_stepRa=self->raPano[0]->widthStep;
  
  m_wPl=self->plPano[0]->width;
  m_hPl=self->plPano[0]->height;
  m_bandPl=self->plPano[0]->nChannels;
  m_stepPl=self->plPano[0]->widthStep;  
  
  
  m_curRaRate=new double[m_wRa*m_hRa];
  memset(m_curRaRate,0,m_wRa*m_hRa*sizeof(double));
  
  m_curImSrc=new uchar[m_wRa*m_hRa];
  memset(m_curImSrc,0xFF,m_wRa*m_hRa);

  m_nPano=self->statesN;
*/
  v3Dpano::v3Dpano(int wIm, int hIm, int bIm, int wRa, int hRa, int bRa, int wPl, int hPl, int bPl, int nPano)
{
  m_bandIm=bIm;
  m_wIm=wIm;
  m_hIm=hIm;
  m_stepIm=wIm*bIm;

  m_bandRa=bRa;
  m_wRa=wRa;
  m_hRa=hRa;
  m_stepRa=wRa*bRa;
  
  m_wPl=wPl;
  m_hPl=hPl;
  m_bandPl=bPl;
  m_stepPl=wPl*m_bandPl;  
  
  
  m_curRaRate=new double[wRa*hRa];
  memset(m_curRaRate,0,wRa*hRa*sizeof(double));
  
  m_curImSrc=new uchar[wRa*hRa];
  memset(m_curImSrc,0xFF,wRa*hRa);
  
   m_nPano=nPano;
  if(m_nPano>10)m_nPano=10;
  for(int i=0; i<m_nPano; i++)
  {
    m_imProj[i]=new int[m_wIm*m_hIm];
    memset(m_imProj[i],0,sizeof(int)*m_wIm*m_hIm);
  }
  
}
v3Dpano::~v3Dpano()
{
  delete []m_curRaRate;
  delete []m_curImSrc;
  for(int i=0; i<m_nPano; i++)delete []m_imProj[i];  

}

bool v3Dpano::creatPanoBySingle(gpano_t &self, int i)
{
   if(self.statesN<=0 || i>=m_nPano)return false;
  
   if(self.raNewPano==NULL && self.statesN>0)
   {
      self.raNewPano=cvCloneImage(self.raPano[0]);
      self.imNewPano=cvCloneImage(self.imPano[0]);
      
      m_raNew=(ushort*)self.raNewPano->imageData;
      memset(m_raNew,0xFF, sizeof(ushort)*m_stepRa*m_hRa);
   
      m_imNew=(uchar *)self.imNewPano->imageData;
      memset(m_imNew,0xFF, sizeof(uchar)*m_stepIm*m_hIm);
   }
   
   CGeopano gp;
   int u,v,po;
   double x,y,xpl, ypl, x1,x2,y1,y2;
   double r,r1,r2, yaw, yaw1,yaw11, yaw2, pitch, pitch2, cos_pitch, sin_pitch, N[4];
   double angle;
   vec3_t vec0,vec1,vec2, P0,P2;
   int xd, yd, u2, v2;
   int poRate, poRa, poIm, poImSrc;
   ushort tmp;
   
   pose3_t pose=pose3_div(self.poseNew, self.states[i].pose);
   
   uchar  *im=(uchar *)self.imPano[i]->imageData;
   ushort *ra=(ushort*)self.raPano[i]->imageData;
   uchar  *pl=(uchar *)self.plPano[i]->imageData;
   double *plane=self.plane[i];
   int planeN=self.planeN[i];
   
   for(v=0; v<m_hIm; v++)
   {
      y=(double)(v*m_hRa)/(double)(m_hIm);
      ypl=(double)(v*m_hPl)/(double)(m_hIm);
      pitch = ( - M_PI * y /(double)(m_hRa) ) + M_PI / 2.0; // + 90 degrees to - 90 degrees
      cos_pitch=cos ( pitch );
      sin_pitch=sin ( pitch );     
      
      for(u=0; u<m_wIm; u++)
      {
	x=(double)(u*m_wRa)/(double)(m_wIm);

	tmp=ra[(int)(x+0.5)*m_bandRa+(int)(y+0.5)*m_stepRa];
	if(tmp<100 || tmp>65500)continue;
	 
	yaw = ( 2.0 * M_PI * x /(double)( m_wRa) ) - M_PI;   // - 180 degrees to + 180 degrees
        yaw1 = -yaw + M_PI / 2.0;

	P0 =vec3_set( tmp * cos_pitch * cos(yaw1),
		      tmp * cos_pitch * sin(yaw1),
		      tmp * sin_pitch);
	
	/// coordinate system change from neighbor PaAno to newPano
	P0.x*=0.01;	P0.y*=0.01;	P0.z*=0.01;///[cm]->[m]	
	//P2=gp.getGlobalPostion(self.states[i].pose,P0);
	P2=gp.getGlobalPostion(pose,P0);
	P2.x*=100.0;	P2.y*=100.0;	P2.z*=100.0;///[m]->[cm]	
	
	r2=vec3_mag(P2);
	if(r2<100.0)continue;
	
	///occlusion test using normal vector informatio
	xpl=(double)(u*m_wPl)/(double)(m_wIm);
	int planeIndex=255-(int)(pl[(int)(xpl+0.5)*m_bandPl+(int)(ypl+0.5)*m_stepPl]);
	if(planeIndex>=planeN || planeIndex==0)continue;
	
	vec0.x=plane[4*planeIndex];
	vec0.y=plane[4*planeIndex+1];
	vec0.z=plane[4*planeIndex+2];
	
	vec1=vec3_rotate(pose.rot, vec0);
	vec2=vec3_unit(P2);
	angle=vec3_dot(vec1, vec2);
	if(angle>=0.0)continue;
	
	///project on the new pano
	pitch2=asin(P2.z/r2);
	yaw11=atan2(P2.y/r2,P2.x/r2);
	yaw2=M_PI/2.0-yaw11;
	x2=(yaw2+M_PI)*m_wRa/(2.0*M_PI);
	y2=-(pitch2-M_PI/2.0)*(double)(m_hRa)/M_PI;
	
	xd=(int)(x2+0.5)%m_wRa;
	yd=(int)(y2+0.5)%m_hRa;
		
	if(xd<0 || xd>=m_wRa || yd<0 || yd>=m_hRa)continue;
	
	tmp=(ushort)(r2);
	u2=(int)(m_wIm*x2/(double)(m_wRa)+0.5);
	v2=(int)(m_hIm*y2/(double)(m_hRa)+0.5);
	
	poRate=xd       +yd*m_wRa;
	poRa  =xd*m_bandRa+yd*m_stepRa;
	poIm   =m_bandIm*u2+m_stepIm*v2;
	poImSrc=m_bandIm*u +m_stepIm*v;
	
	if(m_curRaRate[poRate]==0)
	{
	  m_raNew[poRa]=tmp;
	  m_curRaRate[poRate]=r;
	  m_curImSrc[u2+m_wIm*v2]=i;
	  memcpy(&m_imNew[poIm], &im[poImSrc],m_bandIm);
	  
	  m_imProj[i][u+m_wIm*v]=u2+m_wIm*v2;
	}
	else //already exist
	{
	  if(m_raNew[poRa]>tmp+0.05*m_raNew[poRa])
	  {
	    m_raNew[poRa]=tmp;
	    m_curRaRate[poRate]=r;
	    m_curImSrc[u2+m_wIm*v2]=i;
	    memcpy(&m_imNew[poIm], &im[poImSrc],m_bandIm);
	  
	    m_imProj[i][u+m_wIm*v]=u2+m_wIm*v2;
	  }
	  else if(m_raNew[poRa]>tmp-0.05*m_raNew[poRa] && m_curRaRate[poRate]>r)
	  {
	    m_raNew[poRa]=tmp;
	    m_curRaRate[poRate]=r;
	    m_curImSrc[u2+m_wIm*v2]=i;
	    memcpy(&m_imNew[poIm], &im[poImSrc],m_bandIm);
	  
	    m_imProj[i][u+m_wIm*v]=u2+m_wIm*v2;
	  }
	}
     }
   }
    
  return true;
}

bool v3Dpano::fillPanoHoleBySingle(gpano_t &self, int i)
{
  if(self.statesN<=0 || i>=m_nPano || self.raNewPano==NULL)return false;
  
   int u,v,u0,v0,u1,v1,u2,v2, u3,v3,po, po0, po1, po2;
   int x0,y0,x,y,x1,x2,x3,y1,y2,y3;
   double r,r0,r1,r2,r3;
   uchar color1[6], color2[6], color3[6];
  
   int du[3]={-1, 0, -1};
   int dv[3]={-1,-1,  0};
   int k;
   printf("...filling %d pano's holes...\n",i+1);
   
   uchar  *imNew=(uchar *)self.imNewPano->imageData;
   ushort *raNew=(ushort*)self.raNewPano->imageData;   
   ushort *ra=(ushort*)self.raPano[i]->imageData;
      
   for(v=2; v<m_hIm-2; v++)
   for(u=2; u<m_wIm-2; u++)
   {
     po=u +m_wIm*v;
     if(m_imProj[i][po]==0 || m_curImSrc[m_imProj[i][po]]!=i)continue;
     
     x=(int)((double)(u*m_wRa)/(double)(m_wIm)+0.5);
     y=(int)((double)(v*m_hRa)/(double)(m_hIm)+0.5);
     r=(double)(raNew[x*m_bandRa+y*m_stepRa]);
      
     u1=m_imProj[i][po]%m_wIm;
     v1=m_imProj[i][po]/m_wIm;
     memcpy(color1,&imNew[u1*m_bandIm+v1*m_stepIm],m_bandIm);
     
     x1=(int)((double)(u1*m_wRa)/(double)(m_wIm)+0.5);
     y1=(int)((double)(v1*m_hRa)/(double)(m_hIm)+0.5);
     r1=(double)(raNew[x1*m_bandRa+y1*m_stepRa]);
 	
     for( k=0; k<3; k++)
     {
	u0=u+du[k];
	v0=v+dv[k];
	po0=u0+m_wIm*v0;
	if(m_imProj[i][po0]==0 || m_curImSrc[m_imProj[i][po0]]!=i)continue;
	
	x0=(int)((double)(u0*m_wRa)/(double)(m_wIm)+0.5);
	y0=(int)((double)(v0*m_hRa)/(double)(m_hIm)+0.5);
	r0=(double)(raNew[x0*m_bandRa+y0*m_stepRa]);
     
     	u2=m_imProj[i][po0]%m_wIm;
	v2=m_imProj[i][po0]/m_wIm;
	memcpy(color2,&imNew[u2*m_bandIm+v2*m_stepIm],m_bandIm);
	
	x2=(int)((double)(u1*m_wRa)/(double)(m_wIm)+0.5);
        y2=(int)((double)(v1*m_hRa)/(double)(m_hIm)+0.5);
	r2=(double)(raNew[x2*m_bandRa+y2*m_stepRa]);
	  
	if(fabs(r1-r2)>0.1*r1)continue;// in range space 1[m]->0.1[m]. 10[m]->1[m], 100[m]->10[m]
	  
	int dVu=u2-u1;
	int dVv=v2-v1;
	double dL=sqrt(dVu*dVu+dVv*dVv);
	if(dL<2.0 || dL>100.0)continue; //in image space
	  
	double dVu1=(double)(dVu)/dL;
	double dVv1=(double)(dVv)/dL;  
	
	  
	/// interpolation between po1 and po2
	for(double j=1.0;j<dL;j+=1)
	{
	  u3=(int)(j*dVu1+0.5)+u1;
	  v3=(int)(j*dVv1+0.5)+v1;
	  for(int ii=0;ii<m_bandIm; ii++) color3[ii]=(uchar)(j/dL*(color2[ii]-color1[ii])+0.5)+color1[ii];
	  
	  x3=(int)((double)(u3*m_wRa)/(double)(m_wIm)+0.5);
	  y3=(int)((double)(v3*m_hRa)/(double)(m_hIm)+0.5);
	  r3=(ushort)( j/dL*(r2-r1)+r1+0.5);
	  
	  int poRate=x3         +y3*m_wRa;
	  int poRa  =x3*m_bandRa+y3*m_stepRa;
	  int poIm   =m_bandIm*u3+m_stepIm*v3;
	  int poImSrc=m_bandIm*u +m_stepIm*v;
	
	  if(m_curRaRate[poRate]==0)
	  {
	    raNew[poRa]=r3;
	    m_curRaRate[poRate]=(r+r0)/2.0;
	    m_imProj[i][poImSrc]=poIm;
	    memcpy(&imNew[poIm], color3,m_bandIm);
	  }
	else //already exist
	{
	  if(m_raNew[poRa]>r3+0.05*m_raNew[poRa])
	  {
	    raNew[poRa]=r3;
	    m_curRaRate[poRate]=(r+r0)/2.0;
	    m_imProj[i][poImSrc]=poIm;
	    memcpy(&imNew[poIm], color3,m_bandIm);
	  }
	  else if(m_raNew[poRa]>r3-0.05*m_raNew[poRa] && m_curRaRate[poRate]>(r+r0)/2.0)
	  {
	    raNew[poRa]=r3;
	    m_curRaRate[poRate]=(r+r0)/2.0;
	    m_imProj[i][poImSrc]=poIm;
	    memcpy(&imNew[poIm], color3,m_bandIm);
	  }
	}
     }
   }
}
}

bool v3Dpano::calNormalVector(ushort *ra, int x, int y, vec3_t &P,vec3_t &v, ushort &r)
{
    r=ra[x*m_bandRa+y*m_stepRa];
    
    if(r<100 || r>65500)return false;
    
    /// pick kernel size based on depth
     int kernel_size = 0;
     if(r < 700)kernel_size = 4;
     else if (r < 1500)kernel_size = 3;
     else if (r < 3000)kernel_size = 2;
     else  kernel_size = 1;
     
     double r0, pitch,cos_pitch,sin_pitch,yaw,yaw1, X,Y,Z;
     int pointOffset,numpoints = 0;
     double points[3*81];
     
     for (int kx = x - kernel_size; kx <= x + kernel_size; kx++)
     {
         if (kx < 0 || kx >= m_wRa) continue;
                
         for (int ky = y - kernel_size; ky <= y + kernel_size; ky++)
         {
             if (ky < 0 || ky >= m_hRa) continue;
	     
	     pitch = ( - M_PI * ky /(double)(m_hRa) ) + M_PI / 2.0; // + 90 degrees to - 90 degrees
	     cos_pitch=cos ( pitch );
	     sin_pitch=sin ( pitch );     
	     yaw = ( 2.0 * M_PI * kx /(double)( m_wRa) ) - M_PI;   // - 180 degrees to + 180 degrees
	     yaw1 = -yaw + M_PI / 2.0;

	     r0=ra[kx*m_bandRa+ky*m_stepRa];
	     X=r0 * cos_pitch * cos(yaw1);
	     Y=r0 * cos_pitch * sin(yaw1);
	     Z=r0 * sin_pitch;
			     
             pointOffset = numpoints * 3;
             points[pointOffset]     = X;
             points[pointOffset + 1] = Y;
             points[pointOffset + 2] = Z;
             numpoints++;
          }
      }
      
      if (numpoints < 3) return false;
      
      double eigs[3], N[4];
      getBestFitPlane(numpoints,           // number of input data points
                      points,              // starting address of points array.
                      sizeof(double) * 3,  // stride between input points.
                      NULL,                // *optional point weighting values.
                      0,                   // weight stride for each vertex.
                      N,               // the output best fit plane
                      eigs);

		          
    pitch = ( - M_PI * y /(double)(m_hRa) ) + M_PI / 2.0; // + 90 degrees to - 90 degrees
    cos_pitch=cos ( pitch );
    sin_pitch=sin ( pitch );     
    yaw = ( 2.0 * M_PI * x /(double)( m_wRa) ) - M_PI;   // - 180 degrees to + 180 degrees
    yaw1 = -yaw + M_PI / 2.0;
    
    P.x= cos_pitch * cos(yaw1);
    P.y= cos_pitch * sin(yaw1);
    P.z= sin_pitch;
	// if the plane normal is pointing away from the IMU, flip the sign
    double imu_dot_product = P.x * N[0] + P.y * N[1] + P.z * N[2];
    if (imu_dot_product > 0.0)
    {
             N[0] = -N[0];
             N[1] = -N[1];
             N[2] = -N[2];
             N[3] = -N[3];
    }
    
    v.x=N[0]; v.y=N[1]; v.z=N[2];
    P.x*=r; P.y*=r; P.z*=r;
    
    return true;
}


bool v3Dpano::newRaIm2PTS(gpano_t &self, string filename)
{
  FILE *file = fopen(filename.c_str(), "w");
  if (!file)
   return GPANO_ERROR("unable to open %s : %s", filename.c_str(), strerror(errno));
  
   pose3_t pose=self.poseNew;
   IplImage *raIm=self.raNewPano;
   IplImage *imIm=self.imNewPano;
   if(self.imPanoTarget!=NULL)imIm=self.imPanoTarget;
 
   uchar *im=(uchar*)imIm->imageData;
   int wIm=imIm->width;
   int hIm=imIm->height;
   int bIm=imIm->nChannels;
   int stepIm=imIm->widthStep/sizeof(uchar);
  
   ushort *ra=(ushort*)raIm->imageData;
   int wRa=raIm->width;
   int hRa=raIm->height;
   int bRa=raIm->nChannels;
   int stepRa=raIm->widthStep/sizeof(uchar);
 
   CGeopano gp;
   
   int x,y, u,v;
   double cos_pitch, sin_pitch, pitch, yaw, yaw1, tmp;
   vec3_t P0, P;
   
   for(y=0; y<hRa; y++)
   {
      v=(int)((y*hIm)/(double)(hRa)+0.5);
      if(v<0 || v>=hIm)continue;
      
      pitch = ( - M_PI * y /(double)(hRa) ) + M_PI / 2.0; // + 90 degrees to - 90 degrees
      cos_pitch=cos ( pitch );
      sin_pitch=sin ( pitch );     
      
      
      for(x=0; x<wRa; x++)
      {
	u=(int)((x*wIm)/(double)(wRa)+0.5);
	if(u<0 || u>=wIm)continue;
	
	tmp=(double)(ra[x*bRa+y*m_stepRa]);
	if(tmp<100.0 || tmp>65500.0)continue;
	tmp/=100.0;
	
	yaw = ( 2.0 * M_PI * x /(double)( m_wRa) ) - M_PI;   // - 180 degrees to + 180 degrees
        yaw1 = -yaw + M_PI / 2.0;

	P0 =vec3_set( tmp * cos_pitch * cos(yaw1),
		      tmp * cos_pitch * sin(yaw1),
		      tmp * sin_pitch);
		      
        P=gp.getGlobalPostion(pose,P0);
		      
	if(vec3_mag(P)<200)
	fprintf(file, "%.4f %.4f %.4f 0 %d %d %d\n",
		P.x+self.utm_offset.x,P.y+self.utm_offset.y,P.z+self.utm_offset.z,
		(int)(im[bIm*u+2+v*stepIm]),(int)(im[bIm*u+1+v*stepIm]),(int)(im[bIm*u+v*stepIm]));
      }
   }
   fclose(file);
   return true;
}
   

