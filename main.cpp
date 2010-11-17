///for Debug
///add " SET(CMAKE_BUILD_TYPE Debug) " to CMakeLists.txt
///sudo jaechoon@ubuntu:~$ chmod ug+rw /dev/tty*

#include "global.h"
#include "geopano.h"
#include <unistd.h>     // Header File for sleeping.
#include "v3dpano.h"


/* ASCII code for the escape key. */
#define ESCAPE 27
#define PAGE_UP 73
#define PAGE_DOWN 81
#define HOME 71
#define END 79
#define UP_ARROW 72
#define DOWN_ARROW 80
#define LEFT_ARROW 75
#define RIGHT_ARROW 77

gpano_t gpano; 
void displayAll(GLuint list);
void displayNew(GLuint list);

GLfloat xrot=45.0f,yrot=0.0f,zrot=0.0f;
GLfloat xpos=0.0f, ypos=0.0f, zpos=2.0f;
int mouse[2]={0,0};
GLuint glList;
/* The number of our GLUT window */
int window; 
int curButton=0;


/* A general OpenGL initialization function.  Sets all of the initial parameters. */
void InitGL(int Width, int Height)	        // We call this right after our OpenGL window is created.
{
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);		// This Will Clear The Background Color To Black
  glClearDepth(1.0);				// Enables Clearing Of The Depth Buffer
  glDepthFunc(GL_LESS);			        // The Type Of Depth Test To Do
  glEnable(GL_DEPTH_TEST);		        // Enables Depth Testing
  glShadeModel(GL_SMOOTH);			// Enables Smooth Color Shading

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();				// Reset The Projection Matrix

  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,1.0f,1000.0f);	// Calculate The Aspect Ratio Of The Window

  glMatrixMode(GL_MODELVIEW);
}

/* The function called when our window is resized (which shouldn't happen, because we're fullscreen) */
void ReSizeGLScene(int Width, int Height)
{
  if (Height==0)				// Prevent A Divide By Zero If The Window Is Too Small
    Height=1;

  glViewport(0, 0, Width, Height);		// Reset The Current Viewport And Perspective Transformation

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,1.0f,1000.0f);
  glMatrixMode(GL_MODELVIEW);
}

/* The main drawing function. */
void DrawGLScene()
{    	
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		// Clear The Screen And The Depth Buffer
    glLoadIdentity();

    glRotatef(xrot, 1.0f,0,0);
    glRotatef(yrot, 0, 1.0f, 0);
    glRotatef(zrot, 0, 0, 1.0f);
    glTranslatef(xpos, ypos, zpos); 
    
    glCallList(glList);

  // swap the buffers to display, since double buffering is used.
  glutSwapBuffers();
}

void mouseFunc(int button, int state, int x, int y)
{    
	mouse[0]=x;
	mouse[1]=y;
	curButton=button;
}

void processMouseActiveMotion(int x, int y)
{
 // if(curButton!= GLUT_DOWN )return;
  
  if (curButton == GLUT_RIGHT_BUTTON) {
	yrot+=(y-mouse[1])/2;
	xrot+=(x-mouse[0])/2;
	
	//printf("rot(%f,%f)\n",xrot,zrot);
  }  
  
  if (curButton == GLUT_LEFT_BUTTON ) {
	xpos+=(x-mouse[0])/2;
	ypos+=(y-mouse[1])/2;
	
	//printf("pos(%f,%f)\n",xpos,ypos);
  }

  mouse[0]=x;
  mouse[1]=y;
 
}

void mouseWheel(int button, int dir, int x, int y)
{

    if (dir > 0)
    {
      if(curButton==GLUT_MIDDLE_BUTTON)zrot+=5.0f;
      else zpos+=10.0f;
       
    }
    else
    {
       if(curButton==GLUT_MIDDLE_BUTTON)zrot-=5.0f;
      else zpos-=10.0f;
    }

    DrawGLScene();
    return;
}


/* The function called whenever a key is pressed. */
void keyPressed(unsigned char key, int x, int y) 
{
   /* avoid thrashing this procedure */
    usleep(10);

    switch (key) {    
    case ESCAPE: // kill everything.
	/* exit the program...normal termination. */
	exit(1);                   	
	break; // redundant.

    case 'b': 
     case 'f': 
   
    default:
      printf ("Key %d pressed. No action there yet.\n", key);
      break;
    }	
}


/* The function called whenever a normal key is pressed. */
void specialKeyPressed(int key, int x, int y) 
{
    /* avoid thrashing this procedure */
    usleep(10);

    switch (key) { 
      
    case GLUT_KEY_HOME: // tilt up
	ypos+=1.0f;
	break;
    
    case GLUT_KEY_END: // tilt down
        ypos-=1.0f;
	break;
	
    case GLUT_KEY_PAGE_UP: // tilt up
	zpos+=1.0f;
	break;
    
    case GLUT_KEY_PAGE_DOWN: // tilt down
        zpos-=1.0f;
	break;

    case GLUT_KEY_UP: // walk forward (bob head)
	xrot+=3.0f;
	break;

    case GLUT_KEY_DOWN: // walk back (bob head)
	xrot-=3.0f;
	break;

    case GLUT_KEY_LEFT: // look left
	yrot += 3.0f;
	break;
    
    case GLUT_KEY_RIGHT: // look right
	yrot -= 3.0f;
	break;

    default:
	printf ("Special key %d pressed. No action there yet.\n", key);
	break;
    }	
}

int main(int argc, char **argv) 
{  
  gpano.statesN=0;
  gpano.raNewPano=NULL;
  gpano.imNewPano=NULL;
  gpano.imPanoTarget=NULL;
  
  int i;
  for(i=0; i<GPANO_MAX_FRAMES ; i++)
  { 
    gpano.imPano[i]=NULL;
    gpano.raPano[i]=NULL;
    gpano.plPano[i]=NULL;
    gpano.plane[i]=NULL;
  }
  
  
  /////////////////////////////////////////////////////////
  ///  read meta data and image/range panos
  /////////////////////////////////////////////////////////
  int nNeiPano=2;
  string dir[50], filename;
  dir[0]="/home/jaechoon/Desktop/VPano/Data/1/";
  dir[1]="/home/jaechoon/Desktop/VPano/Data/2/";
  dir[2]="/home/jaechoon/Desktop/pipline/test_data/set2/1000004009194/";
  
  CGeopano gp;
  for(i=0; i<nNeiPano; i++)
  {
    gp.read_state(&gpano,dir[i]+"transform.txt");
    gp.read_plane(&gpano,dir[i]+"plane_palette.txt");
  
    filename=dir[i]+"raster.jpg";
    gpano.imPano[i]= cvLoadImage(filename.c_str(),CV_LOAD_IMAGE_UNCHANGED);
    if(!gpano.imPano[i])
    {
      printf("Could not load image pano file: %s\n",filename.c_str());
      exit(0);
    }
        
    filename=dir[i]+"depth.png";
    gpano.raPano[i]= cvLoadImage(filename.c_str(),CV_LOAD_IMAGE_UNCHANGED);
    if(!gpano.raPano[i])
    {
      printf("Could not load range pano file: %s\n",filename.c_str());
      exit(0);
    }
    
    //OpenCV does not support alpha channel, only masking. 
    //If you want to read in PNG with alpha, use imagemagick first to extract alpha channel:
    //>convert plane_pano.png -channel Alpha -negate -separate plane_pano1.png

    chdir(dir[i].c_str());
    system("convert plane_pano.png -channel Alpha -negate -separate plane_alpha.png");
    
    filename=dir[i]+string("plane_alpha.png");
    gpano.plPano[i]=cvLoadImage(filename.c_str(), CV_LOAD_IMAGE_UNCHANGED);

    if(!gpano.plPano[i])
    {
      printf("Could not load plane pano file: %s\n",filename.c_str());
      exit(0);
    }
    system("rm plane_alpha.png");
  }


  ////////////////////////////////////////////////////////////////////
  /// read new pano's meta data and image data
  ///////////////////////////////////////////////////////////////////
  gp.set_poseNew(&gpano,"/home/jaechoon/Desktop/transform3.txt");
  filename="/home/jaechoon/Desktop/raster3.jpg";
  gpano.imPanoTarget= cvLoadImage(filename.c_str(),CV_LOAD_IMAGE_UNCHANGED);
  if(!gpano.imPanoTarget)
  {
      printf("Could not load Targrt image pano file: %s\n",filename.c_str());
      gpano.imPanoTarget=NULL;
  }
  ///////////////////////////////////////////////////////
  /// build new pano 
  //////////////////////////////////////////////////////
  //v3Dpano v3(&gpano);
  v3Dpano v3(gpano.imPano[0]->width,gpano.imPano[0]->height, gpano.imPano[0]->nChannels,
	     gpano.raPano[0]->width,gpano.raPano[0]->height, gpano.raPano[0]->nChannels,
	     gpano.plPano[0]->width,gpano.plPano[0]->height, gpano.plPano[0]->nChannels,
	     gpano.statesN);
  for(i=0; i<gpano.statesN; i++)
  {
    v3.creatPanoBySingle(gpano,i);
    v3.fillPanoHoleBySingle(gpano,i);
  }
  ///////////////////////////////////////////////////////
  /// save results
  //////////////////////////////////////////////////////
  v3.newRaIm2PTS(gpano, "/home/jaechoon/Desktop/pts3.txt");
  if(!cvSaveImage("/home/jaechoon/Desktop/imNewPano.jpg",gpano.imNewPano))printf("Could not save: %s\n","imNewPano.jpg");
  if(!cvSaveImage("/home/jaechoon/Desktop/raNewPano.png",gpano.raNewPano))printf("Could not save: %s\n","raNewPano.png");
  
  
  
  //cvNamedWindow("Image_New_Pano",0);//1-MAxmized window,0-minimized window
  //cvShowImage("Image_New_Pano",gpano.imNewPano);
  
  //cvNamedWindow("Range_New_Pano",0);//1-MAxmized window,0-minimized window
  //cvShowImage("Range_New_Pano",gpano.raNewPano);
  //cvWaitKey(0);

  //return 0;
  glutInit(&argc, argv);  
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);  
  glutInitWindowSize(640, 480);  
  glutInitWindowPosition(0, 0);  
  window = glutCreateWindow("Virtual 3D Pano Generation");  
  
  glList=glGenLists(1);
  //displayAll(glList);
  displayNew(glList);
  
  glutDisplayFunc(&DrawGLScene);  
  //glutFullScreen();
  glutIdleFunc(&DrawGLScene);
  glutReshapeFunc(&ReSizeGLScene);
  glutMouseFunc(&mouseFunc);
  glutMotionFunc(processMouseActiveMotion);
  glutMouseWheelFunc(&mouseWheel);
  glutKeyboardFunc(&keyPressed);  
  glutSpecialFunc(&specialKeyPressed);  InitGL(640, 480);
   
  glutMainLoop();  

  for(i=0;i<gpano.statesN;i++)
  {
    cvReleaseImage( &gpano.imPano[i]);
    cvReleaseImage( &gpano.raPano[i]);
    cvReleaseImage( &gpano.plPano[i]);
    
    delete []gpano.plane[i];
  }
  cvReleaseImage( &gpano.imNewPano);
  cvReleaseImage( &gpano.raNewPano);
  cvReleaseImage( &gpano.imPanoTarget);
    
  cvDestroyWindow("Image_New_Pano");
  cvDestroyWindow("Range_New_Pano");
  return 0;

}

void displayAll(GLuint	list)
{
   glNewList(list,GL_COMPILE);
   glEnable(GL_DEPTH_TEST);
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glBegin(GL_POINTS);

   CGeopano gp;
   	
   int i;
  for(i=0; i<gpano.statesN; i++)
  {
    uchar  *im=(uchar *)gpano.imPano[i]->imageData;
    int stepIm = gpano.imPano[i]->widthStep/sizeof(uchar);
    int bandIm=gpano.imPano[i]->nChannels;
    int wIm=gpano.imPano[i]->width;
    int hIm=gpano.imPano[i]->height;
    
        
    unsigned short *ra=(unsigned short*)gpano.raPano[i]->imageData;
    int stepRa = gpano.raPano[i]->widthStep/sizeof(uchar)/2;
    int bandRa=gpano.raPano[i]->nChannels;
    int w=gpano.raPano[i]->width;
    int h=gpano.raPano[i]->height;
    
    int u,v,x,y,po,n=0;
    double r, yaw, pitch, cos_pitch, sin_pitch;
    vec3_t P0;
    
    printf("imagePano channal(%d) w(%d) h(%d) wStep(%d)\n",bandIm,gpano.imPano[i]->width,gpano.imPano[i]->height,stepIm);
    printf("rangePano channal(%d) w(%d) h(%d) wStep(%d) depth(%d)\n",bandRa,w,h,stepRa,gpano.raPano[i]->depth);
   
    if(i==1)glColor3f(1.0f,0,0);
    else if(i==2)glColor3f(0,1.0f,0);
    else glColor3f(0,0,1.0f);		 
    
    for(y=0; y<h; y+=10)
    {
      pitch = ( - M_PI * y /(double)(h) ) + M_PI / 2.0; // + 90 degrees to - 90 degrees
      cos_pitch=cos ( pitch );
      sin_pitch=sin ( pitch );
      
      v=hIm*y/h;
      
      for(x=0; x<w; x+=10)
      {
	r=(double)(ra[x*bandRa+y*stepRa])*0.01;
	if(r<=0 || r>200)continue;
      
	u=wIm*x/w;
	
	yaw = ( 2.0 * M_PI * x /(double)( w) ) - M_PI;   // - 180 degrees to + 180 degrees
        yaw = -yaw + M_PI / 2.0;

	P0.x = r * cos_pitch * cos(yaw);
	P0.y = r * cos_pitch * sin(yaw);
	P0.z = r * sin_pitch;
	
	//vec3_t P1=gp.getGlobalPostion(gpano.states[i].pose,P0);
	glColor3b(im[u*bandIm+2+v*stepIm],im[u*bandIm+1+v*stepIm],im[u*bandIm+v*stepIm]);
	glVertex3f(P0.x,P0.y,P0.z);
      
    }
    //printf("\n",r);
    }    
  }  
  
  glEnd();
  glEndList();
}




void displayNew(GLuint	list)
{
   glNewList(list,GL_COMPILE);
   glEnable(GL_DEPTH_TEST);
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glBegin(GL_POINTS);

   CGeopano gp;
   
    IplImage *imIm=gpano.imNewPano;
    if(gpano.imPanoTarget!=NULL)imIm=gpano.imPanoTarget;
   
    uchar  *im=(uchar *)imIm->imageData;
    int stepIm = imIm->widthStep/sizeof(uchar);
    int bandIm=imIm->nChannels;
    int wIm=imIm->width;
    int hIm=imIm->height;
    
        
    unsigned short *ra=(unsigned short*)gpano.raNewPano->imageData;
    int stepRa = gpano.raNewPano->widthStep/sizeof(uchar)/2;
    int bandRa=gpano.raNewPano->nChannels;
    int w=gpano.raNewPano->width;
    int h=gpano.raNewPano->height;
    
    int u,v,x,y,po,n=0;
    double r, yaw, yaw1, pitch, cos_pitch, sin_pitch;
    vec3_t P0;
    
  //  printf("imagePano channal(%d) w(%d) h(%d) wStep(%d)\n",bandIm,gpano.imNewPano->width,gpano.imNewPano->height,stepIm);
  //  printf("rangePano channal(%d) w(%d) h(%d) wStep(%d) depth(%d)\n",bandRa,w,h,stepRa,gpano.raNewPano->depth);
   
    glColor3f(1.0f,0,0);
    
    for(y=0; y<h; y+=2)
    {
      pitch = ( - M_PI * y /(double)(h) ) + M_PI / 2.0; // + 90 degrees to - 90 degrees
      cos_pitch=cos ( pitch );
      sin_pitch=sin ( pitch );
      
      v=hIm*y/h;
      
     
      for(x=0; x<w; x+=2)
      {
	r=(double)(ra[x*bandRa+y*stepRa])*0.01;
	if(r<=0 || r>655)continue;
      
	u=wIm*x/w;
	
	yaw = ( 2.0 * M_PI * x /(double)( w) ) - M_PI;   // - 180 degrees to + 180 degrees
        yaw1 = -yaw + M_PI / 2.0;

	P0.x = r * cos_pitch * cos(yaw1);
	P0.y = r * cos_pitch * sin(yaw1);
	P0.z = r * sin_pitch;
	
	//vec3_t P1=gp.getGlobalPostion(gpano.poseNew,P0);
	//if(P0.z<-1)continue;
	glColor3b(im[u*bandIm+2+v*stepIm],im[u*bandIm+1+v*stepIm],im[u*bandIm+v*stepIm]);
	glVertex3f(P0.x,P0.y,P0.z);
      
    }
    //printf("\n",r);
    }    
  glEnd();
  glEndList();
}





