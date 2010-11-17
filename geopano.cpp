#include "geopano.h"



static int read_scanstr(
    FILE *fp,       /* input file pointer */
    const char *str)    /* input target string */
{
    const char *s;
    int c;

    s = str;
    for (;;) {
    c = getc(fp);
    if (c == EOF)
        return FAILURE;
    if (c == *s) {
        s++;
        if (*s == '\0')
        return SUCCESS;
        }
    else
        s = str;
    }
}

// Read state data
int CGeopano::read_state(gpano_t *self, string filename1)
{
  const char *filename=filename1.c_str();
   FILE *file;
   int i;
   gpano_state_t *state;
   pose3_t pose;
  
   // Initialize state data
   if(self->statesN<=0)
   for (i = 0; i < GPANO_MAX_FRAMES; i++)
   {
     self->statesN=0; 
     state = &self->states[i];
     state->valid = false;
     state->pose = pose3_ident();
    }

    // Prepend the master path
    GPANO_MSG("loading %s", filename);

    file = fopen(filename, "r");
    if (!file)
      return GPANO_ERROR("unable to open %s : %s", filename, strerror(errno));

    double alt, utm_e, utm_n;
    double rx, ry, rz;
    scanstr_(file, "yaw=");
    scan_(file, 1, (file, "%lf \n", &rz));
    scanstr_(file, "pitch=");
    scan_(file, 1, (file, "%lf \n", &rx));
    scanstr_(file, "roll=");
    scan_(file, 1, (file, "%lf \n", &ry));
    scanstr_(file, "x=");
    scan_(file, 1, (file, "%lf \n", &utm_e));
    scanstr_(file, "y=");
    scan_(file, 1, (file, "%lf \n", &utm_n));
    scanstr_(file, "z=");
    scan_(file, 1, (file, "%lf \n", &alt));
    
    fclose(file);
    
    pose.pos = vec3_set(utm_e, utm_n, alt);
    //pose.rot = quat_from_rpy(ry * M_PI/180, rx * M_PI/180, -rz * M_PI/180);
    pose.rot = quat_from_rpy(rx * M_PI/180, ry * M_PI/180, -rz * M_PI/180);

    if(self->statesN==0)
      self->utm_offset = pose.pos;
    
    // Make sure the frame id is sane
    if ( self->statesN >= GPANO_MAX_FRAMES)
    {
      GPANO_MSG("frame %d is out-of-bounds", self->statesN);
      return -1;
    }
    
    state = &self->states[self->statesN];
    state->valid = true;
    state->pose.pos = vec3_sub(pose.pos, self->utm_offset);
    state->pose.rot = pose.rot;

    self->statesN++;
    GPANO_MSG("loaded state data for %d frames", self->statesN);
      
   
  
    return 0;
}


// Write modified state data
int CGeopano::write_state(gpano_t *self, string filename1)
{
  const char *filename=filename1.c_str();
  FILE *file;
  int i;
  gpano_state_t *state;

  GPANO_MSG("saveing  %s", filename);
  
  file = fopen(filename, "w");
  if (!file)
    return GPANO_ERROR("unable to open %s : %s", filename, strerror(errno));

  for (i = 0; i < self->statesN; i++)
  {
    state = &self->states[i];
    if (!state->valid)
      continue;

    double utm_e, utm_n, alt;
    double rx, ry, rz;

    utm_e = state->pose.pos.x + self->utm_offset.x;
    utm_n = state->pose.pos.y + self->utm_offset.y;
    alt = state->pose.pos.z + self->utm_offset.z;

    //quat_to_rpy(state->pose.rot, &ry, &rx, &rz);
    quat_to_rpy(state->pose.rot, &rx, &ry, &rz);
    ry *= +180/M_PI;
    rx *= +180/M_PI;
    rz *= -180/M_PI;

    fprintf(file, "em1_%d_%05d,%.6f,%.6f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
            0, 0, 0.0, 0.0, alt, rz, ry, rx, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, utm_e, utm_n);
  }

  fclose(file);
  
  return 0;
}


void CGeopano::set_poseNew(gpano_t *self, double utm_e, double utm_n, double alt, double roll, double pitch, double yaw )
//[XYZ]=[utm_e utm_n alt],  [rX rY rZ]=[pitch roll yaw]
{ 
    vec3_t pos=vec3_set(utm_e, utm_n, alt);
    self->poseNew.pos= vec3_sub(pos, self->utm_offset);

    double rz=yaw;
    double rx=pitch;
    double ry=roll;
    //self->poseNew.rot = quat_from_rpy(ry * M_PI/180, rx * M_PI/180, -rz * M_PI/180);
    self->poseNew.rot = quat_from_rpy(rx * M_PI/180, ry * M_PI/180, -rz * M_PI/180);
}

bool CGeopano::set_poseNew(gpano_t *self, string filename1)
{ 
    const char *filename=filename1.c_str();
    GPANO_MSG("loading %s", filename);

    
    FILE *file = fopen(filename, "r");
    if (!file)
      return GPANO_ERROR("unable to open %s : %s", filename, strerror(errno));

    double alt, utm_e, utm_n;
    double rx, ry, rz;
    scanstr_(file, "yaw=");
    scan_(file, 1, (file, "%lf \n", &rz));
    scanstr_(file, "pitch=");
    scan_(file, 1, (file, "%lf \n", &rx));
    scanstr_(file, "roll=");
    scan_(file, 1, (file, "%lf \n", &ry));
    scanstr_(file, "x=");
    scan_(file, 1, (file, "%lf \n", &utm_e));
    scanstr_(file, "y=");
    scan_(file, 1, (file, "%lf \n", &utm_n));
    scanstr_(file, "z=");
    scan_(file, 1, (file, "%lf \n", &alt));
    
    fclose(file);
    
    vec3_t pos = vec3_set(utm_e, utm_n, alt);
    self->poseNew.pos=vec3_sub(pos, self->utm_offset);
    //self->poseNew.rot = quat_from_rpy(ry * M_PI/180, rx * M_PI/180, -rz * M_PI/180);
    self->poseNew.rot = quat_from_rpy(rx * M_PI/180, ry * M_PI/180, -rz * M_PI/180);
    return true;
}

// Read plane data
int CGeopano::read_plane(gpano_t *self, string filename1)
{
  // Initialize state data
   if(self->statesN<=0 ||  self->statesN >= GPANO_MAX_FRAMES) return -1;

   const char *filename=filename1.c_str();
   FILE *file;

   // Prepend the master path
   GPANO_MSG("loading %s", filename);

   file = fopen(filename, "r");
   if (!file)return GPANO_ERROR("unable to open %s : %s", filename, strerror(errno));

    
    fscanf(file, "numplanes %d\n",&self->planeN[self->statesN-1]);
    if(self->planeN[self->statesN-1]<=0)
    {
      self->planeN[self->statesN-1]=0;
      return GPANO_ERROR("worng data format of  %s : %s", filename, strerror(errno));
    }
    if(self->planeN[self->statesN-1]>100)self->planeN[self->statesN-1]=100;
    
   
    if(self->plane[self->statesN-1]!=NULL)delete []self->plane[self->statesN-1];
    
    self->plane[self->statesN-1]=new double[self->planeN[self->statesN-1]*4];//Nx Ny Nz D
    
    double *tmp=self->plane[self->statesN-1];
    memset(tmp,0,sizeof(double)*self->planeN[self->statesN-1]*4);
    
    int label, pixeln;
    for(int i=0; i<self->planeN[self->statesN-1]; i++)
      fscanf(file, "%d %lf %lf %lf %lf %d\n",&label, &tmp[i*4],&tmp[i*4+1],&tmp[i*4+2],&tmp[i*4+3],&pixeln);
    fclose(file);
    
    return 0;
}

vec3_t CGeopano::getGlobalPostion( pose3_t &pose, vec3_t in)
{
	return vec3_transform(pose,in);
}











