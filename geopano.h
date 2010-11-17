#ifndef GEOPANO_H
#define GEOPANO_H

#include "global.h"

class CGeopano
{
  public:
  int read_state(gpano_t *self, string filename);
  int write_state(gpano_t *self, string filename);
  bool set_poseNew(gpano_t *self, string filename1);
  void set_poseNew(gpano_t *self, double utm_e, double utm_n, double alt, double roll, double ptich, double yaw );
  int read_plane(gpano_t *self, string filename);
  vec3_t getGlobalPostion( pose3_t &pose, vec3_t in);
};

#endif // GEOPANO_H
