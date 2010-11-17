#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>
#include <cstdio>
#include <cstring>
#include <string>
#include <iostream>

using namespace std;


#include <jplv/vec3.h>
#include <jplv/pose3.h>
#include <jplv/jplv_cmod.h>

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include <vector>

#include <glut.h>
#include <gl.h>
#include <glu.h>	// Header File For The GLu32 Library
#include <freeglut.h>

#ifndef EPSILON
#define EPSILON 1e-15
#endif

#ifndef M_PI
#define M_PI  3.14159265
#endif

#define SUCCESS 0
#define FAILURE (-1)

/// Maximum number of frames
#define GPANO_MAX_FRAMES 128


#include <assert.h>
#include <errno.h>


#define scanstr_(fp, str) \
    if (read_scanstr(fp, str) == FAILURE) { \
        fclose(fp); \
        return FAILURE; \
        }

#define scan_(fp, num, args) \
    if (fscanf args != num) { \
        fclose(fp); \
        return FAILURE; \
        }
	

/// @brief State data.
typedef struct
{
  /// Is the state data valid?
  bool valid;

  /// Original frame id
  int frameid;
  
  /// Vehicle nav pose.
  pose3_t pose;
  
} gpano_state_t;


/// @brief Extrinsic camera calibration.
typedef struct
{
  /// State data list (indexed by frame id).
  gpano_state_t states[GPANO_MAX_FRAMES];
  pose3_t poseNew;

  /// UTM offset for state data
  vec3_t utm_offset;
  
  int statesN; 
  
  IplImage *imPano[GPANO_MAX_FRAMES], *raPano[GPANO_MAX_FRAMES], *plPano[GPANO_MAX_FRAMES], *imNewPano, *raNewPano;
  IplImage *imPanoTarget;
  
  double *plane[GPANO_MAX_FRAMES];
  int planeN[GPANO_MAX_FRAMES];
} gpano_t;

/// Macro for error handling
#define GPANO_ERROR(fmt, ...) \
  (fprintf(stderr, "%s:%d error " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__) ? -1 : -1)

/// Macro for warnings
#define GPANO_WARNING(fmt, ...) \
  fprintf(stderr, "%s:%d warning " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__)

/// Macro for messages
#define GPANO_MSG(fmt, ...) \
  fprintf(stderr, fmt "\n", ##__VA_ARGS__)


#endif