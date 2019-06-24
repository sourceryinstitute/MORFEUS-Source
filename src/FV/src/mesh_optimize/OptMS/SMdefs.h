/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under 
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#ifndef SM_DEFS_H
#define SM_DEFS_H 1

/* added by DPS */
#include <stdlib.h>

#include <fcntl.h>
#if defined(PARCH_sun4) && !defined(__cplusplus) && defined(_Gnu_)
    extern int  open(const char *, int, ...);
    extern int  creat(const char *, unsigned short);
    extern int  write(int, const void *, unsigned int); 
    extern int  close(int);
    extern int  fprintf(FILE*,const char*,...);
#endif
#if defined(solaris) && !defined(__cplusplus) && defined(_Gnu_)
    extern int  open(const char *, int, ...);
    extern int  write(int, const void *, unsigned int); 
    extern int  close(int);
    extern void *malloc(long unsigned int);
    extern void free(void *);
#endif

/* the debugging system
     Level 0 provides no information and the debug macros are empty
     Level 1 provides user function information only
             e.g. threshold set, function used, etc.
     Level 2 provides basic algorithmic information for each local
             submesh
     Level 3 provides more information, data structures, and details
             than most users would want to know about 

     The default is Level 0
*/
#ifdef OPTMS_DBG0
#define OPTMS_DEBUG_LEVEL 0
#endif
#ifdef OPTMS_DBG1
#define OPTMS_DEBUG_LEVEL 1
#endif
#ifdef OPTMS_DBG2
#define OPTMS_DEBUG_LEVEL 2
#endif
#ifdef OPTMS_DBG3
#define OPTMS_DEBUG_LEVEL 3
#endif

#if OPTMS_DEBUG_LEVEL != 0
#define OPTMS_DEBUG_PRINT(level, statement)\
{\
   if ((level <= OPTMS_DEBUG_LEVEL))\
     {\
     fprintf(stdout,statement);\
     fflush(stdout);\
     }\
}
#define OPTMS_DEBUG_ACTION(level, statement)\
{\
   if ((level <= OPTMS_DEBUG_LEVEL))\
     {\
     statement\
     fflush(stdout);\
     }\
}
#else
#define OPTMS_DEBUG_PRINT(level, statement)\
{\
}
#define OPTMS_DEBUG_ACTION(level, statement)\
{\
}
#endif

#define OPTMS_ERROR(msg,action)\
{ \
   fprintf(stderr,msg);\
   action\
}

/* This should be used with debugging level 3.  It writes out a
series of matlab files that illustrate the initial local submesh,
the search direction at each step, and the final local submesh */

#ifdef OPTMS_LOCALTEST
#define OPTMS_MATLAB_ON(a) a
#else
#define OPTMS_MATLAB_ON(a)
#endif
#define OPTMS_MATLAB_OFF

/* Assertions are currently not used very much throughout the code
and need to be expanded and updated */

#ifdef OPTMS_ASSERT
#define OPTMS_ASSERT_ON(a) a
#else
#define OPTMS_ASSERT_ON(a)
#endif
#define OPTMS_ASSERT_OFF

/* The enables statistics gathering;  the user must initialize the
statistics datastructures in the main code and probably wants to print
them out after every smoothing pass 
   The information provided with stats include the number of grid points
   smoothed, how many required optimizaiton, the avg number of optimization
   steps, the a breakdown of the optimization termination criteria
*/

#ifdef OPTMS_STATS
#define OPTMS_STATS_ON(a) a
#else
#define OPTMS_STATS_ON(a)
#endif
#define OPTMS_STATS_OFF

#include "SMuserDefs.h"

#define OPTMS_MAX_NUM_VTX           75
#define OPTMS_MAX_NUM_TRI           150
#define OPTMS_MAX_DIM               3
#define OPTMS_MAX_G_NUM             150
#define OPTMS_DEFAULT_FUNC_PER_TRI  6

#define OPTMS_XDIR 0
#define OPTMS_YDIR 1
#define OPTMS_ZDIR 2

#define OPTMS_TRUE  1
#define OPTMS_FALSE 0

/* various constants to make the code more readable */
#define OPTMS_VALID_MESH      1
#define OPTMS_INVALID_MESH    0
#define OPTMS_MAX_OPT_ITER    20
#define OPTMS_BIG_POS_NMBR    1E300
#define OPTMS_BIG_NEG_NMBR   -1E300
#define OPTMS_MACHINE_EPS     1E-15

/* new initial point constants */
#define OPTMS_NONE      -1 
#define OPTMS_CENTROID   1

/* optimization step constants */
#define OPTMS_STEP_DONE        101
#define OPTMS_STEP_NOT_DONE    102

/* optimization termination constants */
#define OPTMS_STEP_ACCEPTED     100
#define OPTMS_IMP_TOO_SMALL     101
#define OPTMS_FLAT_NO_IMP       102
#define OPTMS_STEP_TOO_SMALL    103
#define OPTMS_EQUILIBRIUM       104
#define OPTMS_ZERO_SEARCH       105
#define OPTMS_MAX_ITER_EXCEEDED 106
#define OPTMS_LAP_ENOUGH        107

/* simple macros */
#define OPTMS_LESS_THAN_MACHINE_EPS(x)   ( ((fabs(x)+1.0) > 1.0) ? 0 : 1 )
#define OPTMS_ISROOT(procinfo) ((procinfo->myid == 0) ? 1 : 0)
#define OPTMS_MAX(a,b) (a > b ? a : b)
#define OPTMS_MIN(a,b) (a < b ? a : b)

#include "SMdata_structs.h"

/* memory allocation macros */
#ifndef OPTMS_MALLOC
#define OPTMS_MALLOC(a,b,c,d) { \
if (c == 0) { \
    a = NULL; \
} else { \
    a = b malloc(c); \
    if (a==NULL) { \
      OPTMS_CHKERR(OPTMS_MEM_ERR); \
    } \
} \
}
#endif

#ifndef OPTMS_FREE
#define OPTMS_FREE(a) \
{ \
    if (a == NULL) { \
       OPTMS_CHKERR(OPTMS_FREE_ERR); \
    } else { \
      free(a); \
    } \
}
#endif

/* printing macros used in debugging */
#define OPTMS_PRINT_ORDERED_PTS(local_mesh) \
{ \
  int i99,j99; \
  for (i99=0;i99<local_mesh->dimension;i99++) \
      printf(" free_vtx[%d] = %f; ",i99,local_mesh->free_vtx[i99]); \
  printf("\n"); \
  for (i99=0;i99<local_mesh->num_incident_vtx;i99++) { \
      for (j99=0;j99<local_mesh->dimension;j99++) \
          printf(" vtx_list[%d][%d]= %f;",i99,j99,local_mesh->incident_vtx[i99][j99]); \
      printf("\n"); \
  } \
}


#define OPTMS_WRITE_ORDERED_PTS(fp,local_mesh) \
{ \
  int i99,j99; \
  fprintf(fp,"%d  %d\n",local_mesh->num_incident_vtx, local_mesh->num_tri);\
  for (i99=0;i99<local_mesh->dimension;i99++) \
      fprintf(fp,"%f  ",local_mesh->free_vtx[i99]); \
  fprintf(fp,"\n"); \
  for (i99=0;i99<local_mesh->num_incident_vtx;i99++) { \
      for (j99=0;j99<local_mesh->dimension;j99++) \
          fprintf(fp,"%f  ",local_mesh->incident_vtx[i99][j99]); \
      fprintf(fp,"\n"); \
  } \
  for (i99=0;i99<local_mesh->num_tri;i99++) { \
      for (j99=0;j99<local_mesh->dimension;j99++) \
          fprintf(fp,"%d  ",local_mesh->vtx_connectivity[i99][j99]); \
      fprintf(fp,"\n"); \
  } \
}

#define OPTMS_WRITE_BINARY_ORDERED_PTS(local_mesh) \
{ \
  int i99,j99; \
  int fd99; \
  char filename99[128]; \
  double temp99; \
  sprintf(filename99,"test.data"); \
  if ((fd99 = creat(filename99, 0666)) == -1) { \
     printf("cannot create filename for writing\n"); \
     exit(0); \
  } \
  temp99 = (double) local_mesh->num_incident_vtx; \
  write(fd99,&temp99,sizeof(double)); \
  temp99 = (double) local_mesh->num_tri;\
  write(fd99,&temp99,sizeof(double)); \
  for (i99=0;i99<local_mesh->dimension;i99++) {\
      temp99 = local_mesh->original_pt[i99]; \
      write(fd99,&temp99,sizeof(double)); \
  } \
  for (i99=0;i99<local_mesh->num_incident_vtx;i99++) { \
      for (j99=0;j99<local_mesh->dimension;j99++){ \
          temp99 = local_mesh->incident_vtx[i99][j99]; \
          write(fd99,&temp99,sizeof(double)); \
      }\
  } \
  if (local_mesh->dimension ==3 ) {\
    for (i99=0;i99<local_mesh->num_tri;i99++) { \
      for (j99=0;j99<local_mesh->dimension;j99++) {\
          temp99 = (double) local_mesh->vtx_connectivity[i99][j99]; \
          write(fd99,&temp99,sizeof(double)); \
      } \
    } \
  } \
  close(fd99); \
}

#define OPTMS_LOCAL_MIN_VOLUME(local_mesh) \
{ \
  int i99,j99; \
  int num_tet99; \
  int num_values99; \
  double vtx99[4][3]; \
  double function99[6]; \
  double min_function99; \
  min_function99=1E300;\
  num_tet99 = local_mesh->num_tri;\
  vtx99[0][0] = local_mesh->free_vtx[0]; \
  vtx99[0][1] = local_mesh->free_vtx[1]; \
  vtx99[0][2] = local_mesh->free_vtx[2]; \
  for (i99=0;i99<num_tet99;i99++) { \
      for (j99=0;j99<3;j99++) {\
         vtx99[j99+1][0]=local_mesh->incident_vtx[local_mesh->vtx_connectivity[i99][j99]][0];\
         vtx99[j99+1][1]=local_mesh->incident_vtx[local_mesh->vtx_connectivity[i99][j99]][1];\
         vtx99[j99+1][2]=local_mesh->incident_vtx[local_mesh->vtx_connectivity[i99][j99]][2];\
      }\
      vComputeTetVolume(vtx99[0],vtx99[1],vtx99[2],vtx99[3],function99,&num_values99); \
      if (function99[0]<min_function99) min_function99=function99[0];\
  }\
  printf("Minimum volume in the local submesh is %f\n",min_function99);\
}

#define OPTMS_PRINT_FUNCTION_VALUES(opt_info) \
{ \
  int i99; \
  for (i99=0;i99<opt_info->num_values;i99++) { \
      printf("Index %d Function Value %f \n",i99,opt_info->function[i99]); \
  } \
}

#define OPTMS_RECORD_ITER_VALUE(opt_info) \
{ \
    opt_info->prev_active_values[opt_info->iter_count] = opt_info->active->true_active_value; \
}

#define OPTMS_DOT(c,a,b,n) {\
  int i99; \
  if (n==2) c = a[0]*b[0] + a[1]*b[1]; \
  else if (n==3) c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];\
  else { \
    for (i99=0;i99<n;i99++) c += a[i99]*b[i99]; \
  } \
}

#define OPTMS_NORMALIZE(v,n) {\
    int i99; \
    double mag99; \
    if (n==2){ \
       mag99 = sqrt(v[0]*v[0] + v[1]*v[1]) ; \
       if (mag99 != 0) { \
          v[0] = v[0]/mag99; \
          v[1] = v[1]/mag99; \
       } \
    } else if (n==3) {\
     mag99 = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) ; \
     if (mag99 != 0) { \
         v[0] = v[0]/mag99; \
         v[1] = v[1]/mag99; \
         v[2] = v[2]/mag99; \
     } \
   } else { \
     for (i99=0;i99<n;i99++) mag99+=v[i99]+v[i99]; \
     if (mag99 != 0) { \
       for (i99=0;i99<n;i99++) v[i99] = v[i99]/mag99;\
     } \
   }\
}

#define OPTMS_COPY_VECTOR(a,b,n) { \
  int i99; \
  if (n==2) { \
     a[0] = b[0];  a[1] = b[1];  \
  } else if (n==3) {\
     a[0] = b[0];  a[1] = b[1];  a[2] = b[2]; \
  } else { \
     for (i99=0;i99<n;i99++) a[i99] = b[i99]; \
  } \
}

#define OPTMS_PRINT_MATRIX(n,m,A) { \
  int i99, j99; \
  for (i99=0;i99<n;i99++) { \
    for (j99=0;j99<m;j99++) { \
      printf("A[%d][%d]=%f",i99,j99,A[i99][j99]);\
    } \
    printf("\n");\
  }\
}

#endif


