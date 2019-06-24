/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#include <stdio.h>
#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "SMsmooth.h"

#undef __FUNC__
#define __FUNC__ "SMerror"
int SMerror(int line,char *func,char* file,char *dir,int n,int p,char *mess)
{
  if (mess) {
    fprintf(stderr,"OPT-MS ERROR:  %s\n",mess);
  }

  switch(n){
  case OPTMS_MEM_ERR:
     fprintf(stderr,"OPT-MS ERROR:  Out of memory. \n");     break;
  case OPTMS_NULL_ERR:
     fprintf(stderr,"OPT-MS ERROR:  Null pointer. \n");     break;
  case OPTMS_INIT_ERR:
     fprintf(stderr,"OPT-MS ERROR:  Data Structure Not Initialized. \n");     break;
  case OPTMS_INPUT_ERR:
     fprintf(stderr,"OPT-MS ERROR:  Incorrect Input \n");  break;
  case OPTMS_FILE_OPEN_ERR:
     fprintf(stderr,"OPT-MS ERROR:  File open error \n");  break;
  case OPTMS_FREE_ERR:
     fprintf(stderr,"OPT-MS ERROR:  Error freeing memory \n");  break;
  case OPTMS_INVALID_MESH_ERR:
     fprintf(stderr,"OPT-MS ERROR:  Invalid Mesh; use SMuntangle to create a valid mesh prior to smoothing \n");  break;
  case OPTMS_DIVIDE_BY_ZERO_ERR:
     fprintf(stderr,"OPT-MS ERROR:  Division by zero \n");  break;
  case OPTMS_DATA_ERR:
     fprintf(stderr,"OPT-MS ERROR:  Incorrect data \n");  break;
  }

  n = OPTMS_PRINT_STACK; /* set it to print the stack */

  fprintf(stderr,"OPT-MS ERROR:  %s()  line %d  in %s/%s\n",func,line,dir,file);

  return(n);
}

#undef __FUNC__
#define __FUNC__ "SMwrite_ordered_points"
int SMwrite_ordered_points(SMlocal_mesh *local_mesh)
{
  FILE *fp;
  int i99,j99;

  if ((fp = fopen("debug.ascii","w")) == NULL) {
      OPTMS_SETERR(OPTMS_FILE_OPEN_ERR,0,"Can't open debug.ascii for writing");
  }

  fprintf(fp,"%d  %d\n",local_mesh->num_incident_vtx, local_mesh->num_tri);
  for (i99=0;i99<local_mesh->dimension;i99++)
      fprintf(fp,"%f  ",local_mesh->original_pt[i99]);
  fprintf(fp,"\n");
  for (i99=0;i99<local_mesh->num_incident_vtx;i99++) {
      for (j99=0;j99<local_mesh->dimension;j99++)
          fprintf(fp,"%f  ",local_mesh->incident_vtx[i99][j99]);
      fprintf(fp,"\n");
  }
  for (i99=0;i99<local_mesh->num_tri;i99++) {
      for (j99=0;j99<local_mesh->dimension;j99++)
          fprintf(fp,"%d  ",local_mesh->vtx_connectivity[i99][j99]);
      fprintf(fp,"\n");
  }
  return(0);
}
