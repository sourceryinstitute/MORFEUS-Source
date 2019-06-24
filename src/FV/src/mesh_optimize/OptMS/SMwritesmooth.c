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
#define __FUNC__ "SMwriteLocalMesh"
int	SMwriteLocalMesh(FILE *fp, SMlocal_mesh *local_mesh)
{
    int ierr;
    OPTMS_CHECK_NULL(fp);
    OPTMS_CHECK_NULL(local_mesh);
    ierr = SMwriteLocalTriangleList(fp,local_mesh);  OPTMS_CHKERR(ierr);
    return(ierr = 0);
}

#undef __FUNC__
#define __FUNC__ "SMwriteLocalAxes"
int	SMwriteLocalAxes(FILE *fp, SMlocal_mesh *mesh)
{
    int ierr;

    OPTMS_CHECK_NULL(fp);
    OPTMS_CHECK_NULL(mesh);

    if (mesh->dimension == 2) {
       if (mesh->free_vtx[OPTMS_XDIR] < mesh->min[OPTMS_XDIR]) 
                 mesh->min[OPTMS_XDIR]=mesh->free_vtx[OPTMS_XDIR];
       if (mesh->free_vtx[OPTMS_YDIR] < mesh->min[OPTMS_YDIR]) 
                 mesh->min[OPTMS_YDIR]=mesh->free_vtx[OPTMS_YDIR];
       if (mesh->free_vtx[OPTMS_XDIR] > mesh->max[OPTMS_XDIR]) 
                 mesh->max[OPTMS_XDIR]=mesh->free_vtx[OPTMS_XDIR];
       if (mesh->free_vtx[OPTMS_YDIR] > mesh->max[OPTMS_YDIR]) 
                mesh->max[OPTMS_YDIR]=mesh->free_vtx[OPTMS_YDIR];

        fprintf(fp,"axes('Xlim',[%16.10f,%16.10f],'Ylim',[%16.10f,%16.10f],'Box','on ',",
                mesh->min[OPTMS_XDIR],mesh->max[OPTMS_XDIR],mesh->min[OPTMS_YDIR],mesh->max[OPTMS_YDIR]);
        fprintf(fp,"'XTick',[%16.10f,%16.10f],'YTick',[%16.10f,%16.10f]);\n",
                mesh->min[OPTMS_XDIR],mesh->max[OPTMS_XDIR],mesh->min[OPTMS_YDIR],mesh->max[OPTMS_YDIR]);
        fprintf(fp,"hold on;\n");
    } else if (mesh->dimension == 3) {
        OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Can't look at three dimensions in matlab");
    } else {
        OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Dimensions must be 2 or 3\n");
    }
    return(ierr=0);
}  

#undef __FUNC__
#define __FUNC__ "SMwriteLocalTriangleList"
int	SMwriteLocalTriangleList(FILE *fp, SMlocal_mesh *local_mesh)
{
    int ierr;
    int i, num_tri;
    double  *vtx1, *vtx2, *vtx3;
    double  x_0,x_1,x_2;
    double  y_0,y_1,y_2;

    OPTMS_CHECK_NULL(fp);
    OPTMS_CHECK_NULL(local_mesh);

    ierr = SMwriteLocalAxes(fp,local_mesh); OPTMS_CHKERR(ierr);

    num_tri = local_mesh->num_tri;
    for (i=0;i<num_tri;i++) {
        vtx1 = local_mesh->free_vtx;
        vtx2 = local_mesh->incident_vtx[local_mesh->vtx_connectivity[i][0]];
        vtx3 = local_mesh->incident_vtx[local_mesh->vtx_connectivity[i][1]];

        x_0 = vtx1[OPTMS_XDIR];    y_0 = vtx1[OPTMS_YDIR];
        x_1 = vtx2[OPTMS_XDIR];    y_1 = vtx2[OPTMS_YDIR];
        x_2 = vtx3[OPTMS_XDIR];    y_2 = vtx3[OPTMS_YDIR];

        fprintf(fp,"line([%f, %f],[%f, %f]);\n",x_0,x_1,y_0,y_1);
        fprintf(fp,"line([%f, %f],[%f, %f]);\n",x_1,x_2,y_1,y_2);
        fprintf(fp,"line([%f, %f],[%f, %f]);\n",x_0,x_2,y_0,y_2);
    }
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMwriteActiveSet"
int SMwriteActiveSet(FILE *fp,SMlocal_mesh *local_mesh)
{
    int ierr;
    int     i, ind;
    int     tri_num;
    int     vert_id, ind1, ind2;
    double  x_pt, y_pt, x_1, x_2, y_1, y_2;
    double  rise, run, slope, b_int;
    double  midx, midy;
    double  xdiff, x_temp, y_temp;

    OPTMS_CHECK_NULL(fp);
    OPTMS_CHECK_NULL(local_mesh);

    /*OPTMS_PRINT_FUNCTION_VALUES(local_mesh->opt_info); */
    /* marks the minimum angle set with cyan * */
    if (local_mesh->opt_info->num_values == local_mesh->num_tri*3) {
      ind1 = ind2 = 0;
      for (i=0;i<local_mesh->opt_info->active->num_active;i++) {
	ind = local_mesh->opt_info->active->active_ind[i];

        tri_num = (int) ind/3;	vert_id = ind % 3;
        ind1 = local_mesh->vtx_connectivity[tri_num][0];
        ind2 = local_mesh->vtx_connectivity[tri_num][1];
        if (vert_id == 0) {
            x_pt = local_mesh->free_vtx[OPTMS_XDIR];
	    y_pt = local_mesh->free_vtx[OPTMS_YDIR];
	    x_1  = local_mesh->incident_vtx[ind1][OPTMS_XDIR];
	    x_2  = local_mesh->incident_vtx[ind2][OPTMS_XDIR];
	    y_1  = local_mesh->incident_vtx[ind1][OPTMS_YDIR];
	    y_2  = local_mesh->incident_vtx[ind2][OPTMS_YDIR];
        } else if (vert_id == 1) {
	    x_pt = local_mesh->incident_vtx[ind1][OPTMS_XDIR];
	    y_pt = local_mesh->incident_vtx[ind1][OPTMS_YDIR];
                  x_1  = local_mesh->free_vtx[OPTMS_XDIR];
	    y_1  = local_mesh->free_vtx[OPTMS_YDIR];
	    x_2  = local_mesh->incident_vtx[ind2][OPTMS_XDIR];
	    y_2  = local_mesh->incident_vtx[ind2][OPTMS_YDIR];
	} else {
	    x_pt = local_mesh->incident_vtx[ind2][OPTMS_XDIR];
	    y_pt = local_mesh->incident_vtx[ind2][OPTMS_YDIR];
                  x_1  = local_mesh->free_vtx[OPTMS_XDIR];
	    y_1  = local_mesh->free_vtx[OPTMS_YDIR];
	    x_2  = local_mesh->incident_vtx[ind1][OPTMS_XDIR];
	    y_2  = local_mesh->incident_vtx[ind1][OPTMS_YDIR];
	}
 
	midx = (x_1 + x_2) / 2;  midy = (y_1 + y_2) / 2;
	xdiff = x_pt - midx;

	rise = y_pt - midy;      run  = x_pt - midx;

        if (fabs(run)<OPTMS_MACHINE_EPS) {
           OPTMS_SETERR(OPTMS_INPUT_ERR,0,
                        "Colocated vertices resulting in division by zero\n");
        }

	slope = rise/run;
	b_int = midy - slope * midx;
	x_temp = x_pt - .1*xdiff;
	y_temp = slope*x_temp + b_int;

	fprintf(fp,"plot( %f , %f, 'c*'); \n",x_temp,y_temp);
       }
    }
    if (local_mesh->opt_info->num_values == local_mesh->num_tri) {
      ind1 = ind2 = 0;
      for (i=0;i<local_mesh->opt_info->active->num_active;i++) {
        ind = local_mesh->opt_info->active->active_ind[i];

        tri_num = (int) ind;	
        ind1 = local_mesh->vtx_connectivity[tri_num][0];
        ind2 = local_mesh->vtx_connectivity[tri_num][1];
        x_pt = local_mesh->free_vtx[OPTMS_XDIR];
        y_pt = local_mesh->free_vtx[OPTMS_YDIR];
        x_1  = local_mesh->incident_vtx[ind1][OPTMS_XDIR];
        x_2  = local_mesh->incident_vtx[ind2][OPTMS_XDIR];
        y_1  = local_mesh->incident_vtx[ind1][OPTMS_YDIR];
        y_2  = local_mesh->incident_vtx[ind2][OPTMS_YDIR];
 
        fprintf(fp,"line([%f, %f],[%f, %f],'Color','c');\n",x_1,x_2,y_1,y_2);
        fprintf(fp,"line([%f, %f],[%f, %f],'Color','c');\n",x_pt,x_2,y_pt,y_2);
        fprintf(fp,"line([%f, %f],[%f, %f],'Color','c');\n",x_1,x_pt,y_1,y_pt);
      }
    }
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMwriteSearch"
int SMwriteSearch(FILE *fp, SMlocal_mesh *local_mesh)
{
    int ierr;
    double x_1, y_1, x_2, y_2;

    OPTMS_CHECK_NULL(fp);
    OPTMS_CHECK_NULL(local_mesh);

    x_1 = local_mesh->free_vtx[OPTMS_XDIR];
    y_1 = local_mesh->free_vtx[OPTMS_YDIR];
    x_2 = x_1 + local_mesh->opt_info->search[OPTMS_XDIR];
    y_2 = y_1 + local_mesh->opt_info->search[OPTMS_YDIR];

    fprintf(fp,"line([%f, %f],[%f, %f],'Color','m');\n",x_1,x_2,y_1,y_2);
    return(ierr=0);

}

#undef __FUNC__
#define __FUNC__ "SMwritePoint"
int SMwritePoint(FILE *fp, double x, double y)
{
    int ierr;
    OPTMS_CHECK_NULL(fp);
    fprintf(fp,"plot(%f,%f,'r*');\n",x,y);
    return(ierr=0);
}
