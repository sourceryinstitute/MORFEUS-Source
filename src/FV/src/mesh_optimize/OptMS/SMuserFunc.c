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

/****************************************************************************************
                                       SMinitSmoothing
*****************************************************************************************/
#undef __FUNC__
#define __FUNC__ "SMinitSMoothing"
/*@ SMinitSmoothing - Initializes the smoothing data structure and sets
values for the smoothing technique, the mesh quality function,
and threshold usage.

   Input Parameters:
+  argc -  the number of input arguments to the program; used
         to initialize logging.  The logging library checkes for the
         command line argument -logfile filename to print the log
         summary to a file rather then to stdout.  If logging is not
         desired, or argc and argv are not available, NULL may be passed
         in for these arguements.

.  argv - input arguments to the program; used to initialize logging

.  dimension - an integer indicating the problem dimension, either 2 or 3

.  technique - a integer argument to set the smoothing technique used to
        adjust grid point location.
.vb
         Input one of:
                OPTMS_LAPLACIAN_ONLY        (or 1)
                OPTMS_SMART_LAPLACIAN_ONLY  (or 2)
                OPTMS_OPTIMIZATION_ONLY     (or 3)
                OPTMS_COMBINED              (or 4)
                OPTMS_COMBINED1             (or 5)
                OPTMS_COMBINED2             (or 6)
                OPTMS_COMBINED3             (or 7)
                OPTMS_FLOATING_THRESHOLD    (or 8)
                OPTMS_TECHNIQUE_DEFAULT     (or 4)
                OPTMS_DEFAULT               (or -1)

           Note that either the COMBINED or FLOATING_THRESHOLD technique
           is recommended.  The COMBINED approach is the default.
.ve

.  FunctionID - an integer argument used to set the mesh quality measure to be
           optimized.

.vb
       In 2D input one of
              OPTMS_MAX_MIN_ANGLE  (or 1): maximize the minimum angle
              OPTMS_MIN_MAX_COSINE (or 2): minimize the maximum cosine of the angle
              OPTMS_MAX_MIN_COSINE (or 3): maximize the minimum cosine of the angle
              OPTMS_MAX_MIN_SINE   (or 4): maximize the minimum sine of the angle
              OPTMS_MIN_MAX_ANGLE  (or 5): minimize the maximum angle
              OPTMS_MIN_MAX_JACOBIAN_DIFF (or 6): minimize the maximum square of the
                  difference of the current jacobian and the jacobian of an equilateral
                  triangle (scaled by the jacobian of an equilateral triangle)
              OPTMS_MAX_MIN_SCALED_JACOBIAN (or 7): maximize the minimum scaled jacobian
                    for each of the three vertices of a triangle (J/(L1*L2)) where L1 and L2 are
                    the lengths of the incident edges.  Same as OPTMS_MAX_MIN_SINE in the feasible
                    region, but returns negative angles for inverted elements
              OPTMS_MAX_MIN_AREA_LENGTH_RATIO (or 8):  Computes the ratio of the the area
                    of the triangle and the sum of the squares of the length of the edges
              OPTMS_MIN_MAX_LENGTH_AREA_RATIO (or 9): Computes the negtive inverse of the
                    OPTMS_MAX_MIN_AREA_LENGTH_RATIO
              OPTMS_FUNCTION2D_DEFAULT (which is OPTMS_MAX_MIN_SINE) (or 4)
              OPTMS_DEFAULT (-1, which will result in a choice of OPTMS_MAX_MIN_SINE)

            In 3D input one of
               OPTMS_MAX_MIN_DIHEDRAL (or 21): maximize the minimum angle
               OPTMS_MIN_MAX_DIHEDRAL (or 22): minimize the maximum angle
               OPTMS_MAX_MIN_COSINE_DIHEDRAL (or 23): maximize the minimum cosine of the angle
               OPTMS_MIN_MAX_COSINE_DIHEDRAL (or 24): minimize the maximum cosine of the angle
               OPTMS_MAX_SINE_DIHEDRAL (or 25): maximize the minimum sine of the angle
               OPTMS_FUNCTION3D_DEFAULT (which is OPTMS_MAX_SINE_DIHEDRAL) (or 25)
               OPTMS_DEFAULT (-1, which will result in a choice of OPTMS_MAX_SINE_DIHEDRAL)

         The default in 2D is OPTMS_MAX_MIN_SINE and in 3D is OPTMS_MAX_SINE_DIHEDRAL.
.ve
-  Threshold - a double argument that sets the degree value of the threshold
        used in either the COMBINED technique, which has a fixed value, or
        the FLOATING_THRESHOLD technique, which allows the threshold to vary.
        The default values for the quality measures that depend on angle measures
        are 30 and 15 degrees for the COMBINED approach for 2D and 3D, respectively,
        and 10 and 15 degrees for the FLOATING_THRESHOLD approach.  For the
        measures that depend on Jacobian ratios the default is .25.

   Output Parameters:
.   smooth_data - a void data structure that contains the context and data structures
               for smoothing

   Return Values:
.   ierr - an integer error code that is equal to zero if the function call is successful.
           If the function call is unsuccessful, an positive integer will be returned
           and error messages displayed to stderr.

.seealso  SMsetSmoothTechnique(), SMsetSmoothFunction(), SMsetSmoothThreshold(),
              SMfinalizeSmoothing()
@*/
int SMinitSmoothing(int argc, char** argv, int dimension,
                     int technique, int FunctionID, double Threshold,
                     void **ext_smooth_data)
{
    SMsmooth_data *smooth_data;
    int ierr;

    /* check for simple input errors */
    if ((dimension != 2) && (dimension != 3)) {
      OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Dimension not equal 2 or 3\n");
    }

    OPTMS_MALLOC(smooth_data,(SMsmooth_data *),sizeof(SMsmooth_data),1);
    OPTMS_MALLOC(smooth_data->smooth_param,(SMparam *),sizeof(SMparam),1);

#ifdef SUMAA_LOG
    SM_LOG_INIT(argc,argv);
    ierr = SMregisterEvents();  OPTMS_CHKERR(ierr);
#endif

    smooth_data->dimension=dimension;
    (*ext_smooth_data) = (SMsmooth_data *)smooth_data;

    ierr = SMinitSmoothParam(technique,FunctionID,Threshold,(*ext_smooth_data));
           OPTMS_CHKERR(ierr);

    ierr = SMmallocLocalMesh(&smooth_data->local_mesh);  OPTMS_CHKERR(ierr);
    smooth_data->local_mesh->dimension = dimension;
    smooth_data->local_mesh->opt_info->dimension = dimension;

    OPTMS_MALLOC(smooth_data->smooth_stats,(SMstats *),sizeof(SMstats),1);
    smooth_data->smooth_stats->stats_initialized = 0;

    OPTMS_MALLOC(smooth_data->smooth_procinfo,(SMprocinfo *),sizeof(SMprocinfo),1);
    ierr = SMinitProcinfo(smooth_data->smooth_procinfo); OPTMS_CHKERR(ierr);

    ierr = SMmallocQualityTable(&smooth_data->quality_table); OPTMS_CHKERR(ierr);

    OPTMS_MALLOC(smooth_data->untangle_param,(SMuntangle_param *),sizeof(SMuntangle_param),1);
    smooth_data->untangle_param->untangle_technique = OPTMS_UNTANGLING_DEFAULT;

    return(ierr = 0);
}

/******************************************************************************************************
                                        SMsetProblemDimension
*******************************************************************************************************/
#undef __FUNC__
#define __FUNC__ "SMsetProblemDimension"
/*@
  SMsetProblemDimension - This function allows the user to set the dimension of the smoothing
       problem.  Opt-MS currently supports 2D planar smoothing (in the x-y plane)
       and 3D smoothing.  This function call be called any time after SMinitSmoothing
       to set or change the dimension of the problem.

   Input Parameters:
+  smooth_data - a void data structure that contains the context and data structures
       for smoothing.  This structure is created in SMinitSmoothing, which must be called
       prior to calling this routine.
-  dimension - an integer argument to set the dimension of the problem, either 2 or 3.

   Return Values:
.   ierr - an integer error code that is equal to zero if the function call is successful.
           If the function call is unsuccessful, an positive integer will be returned
           and error messages displayed to stderr.

.seealso SMinitSmoothing()
@*/
int SMsetProblemDimension(void *smooth_data, int dimension)
{
   int ierr;
   SMsmooth_data *int_smooth_data;

   int_smooth_data = (SMsmooth_data *) smooth_data;
   OPTMS_CHECK_NULL(smooth_data);
   if ((dimension != 2) && (dimension !=3)) {
     OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Dimension must be 2 or 3");
   }

   int_smooth_data->dimension = dimension;
   int_smooth_data->local_mesh->dimension = dimension;
   int_smooth_data->local_mesh->opt_info->dimension = dimension;
   return(ierr=0);
}

/********************************************************************************************************
                                        SMsetSmoothTechnique
*******************************************************************************************************/
#undef __FUNC__
#define __FUNC__ "SMsetSmoothTechnique"
/*@
  SMsetSmoothTechnique - This function allows the user to change the technique used for mesh
       smoothing at any time during the mesh improvement process.

   Input Parameters:
+  smooth_data - a void data structure that contains the context and data structures
       for smoothing.  This structure is created in SMinitSmoothing, which must be called
       prior to calling this routine.
-  technique - an integer argument to set the smoothing technique used to
        adjust grid point location.
.vb
         Input one of:
                OPTMS_LAPLACIAN_ONLY       (or 1)
                OPTMS_SMART_LAPLACIAN_ONLY (or 2)
                OPTMS_OPTIMIZATION_ONLY    (or 3)
                OPTMS_COMBINED             (or 4)
                OPTMS_COMBINED1            (or 5)
                OPTMS_COMBINED2            (or 6)
                OPTMS_COMBINED3            (or 7)
                OPTMS_FLOATING_THRESHOLD   (or 8)
                OPTMS_TECHNIQUE_DEFAULT    (or 4)
                OPTMS_DEFAULT              (or -1)

           Note that either the COMBINED or FLOATING_THRESHOLD technique
           is recommended.  The COMBINED approach is the default.

      Note that either the COMBINED or FLOATING_THRESHOLD technique
      is recommended.  The COMBINED approach is the default.
.ve

   Return Values:
.   ierr - an integer error code that is equal to zero if the function call is successful.
           If the function call is unsuccessful, an positive integer will be returned
           and error messages displayed to stderr.

.seealso SMinitSmoothing()
@*/
int SMsetSmoothTechnique(void *ext_smooth_data, int technique)
{
    int ierr;
    SMsmooth_data *smooth_data;

    smooth_data = (SMsmooth_data *) ext_smooth_data;
    OPTMS_CHECK_NULL(smooth_data);
    if (technique==OPTMS_DEFAULT) technique = OPTMS_COMBINED;
    smooth_data->smooth_param->smooth_technique = technique;
    return(ierr=0);
}


/********************************************************************************************************
                                        SMsetUntangleTechnique
*******************************************************************************************************/
#undef __FUNC__
#define __FUNC__ "SMsetUntangleTechnique"
/*@
  SMsetUntangleTechnique - This function allows the user to change the technique used for mesh
       untangling

   Input Parameters:
+  smooth_data - a void data structure that contains the context and data structures
       for smoothing.  This structure is created in SMinitSmoothing, which must be called
       prior to calling this routine.
-  technique - an integer argument to set the untangling technique
.vb
         Input one of:
                OPTMS_LAPLACIAN_ONLY       (or 1)
                OPTMS_LINEAR_PROGRAM_ONLY  (or 2)
                OPTMS_COMBINED_UNTANGLING  (or 3)
                OPTMS_UNTANGLING_DEFAULT   (or 3)
                OPTMS_DEFAULT              (or -1)

           Note that either the COMBINED_UNTANGLING technique is recommended.
           The COMBINED approach is the default.
.ve

   Return Values:
.   ierr - an integer error code that is equal to zero if the function call is successful.
           If the function call is unsuccessful, an positive integer will be returned
           and error messages displayed to stderr.

.seealso SMinitSmoothing()
@*/
int SMsetUntangleTechnique(void *ext_smooth_data, int technique)
{
    int ierr;
    SMsmooth_data *smooth_data;

    smooth_data = (SMsmooth_data *) ext_smooth_data;
    OPTMS_CHECK_NULL(smooth_data);
    if (technique==OPTMS_DEFAULT) technique = OPTMS_UNTANGLING_DEFAULT;
    smooth_data->untangle_param->untangle_technique = technique;
    return(ierr=0);
}


/********************************************************************************************************
                                              SMsetSmoothFunction
*******************************************************************************************************/
#undef __FUNC__
#define __FUNC__ "SMsetSmoothFunction"
/*@
    SMsetSmoothFunction - This function allows the user to change the mesh quality function that
         is optimized at any time during the smoothing process.

    Input Parameters:
+   smooth_data - a void data structure that contains the context and data structures
       for smoothing.  This structure is created in SMinitSmoothing, which must be called
       prior to calling this routine.
-   FunctionID - an integer argument used to set the mesh quality measure to be
           optimized.
.vb
       In 2D input one of
              OPTMS_MAX_MIN_ANGLE  (or 1): maximize the minimum angle
              OPTMS_MIN_MAX_COSINE (or 2): minimize the maximum cosine of the angle
              OPTMS_MAX_MIN_COSINE (or 3): maximize the minimum cosine of the angle
              OPTMS_MAX_MIN_SINE   (or 4): maximize the minimum sine of the angle
              OPTMS_MIN_MAX_ANGLE  (or 5): minimize the maximum angle
              OPTMS_MIN_MAX_JACOBIAN_DIFF (or 6): minimize the maximum square of the
                  difference of the current jacobian and the jacobian of an equilateral
                  triangle (scaled by the jacobian of an equilateral triangle)
              OPTMS_MAX_MIN_SCALED_JACOBIAN (or 7): maximize the minimum scaled jacobian
                    for each of the three vertices of a triangle (J/(L1*L2)) where L1 and L2 are
                    the lengths of the incident edges.  Same as OPTMS_MAX_MIN_SINE in the feasible
                    region, but returns negative angles for inverted elements
              OPTMS_MAX_MIN_AREA_LENGTH_RATIO (or 8):  Computes the ratio of the the area
                    of the triangle and the sum of the squares of the length of the edges
              OPTMS_MIN_MAX_LENGTH_AREA_RATIO (or 9): Computes the negtive inverse of the
                    OPTMS_MAX_MIN_AREA_LENGTH_RATIO
              OPTMS_FUNCTION2D_DEFAULT (which is OPTMS_MAX_MIN_SINE) (or 4)
              OPTMS_DEFAULT (-1, which will result in a choice of OPTMS_MAX_MIN_SINE)

            In 3D input one of
               OPTMS_MAX_MIN_DIHEDRAL (or 21): maximize the minimum angle
               OPTMS_MIN_MAX_DIHEDRAL (or 22): minimize the maximum angle
               OPTMS_MAX_MIN_COSINE_DIHEDRAL (or 23): maximize the minimum cosine of the angle
               OPTMS_MIN_MAX_COSINE_DIHEDRAL (or 24): minimize the maximum cosine of the angle
               OPTMS_MAX_SINE_DIHEDRAL (or 25): maximize the minimum sine of the angle
               OPTMS_FUNCTION3D_DEFAULT (which is OPTMS_MAX_SINE_DIHEDRAL) (or 25)
               OPTMS_DEFAULT (-1, which will result in a choice of OPTMS_MAX_SINE_DIHEDRAL)

         The default in 2D is OPTMS_MAX_MIN_SINE and in 3D is OPTMS_MAX_SINE_DIHEDRAL.
.ve

   Return Values:
.   ierr - an integer error code that is equal to zero if the function call is successful.
           If the function call is unsuccessful, an positive integer will be returned
           and error messages displayed to stderr.

.seealso SMinitSmoothing()
@*/
int SMsetSmoothFunction(void *ext_smooth_data,int FunctionID)
{
    int ierr;
    SMsmooth_data *smooth_data;
    SMparam *smooth_param;
    int dimension;

    smooth_data =(SMsmooth_data *)ext_smooth_data;
    OPTMS_CHECK_NULL(smooth_data);
    smooth_param = smooth_data->smooth_param;
    smooth_param->function_id = FunctionID;
    dimension = smooth_data->dimension;

    /* set the default values */
    if (dimension ==2) {
          smooth_param->function_values_per_tri=3;
    } else {
          smooth_param->function_values_per_tri=6;
    }
    smooth_param->ComputeFunctionValues2D = SMcomputeTriSines;
    smooth_param->ComputeGradientValues2D = SMcomputeSineGradients;
    smooth_param->ComputeFunctionValues3D = vSineDihedrals;
    smooth_param->ComputeGradientValues3D = vGradSineDihedrals;

    switch(FunctionID) {
      case OPTMS_MAX_MIN_ANGLE:
        smooth_param->function_values_per_tri=3;
        smooth_param->ComputeFunctionValues2D = SMcomputeTriAngles;
        smooth_param->ComputeGradientValues2D = SMcomputeAngGradients;
        OPTMS_DEBUG_PRINT(1,"Setting the 2D function to OPTMS_MAX_MIN_ANGLE\n");
        break;
      case OPTMS_MIN_MAX_ANGLE:
        smooth_param->function_values_per_tri=3;
        smooth_param->ComputeFunctionValues2D = SMcomputeNegTriAngles;
        smooth_param->ComputeGradientValues2D = SMcomputeNegAngGradients;
        OPTMS_DEBUG_PRINT(1,"Setting the 2D function to OPTMS_MIN_MAX_ANGLE\n");
        break;
      case OPTMS_MIN_MAX_COSINE:
        smooth_param->function_values_per_tri=3;
        smooth_param->ComputeFunctionValues2D = SMcomputeNegTriCosines;
        smooth_param->ComputeGradientValues2D = SMcomputeNegCosGradients;
        OPTMS_DEBUG_PRINT(1,"Setting the 2D function to OPTMS_MIN_MAX_COSINE\n");
        break;
      case OPTMS_MAX_MIN_COSINE:
        smooth_param->function_values_per_tri=3;
        smooth_param->ComputeFunctionValues2D = SMcomputeTriCosines;
        smooth_param->ComputeGradientValues2D = SMcomputeCosGradients;
        OPTMS_DEBUG_PRINT(1,"Setting the 2D function to OPTMS_MAX_MIN_COSINE\n");
        break;
      case OPTMS_MAX_MIN_SINE:
        smooth_param->function_values_per_tri=3;
        smooth_param->ComputeFunctionValues2D = SMcomputeTriSines;
        smooth_param->ComputeGradientValues2D = SMcomputeSineGradients;
        OPTMS_DEBUG_PRINT(1,"Setting the 2D function to OPTMS_MAX_MIN_SINE\n");
        break;
      case OPTMS_MIN_MAX_JACOBIAN_DIFF:
        smooth_param->function_values_per_tri=1;
        smooth_param->ComputeFunctionValues2D = SMcomputeTriJacobians;
        smooth_param->ComputeGradientValues2D = SMcomputeJacobianGradients;
        OPTMS_DEBUG_PRINT(1,"Setting the 2D function to OPTMS_MIN_MAX_JACOBIAN_DIFF\n");
        break;
      case OPTMS_MAX_MIN_SCALED_JACOBIAN:
        smooth_param->function_values_per_tri=3;
        smooth_param->ComputeFunctionValues2D = SMcomputeScaledTriJacobians;
        smooth_param->ComputeGradientValues2D = SMcomputeScaledJacobianGradients;
        OPTMS_DEBUG_PRINT(1,"Setting the 2D function to OPTMS_MAX_MIN_SCALED_JACOBIAN\n");
         break;
      case OPTMS_MAX_MIN_AREA_LENGTH_RATIO:
        smooth_param->function_values_per_tri=1;
        smooth_param->ComputeFunctionValues2D = SMcomputeAreaLengthRatio;
        smooth_param->ComputeGradientValues2D = SMcomputeAreaLengthRatioGradients;
        OPTMS_DEBUG_PRINT(1,"Setting the 2D function to OPTMS_MAX_MIN_AREA_LENGTH_RATIO\n");
        break;
      case OPTMS_MIN_MAX_LENGTH_AREA_RATIO:
        smooth_param->function_values_per_tri=1;
        smooth_param->ComputeFunctionValues2D = SMcomputeLengthAreaRatio;
/*        smooth_param->ComputeGradientValues = SMcomputeLengthAreaRatioGradients; */
        OPTMS_DEBUG_PRINT(1,"Setting the 2D function to OPTMS_MIN_MAX_LENGTH_AREA_RATIO\n");
        break;
      case OPTMS_MAX_MIN_INTERIOR_ANGLE:
        smooth_param->function_values_per_tri=1;
        smooth_param->ComputeFunctionValues2D = SMcomputeInteriorTriAngles;
        smooth_param->ComputeGradientValues2D = SMcomputeInteriorAngGradients;
        break;
      case OPTMS_MAX_MIN_INTERIOR_SINE:
        smooth_param->function_values_per_tri=1;
        smooth_param->ComputeFunctionValues2D = SMcomputeInteriorTriSines;
        smooth_param->ComputeGradientValues2D = SMcomputeInteriorSineGradients;
        break;
      case OPTMS_MIN_MAX_INTERIOR_COSINE:
        smooth_param->function_values_per_tri=1;
        smooth_param->ComputeFunctionValues2D = SMcomputeInteriorTriCosines;
        smooth_param->ComputeGradientValues2D = SMcomputeInteriorCosGradients;
        break;
      case OPTMS_MAX_MIN_INTERIOR_SCALED_JACOBIAN:
        smooth_param->function_values_per_tri=1;
        smooth_param->ComputeFunctionValues2D = SMcomputeInteriorScaledTriJacobians;
        smooth_param->ComputeGradientValues2D = SMcomputeInteriorScaledJacobianGradients;
        break;
      case OPTMS_MIN_MAX_NORM_JAC_SQUARED_2D:
        smooth_param->function_values_per_tri=1;
        smooth_param->ComputeFunctionValues2D = SMNormJacSquared2D;
        smooth_param->ComputeGradientValues2D = SMcomputeNormJacSquaredGradients2D;
        break;
      case OPTMS_MIN_MAX_CONDITION_2D:
        smooth_param->function_values_per_tri=3;
        smooth_param->ComputeFunctionValues2D = SMcondition2D;
        smooth_param->ComputeGradientValues2D = SMgradCondition2D;
        break;

      case OPTMS_MAX_MIN_DIHEDRAL:
        smooth_param->function_values_per_tri=6;
        smooth_param->ComputeFunctionValues3D = vDihedrals;
        smooth_param->ComputeGradientValues3D = vGradDihedrals;
        OPTMS_DEBUG_PRINT(1,"Setting the 3D function to OPTMS_MAX_MIN_DIHEDRAL\n");
        break;
      case OPTMS_MIN_MAX_DIHEDRAL:
        smooth_param->function_values_per_tri=6;
        smooth_param->ComputeFunctionValues3D = vNegateDihedrals;
        smooth_param->ComputeGradientValues3D = vNegateGradDihedrals;
        OPTMS_DEBUG_PRINT(1,"Setting the 3D function to OPTMS_MIN_MAX_DIHEDRAL\n");
        break;
      case OPTMS_MIN_MAX_COSINE_DIHEDRAL:
        smooth_param->function_values_per_tri=6;
        smooth_param->ComputeFunctionValues3D = vNegateCosineDihedrals;
        smooth_param->ComputeGradientValues3D = vNegateGradCosineDihedrals;
        OPTMS_DEBUG_PRINT(1,"Setting the 3D function to OPTMS_MIN_MAX_COSINE_DIHEDRAL\n");
        break;
      case OPTMS_MAX_MIN_COSINE_DIHEDRAL:
        smooth_param->function_values_per_tri=6;
        smooth_param->ComputeFunctionValues3D = vCosineDihedrals;
        smooth_param->ComputeGradientValues3D = vGradCosineDihedrals;
        OPTMS_DEBUG_PRINT(1,"Setting the 3D function to OPTMS_MAX_MIN_COSINE_DIHEDRAL\n");
        break;
      case OPTMS_MAX_SINE_DIHEDRAL:
        smooth_param->function_values_per_tri=6;
        smooth_param->ComputeFunctionValues3D = vSineDihedrals;
        smooth_param->ComputeGradientValues3D = vGradSineDihedrals;
        OPTMS_DEBUG_PRINT(1,"Setting the 3D function to OPTMS_MAX_SINE_DIHEDRAL\n");
        break;
      case OPTMS_MAX_MIN_SCALED_JACOBIAN_3D:
        smooth_param->function_values_per_tri=4;
        smooth_param->ComputeFunctionValues3D = vScaledJacobian;
        smooth_param->ComputeGradientValues3D = vGradScaledJacobian;
        OPTMS_DEBUG_PRINT(1,"Setting the 3D function to OPTMS_MAX_MIN_SCALED_JACOBIAN\n");
        break;
      case OPTMS_MIN_MAX_SRMS_VOLUME_RATIO:
        smooth_param->function_values_per_tri=1;
        smooth_param->ComputeFunctionValues3D = vSMRSVolumeRatio;
        smooth_param->ComputeGradientValues3D = vGradSMRSVolumeRatio;
        OPTMS_DEBUG_PRINT(1,"Setting the 3D function to OPTMS_MIN_MAX_SRMS_VOLUME_RATIO\n");
        break;
      case OPTMS_MIN_MAX_NORM_JAC_SQUARED_3D:
        smooth_param->function_values_per_tri=1;
        smooth_param->ComputeFunctionValues3D = vNormJacSquared;
        smooth_param->ComputeGradientValues3D = vGradNormJacSquared;
        OPTMS_DEBUG_PRINT(1,"Setting the 3D function to OPTMS_MIN_MAX_NORM_JAC_SQUARED\n");
        break;
      case OPTMS_MIN_MAX_CONDITION_3D:
        smooth_param->function_values_per_tri=1;
        smooth_param->ComputeFunctionValues3D = vCondition;
        smooth_param->ComputeGradientValues3D = vGradCondition;
        OPTMS_DEBUG_PRINT(1,"Setting the 3D function to OPTMS_MIN_MAX_CONDITION\n");
        break;
      default:
        OPTMS_DEBUG_PRINT(1,"Using the default functions in both 2 and 3D\n");
        break;
   }
   return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMset2DUserQualityFunction"
void SMset2DUserQualityFunction(void *ext_smooth_data,
                              int values_per_tri,
                              SMfunction_ptr2D userQualityFunc,
                              SMgradfunc_ptr2D userQualityGrad)
{
    SMsmooth_data *smooth_data;
    SMparam *smooth_param;

    smooth_data =(SMsmooth_data *)ext_smooth_data;
    smooth_param = smooth_data->smooth_param;
    smooth_param->function_values_per_tri=values_per_tri;
    smooth_param->ComputeFunctionValues2D = userQualityFunc;
    smooth_param->ComputeGradientValues2D = userQualityGrad;
    OPTMS_DEBUG_PRINT(1,"Setting the quality function to user function \n");
}

/********************************************************************************************************
                                              SMsetSmoothThreshold
*******************************************************************************************************/
#undef __FUNC__
#define __FUNC__ "SMsetSmoothThreshold"
/*@
   SMsetSmoothThreshold - This function allows the user to change the value of the threshold
         when using the combined or the FLOATING_THRESHOLD techniques at any
         time in the smoothing process.

    Input Parameters:
+   smooth_data - a void data structure that contains the context and data structures
       for smoothing.  This structure is created in SMinitSmoothing, which must be called
       prior to calling this routine.
-  Threshold - a double argument that sets the degree value of the threshold
        used in either the COMBINED technique, which has a fixed value, or
        the FLOATING_THRESHOLD technique, which allows the threshold to vary.
        The default values for the quality measures that depend on angles are 15 and
        30 degrees for the COMBINED approach for 3D and 2D respectively, and 10 and
        15 degrees for the FLOATING_THRESHOLD approach.  For the measures that
        depend on Jacobian ratios the default is .25.

   Return Values:
.   ierr - an integer error code that is equal to zero if the function call is successful.
           If the function call is unsuccessful, an positive integer will be returned
           and error messages displayed to stderr.

.seealso SMinitSmoothing()
@*/
int SMsetSmoothThreshold(void *ext_smooth_data, double Threshold)
{
    int ierr;
    SMsmooth_data *smooth_data;
    SMparam *smooth_param;
    int    technique;
    int    function_id;
    double radians;
    double add_value;
    int dimension;

    smooth_data =(SMsmooth_data *) ext_smooth_data;
    OPTMS_CHECK_NULL(smooth_data);
    smooth_param = smooth_data->smooth_param;
    dimension = smooth_data->dimension;
    technique = smooth_param->smooth_technique;
    function_id = smooth_param->function_id;

    /***********************************************************
       set the thresholds for the functions that rely on angles
       for their measure of quality
    ***********************************************************/

    if (function_id == OPTMS_MAX_MIN_ANGLE ||
        function_id == OPTMS_MIN_MAX_COSINE ||
        function_id == OPTMS_MAX_MIN_COSINE ||
        function_id == OPTMS_MAX_MIN_SINE ||
        function_id == OPTMS_MIN_MAX_ANGLE ||
        function_id == OPTMS_MAX_MIN_SCALED_JACOBIAN ||
        function_id == OPTMS_MAX_MIN_DIHEDRAL ||
        function_id == OPTMS_MIN_MAX_DIHEDRAL ||
        function_id == OPTMS_MAX_MIN_COSINE_DIHEDRAL ||
        function_id == OPTMS_MIN_MAX_COSINE_DIHEDRAL ||
        function_id == OPTMS_MAX_SINE_DIHEDRAL ) {

        if (dimension ==2) {
           /* set the default values */
           if (Threshold == OPTMS_DEFAULT) {
             if (technique == OPTMS_FLOATING_THRESHOLD) Threshold = 15.;
             else                                 Threshold = 30.;
          }
          add_value = 5.;
        } else {
           if (Threshold == OPTMS_DEFAULT) {
             if (technique == OPTMS_FLOATING_THRESHOLD) Threshold = 10.;
             else                                 Threshold = 15.;
           }
          add_value = 5.;
       }

      /* set the radians to the default or passed-in value */
      radians = Threshold * 3.14159/180;

    /* modify if necessary for the floating technique */
    if ((technique == OPTMS_FLOATING_THRESHOLD) &&
       (smooth_param->global_min_value != OPTMS_BIG_POS_NMBR)) {
      ierr = SMconvertToDegrees(smooth_param->function_id,
                         &(smooth_param->global_min_value)); OPTMS_CHKERR(ierr);
      /* if it's valid */
      if (smooth_param->global_min_value < 180 &&
	  smooth_param->global_min_value > 0 ) {
          radians = (smooth_param->global_min_value+add_value) * 3.14159/180.;
          OPTMS_DEBUG_ACTION(1,{
             fprintf(stdout,"Set floating threshold: min value %f threshold %f\n",
		 smooth_param->global_min_value,
                 smooth_param->global_min_value+5);
          });
      }
    }
  }

    /***********************************************************
       set the thresholds for the functions that rely on jacobians
       for their measure of quality (note that scaled jacobian is really sine)
    ***********************************************************/

    if ( (function_id == OPTMS_MIN_MAX_JACOBIAN_DIFF) ||
         (function_id == OPTMS_MAX_MIN_AREA_LENGTH_RATIO)) {
       if (Threshold==OPTMS_DEFAULT) {
          Threshold = .25;
       }
       if (Threshold > 1.0) {
          Threshold = .25;
          OPTMS_DEBUG_ACTION(1,{
             fprintf(stdout,"User defined threshold must be less than 1, setting to the default value\n");
          });
       }
    }

      switch(smooth_param->function_id) {
      case OPTMS_MAX_MIN_ANGLE:
        OPTMS_DEBUG_ACTION(1,{
          fprintf(stdout,"Setting the threshold value for OPTMS_MAX_MIN_ANGLE: %f\n",radians);
        });
        smooth_param->lap_accept_value = radians; break;
      case OPTMS_MIN_MAX_ANGLE:
        OPTMS_DEBUG_ACTION(1,{
          fprintf(stdout,"Setting the threshold value for OPTMS_MIN_MAX_ANGLE: %f\n",radians);
        });
        smooth_param->lap_accept_value = radians; break;
      case OPTMS_MIN_MAX_COSINE:
        OPTMS_DEBUG_ACTION(1,{
          fprintf(stdout,"Setting the threshold value for OPTMS_MIN_MAX_COSINE: %f\n",
                cos(radians));
        });
        smooth_param->lap_accept_value = cos(radians); break;
      case OPTMS_MAX_MIN_COSINE:
        OPTMS_DEBUG_ACTION(1,{
          fprintf(stdout,"Setting the threshold value for OPTMS_MAX_MIN_COSINE: %f\n",
                cos(radians));
        });
        smooth_param->lap_accept_value = cos(radians); break;
      case OPTMS_MAX_MIN_SINE:
        OPTMS_DEBUG_ACTION(1,{
          fprintf(stdout,"Setting the threshold value for OPTMS_MAX_MIN_SINE: %f\n",
                  sin(radians));
        });
        smooth_param->lap_accept_value = sin(radians); break;
      case OPTMS_MAX_MIN_SCALED_JACOBIAN:
        OPTMS_DEBUG_ACTION(1,{
          fprintf(stdout,"Setting the threshold value for OPTMS_MAX_MIN_SCALED_JACOBIAN: %f\n",
                  sin(radians));
        });
        smooth_param->lap_accept_value = sin(radians); break;
      case OPTMS_MIN_MAX_JACOBIAN_DIFF:
        OPTMS_DEBUG_ACTION(1,{
          fprintf(stdout,"Setting the threshold value for OPTMS_MIN_MAX_JACOBIAN_DIFF: %f\n",
                  Threshold);
        });
        smooth_param->lap_accept_value = Threshold; break;
      case OPTMS_MAX_MIN_AREA_LENGTH_RATIO:
        OPTMS_DEBUG_ACTION(1,{
          fprintf(stdout,"Setting the threshold value for OPTMS_MAX_MIN_AREA_LENGTH_RATIO: %f\n",
                  Threshold);
        });
        smooth_param->lap_accept_value = Threshold; break;
      case OPTMS_MAX_MIN_DIHEDRAL:
        OPTMS_DEBUG_ACTION(1,{
          fprintf(stdout,"Setting the threshold value for OPTMS_MAX_MIN_DIHEDRAL: %f\n",
                  radians);
        });
        smooth_param->lap_accept_value = radians; break;
      case OPTMS_MIN_MAX_DIHEDRAL:
        OPTMS_DEBUG_ACTION(1,{
          printf("Setting the threshold value for OPTMS_MIN_MAX_DIHEDRAL: %f\n",
                radians);
        });
        smooth_param->lap_accept_value = radians; break;
      case OPTMS_MIN_MAX_COSINE_DIHEDRAL:
        OPTMS_DEBUG_ACTION(1,{
          printf("Setting the threshold value for OPTMS_MIN_MAX_COSINE_DIHEDRAL: %f\n",cos(radians));
        });
        smooth_param->lap_accept_value = cos(radians); break;
      case OPTMS_MAX_MIN_COSINE_DIHEDRAL:
        OPTMS_DEBUG_ACTION(1,{
          fprintf(stdout,"Setting the threshold value for OPTMS_MAX_MIN_COSINE_DIHEDRAL: %f\n",cos(radians));
        });
        smooth_param->lap_accept_value = cos(radians); break;
      case OPTMS_MAX_SINE_DIHEDRAL:
        OPTMS_DEBUG_ACTION(1,{
           fprintf(stdout,"Setting the threshold value for OPTMS_MAX_SINE_DIHEDRAL: %f\n",
                sin(radians));
        });
        smooth_param->lap_accept_value = sin(radians); break;
	/*DPS SRMS Volume Ratio is a special case, since it ranges from -inf to -1 */
	case OPTMS_MIN_MAX_SRMS_VOLUME_RATIO:
        OPTMS_DEBUG_ACTION(1,{
           fprintf(stdout,"Setting the threshold value for OPTMS_MIN_MAX_SRMS_VOLUME_RATIO: %f\n",
		   Threshold);
	});
	if (Threshold == -1) Threshold = -3.0;
	smooth_param->lap_accept_value = Threshold; break;
      default:
        if (dimension==2) {
           OPTMS_DEBUG_PRINT(1,"The smoothing function threshold has not been defined\n");
           OPTMS_DEBUG_PRINT(1,"Setting the default acceptable function value to .5\n");
           smooth_param->lap_accept_value = sin(radians);
        } else {
           OPTMS_DEBUG_PRINT(1,"The smoothing function threshold has not been defined\n");
           OPTMS_DEBUG_PRINT(1,"Setting the default acceptable function value to .5\n");
           smooth_param->lap_accept_value = sin(radians);
      }
   }
   return(ierr=0);
}

/********************************************************************************************************
                                              SMinitGlobalMinValue
*******************************************************************************************************/
#undef __FUNC__
#define __FUNC__ "SMinitGlobalMinValue"
/*@

    SMinitGlobalMinValue - This routine initializes the global minimum value of the quality metric
       to a very large value.   This is used in the floating threshold technique in which the
       minimum value of the quality metric is tracked for the next iteration.

    Input Parameter:
.   smooth_data - a void data structure that contains the context and data structures
              for smoothing.  This structure is created in SMinitSmoothing which must be
              called  prior to calling this routine.  SMinitSmoothingStats must also have been
              called prior to calling this routine.


   Return Values:
.   ierr - an integer error code that is equal to zero if the function call is successful.
           If the function call is unsuccessful, an positive integer will be returned
           and error messages displayed to stderr.

.seealso SMinitSmoothing(), SMsetSmoothThreshold()
@*/
int SMinitGlobalMinValue(void *ext_smooth_data)
{
    int ierr;
    SMsmooth_data *smooth_data;
    smooth_data = (SMsmooth_data *) ext_smooth_data;
    OPTMS_CHECK_NULL(smooth_data);
    OPTMS_DEBUG_PRINT(1,"Setting the global minimum value to 1E300\n");
    smooth_data->smooth_param->global_min_value = OPTMS_BIG_POS_NMBR;
    return(ierr=0);
}

/********************************************************************************************************
                                              SMinitSmoothStats
*******************************************************************************************************/
#undef __FUNC__
#define __FUNC__ "SMinitSmoothStats"
/*@
    SMinitSmoothStats - If statistics gathering has been enabled in the configure process,
         then this routine intializes the statistics structure that records the number
         of cells smoothed, the number of equilibrium points found, and the reason for
         algorithm termination.

    Input Parameter:
.   smooth_data - a void data structure that contains the context and data structures
              for smoothing.  This structure is created in SMinitSmoothing, which must be
              called  prior to calling this routine.

   Return Values:
.   ierr - an integer error code that is equal to zero if the function call is successful.
           If the function call is unsuccessful, an positive integer will be returned
           and error messages displayed to stderr.

    Note:
         It is useful to call this routine in conjunction with SMprintSmoothStats during each
         global pass over the mesh so that the incremental improvment can be monitored.

.seealso SMprintSmoothStats()
@*/
int SMinitSmoothStats(void *ext_smooth_data)
{
    int ierr;
    OPTMS_STATS_ON({
      SMsmooth_data  *smooth_data;
      smooth_data = (SMsmooth_data *) ext_smooth_data;
      OPTMS_CHECK_NULL(smooth_data);
      ierr = SMinitStats(smooth_data->smooth_stats); OPTMS_CHKERR(ierr);
    OPTMS_STATS_OFF});
    return(ierr=0);
}

/********************************************************************************************************
                                              SMprintSmoothStats
*******************************************************************************************************/
#undef __FUNC__
#define __FUNC__ "SMprintSmoothStats"
/*@
    SMprintSmoothStats - If statistics gathering has been enabled in the configure process,
         then this routine prints the statistics that have been accumulated by the smoothing
         code since the last call to SMinitSmoothStats.

    Input Parameter:
.   smooth_data - a void data structure that contains the context and data structures
              for smoothing.  This structure is created in SMinitSmoothing which must be
              called  prior to calling this routine.  SMinitSmoothingStats must also have been
              called prior to calling this routine.

    Note:
.vb
   The information printed includes
        - the total number of nodes smoothed,
        - the number for which Laplacian smoothing was used and the number
           of the those that resulted in an invalid mesh and/or no improvement to
           the mesh
        - the number of nodes for which optimization-based smoothing was used
           (including the average iteration count),
        - the number of cells with no improvement,
        - the averate active value and average improvement, and
        - the termination status for the cells that were smoothed.
.ve

   Return Values:
.   ierr - an integer error code that is equal to zero if the function call is successful.
           If the function call is unsuccessful, an positive integer will be returned
           and error messages displayed to stderr.

.seealso SMinitSmoothing(), SMinitSmoothStats()
@*/
int  SMprintSmoothStats(void *ext_smooth_data)
{
    int ierr;
    OPTMS_STATS_ON({
       SMsmooth_data *smooth_data;
       smooth_data = (SMsmooth_data *) ext_smooth_data;
       OPTMS_CHECK_NULL(smooth_data);
       if (!smooth_data->smooth_stats->stats_initialized) {
          OPTMS_SETERR(OPTMS_INIT_ERR,0,"Stats not initialized");
       }

       ierr = SMprintStats(smooth_data); OPTMS_CHKERR(ierr);

    OPTMS_STATS_OFF});
    return(ierr=0);
}

/********************************************************************************************************
                                              SMinitQualityTable
*******************************************************************************************************/
#undef __FUNC__
#define __FUNC__ "SMinitQualityTable"
/*@
    SMinitQualityTable - This function allows the user to take advantage of the
        quality metrics implemented in the OptMS code to analyze the quality of their
        mesh.

    Input Parameter:
.   smooth_data - a void data structure that contains the context and data structures
              for smoothing.  This structure is created in SMinitSmoothing, which must be
              called  prior to calling this routine.   This routine should be called before each
              global pass of measuring quality or information from the previous pass will
              continue to be accumulated.

   Return Values:
.   ierr - an integer error code that is equal to zero if the function call is successful.
           If the function call is unsuccessful, an positive integer will be returned
           and error messages displayed to stderr.

    Note:
.n      In 2D the min, max, and average values of the triangle angles, deviation from
         an equilateral triangle, the scaled jacobians, and the triangle area are printed.
.n
.n      In 3D, the min, max, and average values of the tetrahedral dihedral angles,
         scaled Jacobians, the ratio of sum of the squares of the length of the edges raised
         to the 3/2 power to the volume, and the tetrahdral volume are printed.
.n
.n      If any triangle areas or tetrahedral volumes are negative, the mesh is considered
         to be invalid, and SMuntangle should be called to try to create a valid mesh.

.seealso SMinitSmoothing(), SMuntangle(), SMaccumulateQualityInformation(),
              SMprintQualityInformation()
@*/
int SMinitQualityTable(void *ext_smooth_data)
{
    int i, ierr;
    SMsmooth_data *smooth_data;

    smooth_data = (SMsmooth_data *) ext_smooth_data;
    OPTMS_CHECK_NULL(smooth_data);
    smooth_data->quality_table->initialized=OPTMS_TRUE;
    smooth_data->quality_table->num_tangled_elements=0;;
    smooth_data->quality_table->mesh_validity=1;

    for (i=0;i<smooth_data->quality_table->num_functions;i++) {
         smooth_data->quality_table->measure[i]->min_value = OPTMS_BIG_POS_NMBR;
         smooth_data->quality_table->measure[i]->max_value = OPTMS_BIG_NEG_NMBR;
         smooth_data->quality_table->measure[i]->avg_min_value = 0.0;
         smooth_data->quality_table->measure[i]->avg_max_value = 0.0;
         smooth_data->quality_table->measure[i]->avg_value=0.0;
         smooth_data->quality_table->measure[i]->num_function_values=0;
         smooth_data->quality_table->measure[i]->num_elements=0;
    }
    return(ierr=0);
}

/********************************************************************************************************
                                          SMaccumulateQualityInformation
*******************************************************************************************************/
#undef __FUNC__
#define __FUNC__ "SMaccumulateQualityInformation"
/*@
    SMaccumulateQualityInformation - This function computes the quality information
        for a right-handed triangle or tetrahedra.

    Input Parameters:
+   smooth_data - a void data structure that contains the context and data structures
              for smoothing.  This structure is created in SMinitSmoothing, which must be
              called  prior to calling this routine.
-   vtx - a matrix containing the coordinates of the nodes of the triangle or tetrahedra.
             Of dimension 3 x 2 for triangles, and 4 x 3 for tetrahedra.

   Return Values:
.   ierr - an integer error code that is equal to zero if the function call is successful.
           If the function call is unsuccessful, an positive integer will be returned
           and error messages displayed to stderr.

    Notes:
         In 2D the min, max, and average values of the triangle angles, deviation from
         an equilateral triangle, the scaled jacobians, and the triangle area are computed.
.n
.n      In 3D, the min, max, and average values of the tetrahedral dihedral angles,
         scaled Jacobians, the ratio of sum of the squares of the length of the edges raised
         to the 3/2 power to the volume, and the tetrahdral volume are computed.

.seealso SMinitSmoothing(), SMinitQualityTable(), SMprintQualityInformation()
@*/
int SMaccumulateQualityInformation(void *ext_smooth_data, double **vtx)
{
     int ierr;
     int i,j;
     int num_values;
     int dimension;
     double function[6];
     int inverted=OPTMS_FALSE;

     SMsmooth_data *smooth_data;
     SMquality_table *quality_table;

     smooth_data = (SMsmooth_data *) ext_smooth_data;
     OPTMS_CHECK_NULL(smooth_data);
     dimension = smooth_data->dimension;
     quality_table = smooth_data->quality_table;

     for (i=0;i<quality_table->num_functions;i++) {
          /* 2D functions */
          if (dimension==2) {
          switch(i) {
          case 0:
              ierr = SMcomputeTriArea(vtx[0],vtx[1],vtx[2],function,&num_values); OPTMS_CHKERR(ierr);
              ierr = SMinsertQualityInfo(quality_table,i,function,num_values); OPTMS_CHKERR(ierr);
              if (function[0] <= OPTMS_MACHINE_EPS) {
                   quality_table->num_tangled_elements++;
                   quality_table->mesh_validity=0;
                   inverted=OPTMS_TRUE;
              }
              break;
          case 1:
	      if (inverted) {
                for (j=0;j<3;j++) function[j]=1E6;
              } else {
                ierr = SMcomputeTriAngles(vtx[0],vtx[1],vtx[2],function,&num_values);
                       OPTMS_CHKERR(ierr);
                for (j=0;j<num_values;j++) function[j]=function[j]/3.14159*180;
              }
              ierr = SMinsertQualityInfo(quality_table,i,function,num_values); OPTMS_CHKERR(ierr);
              break;
          case 2:
	      if (inverted) {
                for (j=0;j<1;j++) function[j]=1E6;
              } else {
                ierr = SMcomputeTriJacobians(vtx[0],vtx[1],vtx[2],function,&num_values);
                       OPTMS_CHKERR(ierr);
              }
              ierr = SMinsertQualityInfo(quality_table,i,function,num_values); OPTMS_CHKERR(ierr);
              break;
          case 3:
	      if (inverted) {
                for (j=0;j<3;j++) function[j]=1E6;
              } else {
                ierr = SMcomputeScaledTriJacobians(vtx[0],vtx[1],vtx[2],function,&num_values);
                       OPTMS_CHKERR(ierr);
              }
              ierr = SMinsertQualityInfo(quality_table,i,function,num_values); OPTMS_CHKERR(ierr);
              break;
          case 4:
	      if (inverted) {
                for (j=0;j<1;j++) function[j]=1E6;
              } else {
                ierr = SMcondition2D(vtx[0],vtx[1],vtx[2],function,&num_values); OPTMS_CHKERR(ierr);
                for (j=0;j<num_values;j++) function[j]=-function[j];
              }
              ierr = SMinsertQualityInfo(quality_table,i,function,num_values); OPTMS_CHKERR(ierr);
              break;
           }
          }
          if (dimension==3) {
          switch(i) {
          /* 3D functions */
          case 5:
              ierr = vComputeTetVolume(vtx[0],vtx[1],vtx[2],vtx[3],function,&num_values);
                     OPTMS_CHKERR(ierr);
              ierr = SMinsertQualityInfo(quality_table,i,function,num_values); OPTMS_CHKERR(ierr);
              if (function[0] <= OPTMS_MACHINE_EPS) {
                   quality_table->num_tangled_elements++;
                   quality_table->mesh_validity=0;
                   inverted = OPTMS_TRUE;
              }
              break;
          case 6:
	      if (inverted) {
                for (j=0;j<6;j++) function[j]=1E6;
              } else {
                ierr = vDihedrals(vtx[0],vtx[1],vtx[2],vtx[3],function,&num_values); OPTMS_CHKERR(ierr);
                for (j=0;j<num_values;j++) function[j]=function[j]/3.14159*180;
              }
              ierr = SMinsertQualityInfo(quality_table,i,function,num_values); OPTMS_CHKERR(ierr);
              break;
          case 7:
	      if (inverted) {
                for (j=0;j<6;j++) function[j]=1E6;
              } else {
                ierr = vScaledJacobian(vtx[0],vtx[1],vtx[2],vtx[3],function,&num_values);
                       OPTMS_CHKERR(ierr);
              }
              ierr = SMinsertQualityInfo(quality_table,i,function,num_values); OPTMS_CHKERR(ierr);
              break;
          case 8:
	      if (inverted) {
                for (j=0;j<1;j++) function[j]=1E6;
              } else {
                ierr = vSMRSVolumeRatio(vtx[0],vtx[1],vtx[2],vtx[3],function,&num_values);
                       OPTMS_CHKERR(ierr);
                for (j=0;j<num_values;j++) function[j]=-function[j];
              }
              ierr = SMinsertQualityInfo(quality_table,i,function,num_values); OPTMS_CHKERR(ierr);
              break;
          case 9:
	      if (inverted) {
                for (j=0;j<1;j++) function[j]=1E6;
              } else {
  	        ierr = vCondition(vtx[0],vtx[1],vtx[2],vtx[3],function,&num_values); OPTMS_CHKERR(ierr);
                for (j=0;j<num_values;j++) function[j]=-function[j];
              }
              ierr = SMinsertQualityInfo(quality_table,i,function,num_values); OPTMS_CHKERR(ierr);
              break;
           }
       }
    }
    return(ierr);
}

#undef __FUNC__
#define __FUNC__ "SMinsertQualityInfo"
int SMinsertQualityInfo(SMquality_table *quality_table, int measure_id,
                                      double *function, int num_values)
{
    int ierr;
    int i;
    double element_min = OPTMS_BIG_POS_NMBR;
    double element_max = OPTMS_BIG_NEG_NMBR;

    quality_table->measure[measure_id]->num_function_values+=num_values;
    quality_table->measure[measure_id]->num_elements++;

    for (i=0;i<num_values;i++) {

      if (function[i] < element_min) element_min = function[i];
      if (function[i] > element_max) element_max = function[i];

       quality_table->measure[measure_id]->avg_value+=function[i];
    }
    if (element_min < quality_table->measure[measure_id]->min_value) {
         quality_table->measure[measure_id]->min_value=element_min;
    }
    if (element_max > quality_table->measure[measure_id]->max_value) {
         quality_table->measure[measure_id]->max_value=element_max;
    }
    quality_table->measure[measure_id]->avg_min_value+=element_min;
    quality_table->measure[measure_id]->avg_max_value+=element_max;

    return(ierr=0);
}

/********************************************************************************************************
                                          SMprintQualityInformation
*******************************************************************************************************/
#undef __FUNC__
#define __FUNC__ "SMprintQualityInformation"
/*@
    SMprintQualityInformation - This function prints the quality information accumulated
        since the last call to SMinitQualityInformation.

    Input Parameters:
.   smooth_data - a void data structure that contains the context and data structures
              for smoothing.  This structure is created in SMinitSmoothing, which must be
              called  prior to calling this routine.

   Return Values:
.   ierr - an integer error code that is equal to zero if the function call is successful.
           If the function call is unsuccessful, an positive integer will be returned
           and error messages displayed to stderr.

    Notes:
         In 2D the min, max, and average values of the triangle angles, deviation from
         an equilateral triangle, the scaled jacobians, and the triangle area are printed.
.n
.n      In 3D, the min, max, and average values of the tetrahedral dihedral angles,
         scaled jacobians, the ratio of sum of the squares of the length of the edges raised
         to the 3/2 power to the volume, and the tetrahdral volume are printed.
.n
.n      If any triangle areas or tetrahedral volumes are negative, the mesh is considered
         to be invalid, and SMuntangle should be called to try to create a valid mesh.

.seealso SMinitSmoothing(), SMinitQualityInformation(), SMaccumulateQualityInformation()
@*/
int SMprintQualityInformation(void *ext_smooth_data)
{
     int ierr;
     int i;
     int dimension;
     int measure_ind;
     SMsmooth_data *smooth_data;
     SMquality_table *quality_table;
     SMqualityMeasure *measure;

     smooth_data = (SMsmooth_data *) ext_smooth_data;
     OPTMS_CHECK_NULL(smooth_data);
     dimension = smooth_data->dimension;
     quality_table = smooth_data->quality_table;

     if (dimension == 2) measure_ind=0;
     if (dimension == 3) measure_ind=5;

     printf("\nMesh Quality Information for %d Elements\n",quality_table->measure[measure_ind]->num_elements);
     printf("-------------------------------------------------------------------------------------\n");
     printf("Quality Metric (Target)           Min Value  Max Value  Avg Value  Avg Min    Avg Max\n");
     printf("-------------------------------------------------------------------------------------\n");
     for (i=0;i<quality_table->num_functions;i++) {
       if (quality_table->measure[i]->dimension == dimension) {
         measure = quality_table->measure[i];
         printf("%20s (%3.2e) %10.2e %10.2e %10.2e %10.2e %10.2e (%d)\n",
                     measure->name,
                     measure->target,
                     measure->min_value,
                     measure->max_value,
                     measure->avg_value/measure->num_function_values,
                     measure->avg_min_value/measure->num_elements,
                     measure->avg_max_value/measure->num_elements,
                     measure->num_function_values);
       }
     }
     printf("-----------------------------------------------------------\n");
     printf("There are %d invalid elements in the mesh \n\n",quality_table->num_tangled_elements);
     return(ierr=0);
}

/********************************************************************************************************
                                          SMinvalidMesh
*******************************************************************************************************/
#undef __FUNC__
#define __FUNC__ "SMinvalidMesh"
/*@
    SMisValidMesh - This function returns a Boolean value after testing to see if the
        mesh contains any invalid elements.  The test is based on information contained
        in the quality table and this information must be accumulated before this test can
        be made.  This is useful in determining if the mesh requires untangling.

    Input Parameters:
.   smooth_data - a void data structure that contains the context and data structures
              for smoothing.  This structure is created in SMinitSmoothing, which must be
              called  prior to calling this routine.
.   valid - an pointer to an integer value.  On return, this value will equal 1 if the mesh
            is valid, 0 if it is not.

   Return Values:
.   ierr - an integer error code that is equal to zero if the function call is successful.
           If the function call is unsuccessful, an positive integer will be returned
           and error messages displayed to stderr.

.seealso SMinitSmoothing, SMinitQualityInformation(), SMaccumulateQualityInformtaion,
             SMuntangle()
@*/
int SMisValidMesh(void *ext_smooth_data, int *valid)
{
   int ierr;
   SMsmooth_data *smooth_data;

   OPTMS_CHECK_NULL(ext_smooth_data);
   smooth_data = (SMsmooth_data *) ext_smooth_data;
   *valid=1;
   if (smooth_data->quality_table->mesh_validity == 0) *valid=0;
   return(ierr=0);
}

/********************************************************************************************************
                                          SMsetMeshValidity
*******************************************************************************************************/
#undef __FUNC__
#define __FUNC__ "SMsetMeshValidity"
/*@
    SMsetMeshValidity - This function allows the user to set the mesh validity without using the
        quality functions provided by the Opt-MS package.  This is useful if the user knows the
        mesh requires untangling and doesn't want to evaluate the quality of the entire mesh to
        determine if there are inverted elements.

    Input Parameters:
+   mesh_validity - a Boolean value of 1 if the mesh is valid and 0 if the mesh contains
              elements with negative area
-   smooth_data - a void data structure that contains the context and data structures
              for smoothing.  This structure is created in SMinitSmoothing, which must be
              called  prior to calling this routine.

   Return Values:
.   ierr - an integer error code that is equal to zero if the function call is successful.
           If the function call is unsuccessful, an positive integer will be returned
           and error messages displayed to stderr.

.seealso SMinitSmoothing, SMinitQualityInformation(), SMaccumulateQualityInformtaion,
             SMuntangle()
@*/
int SMsetMeshValidity(int mesh_validity, void *ext_smooth_data)
{
   int ierr;
   SMsmooth_data *smooth_data;
   OPTMS_CHECK_NULL(ext_smooth_data);
   smooth_data = (SMsmooth_data *) ext_smooth_data;
   smooth_data->quality_table->mesh_validity=mesh_validity;
   return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMfinalizeSmoothing"
/*@
SMfinalizeSmoothing - This routine frees all the memory allocated in SMinitSmoothing
including the data structure smooth_data.  This routine should be called when mesh
optimization is complete.

   Input Parameters:
.  smooth_data - a void data structure that contains the context and data structures
       for smoothing.  This structure is created in SMinitSmoothing, which must be called
       prior to calling this routine.

    Return Values:
.   ierr - an integer error code that is equal to zero if the function call is successful.
           If the function call is unsuccessful, an positive integer will be returned
           and error messages displayed to stderr.

.seealso SMinitSmoothing()
@*/
int SMfinalizeSmoothing(void *ext_smooth_data)
{
    int ierr;
    SMsmooth_data *smooth_data;
    smooth_data = (SMsmooth_data *) ext_smooth_data;

#ifdef PARALLEL_LOG
    SM_LOG_PRINT(smooth_data->smooth_procinfo->nprocs,
                 smooth_data->smooth_procinfo->myid,
                 smooth_data->smooth_procinfo->procset);
#else
    SM_LOG_PRINT();
#endif
    SM_LOG_FREE;

    ierr = SMfreeLocalMesh(smooth_data->local_mesh); OPTMS_CHKERR(ierr);
    ierr = SMfreeParam(smooth_data->smooth_param); OPTMS_CHKERR(ierr);
    ierr = SMfreeProcinfo(smooth_data->smooth_procinfo); OPTMS_CHKERR(ierr);
    ierr = SMfreeQualityTable(smooth_data->quality_table); OPTMS_CHKERR(ierr);
    OPTMS_FREE(smooth_data->smooth_stats);
    /*DPS added */
    OPTMS_FREE(smooth_data->untangle_param);
    /* DPS end */
    OPTMS_FREE(smooth_data);
    return(ierr=0);
}
