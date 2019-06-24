/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under 
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
/* intrinsic.c */
/* 
   Andrew Mauer
   Options:

   SILENT_EXCEPTIONS: No message is printed out when an exception 
                      occurs if this is defined.

   BIG_FABS_PARTIAL: If an exception occurs for in g_ad_fabs, the
                     partial dz/dx is 10e+40 rather than the default 0.

   Other curiosities:

   g_ad_fabs: Contains a call to user_adprint that never occurs. This
           ensures that user_adprint will be linked into the resulting
           executable.

*/

#include "SMderiv.h"
#include <stdio.h>
#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#define IN_AD_INTRINSIC_C
#include "SMintrinsic.h"
#include "SMsmooth.h"

#ifdef __cplusplus
/* To use C++ you will need to compile the axpy library with C++ as well. */
#error This file should be compiled with a standard (ANSI) C compiler.
#endif

#ifdef AD_DERIV_TYPE_ARRAY

/* This is a hack --- this function should be provided somewhere sensible. */

#undef __FUNC__
#define __FUNC__ "is_zero_vector" 
static int
is_zero_vector (DERIV_TYPE x)
{
     int i;
     int result = 1;
     
     for (i=0; i<PMAX; i++)
     {
	  if ( DERIV_GRAD(x)[i] != 0 )
	  {
	       result = 0;
	       break;
	  }
     }

     return result;
}

#define HAS_NONZERO_DERIV(x) ( ! is_zero_vector (x) )

#else /* ndef AD_DERIV_TYPE_ARRAY == AD_DERIV_TYPE_POINTER */

extern int is_zero_vector Proto((void *));
#define HAS_NONZERO_DERIV(x) (! is_zero_vector(DERIV_GRAD_FOR_SAXPY(x)))

#endif /* ndef AD_DERIV_TYPE_ARRAY */

#define ZERO ((double) 0.0)
#define ONE_HALF ((double) 0.5)
#define ONE ((double) 1.0)
#define TWO ((double) 2.0)

#ifndef SILENT_EXCEPTIONS
#define PRINT_ERR_1(string,x) fprintf (stderr,string,x)
#define PRINT_ERR_2(string,x,y) fprintf (stderr,string,x,y)
#else
#define PRINT_ERR_1(string,x)
#define PRINT_ERR_2(string,x,y)
#endif


#ifdef AD_DERIV_TYPE_POINTER

#ifdef AD_SAXPY_INDIRECT

#define AXPY(outvec,alpha,invec)					      \
{									      \
   double aLpHa = alpha;						      \
   ADXPY1_FN ((void*)0,							      \
	   DERIV_GRAD_FOR_SAXPY(outvec),				      \
	   &aLpHa,							      \
	   DERIV_GRAD_FOR_SAXPY(invec));				      \
}

#define AXPY2(outvec,alpha,invec,alpha2,invec2)	\
{						\
   double aLpHa = alpha;			\
   double aLpHa2 = alpha2;			\
   ADXPY2_FN ((void*)0,				\
	   DERIV_GRAD_FOR_SAXPY(outvec),	\
	   &aLpHa,				\
	   DERIV_GRAD_FOR_SAXPY(invec),		\
	   &aLpHa2,				\
	   DERIV_GRAD_FOR_SAXPY(invec2));	\
}
 
#endif /* AD_SAXPY_INDIRECT */

#ifdef AD_SAXPY_DIRECT

#define AXPY(outvec,alpha,invec)					      \
{									      \
   double aLpHa = alpha;						      \
   ADXPY1_FN ((void*)0,							      \
	   DERIV_GRAD_FOR_SAXPY(outvec),				      \
	   aLpHa,							      \
	   DERIV_GRAD(invec));				                      \
}

#define AXPY2(outvec,alpha,invec,alpha2,invec2)	\
{						\
   double aLpHa = alpha;			\
   double aLpHa2 = alpha2;			\
   ADXPY2_FN ((void*)0,				\
	   DERIV_GRAD_FOR_SAXPY(outvec),	\
	   aLpHa,				\
	   DERIV_GRAD(invec),		        \
	   aLpHa2,				\
	   DERIV_GRAD(invec2));	                \
}

#endif /* AD_SAXPY_DIRECT axpy calls */

#ifdef AD_SAXPY_INLINE

#define AXPY(outvec,alpha,invec)					      \
{									      \
   double aLpHa = alpha;						      \
   int g_ad_i_;								      \
									      \
   if (DERIV_GRAD(invec) == 0 ) {					      \
	DERIV_GRAD(invec) = (void *) calloc (PMAX, sizeof (double));	      \
   }									      \
   if (DERIV_GRAD(outvec) == 0 ) {					      \
	DERIV_GRAD(outvec) = (void *) calloc (PMAX, sizeof (double));	      \
   }									      \
									      \
   for (g_ad_i_ = 0 ; g_ad_i_ < PMAX ; g_ad_i_ ++ )				      \
   {									      \
	((double *)DERIV_GRAD(outvec))[g_ad_i_]				      \
	     = aLpHa * ((double *)DERIV_GRAD(invec))[g_ad_i_];		      \
   }									      \
}

#define AXPY2(outvec,alpha,invec,alpha2,invec2)				      \
{									      \
   double aLpHa = alpha;						      \
   double aLpHa2 = alpha2;						      \
   int g_ad_i_;								      \
									      \
   if (DERIV_GRAD(invec) == 0 ) {					      \
	DERIV_GRAD(invec) = (void *) calloc (PMAX, sizeof (double));	      \
   }									      \
   if (DERIV_GRAD(invec2) == 0 ) {					      \
	DERIV_GRAD(invec2) = (void *) calloc (PMAX, sizeof (double));	      \
   }									      \
   if (DERIV_GRAD(outvec) == 0 ) {					      \
	DERIV_GRAD(outvec) = (void *) calloc (PMAX, sizeof (double));	      \
   }									      \
									      \
   for (g_ad_i_ = 0 ; g_ad_i_ < PMAX ; g_ad_i_ ++ )				      \
   {									      \
	((double *)DERIV_GRAD(outvec))[g_ad_i_]				      \
	     = aLpHa * ((double *)DERIV_GRAD(invec))[g_ad_i_]		      \
	     + aLpHa2 * ((double *)DERIV_GRAD(invec2))[g_ad_i_];		      \
   }									      \
}

#endif /* AD_SAXPY_INLINE */


#else /* not AD_DERIV_TYPE_POINTER == AD_DERIV_TYPE_ARRAY axpy calls */


#define AXPY(outvec,alpha,invec)					      \
{									      \
   double aLpHa = alpha;						      \
   int g_ad_i_;								      \
   for (g_ad_i_ = 0 ; g_ad_i_ < PMAX ; g_ad_i_ ++ )				      \
   {									      \
	DERIV_GRAD(outvec)[g_ad_i_] = aLpHa * DERIV_GRAD(invec)[g_ad_i_];	      \
   }									      \
}

#define AXPY2(outvec,alpha,invec,alpha2,invec2)				      \
do {									      \
   double aLpHa = alpha;						      \
   double aLpHa2 = alpha2;						      \
   int g_ad_i_;								      \
   for (g_ad_i_ = 0 ; g_ad_i_ < PMAX ; g_ad_i_ ++ )				      \
   {									      \
	DERIV_GRAD(outvec)[g_ad_i_] = aLpHa * DERIV_GRAD(invec)[g_ad_i_]	      \
	                            + aLpHa2 * DERIV_GRAD(invec2)[g_ad_i_];      \
   }									      \
} while (0)

#endif /* AD_DERIV_TYPE_ARRAY axpy calls */


#undef __FUNC__
#define __FUNC__ "g_ad_log" 
void 
g_ad_log (DERIV_TYPE *result, DERIV_TYPE arg1)
{
     DERIV_VAL(*result) = log (DERIV_VAL(arg1));

     AXPY(*result, ONE / DERIV_VAL(arg1), arg1);
}


#undef __FUNC__
#define __FUNC__ "g_ad_sqrt" 
void 
g_ad_sqrt (DERIV_TYPE *result, DERIV_TYPE arg1)
{
     double x = DERIV_VAL(arg1);
     double fx;
     double z = sqrt(x);
     
     DERIV_VAL(*result) = z;

     if ( x > ZERO )
     {
	  fx = ONE / (TWO * z);
     }
     else
     {
	  if ( HAS_NONZERO_DERIV (arg1) )
	  {
	       PRINT_ERR_1 ("g_ad_sqrt: Exception at sqrt(%f)\n",x);
	  }

	  fx = ZERO;
     }

     AXPY(*result, fx, arg1);
}

#undef __FUNC__
#define __FUNC__ "g_ad_cos" 
void
g_ad_cos (DERIV_TYPE *result, DERIV_TYPE arg1)
{
     double x = DERIV_VAL(arg1);
     double z = cos(x);
     double fx = -sin(x);

     DERIV_VAL(*result) = z;

     AXPY (*result, fx,arg1);
}


#undef __FUNC__
#define __FUNC__ "g_ad_sin" 
void
g_ad_sin (DERIV_TYPE *result, DERIV_TYPE arg1)
{
     double x = DERIV_VAL(arg1);
     double z = sin(x);
     double fx = cos(x);

     DERIV_VAL(*result) = z;

     AXPY (*result, fx,arg1);
}

#undef __FUNC__
#define __FUNC__ "g_ad_exp" 
void
g_ad_exp (DERIV_TYPE *result, DERIV_TYPE arg1)
{
     double x = DERIV_VAL(arg1);
     double z = exp(x);
     double fx = z;

     DERIV_VAL(*result) = z;

     AXPY (*result, fx,arg1);
}


#undef __FUNC__
#define __FUNC__ "g_ad_tan" 
void
g_ad_tan (DERIV_TYPE *result, DERIV_TYPE arg1)
{
     double x = DERIV_VAL(arg1);
     double z = tan(x);
     double fx = ONE + z*z;

     DERIV_VAL(*result) = z;

     AXPY (*result, fx,arg1);
}


/*
  This is not super-efficient, but I'm not worried about that right now.
  */

#undef __FUNC__
#define __FUNC__ "g_ad_pow" 
void
g_ad_pow (DERIV_TYPE *result, DERIV_TYPE arg1, DERIV_TYPE arg2)
{
     double x = DERIV_VAL(arg1);
     double y = DERIV_VAL(arg2);
     double z;
     double fx, fy;

     z = pow(x,y);
     
     DERIV_VAL(*result) = z;

     if ( x != ZERO )
     {
	  fx = y * pow(x,y-1);
     }
     else if ( x == ZERO && (ZERO < y) && (y < ONE) )
     {
	  fx = ZERO;
     }
     else
     {
	  if ( HAS_NONZERO_DERIV (arg1) 
	      || HAS_NONZERO_DERIV(arg2) )
	  {
	       PRINT_ERR_2 ("g_ad_power: Exception at pow(%f,%f)\n",x,y);
	  }

	  fx = ZERO;
     }

     if ( x > 0 )
     {
	  fy = log(x) * z;
     }
     else if ( x == ZERO  && y != ZERO )
     {
	  fy = ZERO;
     }
     else
     {
	  /* No need for an exception if y is an integer, even if x <= 0 */
	  if ( y != (int) y)
	  {
	       /* Maybe this should just check `y', but I do not
		  want to think about it right now. */

	       if ( HAS_NONZERO_DERIV(arg1)
		   || HAS_NONZERO_DERIV(arg2))
	       {
		    PRINT_ERR_2 ("g_ad_power: Exception at pow(%f,%f)\n",x,y);
	       }
	  }

	  /* But if y is an integer, the partial with respect to the
	     exponent had better be zero! */
	  fy = ZERO;
     }

     AXPY2 (*result, fx,arg1,fy,arg2 );
}


#undef __FUNC__
#define __FUNC__ "g_ad_fmin" 
void
g_ad_fmin (DERIV_TYPE *result, DERIV_TYPE arg1, DERIV_TYPE arg2)
{
     double x = DERIV_VAL(arg1);
     double y = DERIV_VAL(arg2);
     double z;
     double fx, fy;

     if ( x < y )
     {
	  fx = ONE;
	  fy = ZERO;
	  z = x;
     }
     else if ( x > y )
     {
	  fx = ZERO;
	  fy = ONE;
	  z = y;
     }
     else
     {
	  if ( HAS_NONZERO_DERIV(arg1)
	      || HAS_NONZERO_DERIV(arg2))
	  {
	       PRINT_ERR_2 ("g_ad_fmin: Exception at min(%f,%f)\n",x,y);
	  }

	  fx = ONE_HALF;
	  fy = ONE_HALF;
	  z = x;
     }

     DERIV_VAL(*result) = z;
     AXPY2 (*result, fx, arg1, fy, arg2);
}


#undef __FUNC__
#define __FUNC__ "g_ad_fmax" 
void
g_ad_fmax (DERIV_TYPE *result, DERIV_TYPE arg1, DERIV_TYPE arg2)
{
     double x = DERIV_VAL(arg1);
     double y = DERIV_VAL(arg2);
     double z;
     double fx, fy;

     if ( x > y )
     {
	  fx = ONE;
	  fy = ZERO;
	  z = x;
     }
     else if ( x < y )
     {
	  fx = ZERO;
	  fy = ONE;
	  z = y;
     }
     else
     {
	  if ( HAS_NONZERO_DERIV(arg1)
	      || HAS_NONZERO_DERIV(arg2))
	  {
	       PRINT_ERR_2 ("g_ad_fmax: Exception at min(%f,%f)\n",x,y);
	  }
	  fx = ONE_HALF;
	  fy = ONE_HALF;
	  z = x;
     }

     DERIV_VAL(*result) = z; 
     AXPY2 (*result, fx, arg1, fy, arg2);
}


#undef __FUNC__
#define __FUNC__ "g_ad_fabs" 
void
g_ad_fabs (DERIV_TYPE *result, DERIV_TYPE arg1)
{
     double x = DERIV_VAL(arg1);
     double z;
     double fx;
     
#ifdef DEBUG
     int never = 0; 
     if ( never ) 
     {
	  user_adprint(arg1);
     }
#endif

     if ( x < ZERO )
     {
	  fx = -ONE;
	  z = -x;
     }
     else if ( x > ZERO )
     {
	  fx = ONE;
	  z = x;
     }
     else
     {
	  if ( HAS_NONZERO_DERIV(arg1) )
	  {
	       PRINT_ERR_1 ("g_ad_fabs: Exception at abs(%f)\n",x);
	  }
#ifndef BIG_FABS_PARTIAL
	  fx = ZERO;
#else
	  fx = 10e+40;
#endif
	  z = x;
     }

     DERIV_VAL(*result) = z;
     AXPY (*result, fx, arg1);
}


#undef __FUNC__
#define __FUNC__ "g_ad_sinh" 
void
g_ad_sinh (DERIV_TYPE *result, DERIV_TYPE arg1)
{
     double x = DERIV_VAL(arg1);
     double z = sinh(x);
     double fx = cosh(x);

     DERIV_VAL(*result) = z;

     AXPY (*result, fx,arg1);
}

#undef __FUNC__
#define __FUNC__ "g_ad_cosh" 
void
g_ad_cosh (DERIV_TYPE *result, DERIV_TYPE arg1)
{
     double x = DERIV_VAL(arg1);
     double z = cosh(x);
     double fx = sinh(x);

     DERIV_VAL(*result) = z;

     AXPY (*result, fx,arg1);
}


#undef __FUNC__
#define __FUNC__ "g_ad_tanh" 
void
g_ad_tanh (DERIV_TYPE *result, DERIV_TYPE arg1)
{
     double x = DERIV_VAL(arg1);
     double z = tanh(x);
     double fx = (ONE - z)*(ONE + z);

     DERIV_VAL(*result) = z;

     AXPY (*result, fx,arg1);
}


#undef __FUNC__
#define __FUNC__ "g_ad_acos" 
void
g_ad_acos (DERIV_TYPE *result, DERIV_TYPE arg1)
{
     double x = DERIV_VAL(arg1);
     double z = acos(x);
     double fx;

     if ( fabs (x) <= ONE ) {
	  fx = - ONE / sqrt ( (ONE + x)*(ONE - x) );
     }
     else {
	  if ( HAS_NONZERO_DERIV(arg1) )
	  {
	       PRINT_ERR_1 ("g_ad_acos: Exception at acos(%f)\n",x);
	  }
	  fx = ZERO;
     }

     DERIV_VAL(*result) = z;

     AXPY (*result, fx,arg1);
}


#undef __FUNC__
#define __FUNC__ "g_ad_asin" 
void
g_ad_asin (DERIV_TYPE *result, DERIV_TYPE arg1)
{
     double x = DERIV_VAL(arg1);
     double z = asin(x);
     double fx;

     if ( fabs (x) <= ONE ) {
	  fx = ONE / sqrt ( (ONE + x)*(ONE - x) );
     }
     else {
	  if (HAS_NONZERO_DERIV(arg1)) 
	  {
	       PRINT_ERR_1 ("g_ad_asin: Exception at asin(%f)\n",x);
	  }
	  fx = ZERO;
     }

     DERIV_VAL(*result) = z;

     AXPY (*result, fx,arg1);
}


#undef __FUNC__
#define __FUNC__ "g_ad_atan" 
void
g_ad_atan (DERIV_TYPE *result, DERIV_TYPE arg1)
{
     double x = DERIV_VAL(arg1);
     double z = atan(x);
     double fx = ONE / ( ONE + x * x);

     DERIV_VAL(*result) = z;

     AXPY (*result, fx,arg1);
}


#undef __FUNC__
#define __FUNC__ "g_ad_log10" 
void
g_ad_log10 (DERIV_TYPE *result, DERIV_TYPE arg1)
{
     const double LOGTEN = log ((double) 10.0);

     double x = DERIV_VAL(arg1);
     /*     double z = log10(x); */
     double fx = ONE / ( x * LOGTEN);

     AXPY (*result, fx,arg1);
}

#undef __FUNC__
#define __FUNC__ "g_ad_atan2" 
void
g_ad_atan2 (DERIV_TYPE *result, DERIV_TYPE arg1, DERIV_TYPE arg2)
{
     double x = DERIV_VAL(arg1);
     double y = DERIV_VAL(arg2);
     /*     double z = atan2(x,y); */

     double fx, fy, scratch;

     scratch = ONE / ( x*x + y*y );
     fx = y * scratch;
     fy = - ( x * scratch );
     AXPY2 (*result, fx, arg1, fy, arg2);
}

#undef __FUNC__
#define __FUNC__ "g_ad_fmod" 
void
g_ad_fmod (DERIV_TYPE *result, DERIV_TYPE arg1, DERIV_TYPE arg2)
{
     double x = DERIV_VAL(arg1);
     double y = DERIV_VAL(arg2);
     double z = fmod(x,y);

     double fx;
     double fy;
     fx = 1.0;
     if ( z == 0.0 ) {
          PRINT_ERR_2 ("g_ad_fmod: Exception at fmod(%f,%f)\n",x,y);
           fy = 0.;
     }
     else {
          fy = - (int) (x / y);
     }

     AXPY2 (*result, fx, arg1, fy, arg2);
}
