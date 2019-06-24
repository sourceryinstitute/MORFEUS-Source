/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#ifndef SM_DERIV_H_
#define SM_DERIV_H_ 1

/*added by DPS */
#define AD_DERIV_TYPE_ARRAY

#define PMAX 3

#define AD_COMPUTATION 1

#ifdef AD_DERIV_TYPE_ARRAY
#include <stdio.h>
#endif

#ifdef __sun
#ifdef USE_MACHINE_SPECIFIC_NAMES /* Use SparseLinC */
#define AD_APPEND_UNDERSCORE 1
#endif /* USE_MACHINE_SPECIFIC_NAMES */
#endif /* is sun */

#ifdef AD_APPEND_UNDERSCORE

#define ADINIT_FN adinit_
#define ADEXTRACT_FN adextract_
#define ADALLOCATE_FN adallocate_
#define ADCOPY_FN adcopy_
#define ADFREE_FN adfree_vectors_buckets_
#define ADZERO_FN adzero_

/* These are used in intrinsic.c */
#define ADXPY1_FN adxpy1_
#define ADXPY2_FN adxpy2_

#else /* ndef AD_APPEND_UNDERSCORE == No Underscore */

#define ADINIT_FN adinit
#define ADEXTRACT_FN adextract
#define ADALLOCATE_FN adallocate
#define ADCOPY_FN adcopy
#define ADFREE_FN adfree_vectors_buckets
#define ADZERO_FN adzero

/* These are used in intrinsic.c */
#define ADXPY1_FN adxpy1
#define ADXPY2_FN adxpy2

#endif /* ndef AD_APPEND_UNDERSCORE == No Underscore */

#ifdef AD_DERIV_TYPE_ARRAY

typedef struct DERIV_TYPE {
     double val;
     double grad[PMAX];
} DERIV_TYPE;

#endif

#ifdef AD_DERIV_TYPE_POINTER

#ifndef AD_SAXPY_INLINE

typedef struct DERIV_TYPE {
     double val;
     void *grad;
} DERIV_TYPE;

#else /* AD_SAXPY_INLINE */

typedef struct DERIV_TYPE {
     double val;
     double *grad;
} DERIV_TYPE;

#endif /* AD_SAXPY_INLINE */

#endif /* AD_DERIV_TYPE_POINTER */


#define DERIV_VAL(x) ((x).val)
#define DERIV_GRAD(x) ((x).grad)
#define DERIV_IN_FN_CALL(x) x

#define g_ad__FLOAT_INITIALIZER_(x) { x, 0.0 }
#define _FLOAT_INITIALIZER_(x) { x, 0.0 }

#define COPY_VAL_TO_DERIV_ARRAY(x, y, n) \
    { int i; for (i = 0; i < n; i++) (x)[i].val = y[i]; }

#ifdef AD_DERIV_TYPE_POINTER
#define DERIV_GRAD_FOR_SAXPY(x) (&(DERIV_GRAD(x)))
#endif

#ifdef AD_DERIV_TYPE_ARRAY
#define ADINIT(value,pos,vector)                                              \
{                                                                          \
     double vAlUe = value;                                                    \
     int pOsItIoN = pos;                                                      \
     if ( (pOsItIoN <= 0) || (pOsItIoN > PMAX) )                              \
     {                                                                        \
          fprintf (stderr,"ADINIT: Index %d required to be positive, less than PMAX = %d\n",pOsItIoN,PMAX); \
          abort();                                                            \
     }                                                                        \
     DERIV_GRAD(vector)[pOsItIoN - 1] = vAlUe;                                \
}

#define ADEXTRACT(src,dest,max_indep)                                         \
do {                                                                          \
     int adIloop;                                                             \
     int mAX_inDEp = (max_indep);                                             \
     if ( (mAX_inDEp) < PMAX )                                                \
     {                                                                        \
          fprintf(stderr,"ADEXTRACT: max_indep (%d) < PMAX (%d)\n",           \
                  mAX_inDEp, PMAX);                                           \
          abort();                                                            \
     }                                                                        \
     for (adIloop = 0 ; adIloop < PMAX ; adIloop ++ )                         \
     {                                                                        \
          (dest)[adIloop] = DERIV_GRAD(src)[adIloop];                         \
     }                                                                        \
} while (0)

#define ADALLOCATE(size)                                                      \
do {                                                                          \
     int sIZe = (size);                                                       \
                                                                              \
     if ((sIZe) != PMAX)                                                      \
     {                                                                        \
          fprintf (stderr,"ADALLOCATE: Tried to allocate PMAX = %d rather than the #defined PMAX = %d in array mode.\n",sIZe,PMAX); \
          abort();                                                            \
     }                                                                        \
} while (0)

#define ADCOPY(dest,src)                                                      \
do {                                                                          \
     int adIloop;                                                             \
     for (adIloop = 0 ; adIloop < PMAX ; adIloop ++ )                         \
     {                                                                        \
          DERIV_GRAD(dest)[adIloop] = DERIV_GRAD(src)[adIloop];               \
     }                                                                        \
} while (0)

#define ADFREE(vector) /* No-op */
#define ADZERO(vector)                                                        \
do {                                                                          \
     int adIloop;                                                             \
     for (adIloop = 0 ; adIloop < PMAX ; adIloop ++ )                         \
     {                                                                        \
          DERIV_GRAD(vector)[adIloop] = 0.0;                                  \
     }                                                                        \
} while(0)

#define AD_AUTHORIZE_ONE(deriv_object) /* Do nothing */

#define AD_DEAUTHORIZE_ONE(deriv_object) /* Do nothing */

#define AD_LOCALIZE_ONE(deriv_object) /* Do nothing */
#elif AD_DERIV_TYPE_POINTER
#define ADINIT(value,pos,vector)                                              \
do {                                                                          \
     double vAlUe = value;                                                    \
     int pOsItIoN = pos;                                                      \
     int jUsToNe = 1;                                                         \
     ADINIT_FN ( &vAlUe, &pOsItIoN, &jUsToNe, DERIV_GRAD_FOR_SAXPY(vector));  \
} while(0)

#define ADEXTRACT(source,dest,max_indep)                                      \
do {                                                                          \
     int mAx_InDeP = max_indep;                                               \
                                                                              \
    ADEXTRACT_FN (DERIV_GRAD_FOR_SAXPY(source), dest, &mAx_InDeP);            \
} while (0)

#define ADALLOCATE(size)                                                      \
do {                                                                          \
    int sIze = (size) ;                                                       \
    ADALLOCATE_FN (&sIze);                                                    \
} while (0)

#define ADCOPY(dest,src) \
    ADCOPY_FN ((int *)0,DERIV_GRAD_FOR_SAXPY(dest),DERIV_GRAD_FOR_SAXPY(src))

#define ADFREE(vector)                                                        \
do {                                                                          \
     if ( DERIV_GRAD(vector) != 0 )                                           \
     {                                                                        \
          ADFREE_FN(DERIV_GRAD_FOR_SAXPY(vector));                            \
          DERIV_GRAD(vector) = 0;                                             \
     }                                                                        \
} while (0)

#ifdef AD_SAXPY_INDIRECT
#define ADZERO(vector) ADZERO_FN ((int *)0, DERIV_GRAD_FOR_SAXPY(vector))
#endif

#ifdef AD_SAXPY_DIRECT /* -- inlined stuff for more speed */

#define ADZERO(vector)                                                        \
do {                                                                          \
     if ( DERIV_GRAD(vector) != 0 )                                           \
     {                                                                        \
          ADZERO_FN((int *)0, DERIV_GRAD_FOR_SAXPY(vector));                  \
          DERIV_GRAD(vector) = 0;                                             \
     }                                                                        \
} while (0)

#endif /* AD_SAXPY_DIRECT */

#ifdef AD_SAXPY_INLINE

#define ADZERO(vector)                                                        \
do {                                                                          \
     int adIloop;                                                             \
     if (DERIV_GRAD(vector) == 0) break;                                      \
                                                                              \
     for (adIloop = 0 ; adIloop < PMAX ; adIloop ++ )                         \
     {                                                                        \
          ((double *)DERIV_GRAD(vector))[adIloop] = 0.0;                      \
     }                                                                        \
} while(0)

#endif /* AD_SAXPY_INLINE */

#ifdef AD_SAXPY_INLINE

#define AD_AUTHORIZE_ONE(deriv_object) \
     DERIV_GRAD(deriv_object) = calloc (PMAX, sizeof(struct DERIV_TYPE))

#else /* Normal behaviour */

#define AD_AUTHORIZE_ONE(deriv_object) \
     DERIV_GRAD(deriv_object) = 0

#endif /* not AD_SAXPY_INLINE */

#define AD_DEAUTHORIZE_ONE(deriv_object) \
     ADFREE(deriv_object)

#define AD_LOCALIZE_ONE(deriv_object)                                         \
do {                                                                          \
     DERIV_TYPE temPVAR = deriv_object;                                       \
     AD_AUTHORIZE_ONE(temPVAR);                                               \
     ADCOPY(temPVAR,(deriv_object));                                          \
     deriv_object = temPVAR;                                                  \
} while(0)

#endif
#ifdef AD_DERIV_TYPE_ARRAY
#undef AD_DO_LOCALIZE_GRADIENTS
#endif

#ifdef AD_DERIV_TYPE_POINTER
#define AD_DO_LOCALIZE_GRADIENTS 1
#endif
#ifdef AD_DERIV_TYPE_ARRAY
#undef AD_DO_AUTHORIZE_ARRAYS
#endif

#ifdef AD_DERIV_TYPE_POINTER
#define AD_DO_AUTHORIZE_ARRAYS 1
#endif

#define ADZERO_ARRAY(x, n) \
        { int i; for (i = 0; i < n; i++) ADZERO((x)[i]); }

#define ADINIT_ARRAY(V, POS, X, N) \
        { int i; for (i = 0; i < N; i++) ADINIT(V, POS+i, X[i]); }

#endif /* ndef _DERIV_H_ */
