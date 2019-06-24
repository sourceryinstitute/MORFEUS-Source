/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#if !defined(__SM_ERROR_H)
#define __SM_ERROR_H

/*
   Defines the directory where the compiled source is located; used
   in printing error messages. Each makefile has an entry
   LOCDIR     =  thedirectory
   and bmake/common includes in CFLAGS -D__SDIR__='"${LOCDIR}"'
   which is a flag passed to the compilers.
*/
#if !defined(__SDIR__)
#define __SDIR__ "unknowndirectory/"
#endif

/*
   Defines the function where the compiled source is located; used
   in printing error messages.
*/
#if !defined(__FUNC__)
#define __FUNC__ "unknownfunction"
#endif

/* Error Routines - Accessed through OPTMS_CHKERR(ierr) */
int  SMerror(int line,char *func,char* file,char *dir,int n,int p,char *mess);

/*
     These are the generic error codes. These error codes are used
     many different places in the PETSc source code.

*/
#define OPTMS_PRINT_STACK -1

#define OPTMS_MEM_ERR             55   /* unable to allocate the requested memory */
#define OPTMS_NULL_ERR            56   /* null data pointer */
#define OPTMS_INPUT_ERR           57   /* something is wrong with input to function */
#define OPTMS_INIT_ERR            58   /* data structure not initialized */
#define OPTMS_FILE_OPEN_ERR       59   /* unable to open file */
#define OPTMS_FREE_ERR            60   /* unable to free memory */
#define OPTMS_INVALID_MESH_ERR    61   /* unable to free memory */
#define OPTMS_DIVIDE_BY_ZERO_ERR  62   /* division by zero */
#define OPTMS_DATA_ERR            63   /* incorrect data */

#define OPTMS_SETERR(n,p,s) {return SMerror(__LINE__,__FUNC__,__FILE__,__SDIR__,n,p,s);}
#define OPTMS_CHKERR(n)     {if (n) OPTMS_SETERR(n,0,(char *)0);}

#ifndef OPTMS_CHECK_NULL
#define OPTMS_CHECK_NULL(a) \
{ \
   if (a==NULL) OPTMS_CHKERR(OPTMS_NULL_ERR); \
}
#endif

#endif
