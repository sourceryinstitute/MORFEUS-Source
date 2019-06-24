!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under 
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
!
!    NEMO - Numerical Engine (for) Multiphysics Operators
! Copyright (c) 2007, Stefano Toninel
!                     Gian Marco Bianchi  University of Bologna
!              David P. Schmidt    University of Massachusetts - Amherst
!              Salvatore Filippone University of Rome Tor Vergata
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification,
! are permitted provided that the following conditions are met:
!
!     1. Redistributions of source code must retain the above copyright notice,
!        this list of conditions and the following disclaimer.
!     2. Redistributions in binary form must reproduce the above copyright notice,
!        this list of conditions and the following disclaimer in the documentation
!        and/or other materials provided with the distribution.
!     3. Neither the name of the NEMO project nor the names of its contributors
!        may be used to endorse or promote products derived from this software
!        without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!---------------------------------------------------------------------------------
!
! $Id: class_connectivity.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
!****c* NEMO/Connectivity
! NAME
!    Connectivity
! DESCRIPTION
!    Handles connectivity between general elements.  The general form is A2B, where
!    the data indicates how a given element of set B is connected to zero or more
!    elements of A.
!
! METHODS
!    ALLOC_CONN       constructor for CONNECTIVITY class.
!    FREE_CONN        destructor for CONNECTIVITY class.
!    BCAST_CONN       broadcasts CONNECTIVITY from P0 to other processes.
!    G2L_CONN         global to local reallocation for a A2B connectivity.
!    L2G_CONN         local to global reallocation for a A2B connectivity.
!    GET_ITH_CONN     given an element and a connectivity structure,
!                     returns a pointer for the connected elements.
!    GET_CONN_CSR     returns a copy of a connectivity structure in CSR format.
!    GET_DUAL_CONN    given "A to B" connectivity, returns "B to A" connectivity.
!    NEL_             returns the number of global or local (after G2L) items
!                     "A" total in a connectivity structure A2B.
!    NCONN_           returns the number of global or local (after G2L)
!                     connections "B" total in a connectivity structure A2B.
!    MAX_CONN         gets the maximum connectivity degree.
!
!    UNUSED_ELEMENTS  are there any unreferenced elements in the connectivity?
!    COUNT_REFERENCES counts references to elements in connectivity structures.
!*******
MODULE class_connectivity

    USE class_psblas

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: connectivity                   ! Class
    PUBLIC :: alloc_conn, free_conn          ! Constructor/Destructor
    PUBLIC :: bcast_conn, g2l_conn, l2g_conn ! Parallel ops.
    PUBLIC :: ASSIGNMENT(=)    ! Setters
    PUBLIC :: print_conn

    ! `public :: assignment(=)'  ==  `public :: copy_conn'

    ! Abstract data type for handling general connectivity (vertex to cell,
    ! cell to vertex, etc.) In general, A2B. The index of the LOOKUP array
    ! corresponds to the numbering of A and the contents of LOOKUP indicate
    ! where in the connectivity list (CONN) to find the required information.
    ! The contents of CONN correspond to the numbering of B. The last entry
    ! of LOOKUP points to one past the last CONN entry, for convenience.

    !****it* Connectivity/connectivity
    !  NAME
    !    connectivity
    !
    !  SOURCE
    TYPE connectivity
        PRIVATE
        INTEGER, ALLOCATABLE :: lookup(:)
        INTEGER, ALLOCATABLE :: conn(:)
        INTEGER :: nel_glob
        INTEGER :: nconn_glob
    CONTAINS
        PROCEDURE :: get_ith_conn, get_conn_csr, get_dual_conn    ! Getters
        PROCEDURE :: nel_, nconn_                                 ! Getters cont.
        PROCEDURE :: set_ith_conn                   ! Setter
        PROCEDURE :: max_conn                       ! Other
        PROCEDURE :: unused_elements
        PROCEDURE, PRIVATE :: nemo_connectivity_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_connectivity_sizeof
        PROCEDURE, PRIVATE :: g2l_conn_core
    END TYPE connectivity
    !*******

    ! Description of class members.
    ! LOOKUP:     each entry contains an index for lookups in CONN(:).
    ! CONN:       connectivity array, packed with indices of elements.
    ! NEL_GLOB:   the global number of elements (e.g. in vertex 2 cell,
    !             this is the global number of cells)
    ! NCONN_GLOB: global number of connectivities stored


    ! ----- Generic Interfaces -----

    ! Setter
    INTERFACE ASSIGNMENT(=)
        MODULE SUBROUTINE copy_conn(c2d,a2b)
            IMPLICIT NONE
            TYPE(connectivity), INTENT(OUT) :: c2d
            TYPE(connectivity), INTENT(IN) :: a2b
        END SUBROUTINE copy_conn
    END INTERFACE ASSIGNMENT(=)

    ! REMARK: the initialization of CONNECTIVITY objects is supposed to
    ! happen sequentially, involving all elements.

    ! Global 2 Local
    INTERFACE g2l_conn
         MODULE SUBROUTINE g2l_conn_1d(a2b,desc_a,is_a2a)
            IMPLICIT NONE
            TYPE(connectivity), INTENT(INOUT) :: a2b
            TYPE(psb_desc_type), INTENT(IN) :: desc_a
            LOGICAL, INTENT(IN), OPTIONAL :: is_a2a
        END SUBROUTINE g2l_conn_1d

        MODULE SUBROUTINE g2l_conn_2d(a2b,desc_a,desc_b)
            IMPLICIT NONE
            TYPE(connectivity), INTENT(INOUT) :: a2b
            TYPE(psb_desc_type), INTENT(IN) :: desc_a
            TYPE(psb_desc_type), INTENT(IN) :: desc_b
        END SUBROUTINE g2l_conn_2d
    END INTERFACE g2l_conn

    INTERFACE

        ELEMENTAL MODULE FUNCTION nemo_connectivity_sizeof(conn)
            USE psb_base_mod
            IMPLICIT NONE
            CLASS(connectivity), INTENT(IN) :: conn
            INTEGER(kind=nemo_int_long_)   :: nemo_connectivity_sizeof
        END FUNCTION nemo_connectivity_sizeof

        MODULE SUBROUTINE alloc_conn(a2b,nel,nconn,lb)
            IMPLICIT NONE
            TYPE(connectivity), INTENT(INOUT) :: a2b
            INTEGER, INTENT(IN) :: nel, nconn
            INTEGER, INTENT(IN), OPTIONAL :: lb
        END SUBROUTINE alloc_conn

    ! ----- Destructor -----

    !****m* Connectivity/free_conn
    !
    !  NAME
    !     free_conn
    !  USAGE
    !     call free_conn(a2b)
    !
    !  DESCRIPTION
    !     Destructor
    !
    !  INPUTS
    !     a2b   ::  type(connectivity)  The connectivity to be destroyed;
    !
    !
    !*******

        MODULE SUBROUTINE free_conn(a2b)
            IMPLICIT NONE
            TYPE(connectivity), INTENT(INOUT) :: a2b
        END SUBROUTINE free_conn

        MODULE SUBROUTINE print_conn(iout,a2b,head)
            IMPLICIT NONE
            INTEGER, INTENT(IN)            :: iout
            TYPE(connectivity), INTENT(IN) :: a2b
            CHARACTER(len=*), INTENT(IN), OPTIONAL :: head
        END SUBROUTINE print_conn

    ! ----- Parallel Operations -----

    !****f* Connectivity/bcast_conn
    !
    !  NAME
    !     bcast_conn
    !  USAGE
    !
    !     call bcast_conn(a2b)
    !
    !
    !  DESCRIPTION
    !     Distribute a copy of the connectivity from process 0 to everybody else.
    !
    !  INPUTS
    !     a2b   ::  type(connectivity)  The connectivity to be broadcast;
    !
    !
    !*******

        MODULE SUBROUTINE bcast_conn(a2b)
            IMPLICIT NONE
            TYPE(connectivity) :: a2b
        END SUBROUTINE bcast_conn

    ! ----- Global 2 Local -----

        MODULE SUBROUTINE g2l_conn_core(a2b,iglob_to_loc_a,iloc_to_glob_b)
            IMPLICIT NONE
            CLASS(connectivity), INTENT(INOUT) :: a2b
            INTEGER, INTENT(IN) :: iglob_to_loc_a(:)
            !INTEGER, INTENT(IN) :: iloc_to_glob_b(LBOUND(a2b%lookup,DIM=1):)
            INTEGER, INTENT(IN) :: iloc_to_glob_b(:)
            !! IP - 5/27/2019: This work-around fixes an ICE with Intel 18.0.5.274
            !! Note: The ICE only shows up when building the submodule that implements this procedure
        END SUBROUTINE g2l_conn_core

        MODULE SUBROUTINE l2g_conn(a2b_loc,a2b_glob,desc_a,desc_b)
            IMPLICIT NONE
            TYPE(connectivity), INTENT(IN) :: a2b_loc
            TYPE(connectivity), INTENT(OUT) :: a2b_glob
            TYPE(psb_desc_type), INTENT(IN) :: desc_a
            TYPE(psb_desc_type), INTENT(IN) :: desc_b
        END SUBROUTINE  l2g_conn

    ! ----- Getters -----

    !****m* Connectivity/get_ith_conn
    !
    !  NAME
    !     get_ith_conn
    !  USAGE
    !     call get_ith_conn(ith_conn,a2b,i)
    !
    !  DESCRIPTION
    !     Getter. Access the I-th row in the connectivity table.
    !     Returns a pointer to the requested row. Note: returning
    !     a pointer and not a copy, because it is used in tight loops;
    !     user should NEVER modify nor deallocate the returned pointer.
    !
    !  INPUTS
    !     a2b   ::  type(connectivity), intent(in), target
    !                                   The connectivity to be accessed;
    !     i     ::  integer             Request access to the I-th row;
    !  OUTPUT
    !     ith_conn(:) :: integer, pointer
    !                                   A pointer to the i-th row;
    !
    !*******

        MODULE SUBROUTINE get_ith_conn(a2b,ith_conn,i)
        ! Retrieves connectivity data for the B(i) element.
        !
        ! NOTE: this getter has been designed for high frequency access,
        ! therefore the class content is POINTED instead of being COPIED.
        ! WARNING! In the calling unit, don't use ITH_CONN on the LHS of
        ! assignment expressions, if you want to keep intact the original
        ! connectivity data.
            IMPLICIT NONE
            INTEGER, POINTER :: ith_conn(:)
            CLASS(connectivity), INTENT(IN), TARGET :: a2b
            INTEGER, INTENT(IN) :: i
        END SUBROUTINE get_ith_conn

    !****m* Connectivity/get_conn_csr
    !
    !  NAME
    !     get_conn_csr
    !  USAGE
    !     call get_conn_csr(a2b,lookup,conn)
    !
    !  DESCRIPTION
    !     Get a full copy of the connectivity table in CSR format.
    !     Being a copy the user may do as he/she pleases with it,
    !     and should remember to deallocate it.
    !
    !  INPUTS
    !     a2b   ::  type(connectivity), intent(in)
    !                                   The connectivity to be accessed;
    !  OUTPUT
    !   lookup(:) :: integer, allocatable
    !     conn(:) :: integer, allocatable
    !                                   A copy of the lookup table.
    !                                   The neighbours of item I are listed
    !                                   at CONN(lookup(i):lookup(i+1)-1)
    !
    !*******

        MODULE SUBROUTINE get_conn_csr(a2b,lookup,conn)
        ! Retrieves all connectivity data in CSR format.
        !
        ! NOTE: this getter has been designed in order to get all connectivity
        ! data in a single shot. The pointers LOOKUP and CONN receive a COPY of
        ! A2B%LOOKUP and A2B%CONN respectively.
        ! They can be manipulated in the calling unit without affecting the
        ! original class content.
            IMPLICIT NONE
            CLASS(connectivity), INTENT(IN) :: a2b
            INTEGER, ALLOCATABLE :: lookup(:), conn(:)
        END SUBROUTINE get_conn_csr

        MODULE SUBROUTINE get_dual_conn(a2b,b2a)
        ! The B2A dual CONNECTIVITY can be thought as a set of stacks which
        ! are eventually stuck together sequentially.
        ! There is one stack for every element of B which appears at least
        ! once in A2B%CONN. The I-th stack is built with the indices of
        ! the elements of A (A2B%LOOKUP superscripts) which are connected to
        ! the I-th element of B.
        ! A practical application is extracting the "face to vertex" conn.,
        ! starting from the "vertex to face" one.
        ! WARNING: we are not referring to a STACK according to the common
        ! vocabulary of the programming science, i.e. as a special container
        ! implemented with pointers etc.
            IMPLICIT NONE
            CLASS(connectivity), INTENT(IN) :: a2b
            TYPE(connectivity), INTENT(OUT) :: b2a
        END SUBROUTINE get_dual_conn

        MODULE FUNCTION nel_(a2b,gl)
        ! get number of elements (either local or global)
            IMPLICIT NONE 
            INTEGER :: nel_
            CLASS(connectivity), INTENT(IN) :: a2b
            CHARACTER(len=1), INTENT(IN), OPTIONAL :: gl
        END FUNCTION nel_

        MODULE FUNCTION nconn_(a2b,gl)
            IMPLICIT NONE
            INTEGER :: nconn_
            CLASS(connectivity), INTENT(IN) :: a2b
            CHARACTER(len=1), INTENT(IN), OPTIONAL :: gl
        END FUNCTION nconn_

    ! ----- Setters -----

        MODULE SUBROUTINE set_ith_conn(a2b,i,ith_conn)
            IMPLICIT NONE
            CLASS(connectivity), INTENT(INOUT) :: a2b
            INTEGER, INTENT(IN) :: i
            INTEGER, INTENT(IN) :: ith_conn(:)
        ! WARNING: if A2B is declared only as INTENT(OUT), it becomes UNDEFINED
        ! on entry to the procedure (See Metcalf). Not detected by Gfortran.
        END SUBROUTINE set_ith_conn

    ! ----- Other -----

        MODULE FUNCTION max_conn(a2b,lb)
        ! Computes the maximum connectivity degree
        IMPLICIT NONE
            INTEGER :: max_conn
            CLASS(connectivity), INTENT(IN) :: a2b
            INTEGER, INTENT(IN), OPTIONAL :: lb
        END FUNCTION max_conn

    ! ----- Integrity Checks -----

        MODULE FUNCTION unused_elements(a2b)
        ! For debugging connectivity changes
        ! Indicates if unused elements, B, occur in the A2B structure.
        ! This function will probably be deleted in the future, but for
        ! now it gives me a place for the debugger to break when something
        ! is wrong with the mesh.
            IMPLICIT NONE
            INTEGER :: unused_elements             ! > 0 indicates a problem
            CLASS(connectivity), INTENT(IN) :: a2b  ! connectivity structure
        END FUNCTION unused_elements

        MODULE FUNCTION count_references(a2b,nb)
        !! Counts the number of references to each B in the A2B structure
            IMPLICIT NONE
            TYPE(connectivity), INTENT(IN) :: a2b ! connectivity structure
            INTEGER, INTENT(IN)            :: nb  ! number of unique b elements
            INTEGER :: count_references(nb)
        END FUNCTION count_references

    END INTERFACE

END MODULE class_connectivity
