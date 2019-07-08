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
SUBMODULE(class_connectivity) class_connectivity_procedures

    USE class_psblas

    IMPLICIT NONE

CONTAINS

    ! ----- Constructor -----

    !****m* Connectivity/alloc_conn
    !
    !  NAME
    !     alloc_conn
    !  USAGE
    !     call alloc_conn(a2b,nel=nv,nconn=nedg)
    !
    !  DESCRIPTION
    !     Constructor.
    !     (re)Allocates a connectivity table.
    !
    !  INPUTS
    !     a2b   ::  type(connectivity)  The connectivity to be allocated;
    !                                   reused if already available.
    !     nel   ::  integer             Number of items (vertices) to be connected
    !     nconn ::  integer             Number of connections (edges)
    !
    !
    !*******
    MODULE PROCEDURE alloc_conn
        !
        INTEGER :: info(2)
        INTEGER :: lb_

        ! Sets lower bound LB_ of A2B%LOOKUP
        IF(PRESENT(lb)) THEN
            lb_ = lb  ! Optional
        ELSE
            lb_ = 1   ! Default
        END IF

        info = 0
        CALL psb_realloc(nel+1,a2b%lookup,info(1),lb=lb_)
        CALL psb_realloc(nconn,a2b%conn,info(2))
        IF(ANY(info /= 0)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! Initialization
        a2b%lookup(lb_) = 1
        a2b%nel_glob   = nel
        a2b%nconn_glob = nconn

100     FORMAT(' ERROR! Memory allocation failure in ALLOC_CONN')

    END PROCEDURE alloc_conn


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
    MODULE PROCEDURE free_conn
        !
        INTEGER :: info(2)

        info = 0
        IF (ALLOCATED(a2b%lookup)) DEALLOCATE(a2b%lookup,stat=info(1))
        IF (ALLOCATED(a2b%conn))   DEALLOCATE(a2b%conn,stat=info(2))
        IF(ANY(info /= 0)) THEN
            WRITE(*,100)
            WRITE(*,*) info
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Memory deallocation failure in FREE_CONN')

    END PROCEDURE free_conn


    MODULE PROCEDURE print_conn
        !
        INTEGER :: i

        IF (.NOT.ALLOCATED(a2b%lookup)) THEN
            WRITE(iout,*) 'Unallocated lookup entry '
            RETURN
        END IF
        IF (.NOT.ALLOCATED(a2b%conn)) THEN
            WRITE(iout,*) 'Unallocated conn entry '
            RETURN
        END IF
        WRITE(iout,*) 'Printing connectivity'
        IF (PRESENT(head)) WRITE(iout,*) 'Connectivity name: ',head
        DO i=LBOUND(a2b%lookup,1),UBOUND(a2b%lookup,1)-1
            WRITE(iout,'(a,3(i8,1x))')  'row    : ', i, a2b%lookup(i), a2b%lookup(i+1)
            WRITE(iout,'(a,50(i8,1x))') 'entries: ',a2b%conn(a2b%lookup(i):a2b%lookup(i+1)-1)
        END DO
    END PROCEDURE print_conn


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

    MODULE PROCEDURE bcast_conn
        !
        INTEGER :: mypnum, icontxt
        INTEGER :: lb, nel, nconn

        icontxt = icontxt_()
        mypnum  = mypnum_()

        ! A2B is supposed to be associated only in P0
        ! => Check and allocation on P > 0
        IF(mypnum == 0) THEN
            IF(    .NOT.ALLOCATED(a2b%lookup) .OR. &
                & .NOT.ALLOCATED(a2b%conn)) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            lb    = LBOUND(a2b%lookup,1)
            nel   = SIZE(a2b%lookup) - 1
            nconn = SIZE(a2b%conn)

            CALL psb_bcast(icontxt,lb)
            CALL psb_bcast(icontxt,nel)
            CALL psb_bcast(icontxt,nconn)

            CALL psb_bcast(icontxt,a2b%lookup)
            CALL psb_bcast(icontxt,a2b%conn)

        ELSE
            IF(    ALLOCATED(a2b%lookup) .OR. &
                & ALLOCATED(a2b%conn)) THEN
                WRITE(*,200)
                CALL abort_psblas
            END IF

            CALL psb_bcast(icontxt,lb)
            CALL psb_bcast(icontxt,nel)
            CALL psb_bcast(icontxt,nconn)

            CALL alloc_conn(a2b,nel,nconn,lb=lb)

            CALL psb_bcast(icontxt,a2b%lookup)
            CALL psb_bcast(icontxt,a2b%conn)

        END IF

100     FORMAT(' ERROR! Illegal Send: CONNECTIVITY not allocated')
200     FORMAT(' ERROR! Illegal Recv: CONNECTIVITY already allocated')

    END PROCEDURE bcast_conn

    MODULE PROCEDURE nemo_connectivity_sizeof
        USE psb_base_mod

        INTEGER(kind=nemo_int_long_)   :: val

        val = 2 * nemo_sizeof_int
        IF (ALLOCATED(conn%lookup))&
            & val = val + nemo_sizeof_int * SIZE(conn%lookup)
        IF (ALLOCATED(conn%conn))&
            &  val = val + nemo_sizeof_int * SIZE(conn%conn)
        nemo_connectivity_sizeof = val

    END PROCEDURE nemo_connectivity_sizeof

    ! ----- Global 2 Local -----

    MODULE PROCEDURE g2l_conn_1d
        ! 1D -> means "one decriptor"
        !
        LOGICAL :: is_a2a_
        INTEGER :: ib, info, lb, ub
        INTEGER, ALLOCATABLE :: iglob_to_loc_a(:)
        INTEGER, ALLOCATABLE :: iloc_to_glob_a(:)
        INTEGER, ALLOCATABLE :: iloc_to_glob_b(:)

        IF(PRESENT(is_a2a)) THEN
            is_a2a_ = is_a2a
        ELSE
            is_a2a_ = .FALSE. ! Default
        END IF

        CALL psb_get_glob_to_loc(desc_a,iglob_to_loc_a)

        IF(is_a2a_) THEN
            CALL psb_get_loc_to_glob(desc_a,iloc_to_glob_a)

            CALL a2b%g2l_conn_core(iglob_to_loc_a,iloc_to_glob_a)

            DEALLOCATE(iloc_to_glob_a)
        ELSE
            lb = LBOUND(a2b%lookup,1)     ! Lower bound
            ub = UBOUND(a2b%lookup,1) - 1 ! Upper bound
            ALLOCATE(iloc_to_glob_b(lb:ub),stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            iloc_to_glob_b = [ (ib, ib = lb, ub) ]

            CALL a2b%g2l_conn_core(iglob_to_loc_a,iloc_to_glob_b)

            DEALLOCATE(iloc_to_glob_b)
        END IF

        DEALLOCATE(iglob_to_loc_a)

100     FORMAT(' ERROR! Memory allocation failure in G2L_CONN_1D')

    END PROCEDURE g2l_conn_1d


    MODULE PROCEDURE g2l_conn_2d
        ! 2D -> means "two decriptors"
        !
        INTEGER, ALLOCATABLE :: iglob_to_loc_a(:)
        INTEGER, ALLOCATABLE :: iloc_to_glob_b(:)

        CALL psb_get_glob_to_loc(desc_a,iglob_to_loc_a)
        CALL psb_get_loc_to_glob(desc_b,iloc_to_glob_b)

        CALL a2b%g2l_conn_core(iglob_to_loc_a,iloc_to_glob_b)

        DEALLOCATE(iglob_to_loc_a,iloc_to_glob_b)

    END PROCEDURE g2l_conn_2d


    MODULE PROCEDURE g2l_conn_core
        !
        INTEGER, PARAMETER :: maxconn = 8
        INTEGER :: i, i1, i2, j
        INTEGER :: ifirst, ilast
        INTEGER :: ia_glob, ia_loc, ib_glob, ib_loc
        INTEGER :: nel_glob, nel_loc, nconn_glob, nconn_loc
        TYPE(connectivity) :: work
        !
        ! WARNING: MAXCONN = maximum number of a2b connectivities. Set = 8.
        ! It must be increased in case of POLYHEDRAL cells.

        ! Saves NEL_GLOB and NCONN_GLOB
        nel_glob   = a2b%nel_glob
        nconn_glob = a2b%nconn_glob

        nel_loc  = SIZE(iloc_to_glob_b)
        ifirst = LBOUND(a2b%lookup,1)
        ilast  = ifirst + nel_loc - 1

        CALL alloc_conn(work,nel=nel_loc,nconn=nconn_glob,lb=ifirst)

        j = 1; work%lookup(ifirst) = 1
        DO ib_loc = ifirst, ilast
            ASSOCIATE (ib_loc_tmp => ib_loc + (LBOUND(iloc_to_glob_b,DIM=1) - ifirst))
                ib_glob = iloc_to_glob_b(ib_loc_tmp)
            END ASSOCIATE
            i1 = a2b%lookup(ib_glob)
            i2 = a2b%lookup(ib_glob+1) - 1
            DO i = i1, i2
                ia_glob = a2b%conn(i)
                ia_loc  = iglob_to_loc_a(ia_glob)
                IF(ia_loc > 0) THEN
                    work%conn(j) = ia_loc
                    j = j + 1
                END IF
            END DO
            work%lookup(ib_loc+1) = j
        END DO

        CALL free_conn(a2b)

        nconn_loc = work%lookup(ilast+1) - 1
        CALL alloc_conn(a2b,nel=nel_loc,nconn=nconn_loc,lb=ifirst)

        a2b%lookup = work%lookup(ifirst:ilast+1)
        a2b%conn = work%conn(1:nconn_loc)

        ! Restores original NEL_GLOB and NCONN_GLOB
        a2b%nel_glob   = nel_glob
        a2b%nconn_glob = nconn_glob

        CALL free_conn(work)

    END PROCEDURE g2l_conn_core


    MODULE PROCEDURE l2g_conn
        ! WARNING! The global results is allocated only on P0. After its usage
        ! it must be deallocated in the calling unit by means of the statement:
        ! "if(associated(glob_res)) deallocate(glob_res)"

        !
        INTEGER :: err_act, icontxt, info
        INTEGER :: maxconn, nconn_glob, nel_glob, nel_loc
        INTEGER :: i, j, k, k1, k2, n
        INTEGER, ALLOCATABLE :: ia_loc_to_glob(:)
        INTEGER, ALLOCATABLE :: ia2b(:)
        INTEGER, ALLOCATABLE :: itab_glob(:,:)
        INTEGER, ALLOCATABLE :: itab_loc(:,:)
        INTEGER, ALLOCATABLE :: kconn_glob(:)
        INTEGER, ALLOCATABLE :: kconn_loc(:)


        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)

        icontxt = icontxt_()

        nconn_glob = a2b_loc%nconn_glob
        nel_glob = a2b_loc%nel_glob
        nel_loc = SIZE(a2b_loc%lookup) - 1

        maxconn = a2b_loc%max_conn()            ! Local maximum
        CALL psb_amx(icontxt,maxconn)          ! Global maximum


        ! Allocation of local objects
        CALL psb_geall(itab_loc,desc_b,info,maxconn)
        CALL psb_check_error(info,'l2g_conn','psb_geall',icontxt)

        CALL psb_geall(kconn_loc,desc_b,info,maxconn)
        CALL psb_check_error(info,'l2g_conn','psb_geall',icontxt)

        ! Allocation of global objects
        ALLOCATE(itab_glob(nel_glob,maxconn),kconn_glob(nel_glob),&
            &   ia2b(maxconn),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        CALL psb_get_loc_to_glob(desc_a,ia_loc_to_glob)


        DO i = 1, nel_loc
            k1 = a2b_loc%lookup(i)
            k2 = a2b_loc%lookup(i+1)

            n = k2 - k1
            kconn_loc(i) = n

            k = k1 - 1
            DO j = 1, n
                k = k + 1
                itab_loc(i,j) = ia_loc_to_glob(a2b_loc%conn(k))
            END DO
        END DO



        CALL psb_gather(kconn_glob,kconn_loc,desc_b,info,root=0)
        CALL psb_check_error(info,'l2g_conn','psb_igatherv',icontxt)

        CALL psb_gather(itab_glob,itab_loc,desc_b,info,root=0)
        CALL psb_check_error(info,'l2g_conn','psb_igatherm',icontxt)


        IF(mypnum_() == 0) THEN
            CALL free_conn(a2b_glob)
            CALL alloc_conn(a2b_glob,nel=nel_glob,nconn=nconn_glob)

            DO i = 1, nel_glob
                n = kconn_glob(i)
                DO j = 1, n
                    ia2b(j) = itab_glob(i,j)
                END DO
                CALL a2b_glob%set_ith_conn(i,ia2b(1:n))
            END DO
        END IF


        ! Frees memory
        DEALLOCATE(ia2b)
        DEALLOCATE(ia_loc_to_glob)
        DEALLOCATE(itab_glob,kconn_glob)

        CALL psb_gefree(kconn_loc,desc_b,info)
        CALL psb_check_error(info,'l2g_conn','psb_gefree',icontxt)

        CALL psb_gefree(itab_loc,desc_b,info)
        CALL psb_check_error(info,'l2g_conn','psb_gefree',icontxt)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Memory allocation failure in L2G_CONN')

    END PROCEDURE l2g_conn


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

    MODULE PROCEDURE get_ith_conn
        ! Retrieves connectivity data for the B(i) element.
        !
        ! NOTE: this getter has been designed for high frequency access,
        ! therefore the class content is POINTED instead of being COPIED.
        ! WARNING! In the calling unit, don't use ITH_CONN on the LHS of
        ! assignment expressions, if you want to keep intact the original
        ! connectivity data.

        !
        INTEGER :: i1, i2

        ! REMARK: preliminary checks about pointers status have been omitted
        ! in order to avoid repeated IF cycles when the getter is called in a
        ! do loop.

        i1 = a2b%lookup(i)
        i2 = a2b%lookup(i+1) - 1

        ith_conn => a2b%conn(i1:i2) ! That's a pointing.

    END PROCEDURE get_ith_conn



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
    MODULE PROCEDURE get_conn_csr
        ! Retrieves all connectivity data in CSR format.
        !
        ! NOTE: this getter has been designed in order to get all connectivity
        ! data in a single shot. The pointers LOOKUP and CONN receive a COPY of
        ! A2B%LOOKUP and A2B%CONN respectively.
        ! They can be manipulated in the calling unit without affecting the
        ! original class content.

        !
        INTEGER :: lb,ub
        INTEGER :: info,  nconn

        IF(    ALLOCATED(lookup) .OR. &
            & ALLOCATED(conn)) THEN
            DEALLOCATE(lookup,conn,stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF
        END IF

        lb = LBOUND(a2b%lookup,1)
        ub = UBOUND(a2b%lookup,1)

        nconn = SIZE(a2b%conn)

        ALLOCATE(lookup(lb:ub),conn(nconn),stat=info)
        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

        lookup = a2b%lookup
        conn = a2b%conn

100     FORMAT(' ERROR! Memory deallocation failure in GET_CONN_CSR')
200     FORMAT(' ERROR! Memory allocation failure in GET_CONN_CSR')

    END PROCEDURE get_conn_csr


    MODULE PROCEDURE get_dual_conn
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

        !
        INTEGER :: i, lookup, j, info, ir, is
        INTEGER :: nel, nrd, nrt, ntot
        INTEGER, ALLOCATABLE :: level(:)
        INTEGER :: clb, cub ! lower & upper bound of a2b connectivity (lb can be negative)
        INTEGER :: llb, lub ! lower & upper bound of a2b lookup (lb can be negative)

        ! Number of rows in the original A2B CONNECTIVITY
        nrt = SIZE(a2b%lookup) - 1

        ! Note bounds of lookup

        llb = LBOUND(a2b%lookup,1)
        lub = UBOUND(a2b%lookup,1)

        ! Note bounds of connectivity
        clb = MINVAL(a2b%conn)

        ! if the numbering is negative, take note.  Otherwise, start at 1
        IF (clb > 1) THEN
            WRITE(0,*) 'Warning: get_dual_connectivity: dual starts at ',clb
        END IF
        clb = MIN (clb,1)

        cub = MAXVAL(a2b%conn)

        ! Number of rows of the B2A CONNECTIVITY, i.e. number of stacks
        nrd = cub - clb + 1

        ! The total amount of elements in the stacks is equal to:
        ! - the number of elements in A2B%CONN
        ! - the number of elements in B2A%CONN
        ! I.e. the number of links expressed by a CONNECTIVITY is equal
        ! in the DUAL one.
        ntot = SIZE(a2b%conn)

        ! LEVEL(I) = level of I-th stack container in B2A dual CONNECTIVITY.
        ALLOCATE(level(clb:cub),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! 1st swap of A2B%CONN -> Evaluates the level of each stack
        level(:) = 0
        DO i=1,SIZE(a2b%conn)

            is = a2b%conn(i)

            ! Increase the filling level of the IS-th stack
            level(is) = level(is) + 1

        END DO

        CALL alloc_conn(b2a,nel=nrd,nconn=ntot,lb=clb)

        ! The following operation is equivalent to size the empty stacks
        b2a%lookup(clb) = 1
        DO is = clb, cub
            b2a%lookup(is+1) = b2a%lookup(is) + level(is)
        END DO

        ! 2nd swap of A2B%CONN --> Building of the stacks set B2A%conn
        !
        ! Resets to zero the level counter of the stacks set
        level(:) = 0

        DO ir = llb,lub-1

            lookup = a2b%lookup(ir) - 1 !where (-1) this element ir exists in a2b%conn

            nel = a2b%lookup(ir+1) - a2b%lookup(ir) ! number of elements to scan

            DO i = 1, nel

                is = a2b%conn(lookup+i) !which b element is referencing a

                ! Increase the filling level of the IS-th stack
                level(is) = level(is) + 1
                !
                ! Sets the position on the B2A CONNECTIVITY of the new element
                ! found to belong to the IS-th stack
                j = b2a%lookup(is) -1 + level(is) ! the -1 because we pre-incremented
                !
                ! Stores index IR in the IS-th stack
                b2a%conn(j) = ir

            END DO

        END DO

        DEALLOCATE(level)

100     FORMAT(' ERROR! Memory allocation failure in GET_DUAL_CONN')

    END PROCEDURE get_dual_conn


    MODULE PROCEDURE nel_
        ! get number of elements (either local or global)

        !
        INTEGER :: nglob, nloc  ! global and local number of elements

        nglob = a2b%nel_glob
        nloc  = SIZE(a2b%lookup) - 1

        ! Default is local
        nel_ = nloc

        ! Optional
        IF(PRESENT(gl)) THEN
            SELECT CASE(gl)
            CASE('g')
                nel_ = nglob
            CASE('l')
                nel_ = nloc
            CASE default
                WRITE(*,100)
                CALL abort_psblas
            END SELECT
        END IF

100     FORMAT(' ERROR! Optional arg GL in NEL_ must be ''g'' or ''l''')

    END PROCEDURE nel_


    MODULE PROCEDURE nconn_
        !
        INTEGER :: nglob, nloc

        nglob = a2b%nconn_glob
        nloc  = SIZE(a2b%conn)

        ! Default is local
        nconn_ = nloc

        ! Optional
        IF(PRESENT(gl)) THEN
            SELECT CASE(gl)
            CASE('g')
                nconn_ = nglob
            CASE('l')
                nconn_ = nloc
            CASE default
                WRITE(*,100)
                CALL abort_psblas
            END SELECT
        END IF

100     FORMAT(' ERROR! Optional arg GL in NCONN_ must be ''g'' or ''l''')

    END PROCEDURE nconn_


    ! ----- Setters -----

    MODULE PROCEDURE set_ith_conn
        !
        INTEGER :: i1, i2

        ! WARNING: if A2B is declared only as INTENT(OUT), it becomes UNDEFINED
        ! on entry to the procedure (See Metcalf). Not detected by Gfortran.

        ! WARNING: preliminary checks about pointers status have been omitted
        ! in order to avoid repeated IF cycles when the setter is called in a
        ! do loop.

        i1 = a2b%lookup(i)
        i2 = i1 + SIZE(ith_conn) - 1

        a2b%conn(i1:i2) = ith_conn(:) ! That's a copy.

        a2b%lookup(i+1) = i2 + 1

    END PROCEDURE set_ith_conn


    MODULE PROCEDURE copy_conn
        !
        INTEGER :: lb, nel, nconn

        lb = LBOUND(a2b%lookup,1)
        nel = SIZE(a2b%lookup) - 1
        nconn = SIZE(a2b%conn)

        CALL alloc_conn(c2d,nel=nel,nconn=nconn,lb=lb)

        c2d%lookup(:) = a2b%lookup(:)
        c2d%conn(:)   = a2b%conn(:)

    END PROCEDURE copy_conn


    ! ----- Other -----

    MODULE PROCEDURE max_conn
        ! Computes the maximum connectivity degree

        !
        INTEGER :: i, j, ideg
        INTEGER :: lb_, ub_

        lb_ = LBOUND(a2b%lookup,1)
        ub_ = UBOUND(a2b%lookup,1) - 1

        IF(PRESENT(lb)) THEN
            IF(lb < lb_ .OR. lb > ub_) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF
            lb_ = lb
        END IF

        ideg = 0
        DO i = lb_, ub_
            j = a2b%lookup(i+1) - a2b%lookup(i)
            IF(j > ideg) ideg = j
        END DO

        max_conn = ideg

100     FORMAT(' ERROR! LB argument in MAX_CONN is out of bounds')

    END PROCEDURE max_conn


    ! ----- Integrity Checks -----

    MODULE PROCEDURE unused_elements
        ! For debugging connectivity changes
        ! Indicates if unused elements, B, occur in the A2B structure.
        ! This function will probably be deleted in the future, but for
        ! now it gives me a place for the debugger to break when something
        ! is wrong with the mesh.

        !
        INTEGER,ALLOCATABLE :: occurrence(:)
        INTEGER :: ielement, icount
        INTEGER :: info                       ! indicates success/failure
        INTEGER :: nb                         ! number of elements referred to

        nb = MAXVAL(a2b%conn)

        ALLOCATE(occurrence(nb),stat=info) ! assumes first element is index 1
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        occurrence(:) = 0  ! initialize the simple way

        occurrence = count_references(a2b,nb)

        icount = 0
        DO ielement = 1, SIZE(occurrence)
            IF (occurrence(ielement) == 0) THEN
                icount = icount + 1
                WRITE(*,200)
            ENDIF
        ENDDO

        unused_elements = icount

        DEALLOCATE(occurrence,stat=info)
        IF(info /= 0) THEN
            WRITE(*,300)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Memory allocation failure in UNUSED_ELEMENTS function')
200     FORMAT(' WARNING! Unused element in connectivity structure')
300     FORMAT(' ERROR! Memory deallocation failure in UNUSED_ELEMENTS function')

    END PROCEDURE unused_elements


    MODULE PROCEDURE count_references
        ! Counts the number of references to each B in the A2B structure

        INTEGER :: occurrence(nb)             ! how many times each B is referenced
        INTEGER :: ib, ielement               ! loop indices

        occurrence(:) = 0  ! initialize the simple way

        DO ib = 1, SIZE(a2b%conn)
            ielement = a2b%conn(ib)
            occurrence(ielement) = occurrence(ielement) + 1
        ENDDO

        count_references = occurrence

    END PROCEDURE count_references


END SUBMODULE class_connectivity_procedures
