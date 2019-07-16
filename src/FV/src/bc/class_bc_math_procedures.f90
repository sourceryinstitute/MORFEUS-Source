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
! $Id: class_bc_math.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    MATHematical boundary condition class.
!
SUBMODULE(class_bc_math) class_bc_math_procedures
    USE class_face

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_bc_math_sizeof
        USE psb_base_mod

        INTEGER(kind=nemo_int_long_)   :: val

        val = 2 * nemo_sizeof_int
        IF (ALLOCATED(bc%a)) &
            & val = val + nemo_sizeof_dp * SIZE(bc%a)
        IF (ALLOCATED(bc%b)) &
            & val = val + nemo_sizeof_dp * SIZE(bc%b)
        IF (ALLOCATED(bc%c)) &
            & val = val + nemo_sizeof_dp * SIZE(bc%c)
        nemo_bc_math_sizeof = val

    END PROCEDURE nemo_bc_math_sizeof


    ! REMARK: the implementation of run-time polymorphism requires
    ! specific BC object as POINTERS!

    ! ----- Constructor -----

    MODULE PROCEDURE create_bc_math
        USE class_mesh
        USE tools_bc

        INTEGER :: info

        ! Allocates bc_math target on every process
        IF(ASSOCIATED(bc)) THEN
            WRITE(*,100)
            CALL abort_psblas
        ELSE
            ALLOCATE(bc,stat=info)
            IF(info /= 0) THEN
                WRITE(*,200)
                CALL abort_psblas
            END IF
        END IF

        ! Allocates and sets BC class members, according to the parameters
        ! get from the input file.
        CALL rd_inp_bc_math(input_file,sec,nbf,bc%id,bc%a,bc%b,bc%c)

        ! Set number of boundary faces
        bc%nbf = nbf

100     FORMAT(' ERROR! Illegal call to CREATE_BC_MATH: pointer already associated')
200     FORMAT(' ERROR! Memory allocation failure in CREATE_BC_MATH')

    END PROCEDURE create_bc_math


    MODULE PROCEDURE alloc_bc_math
        USE tools_bc

        INTEGER :: na, nb, nc
        INTEGER :: info

        IF(    ALLOCATED(bc%a) .OR. &
            & ALLOCATED(bc%b) .OR. &
            & ALLOCATED(bc%c)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        SELECT CASE(id)
        CASE(  bc_dirichlet_, &
            & bc_neumann_, bc_neumann_flux_, &
            & bc_robin_, bc_robin_convection_ )
            na = 1; nb = 1; nc = 1
        CASE(  bc_dirichlet_map_, &
            & bc_neumann_map_)
            na = 1; nb = 1; nc = nbf
        CASE(  bc_robin_map_)
            na = nbf; nb = nbf; nc = nbf
        CASE default
            WRITE(*,200)
            CALL abort_psblas
        END SELECT

        IF(    na /= SIZE(a) .OR. &
            & nb /= SIZE(b) .OR. &
            & nc /= SIZE(c) ) THEN
            WRITE(*,300)
            CALL abort_psblas
        END IF

        ALLOCATE(bc%a(na),bc%b(nb),bc%c(nc),stat=info)
        IF(info /= 0) THEN
            WRITE(*,400)
            CALL abort_psblas
        END IF

        bc%id = id
        bc%nbf = nbf
        bc%a = a
        bc%b = b
        bc%c = c

100     FORMAT(' ERROR! Illegal call to ALLOC_BC_MATH.',/,&
            & ' One among A, B, C pointers is already allocated.')
200     FORMAT(' ERROR! Unsupported BC type in ALLOC_BC_MATH')
300     FORMAT(' ERROR! Size mismatch in ALLOC_BC_MATH')
400     FORMAT(' ERROR! Memory allocation failure in ALLOC_BC_MATH')

    END PROCEDURE alloc_bc_math


    ! ----- Destructor -----

    MODULE PROCEDURE free_bc_math
        !! To be invoked when BC_MATH is used as high-level BC.
        !
        INTEGER :: info

        CALL dealloc_bc_math(bc)

        DEALLOCATE(bc,stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Memory allocation failure in FREE_BC_MATH')

    END PROCEDURE free_bc_math


    MODULE PROCEDURE dealloc_bc_math
        !! To be invoked when BC_MATH is a member of another BC class.
        !
        INTEGER :: info

        DEALLOCATE(bc%a,bc%b,bc%c,stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Memory allocation failure in DEALLOC_BC_MATH')

    END PROCEDURE dealloc_bc_math


    ! ----- Getter -----

    MODULE PROCEDURE get_abc_math_s
        !USE class_connectivity
        USE class_mesh
        USE tools_bc
        !
        INTEGER :: IF, n

        ! WARNING! Use intent(inout) for A, B and C.

        n = bc%nbf ! Number of Boundary Faces

        id = bc%id
        IF(  SIZE(a) /= n .OR. &
            SIZE(b) /= n .OR. &
            SIZE(c) /= n) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        SELECT CASE(bc%id)
        CASE(  bc_dirichlet_,    &
            & bc_neumann_  ,    &
            & bc_robin_    ,    &
            & bc_neumann_flux_, &
            & bc_robin_convection_)
            DO IF = 1, n
                a(IF) = bc%a(1)
                b(IF) = bc%b(1)
                c(IF) = bc%c(1)
            END DO
        CASE(  bc_dirichlet_map_, &
            & bc_neumann_map_)
            DO IF = 1, n
                a(IF) = bc%a(1)
                b(IF) = bc%b(1)
                c(IF) = bc%c(IF)
            END DO
        CASE(  bc_robin_map_)
            DO IF = 1, n
                a(IF) = bc%a(IF)
                b(IF) = bc%b(IF)
                c(IF) = bc%c(IF)
            END DO
        CASE default
            WRITE(*,200)
            CALL abort_psblas
        END SELECT


100     FORMAT(' ERROR! Size mismatch in GET_ABC_MATH_S')
200     FORMAT(' ERROR! Unsupported BC in GET_ABC_MATH_S')

    END PROCEDURE get_abc_math_s


    MODULE PROCEDURE get_abc_math_v
    !! WARNING! Use intent(inout) for A, B and C.
    !! REMARK: BC(:) elements are supposed to differ only in "C" term

        !USE class_connectivity
        USE class_mesh
        USE class_vector
        USE tools_bc

        INTEGER :: IF, n

        IF(SIZE(bc) /= 3) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        n = bc(1)%nbf ! Number of Boundary Faces

        id = bc(1)%id
        IF(  SIZE(a) /= n .OR. &
            SIZE(b) /= n .OR. &
            SIZE(c) /= n) THEN
            PRINT *, n, SIZE(a), SIZE(b), SIZE(c)
            WRITE(*,200)
            CALL abort_psblas
        END IF

        SELECT CASE(bc(1)%id)
        CASE(  bc_dirichlet_,    &
            & bc_neumann_)
            DO IF = 1, n
                a(IF) = bc(1)%a(1)
                b(IF) = bc(1)%b(1)
                c(IF) = vector_(bc(1)%c(1),bc(2)%c(1),bc(3)%c(1))
            END DO
        CASE(  bc_dirichlet_map_, &
            & bc_neumann_map_)
            DO IF = 1, n
                a(IF) = bc(1)%a(1)
                b(IF) = bc(1)%b(1)
                c(IF) = vector_(bc(1)%c(IF),bc(2)%c(IF),bc(3)%c(IF))
            END DO
        CASE default
            WRITE(*,300)
            CALL abort_psblas
        END SELECT

100     FORMAT(' ERROR! BC argument in GET_ABC_MATH_V is not 3D')
200     FORMAT(' ERROR! Size mismatch in GET_ABC_MATH_V')
300     FORMAT(' ERROR! Unsupported BC in GET_ABC_MATH_V')

    END PROCEDURE get_abc_math_v


    MODULE PROCEDURE is_allocated

        IF(    ALLOCATED(bc%a) .OR. &
            & ALLOCATED(bc%b) .OR. &
            & ALLOCATED(bc%c)) THEN
            is_allocated = .TRUE.
        ELSE
            is_allocated = .FALSE.
        END IF

    END PROCEDURE is_allocated


    ! ----- Setter -----

    MODULE PROCEDURE set_bc_math_map
        USE tools_bc

        SELECT CASE(bc%id)
        CASE(bc_dirichlet_map_, bc_neumann_map_)
            bc%a = a
            bc%b = b
            bc%c(i) = c
        CASE(bc_robin_map_)
            bc%a(i) = a
            bc%b(i) = b
            bc%c(i) = c
        CASE default
            WRITE(*,100)
            CALL abort_psblas
        END SELECT

100     FORMAT(' ERROR! Unsupported BC type in SET_BC_MATH_MAP')

    END PROCEDURE set_bc_math_map


    ! ----- Boundary Values Updating -----

    MODULE PROCEDURE apply_abc_to_boundary_s
        !! WARNING!
        !! - Use intent(inout) for BX(:)
        !! - Do loop on the faces subset corresponding to IB bc.
        !! - Only this section of BX(:) is going to be modified.
        !! - BX(:) indexing starts from 1: when BX(:) is referenced an offset
        !!   equal to the # of boundary faces with flag > 0 and < IB must be
        !!   added to the I counter.
        USE class_connectivity
        USE class_mesh
        USE tools_bc

        INTEGER, POINTER :: if2b(:) => NULL()
        INTEGER :: i, ibf, IF, im, ib_offset, n
        REAL(psb_dpk_) :: d, r
        TYPE(connectivity) :: conn_temp


        ! Boundary faces with flag < IB
        ib_offset = COUNT(msh%faces%flag_() > 0 .AND. msh%faces%flag_() < ib)

        !CALL msh%f2b%get_ith_conn(if2b,ib)
        conn_temp=msh%f2b
        CALL conn_temp%get_ith_conn(if2b, ib)
        n = SIZE(if2b)

        SELECT CASE(id)
        CASE(bc_dirichlet_)
            DO i = 1, n
                ibf = ib_offset + i
                bx(ibf) = c(1)
            END DO
        CASE(bc_dirichlet_map_)
            DO i = 1, n
                ibf = ib_offset + i
                bx(ibf) = c(i)
            END DO
        CASE(bc_neumann_)
            DO i = 1, n
                IF = if2b(i)
                im = msh%faces(IF)%master_()

                d = msh%dist(IF)

                ibf = ib_offset + i
                bx(ibf) = x(im) + c(1) * d
            END DO
        CASE(bc_robin_)
            DO i = 1, n
                IF = if2b(i)
                im = msh%faces(IF)%master_()

                d = msh%dist(IF)

                ibf = ib_offset + i
                r = a(1) * d + b(1)

                bx(ibf) = (b(1) * x(im) + c(1) * d) / r
            END DO
        CASE(bc_neumann_flux_)
            DO i = 1, n
                IF = if2b(i)
                im = msh%faces(IF)%master_()

                d = msh%dist(IF)

                ibf = ib_offset + i
                bx(ibf) = x(im) + c(1) * d / b(i)
            END DO
        CASE(bc_robin_convection_)
            DO i = 1, n
                IF = if2b(i)
                im = msh%faces(IF)%master_()

                d = msh%dist(IF)

                ibf = ib_offset + i
                r = a(1) * d + b(i)

                bx(ibf) = (b(i) * x(im) + c(1) * d) / r
            END DO
        CASE(bc_robin_map_)
            DO i = 1, n
                IF = if2b(i)
                im = msh%faces(IF)%master_()

                d = msh%dist(IF)

                ibf = ib_offset + i
                r = a(i) * d + b(i)

                bx(ibf) = (b(i) * x(im) + c(i) * d) / r
            END DO
        CASE default
            WRITE(*,100)
            CALL abort_psblas
        END SELECT

        NULLIFY(if2b)

100     FORMAT(' ERROR! Unsupported BC type in APPLY_BC_TO_BOUNDARY_S')

    END PROCEDURE apply_abc_to_boundary_s


    MODULE PROCEDURE apply_abc_to_boundary_v
        !! WARNING!
        !! - Use intent(inout) for BX(:)
        !! - Do loop on the faces subset corresponding to IB bc.
        !! - Only this section of BX(:) is going to be modified.
        !! - BX(:) idexing starts from 1: when BX(:) is referenced an offset
        !!   equal to the # of boundary faces with flag > 0 and < IB must be
        !!   added to the I counter.

        USE class_connectivity
        USE class_mesh
        USE class_vector
        USE tools_bc

        INTEGER, POINTER :: if2b(:) => NULL()
        INTEGER :: i, ibf, IF, im, ib_offset, n
        REAL(psb_dpk_) :: d
        TYPE(connectivity) :: conn_temp

        ! Dummy usage of A and B
        i = SIZE(a)
        i = SIZE(b)

        ! Boundary faces with flag < IB
        ib_offset = COUNT(msh%faces%flag_() > 0 .AND. msh%faces%flag_() < ib)

        !CALL msh%f2b%get_ith_conn(if2b,ib)
        conn_temp=msh%f2b
        CALL conn_temp%get_ith_conn(if2b,ib)

        n = SIZE(if2b)

        SELECT CASE(id)
        CASE(bc_dirichlet_)
!!$      write(0,*) 'apply_abc_to_boundary_v: Debug: dirichlet',n
            DO i = 1, n
                ibf = ib_offset + i
                bx(ibf) = c(1)
            END DO
        CASE(bc_dirichlet_map_)
!!$      write(0,*) 'apply_abc_to_boundary_v: Debug: dirichlet_map',n
            DO i = 1, n
                ibf = ib_offset + i
                bx(ibf) = c(i)
            END DO
        CASE(bc_neumann_)
            DO i = 1, n
                IF = if2b(i)
                im = msh%faces(IF)%master_()

                d = msh%dist(IF)

                ibf = ib_offset + i
                bx(ibf) = x(im) + d * c(1)
            END DO
        CASE default
            WRITE(*,100)
            CALL abort_psblas
        END SELECT

        NULLIFY(if2b)

100     FORMAT(' ERROR! Unsupported BC type in APPLY_BC_TO_BOUNDARY_V')

    END PROCEDURE apply_abc_to_boundary_v


    MODULE PROCEDURE update_boundary_math
        !! WARNING! Use intent(inout) for BX(:)
        USE class_mesh

        CALL apply_abc_to_boundary(bc%id,bc%a,bc%b,bc%c,ib,msh,x,bx)

    END PROCEDURE update_boundary_math


    ! ----- Debug -----

    MODULE PROCEDURE debug_bc_math
        USE tools_bc

        CHARACTER(len=30) :: bc_type

        SELECT CASE(bc%id)
        CASE(bc_dirichlet_)
            bc_type = 'Homogeneous Dirichlet'
        CASE(bc_neumann_, bc_neumann_flux_)
            bc_type = 'Homogeneous Neumann'
        CASE(bc_robin_, bc_robin_convection_)
            bc_type = 'Homogeneous Robin'
        CASE default
            bc_type = 'Unknown '
            WRITE(0,*) bc_type, ' ',bc%id, 'Bailing out'
            RETURN
        END SELECT

        WRITE(*,100) '  Type of BC_MATH = ', TRIM(bc_type)
        WRITE(*,200) '  BC_MATH%a = ', bc%a
        WRITE(*,200) '  BC_MATH%b = ', bc%b
        WRITE(*,200) '  BC_MATH%c = ', bc%c

100     FORMAT(1x,a,a)
200     FORMAT(1x,a,es10.3)

    END PROCEDURE debug_bc_math


END SUBMODULE class_bc_math_procedures
