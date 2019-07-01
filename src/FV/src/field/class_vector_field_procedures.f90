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
! $Id: class_vector_field.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(class_vector_field) class_vector_field_procedures

    USE class_psblas, ONLY : psb_dpk_, nemo_int_long_, nemo_sizeof_dp, nemo_sizeof_int,&
        & icontxt_, psb_erractionsave, abort_psblas, psb_check_error, psb_erractionrestore,&
        & psb_halo
    USE class_field
    USE class_vector

    IMPLICIT NONE

CONTAINS


    MODULE PROCEDURE nemo_vector_field_sizeof
        IMPLICIT NONE

        INTEGER(kind=nemo_int_long_)   :: val

        val = fld%base%nemo_sizeof()
        IF (ALLOCATED(fld%x)) &
            & val = val + nemo_sizeof_dp * SIZE(fld%x)
        IF (ALLOCATED(fld%xp)) &
            & val = val + nemo_sizeof_dp * SIZE(fld%xp)
        IF (ALLOCATED(fld%bx)) &
            & val = val + nemo_sizeof_dp * SIZE(fld%bx)

        nemo_vector_field_sizeof = val

    END PROCEDURE nemo_vector_field_sizeof


    ! ----- Constructor -----

    MODULE PROCEDURE vector_field_
        IMPLICIT NONE
        !! Default public constructor, necessary with ifort

        !vector_field_ = vector_field(base,x,bx)
        !! Workaround for Intel 18 error #6053: Structure constructor may not have fields with a PRIVATE attribute
        vector_field_%base = base
        vector_field_%x    = x
        vector_field_%bx   = bx

        ! Allocation check
        IF(  .NOT.ALLOCATED(vector_field_%x) .OR. &
            .NOT.ALLOCATED(vector_field_%x)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Memory allocation failure in VECTOR_FIELD_')

    END PROCEDURE vector_field_


    MODULE PROCEDURE create_vector_field
        USE class_bc
        USE class_connectivity
        USE class_dimensions
        USE class_material
        USE class_mesh
        USE class_vector
        USE tools_material
        IMPLICIT NONE

        INTEGER          :: info, nel, nbf, isize(2), i, im, ic, ncg, n_mats
        INTEGER, POINTER :: ic2g(:) => NULL()
        TYPE(vector)     :: x0_
        TYPE(dimensions) :: fdim


        ! Creates the base-class member
        CALL fld%base%create_field(msh,dim,bc,mats,on_faces)

        ! Gets field dimensions
        fdim = fld%base%dim_()

        ! Gets field size
        isize = fld%base%fld_size()
        nel   = isize(fld_internal_)
        nbf   = isize(fld_boundary_)

        ! Allocates arrays for inner and boundary elements
        CALL alloc_vector(fld%x,nel)
        CALL alloc_vector(fld%xp,nel)
        CALL alloc_vector(fld%bx,nbf)
        ALLOCATE(fld%mat(nel),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! Material allocation
        n_mats = 0
        IF(PRESENT(mats)) n_mats = SIZE(mats)
        IF (n_mats == 1) THEN
            fld%mat(:) = 1
        ELSE IF(PRESENT(mats)) THEN
            DO im = 1, SIZE(mats)
                CALL msh%c2g%get_ith_conn(ic2g,im)
                ncg = SIZE(ic2g)
                DO i = 1, ncg
                    ic = ic2g(i)
                    fld%mat(ic) = im
                END DO
            END DO
        END IF

        ! Field initialization
        IF(PRESENT(x0)) THEN
            x0_ = x0
        ELSE
            x0_ = vector_(0.d0,0.d0,0.d0)
        END IF

        fld%x(:)  = x0_
        fld%xp(:)  = x0_
        fld%bx(:) = x0_

100     FORMAT(' ERROR! Memory allocation failure in CREATE_VECTOR_FIELD')

    END PROCEDURE create_vector_field


    ! ----- Destructor -----

    MODULE PROCEDURE free_vector_field
        IMPLICIT NONE

        CALL fld%base%free_field()
        CALL free_vector(fld%x)
        CALL free_vector(fld%bx)

    END PROCEDURE free_vector_field


    ! ----- Getters for Inherited Members -----

    MODULE PROCEDURE get_vector_field_name
        IMPLICIT NONE

        get_vector_field_name = fld%base%name_()

    END PROCEDURE get_vector_field_name


    MODULE PROCEDURE get_vector_field_dim
        USE class_dimensions
        IMPLICIT NONE

        get_vector_field_dim = fld%base%dim_()

    END PROCEDURE get_vector_field_dim


    MODULE PROCEDURE get_vector_field_msh_fun
        USE class_mesh
        IMPLICIT NONE

        get_vector_field_msh_fun => fld%base%msh_()

    END PROCEDURE get_vector_field_msh_fun


    MODULE PROCEDURE get_vector_field_msh_sub
        IMPLICIT NONE

        CALL fld%base%get_mesh(msh)

    END PROCEDURE get_vector_field_msh_sub


    MODULE PROCEDURE get_vector_field_on_faces
        IMPLICIT NONE

        get_vector_field_on_faces = fld%base%on_faces_()

    END PROCEDURE get_vector_field_on_faces


    MODULE PROCEDURE get_vector_field_bc
        USE class_bc
        IMPLICIT NONE

        get_vector_field_bc => fld%base%bc_()

    END PROCEDURE get_vector_field_bc


    MODULE PROCEDURE get_vector_field_mat
        USE class_material
        IMPLICIT NONE

        get_vector_field_mat => fld%base%mat_()

    END PROCEDURE get_vector_field_mat


    MODULE PROCEDURE get_vector_field_mat_sub
        USE class_material
        IMPLICIT NONE

        IF (PRESENT(i)) THEN
            CALL fld%base%get_material(i,mat)
        ELSE
            CALL fld%base%get_material(1,mat)
        END IF
    END PROCEDURE get_vector_field_mat_sub


    MODULE PROCEDURE get_vector_field_base
        IMPLICIT NONE

        base = fld%base

    END PROCEDURE get_vector_field_base


    ! ----- Getters for Additional Members -----

    MODULE PROCEDURE get_vector_field_x_r
        IMPLICIT NONE
        INTEGER :: i, info, n

        n = SIZE(fld%x)

        ALLOCATE(x(n,3),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        DO i = 1, n
            x(i,1) = fld%x(i)%x_()
            x(i,2) = fld%x(i)%y_()
            x(i,3) = fld%x(i)%z_()
        END DO

100     FORMAT(' ERROR! Memory allocation failure in GET_VECTOR_FIELD_X_R')

    END PROCEDURE get_vector_field_x_r


    MODULE PROCEDURE get_vector_field_x_v
        IMPLICIT NONE

        CALL alloc_vector(x,SIZE(fld%x))

        x = fld%x

    END PROCEDURE get_vector_field_x_v

    MODULE PROCEDURE get_vector_field_xp_r
        IMPLICIT NONE
        INTEGER :: i, info, n

        n = SIZE(fld%xp)

        ALLOCATE(xp(n,3),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        DO i = 1, n
            xp(i,1) = fld%xp(i)%x_()
            xp(i,2) = fld%xp(i)%y_()
            xp(i,3) = fld%xp(i)%z_()
        END DO

100     FORMAT(' ERROR! Memory allocation failure in GET_VECTOR_FIELD_XP_R')

    END PROCEDURE get_vector_field_xp_r


    MODULE PROCEDURE get_vector_field_xp_v
        IMPLICIT NONE

        CALL alloc_vector(xp,SIZE(fld%xp))

        xp = fld%xp

    END PROCEDURE get_vector_field_xp_v


    MODULE PROCEDURE get_vector_field_bx_r
        IMPLICIT NONE
        INTEGER :: info

        ALLOCATE(bx(SIZE(fld%bx),3),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        bx(:,1) = fld%bx%x_()
        bx(:,2) = fld%bx%y_()
        bx(:,3) = fld%bx%z_()

100     FORMAT(' ERROR! Memory allocation failure in GET_VECTOR_FIELD_BX_R')

    END PROCEDURE get_vector_field_bx_r


    MODULE PROCEDURE get_vector_field_bx_v
        IMPLICIT NONE

        CALL alloc_vector(bx,SIZE(fld%bx))

        bx = fld%bx

    END PROCEDURE get_vector_field_bx_v


    ! ----- Setters -----


    MODULE PROCEDURE set_vector_field_x
        IMPLICIT NONE
        INTEGER :: i, info, n

        n = SIZE(fld%x)

        DO i = 1, n
            CALL fld%x(i)%set_x_(x(i,1))
            CALL fld%x(i)%set_y_(x(i,2))
            CALL fld%x(i)%set_z_(x(i,3))
        END DO

    END PROCEDURE set_vector_field_x


    MODULE PROCEDURE update_vector_field
        USE class_bc
        USE class_dimensions
        USE class_material
        USE class_mesh
        IMPLICIT NONE

        INTEGER :: ib, info, err_act
        TYPE(mesh),     POINTER :: msh   => NULL()
        TYPE(bc_poly),  POINTER :: bc(:) => NULL()
        TYPE(material), POINTER :: mat   => NULL()

        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)


        ! Gets pointer base-class members
!!$    msh => msh_(fld%base)
        CALL fld%base%get_mesh(msh)
        bc  => fld%base%bc_()


        ! Face-centered unknown
        IF(fld%on_faces_()) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF


        ! Cell-centered unknown
!!$      mat => mat_(fld%base)
        CALL fld%base%get_material(1,mat)

        CALL update_vector_halo(fld%x,msh%desc_c)

        DO ib = 1, msh%nbc
            CALL update_boundary(ib,bc(ib),fld%dim_(),msh,fld%x,fld%bx)
        END DO

        NULLIFY(mat)
        NULLIFY(bc)
        NULLIFY(msh)


        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)


100     FORMAT(' ERROR! Vector unknown field is face-centered.',/,&
            & ' Illegal call to UPDATE_VECTOR_FIELD. Use INTERP_FIELD instead!')

    END PROCEDURE update_vector_field


    MODULE PROCEDURE assign_vector_field_s
        USE class_vector
        IMPLICIT NONE

        f%x(:) = x

    END PROCEDURE assign_vector_field_s


    MODULE PROCEDURE assign_vector_field_v
        USE class_vector
        IMPLICIT NONE

        IF(SIZE(f%x) /= SIZE(x)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        f%x = x

100     FORMAT(' ERROR! Size mismatch in ASSIGN_VECTOR_FIELD_V')

    END PROCEDURE assign_vector_field_v


    MODULE PROCEDURE set_vector_field_element
        USE class_vector
        IMPLICIT NONE

        IF(i < 1 .OR. i > SIZE(f%x,1)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        f%x(i) = x

100     FORMAT(' ERROR! Out of bounds index in SET_VECTOR_FIELD_ELEMENT')

    END PROCEDURE set_vector_field_element

    MODULE PROCEDURE set_vector_field_bound_element
        USE class_vector
        IMPLICIT NONE

        IF(i < 1 .OR. i > SIZE(f%bx,1)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        f%bx(i) = x

100     FORMAT(' ERROR! Out of bounds index in SET_VECTOR_FIELD_ELEMENT')

    END PROCEDURE set_vector_field_bound_element



    MODULE PROCEDURE set_vector_field_group
        USE class_connectivity
        USE class_mesh
        USE class_vector
        IMPLICIT NONE

        INTEGER :: i, ic, ncg
        INTEGER, POINTER :: ic2g(:) => NULL()
        TYPE(mesh), POINTER :: msh  => NULL()

        ! MESH pointer dereferencing
!!$    msh => msh_(f%base)
        CALL f%base%get_mesh(msh)

        IF(f%on_faces_()) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        IF(ig < 0 .OR. ig > msh%ngp) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

        CALL msh%c2g%get_ith_conn(ic2g,ig)
        ncg = SIZE(ic2g)

        DO i = 1, ncg
            ic = ic2g(i)
            f%x(ic) = x
        END DO

        NULLIFY(ic2g)
        NULLIFY(msh)

200     FORMAT(' ERROR! F field in SET_VECTOR_FIELD_GROUP is face-centered')
100     FORMAT(' ERROR! Group ',i1,' in SET_VECTOR_FIELD_GROUP',&
            & ' is out of bound')

    END PROCEDURE set_vector_field_group


    ! ----- Algebra Operations -----

    MODULE PROCEDURE vector_field_sum
        USE class_dimensions
        IMPLICIT NONE

        INTEGER :: isize(2), nel, nbf, info

        ! Check consistency of operands
        CALL f1%base%check_field_operands(f2%base,'VECTOR_FIELD_SUM')

        ! Check consistency of operand dimensions
        IF(f1%dim_() /= f2%dim_()) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF
        r%base = f1%base

        isize = f1%base%fld_size()
        nel   = isize(fld_internal_)
        nbf   = isize(fld_boundary_)

        ! Allocates arrays for inner and boundary elements
        ALLOCATE(r%x(nel),r%bx(nbf),stat=info)
        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

        r%x  = f1%x + f2%x
        r%bx = f1%bx + f2%bx

100     FORMAT(' ERROR! The dimensions of the operands in VECTOR_FIELD_SUM',&
            & ' are not consistent')
200     FORMAT(' ERROR! Memory allocation failure in VECTOR_FIELD_SUM')

    END PROCEDURE vector_field_sum

    MODULE PROCEDURE vector_field_scal
        USE class_dimensions
        IMPLICIT NONE

        INTEGER :: isize(2), nel, nbf, info

        r%base = f2%base

        isize = f2%base%fld_size()
        nel   = isize(fld_internal_)
        nbf   = isize(fld_boundary_)

        ! Allocates arrays for inner and boundary elements
        ALLOCATE(r%x(nel),r%bx(nbf),stat=info)
        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

        r%x  = a*f2%x
        r%bx = a*f2%bx

100     FORMAT(' ERROR! The dimensions of the operands in VECTOR_FIELD_SUM',&
            & ' are not consistent')
200     FORMAT(' ERROR! Memory allocation failure in VECTOR_FIELD_SUM')

    END PROCEDURE vector_field_scal


    MODULE PROCEDURE vector_field_dif
        USE class_dimensions
        IMPLICIT NONE
        INTEGER :: isize(2), nel, nbf, info

        ! Check consistency of operands
        CALL f1%base%check_field_operands(f2%base,'VECTOR_FIELD_SUM')

        ! Check consistency of operand dimensions
        IF(f1%dim_() /= f2%dim_()) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF
        r%base = f1%base

        isize = f1%base%fld_size()
        nel   = isize(fld_internal_)
        nbf   = isize(fld_boundary_)

        ! Allocates arrays for inner and boundary elements
        ALLOCATE(r%x(nel),r%bx(nbf),stat=info)
        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

        r%x  = f1%x - f2%x
        r%bx = f1%bx - f2%bx

100     FORMAT(' ERROR! The dimensions of the operands in VECTOR_FIELD_DIF',&
            & ' are not consistent')
200     FORMAT(' ERROR! Memory allocation failure in VECTOR_FIELD_DIF')

    END PROCEDURE vector_field_dif


    MODULE PROCEDURE interp_on_faces_v
        USE class_connectivity
        USE class_face
        USE class_mesh
        USE tools_math
        IMPLICIT NONE

        INTEGER :: i, IF, im, inb, is, ioffset
        INTEGER :: k, n, nel, nbf
        INTEGER, POINTER :: if2b(:) => NULL()
        REAL(psb_dpk_) :: w
        TYPE(mesh), POINTER :: msh => NULL()


        IF(fld%on_faces_()) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! Sets base member
        r%base = fld%base

        CALL fld%base%get_mesh(msh)
        nel = COUNT(msh%faces%flag_() <= 0)
        nbf = COUNT(msh%faces%flag_() > 0)

        ! Allocates arrays for inner and boundary elements
        CALL alloc_vector(r%x,nel)
        CALL alloc_vector(r%bx,nbf)


        ! 1) flag = 0 => interpolates master/slave cell values.
        CALL msh%f2b%get_ith_conn(if2b,0)
        n = SIZE(if2b) ! # of internal fluid faces

        DO i = 1, n
            IF = if2b(i)
            im = msh%faces(IF)%master_()
            is = msh%faces(IF)%slave_()
            w  = msh%interp(IF)
            r%x(i) = lin_interp(fld%x(im),fld%x(is),w)
        END DO

        ioffset = n

        ! 2) flag = -1 => uses cell value (same process)
        CALL msh%f2b%get_ith_conn(if2b,-1)
        n = SIZE(if2b)

        DO i = 1, n
            IF = if2b(i)
            im = msh%faces(IF)%master_()
            is = msh%faces(IF)%slave_()
            k = ioffset + i
            inb = MAX(im,is)
            r%x(k) = fld%x(inb)
        END DO


        ! 3) Copies boundary values
        r%bx = fld%bx
        CALL r%base%set_field_on_faces(.TRUE.)
        NULLIFY(if2b)
        NULLIFY(msh)

100     FORMAT(' ERROR! Argument of INTERP_ON_FACES_V is already on faces')

    END PROCEDURE interp_on_faces_v


    ! ----- Check Procedures -----

    MODULE PROCEDURE check_mesh_consistency_vf
        IMPLICIT NONE

        CALL f1%base%check_mesh_consistency(f2%base,WHERE)

    END PROCEDURE check_mesh_consistency_vf

END SUBMODULE class_vector_field_procedures
