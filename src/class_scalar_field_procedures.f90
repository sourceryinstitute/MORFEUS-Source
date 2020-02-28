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
! $Id: class_scalar_field.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(class_scalar_field) class_scalar_field_procedures
    USE class_field

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_sizeof
        !use psb_base_mod

        INTEGER(kind=nemo_int_long_)   :: val

        val = fld%field%nemo_sizeof()
        IF (ALLOCATED(fld%x)) &
            & val = val + nemo_sizeof_dp * SIZE(fld%x)
        IF (ALLOCATED(fld%xp)) &
            & val = val + nemo_sizeof_dp * SIZE(fld%xp)
        IF (ALLOCATED(fld%bx)) &
            & val = val + nemo_sizeof_dp * SIZE(fld%bx)

        nemo_sizeof = val

    END PROCEDURE nemo_sizeof

    ! ----- Constructor -----

    ! Default public constructor, necessary with ifort
    MODULE PROCEDURE scalar_field_

        !scalar_field_ = scalar_field(base,x,bx)
        !! Workaround for Intel 18 error #6053: Structure constructor may not have fields with a PRIVATE attribute
        scalar_field_%field = base
        scalar_field_%x    = x
        scalar_field_%bx   = bx

        ! Allocation check
        IF(  .NOT.ALLOCATED(scalar_field_%x) .OR. &
            .NOT.ALLOCATED(scalar_field_%bx)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Memory allocation failure in SCALAR_FIELD_')

    END PROCEDURE scalar_field_


    MODULE PROCEDURE create_scalar_field
        USE class_bc
        USE class_connectivity
        USE class_dimensions
        USE class_material
        USE class_mesh
        USE class_cell
        USE class_face
        USE tools_material

        INTEGER :: info, nel, nbf, isize(2), i, j, ii, &
            im, ic, ncg, ib, ib_offset, IF, ibf, n, n_mats
        INTEGER, POINTER :: ic2g(:) => NULL(), if2c(:) => NULL(), if2b(:) => NULL()
        REAL(psb_dpk_), ALLOCATABLE :: x0_(:), xb0_(:), t0_(:), tb0_(:)
        INTEGER, ALLOCATABLE :: im_(:), imb_(:)
        TYPE(dimensions) :: fdim

        ! Creates the base-class member
        CALL fld%field%create_field(msh,dim,bc,mats,on_faces,name)

        ! Gets field dimensions
        fdim = fld%field%dim_()

        ! Gets field size
        isize = fld%field%fld_size()
        nel   = isize(fld_internal_)
        nbf   = isize(fld_boundary_)

        ! Allocates arrays for inner and boundary elements
        ALLOCATE(fld%x(nel),fld%xp(nel),fld%bx(nbf),fld%mat(nel),fld%bmat(nbf),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! Material allocation
        ! for the cells
        n_mats = 0
        IF(PRESENT(mats)) n_mats = SIZE(mats)
        IF (n_mats == 1) THEN
            fld%mat(:) = 1
            fld%bmat(:) = 1
        ELSE
            IF (.NOT.fld%field%on_faces_()) THEN
                DO ii = 1, SIZE(fld%mat)
                    fld%mat(ii) = msh%cells(ii)%group_()
                END DO
            ELSE
                DO ii = 1, SIZE(fld%mat)
                    ic = msh%faces(ii)%master_()
                    IF (ic < 0) THEN
                        ic = msh%faces(ii)%slave_()
                    END IF
                    fld%mat(ii) = msh%cells(ic)%group_()
                END DO
            END IF

            ! for the faces that are on the boundary
            DO ib = 1, msh%nbc
                CALL msh%f2b%get_ith_conn(if2b,ib)
                n = SIZE(if2b)
                ib_offset = COUNT(msh%faces%flag_() > 0 .AND. msh%faces%flag_() < ib)
                DO i = 1, n
                    IF = if2b(i)
                    ibf = ib_offset + i
                    ic = msh%faces(IF)%master_()
                    fld%bmat(ibf) = msh%cells(ic)%group_()
                END DO
            END DO
        END IF

        ! Field initialization

        ALLOCATE(t0_(nel),x0_(nel),im_(nel),tb0_(nbf),xb0_(nbf),imb_(nbf),stat=info)
        im_(:) = fld%mat(:)
        imb_(:) = fld%bmat(:)
        t0_(:) = std_temp_
        tb0_(:) = std_temp_

        IF(PRESENT(x0)) THEN
            x0_ = x0
            xb0_ = x0
        ELSE
            IF(    fdim == temperature_)      THEN
                x0_  = std_temp_
                xb0_ = std_temp_
            ELSEIF(fdim == density_)       THEN
                CALL matlaw(mats,im_,t0_,dim,x0_)
                CALL matlaw(mats,imb_,tb0_,dim,xb0_)
            ELSEIF(fdim == conductivity_)  THEN
                CALL matlaw(mats,im_,t0_,dim,x0_)
                CALL matlaw(mats,imb_,tb0_,dim,xb0_)
            ELSEIF(fdim == surface_/time_)  THEN
                CALL matlaw(mats,im_,t0_,dim,x0_)
                CALL matlaw(mats,imb_,tb0_,dim,xb0_)
            ELSEIF(fdim == specific_heat_) THEN
                CALL matlaw(mats,im_,t0_,dim,x0_)
                CALL matlaw(mats,imb_,tb0_,dim,xb0_)
            ELSEIF(fdim == youngs_modulus_)  THEN
                CALL matlaw(mats,im_,t0_,dim,x0_)
                CALL matlaw(mats,imb_,tb0_,dim,xb0_)
            ELSEIF(fdim == therm_exp_coeff_) THEN
                CALL matlaw(mats,im_,t0_,dim,x0_)
                CALL matlaw(mats,imb_,tb0_,dim,xb0_)
            ELSEIF(fdim == viscosity_) THEN
                CALL matlaw(mats,im_,t0_,dim,x0_)
                CALL matlaw(mats,imb_,tb0_,dim,xb0_)
            ELSE
                x0_ = 0.d0
                xb0_ = 0.d0
            END IF
        END IF

        fld%x(:)  = x0_
        fld%bx(:) = xb0_
        fld%xp(:) = x0_

        ! WARNING!
        ! Cell-centered field: fld%x is numbered according to msh%desc_c
        ! Face-centered field: fld%x contains first flag = 0 then elements,
        ! then flag = -1 elements, as returned by msh%f2b.

100     FORMAT(' ERROR! Memory allocation failure in CREATE_SCALAR_FIELD')

    END PROCEDURE create_scalar_field

    ! ----- Destructor -----

    MODULE PROCEDURE free_field
        !
        INTEGER :: info

        CALL fld%field%free_field()

        DEALLOCATE(fld%x,fld%xp,fld%bx,stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Memory allocation failure in FREE_SCALAR_FIELD')

    END PROCEDURE free_field


    ! ----- Getters for Inherited Members -----

    MODULE PROCEDURE get_scalar_field_mat_id

        get_scalar_field_mat_id = fld%mat(i)

    END PROCEDURE get_scalar_field_mat_id

    MODULE PROCEDURE get_base

        base = fld%field

    END PROCEDURE get_base


    ! ----- Getters for Additional Members -----

    MODULE PROCEDURE get_scalar_field_x
        INTEGER :: info

        ALLOCATE(x(SIZE(fld%x)),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        x = fld%x

100     FORMAT(' ERROR! Memory allocation failure in GET_SCALAR_FIELD_X')

    END PROCEDURE get_scalar_field_x

    MODULE PROCEDURE get_scalar_field_element

        x=fld%x(i)
    END PROCEDURE get_scalar_field_element

    MODULE PROCEDURE get_scalar_field_xp

        INTEGER :: info

        ALLOCATE(xp(SIZE(fld%xp)),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        xp = fld%xp

100     FORMAT(' ERROR! Memory allocation failure in GET_SCALAR_FIELD_XP')

    END PROCEDURE get_scalar_field_xp

    MODULE PROCEDURE get_scalar_field_element_prev

        xp=fld%xp(i)

    END PROCEDURE get_scalar_field_element_prev


    MODULE PROCEDURE get_scalar_field_bx

        INTEGER :: info

        ALLOCATE(bx(SIZE(fld%bx)),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        bx = fld%bx

100     FORMAT(' ERROR! Memory allocation failure in GET_SCALAR_FIELD_BX')

    END PROCEDURE get_scalar_field_bx

    ! ----- Setters -----

    MODULE PROCEDURE update_scalar_field
        USE class_bc
        USE class_connectivity
        USE class_dimensions
        USE class_face
        USE class_material
        USE class_mesh
        USE tools_math

        INTEGER :: i, ib, IF, im, is, info, err_act, n
        INTEGER, POINTER :: if2b(:) => NULL()
        REAL(psb_dpk_) :: w
        REAL(psb_dpk_), ALLOCATABLE :: fld_c(:)
        TYPE(mesh),     POINTER :: msh   => NULL()
        TYPE(bc_poly),  POINTER :: bc(:) => NULL()
        TYPE(material), POINTER :: mat   => NULL()

        ! The temperature TEMP acts as a possible independent variable

        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)

        ! Copy fld%x to previous timestep value
        fld%xp = fld%x

        ! Gets pointer base-class members
!!$    msh => msh_(fld%field)
        CALL fld%field%get_mesh(msh)
        bc  => fld%field%bc_()


        ! Preliminary checks based on TEMP
        IF(PRESENT(temp)) THEN
            IF(temp%dim_() /= temperature_) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            ALLOCATE(fld_c(SIZE(msh%cells)),stat=info)
            IF(info /= 0) THEN
                WRITE(*,300)
                CALL abort_psblas
            END IF

            ! TEMP must be FACE-CENTERED
            IF(temp%on_faces_()) THEN
                WRITE(*,400)
                CALL abort_psblas
            END IF
        END IF


        ! 4-ways router
        IF(.NOT.fld%on_faces_() .AND. .NOT.PRESENT(temp)) THEN
            ! Cell-centered unknown

            CALL psb_halo(fld%x,msh%desc_c,info)
            CALL psb_check_error(info,'updscalar_field','psb_halo',icontxt_())

            DO ib = 1, msh%nbc
                CALL update_boundary(ib,bc(ib),fld%dim_(),msh,mats,fld%mat,fld%x,fld%bx)
            END DO

        ELSEIF(.NOT.fld%on_faces_() .AND. PRESENT(temp)) THEN
            ! Cell-centered phys. prop.

            CALL matlaw(mats,temp%mat,temp%x,fld%dim_(),fld%x)
            CALL matlaw(mats,temp%bmat,temp%bx,fld%dim_(),fld%bx)

        ELSEIF(fld%on_faces_() .AND. .NOT.PRESENT(temp)) THEN
            ! Face-centered unknown

            WRITE(*,200)
            CALL abort_psblas

        ELSEIF(fld%on_faces_() .AND. PRESENT(temp)) THEN
            ! Face-centered phys. prop.

            ! First compute cell-centered field...
            CALL matlaw(mats,temp%mat,temp%x,fld%dim_(),fld_c)

            ! ... then Interpolate on faces
            CALL msh%f2b%get_ith_conn(if2b,0)
            n = SIZE(if2b) ! # of internal fluid faces

            DO i = 1, n
                IF = if2b(i)
                im = msh%faces(IF)%master_()
                is = msh%faces(IF)%slave_()
                w  = msh%interp(IF)
                fld%x(i) = lin_interp(fld_c(im),fld_c(is),w)
            END DO

            CALL matlaw(mats,temp%bmat,temp%bx,fld%dim_(),fld%bx)
        END IF

        IF(ALLOCATED(fld_c)) DEALLOCATE(fld_c)
        NULLIFY(if2b)
        NULLIFY(mat)
        NULLIFY(bc)
        NULLIFY(msh)


        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Illegal dimensions of TEMP field in UPDATE_SCALAR_FIELD')
200     FORMAT(' ERROR! Scalar unknown field is face-centered.',/,&
            & ' Illegal call to UPDATE_SCALAR_FIELD. Use INTERP_FIELD instead!')
300     FORMAT(' ERROR! Memory allocation failure in UPDATE_SCALAR_FIELD')
400     FORMAT(' ERROR! Face-centered temperature field in UPDATE_SCALAR_FIELD')

    END PROCEDURE update_scalar_field


    MODULE PROCEDURE assign_scalar_field_s

        f%x = x

    END PROCEDURE assign_scalar_field_s


    MODULE PROCEDURE assign_scalar_field_v

        IF(SIZE(f%x) /= SIZE(x)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        f%x = x

100     FORMAT(' ERROR! Size mismatch in ASSIGN_SCALAR_FIELD_V')

    END PROCEDURE assign_scalar_field_v

    MODULE PROCEDURE set_scalar_field_element

        IF(i < 1 .OR. i > SIZE(f%x)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        f%x(i) = x

100     FORMAT(' ERROR! Out of bounds index in SET_SCALAR_FIELD_ELEMENT')

    END PROCEDURE set_scalar_field_element


    MODULE PROCEDURE set_scalar_field_group
        USE class_connectivity
        USE class_mesh

        INTEGER :: i, ic, ncg
        INTEGER, POINTER :: ic2g(:) => NULL()
        TYPE(mesh), POINTER :: msh  => NULL()

        ! MESH pointer dereferencing
!!$    msh => msh_(f%field)
        CALL f%field%get_mesh(msh)

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

200     FORMAT(' ERROR! F field in SET_SCALAR_FIELD_GROUP is face-centered')
100     FORMAT(' ERROR! Group ',i1,' in SET_SCALAR_FIELD_GROUP',&
            & ' is out of bound')

    END PROCEDURE set_scalar_field_group


    ! ----- Algebra Operations -----


    MODULE PROCEDURE nemo_scalar_field_normi
        USE class_mesh
        USE tools_psblas
        INTEGER :: i, ib, IF, im, is, info, err_act, n
        TYPE(mesh),     POINTER :: msh

        ! The temperature TEMP acts as a possible independent variable

        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)


        ! Gets pointer base-class members
!!$    msh => msh_(fld%field)
        CALL fld%field%get_mesh(msh)

        norm = psb_geamax(fld%x,msh%desc_c,info)


        NULLIFY(msh)


        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)


    END PROCEDURE nemo_scalar_field_normi

    MODULE PROCEDURE nemo_scalar_field_norm1
        USE class_mesh
        USE tools_psblas

        INTEGER :: i, ib, IF, im, is, info, err_act, n
        TYPE(mesh),     POINTER :: msh

        ! The temperature TEMP acts as a possible independent variable

        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)


        ! Gets pointer base-class members
!!$    msh => msh_(fld%field)
        CALL fld%field%get_mesh(msh)

        norm = psb_geasum(fld%x,msh%desc_c,info)


        NULLIFY(msh)


        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)


    END PROCEDURE nemo_scalar_field_norm1

    MODULE PROCEDURE scalar_field_sum
        USE class_dimensions

        INTEGER :: isize(2), nel, nbf, info


        ! Check consistency of operands
        CALL f1%field%check_field_operands(f2%field,'SCALAR_FIELD_SUM')

        ! Check consistency of operand dimensions
        IF(f1%dim_() /= f2%dim_()) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        r%field = f1%field

        isize = f1%field%fld_size()
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

100     FORMAT(' ERROR! The dimensions of the operands in SCALAR_FIELD_SUM',&
            & ' are not consistent')
200     FORMAT(' ERROR! Memory allocation failure in SCALAR_FIELD_SUM')

    END PROCEDURE scalar_field_sum


    MODULE PROCEDURE scalar_field_dif
        USE class_dimensions

        INTEGER :: isize(2), nel, nbf, info


        ! Check consistency of operands
        CALL f1%field%check_field_operands(f2%field,'SCALAR_FIELD_DIF')

        ! Check consistency of operand dimensions
        IF(f1%dim_() /= f2%dim_()) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        r%field = f1%field

        isize = f1%field%fld_size()
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

100     FORMAT(' ERROR! The dimensions of the operands in SCALAR_FIELD_DIF',&
            & ' are not consistent')
200     FORMAT(' ERROR! Memory allocation failure in SCALAR_FIELD_DIF')

    END PROCEDURE scalar_field_dif

    MODULE PROCEDURE scalar_field_dif_s
        INTEGER :: isize(2), nel, nbf, info

        r%field = f1%field

        isize = f1%field%fld_size()
        nel   = isize(fld_internal_)
        nbf   = isize(fld_boundary_)

        ! Allocates arrays for inner and boundary elements
        ALLOCATE(r%x(nel),r%bx(nbf),stat=info)
        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF
        r%x  = f1%x - f2
        r%bx = f1%bx - f2

200     FORMAT(' ERROR! Memory allocation failure in SCALAR_FIELD_DIF_S')

    END PROCEDURE scalar_field_dif_s


    MODULE PROCEDURE scalar_field_scal
        !
        INTEGER :: isize(2), nel, nbf, info

        r%field = f%field

        isize = f%field%fld_size()
        nel   = isize(fld_internal_)
        nbf   = isize(fld_boundary_)

        ! Allocates arrays for inner and boundary elements
        ALLOCATE(r%x(nel),r%bx(nbf),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        r%x  = a * f%x
        r%bx = a * f%bx

100     FORMAT(' ERROR! Memory allocation failure in SCALAR_FIELD_SCAL')

    END PROCEDURE scalar_field_scal


    MODULE PROCEDURE scalar_field_mul
        USE class_dimensions

        INTEGER :: isize(2), nel, nbf, info
        TYPE(dimensions) :: dim


        ! Check consistency of operands
        CALL f1%field%check_field_operands(f2%field,'SCALAR_FIELD_MUL')

        r%field = f1%field

        ! Computes result dimensions
        dim = f1%dim_() * f2%dim_()

        ! Sets DIM member in the base field object
        CALL r%field%set_field_dim(dim)

        isize = f1%field%fld_size()
        nel   = isize(fld_internal_)
        nbf   = isize(fld_boundary_)

        ! Allocates arrays for inner and boundary elements
        ALLOCATE(r%x(nel),r%bx(nbf),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        r%x  = f1%x * f2%x
        r%bx = f1%bx * f2%bx

100     FORMAT(' ERROR! Memory allocation failure in SCALAR_FIELD_MUL')

    END PROCEDURE scalar_field_mul


    MODULE PROCEDURE scalar_field_div
        USE class_dimensions

        INTEGER :: isize(2), nel, nbf, info
        TYPE(dimensions) :: dim


        ! Check consistency of operands
        CALL f1%field%check_field_operands(f2%field,'SCALAR_FIELD_DIV')

        r%field = f1%field

        ! Computes result dimensions
        dim = f1%dim_() / f2%dim_()

        ! Sets DIM member in the base field object
        CALL r%field%set_field_dim(dim)

        isize = f1%field%fld_size()
        nel   = isize(fld_internal_)
        nbf   = isize(fld_boundary_)

        ! Allocates arrays for inner and boundary elements
        ALLOCATE(r%x(nel),r%bx(nbf),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        r%x  = f1%x / f2%x
        r%bx = f1%bx / f2%bx

100     FORMAT(' ERROR! Memory allocation failure in SCALAR_FIELD_DIV')

    END PROCEDURE scalar_field_div


    MODULE PROCEDURE interp_on_faces_s
        USE class_connectivity
        USE class_face
        USE class_mesh
        USE tools_math

        INTEGER :: i, IF, im, info, is, ioffset
        INTEGER :: k, n, nel, nbf
        INTEGER, POINTER :: if2b(:) => NULL()
        REAL(psb_dpk_) :: w
        TYPE(mesh), POINTER :: msh => NULL()


        IF(fld%on_faces_()) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! Sets base member
        r%field = fld%field

        CALL fld%field%get_mesh(msh)
        nel = COUNT(msh%faces%flag_() <= 0)
        nbf = COUNT(msh%faces%flag_() > 0)

        ! Allocates arrays for inner and boundary elements
        ALLOCATE(r%x(nel),r%bx(nbf),stat=info)
        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

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
            r%x(k) = fld%x(MAX(im,is))
        END DO


        ! 3) Copies boundary values
        r%bx = fld%bx
        CALL r%field%set_field_on_faces(.TRUE.)
        NULLIFY(if2b)
        NULLIFY(msh)

100     FORMAT(' ERROR! Argument of INTERP_ON_FACES_S is already on faces')
200     FORMAT(' ERROR! Memory allocation failure in INTERP_ON_FACES_S')

    END PROCEDURE interp_on_faces_s


    ! ----- Check Procedures -----

    MODULE PROCEDURE check_mesh_consistency_sf

        CALL f1%field%check_mesh_consistency(f2%field,WHERE)

    END PROCEDURE check_mesh_consistency_sf

END SUBMODULE class_scalar_field_procedures
