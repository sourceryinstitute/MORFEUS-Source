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
! $Id: class_material.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(class_material) class_material_procedures

    USE class_psblas
    USE MatLib, ONLY : MatProp, Set_Material_ids
    IMPLICIT NONE

CONTAINS

    ! ----- Constructors -----

    MODULE PROCEDURE create_material
      !! Global Constructor
        USE tools_material

        ! Reads and bcast_ input parameters
        CALL rd_inp_material(input_file,sec,mat%name,mat%ilaw,mat%mat_type,mat%mat_id)

        ! Load material physical properties
        CALL load_material(mat%name,mat%state,mat%dtemp,mat%tmin,mat%tmax,&
            &             mat%rho,mat%mu,mat%lambda,mat%sh)

        CALL Set_Material_ids (mat%mat_id)

        IF(debug) CALL debug_material(mat)

    END PROCEDURE create_material


    ! ----- Destructor -----

    MODULE PROCEDURE free_material
        !
        INTEGER :: info

        IF(mat%state == 'l' .OR. mat%state == 'g') THEN
            ! Liquid or gas
            DEALLOCATE(mat%rho,mat%mu,mat%lambda,mat%sh,stat=info)
        ELSE
            ! Solid => no viscosity mat%mu
            DEALLOCATE(mat%rho,mat%lambda,mat%sh,stat=info)
        END IF
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Memory deallocation failure in FREE_MATERIAL')

    END PROCEDURE free_material


    ! ----- Getters -----

    MODULE PROCEDURE get_material_name

        get_material_name = mat%name

    END PROCEDURE get_material_name

    MODULE PROCEDURE get_material_id

        get_material_id = mat%mat_id

    END PROCEDURE get_material_id


    ! ----- Physical Properties Laws -----

    MODULE PROCEDURE matlaw_v
        USE class_dimensions
        USE tools_material
        USE tools_math

        !
        INTEGER :: ilaw, i
        REAL(psb_dpk_), POINTER :: DATA(:) => NULL()
        REAL(psb_dpk_) :: f_s

        IF(SIZE(t) /= SIZE(f) .OR. SIZE(im) /= SIZE(f)) THEN
            WRITE(*,100)
            WRITE(*,*) SIZE(f), SIZE(im), SIZE(t)
            CALL abort_psblas
        END IF

        IF (mat_id_(mats(1)%mat) < 0) THEN
            IF(dim == density_) THEN
                ilaw = mats(1)%mat%ilaw(irho)
                DATA => mats(1)%mat%rho
            ELSEIF(dim == viscosity_) THEN
                ilaw = mats(1)%mat%ilaw(imu)
                DATA => mats(1)%mat%mu
            ELSEIF(dim == conductivity_) THEN
                ilaw = mats(1)%mat%ilaw(ilambda)
                DATA => mats(1)%mat%lambda
            ELSEIF(dim == specific_heat_) THEN
                ilaw = mats(1)%mat%ilaw(ish)
                DATA => mats(1)%mat%sh
            ELSE
                WRITE(*,200)
                CALL abort_psblas
            END IF

            CALL check_temp(t,mats(1)%mat)

            SELECT CASE(ilaw)
            CASE(1)
                IF(SIZE(DATA) == 1) THEN
                    ! Just one value is available
                    f(:) = DATA(1)
                ELSE
                    ! Computes pwl_interp for Tref
                    CALL pwl_interp(f_s,std_temp_,DATA,mats(1)%mat%tmin,mats(1)%mat%dtemp)
                    f(:) = f_s
                END IF
            CASE(2)
                ! Computes pwl_interp for T
                CALL pwl_interp(f,t,DATA,mats(1)%mat%tmin,mats(1)%mat%dtemp)
            END SELECT

            NULLIFY(DATA)

        ELSE
            IF(dim == density_)       THEN
                DO i = 1, SIZE(im)
                    IF (mat_id_(mats(im(i))%mat) < 400) THEN
                        CALL matlaw_fast_s(mats(im(i))%mat,t(i),'ASFABDENSITY',f(i))
                    ELSE
                        CALL matlaw_fast_s(mats(im(i))%mat,t(i),'DENSITY',f(i))
                    END IF
                END DO
            ELSEIF(dim == conductivity_)  THEN
                DO i = 1, SIZE(im)
                    CALL matlaw_fast_s(mats(im(i))%mat,t(i),'THERMCOND',f(i))
                END DO
            ELSEIF(dim == surface_/time_)  THEN
                DO i = 1, SIZE(im)
                    CALL matlaw_fast_s(mats(im(i))%mat,t(i),'THERMCOND',f(i))
                END DO
            ELSEIF(dim == specific_heat_) THEN
                DO i = 1, SIZE(im)
                    CALL matlaw_fast_s(mats(im(i))%mat,t(i),'SPECHEAT',f(i))
                END DO
            ELSEIF(dim == youngs_modulus_) THEN
                DO i = 1, SIZE(im)
                    CALL matlaw_fast_s(mats(im(i))%mat,t(i),'YOUNG_MOD',f(i))
                END DO
            ELSEIF(dim == therm_exp_coeff_) THEN
                DO i = 1, SIZE(im)
                    IF (mat_id_(mats(im(i))%mat) < 400) THEN
                        CALL matlaw_fast_s(mats(im(i))%mat,t(i),'THEXP_COEF',f(i))
                    ELSE
                        CALL matlaw_fast_s(mats(im(i))%mat,t(i),'THEXP',f(i))
                    END IF
                END DO
            ELSE
                WRITE(*,300)
                CALL abort_psblas
            END IF
        END IF

        ! REMARK. DATA buffer: alternatively you could use an ALLOCATABLE
        ! array (+ copy) instead of a POINTER (+ dereferencing).

100     FORMAT(' ERROR! Size mismatch in MATLAW_V')
200     FORMAT(' ERROR! Unsupported physical property in MATLAW_V')
300     FORMAT(' ERROR! Unsupported physical property in MATLAW_V for FAST')

    END PROCEDURE matlaw_v


    MODULE PROCEDURE matlaw_s
        USE class_dimensions
        USE tools_math

        INTEGER :: i(1)
        !
        REAL(psb_dpk_) :: t_s(1), f_s(1)

        t_s(1) = t
        i(1) = im
        CALL matlaw_v(mats,i,t_s,dim,f_s)
        f = f_s(1)

    END PROCEDURE matlaw_s


    MODULE PROCEDURE matlaw_fast_s

        f = MatProp(Mat_id=mat%mat_id, Property=property, Temperature=t)

    END PROCEDURE matlaw_fast_s

    ! ----- Check Procedures -----

    MODULE PROCEDURE check_temp_s

        IF(t < mat%tmin .OR. t > mat%tmax) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Out of range in temperature scalar')

    END PROCEDURE check_temp_s


    MODULE PROCEDURE check_temp_v

        IF(ANY(t < mat%tmin) .OR. ANY(t > mat%tmax)) THEN
            PRINT *, t, mat%tmin, mat%tmax
            WRITE(*,100)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Out of range in temperature array')

    END PROCEDURE check_temp_v


    MODULE PROCEDURE check_material_consistency

        IF(.NOT.ASSOCIATED(mat1,mat2)) THEN
            WRITE(*,100) TRIM(WHERE)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! In ', a,': material pointers are not consistent')

    END PROCEDURE check_material_consistency


    ! ----- Debug -----

    MODULE PROCEDURE debug_material
        USE tools_material

        WRITE(*,100) mypnum_()
        WRITE(*,200) 'material%name          = ', mat%name
        WRITE(*,300) 'material%ilaw(irho)    = ', mat%ilaw(irho)
        WRITE(*,300) 'material%ilaw(imu)     = ', mat%ilaw(imu)
        WRITE(*,300) 'material%ilaw(ilambda) = ', mat%ilaw(ilambda)
        WRITE(*,300) 'material%ilaw(ish)     = ', mat%ilaw(ish)
        WRITE(*,200) 'material%state         = ', mat%state

100     FORMAT(' ----- Process ID = ',i2,' -----')
200     FORMAT(1x,a,a)
300     FORMAT(1x,a,i2)

    END PROCEDURE debug_material


    MODULE PROCEDURE nemo_material_sizeof
        USE psb_base_mod

        INTEGER(kind=nemo_int_long_)   :: val

        val = LEN(mat%name) + nemo_sizeof_int * SIZE(mat%ilaw)
        val = val + 3 * nemo_sizeof_dp
        IF (ALLOCATED(mat%rho)) &
            & val = val + nemo_sizeof_dp * SIZE(mat%rho)
        IF (ALLOCATED(mat%mu)) &
            & val = val + nemo_sizeof_dp * SIZE(mat%mu)
        IF (ALLOCATED(mat%lambda)) &
            & val = val + nemo_sizeof_dp * SIZE(mat%lambda)
        IF (ALLOCATED(mat%sh)) &
            & val = val + nemo_sizeof_dp * SIZE(mat%sh)
        nemo_material_sizeof = val

    END PROCEDURE nemo_material_sizeof

END SUBMODULE class_material_procedures
