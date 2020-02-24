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
! $Id: class_bc_wall.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    WALL boundary condition class.
!
SUBMODULE(class_bc_wall) class_bc_wall_procedures

    USE class_psblas
    USE class_bc_math
    USE class_material
    USE class_vector
    USE class_face

    IMPLICIT NONE

CONTAINS

    ! REMARK: the implementation of run-time polymorphism requires
    ! specific BC object as POINTERS!

    MODULE PROCEDURE nemo_bc_wall_sizeof
        USE psb_base_mod

        INTEGER(kind=nemo_int_long_)   :: val

        val = 2 * nemo_sizeof_int
        val = val + bc%temp%nemo_sizeof() + SUM(bc%vel%nemo_sizeof())
        nemo_bc_wall_sizeof = val

    END PROCEDURE nemo_bc_wall_sizeof

    ! ----- Constructor -----

    MODULE PROCEDURE create_bc_wall
        IMPLICIT NONE
        !
        INTEGER :: info

        ! Alloc bc_wall target on every process

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
        CALL rd_inp_bc_wall(input_file,sec,nbf,bc%id,bc%temp,bc%conc,bc%vel,bc%stress)

100     FORMAT(' ERROR! Illegal call to CREATE_BC_WALL: pointer already associated')
200     FORMAT(' ERROR! Memory allocation failure in CREATE_BC_WALL')

    END PROCEDURE create_bc_wall


    ! ----- Destructor -----

    MODULE PROCEDURE free_bc_wall
        IMPLICIT NONE
        !
        INTEGER :: info

        IF(bc%temp%is_allocated()) CALL bc%temp%dealloc_bc_math()
        IF(ANY(bc%vel%is_allocated()))  THEN
            CALL bc%vel(1)%dealloc_bc_math()
            CALL bc%vel(2)%dealloc_bc_math()
            CALL bc%vel(3)%dealloc_bc_math()
        END IF

        DEALLOCATE(bc,stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Memory allocation failure in FREE_BC_WALL')

    END PROCEDURE free_bc_wall


    ! ----- Getter -----

    MODULE PROCEDURE get_abc_wall_s
        USE class_dimensions
        IMPLICIT NONE

        IF(dim == temperature_) THEN
            CALL get_abc_math(bc%temp,id,a,b,c)
        ELSEIF(dim == density_) THEN
            CALL get_abc_math(bc%temp,id,a,b,c)
        ELSE

            WRITE(*,100)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Unsupported field dimensions in GET_ABC_WALL_S')

    END PROCEDURE get_abc_wall_s


    MODULE PROCEDURE get_abc_wall_v
        USE class_dimensions
        IMPLICIT NONE

        IF(dim == velocity_) THEN
            CALL get_abc_math(bc%vel,id,a,b,c)
        ELSE IF(dim == length_) THEN
            CALL get_abc_math(bc%vel,id,a,b,c)
        ELSE IF(dim == pressure_) THEN
            CALL get_abc_math(bc%stress,id,a,b,c)
        ELSE
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! REMARK: BC(:) elements are supposed to differ only in "C" term

100     FORMAT(' ERROR! Unsupported field dimensions in GET_ABC_WALL_V')

    END PROCEDURE get_abc_wall_v

    ! ----- Setter -----

    MODULE PROCEDURE set_bc_wall_map_s
        USE tools_bc
        USE class_bc_math
        IMPLICIT NONE

        SELECT CASE(bc%id(1))
        CASE(bc_temp_convection_map_)
            CALL bc%temp%set_bc_math_map(i,a,b,c)
        CASE DEFAULT
            WRITE(*,100)
            CALL abort_psblas
        END SELECT

100     FORMAT(' ERROR! Unsupported BC type in SET_BC_WALL_MAP')

    END PROCEDURE set_bc_wall_map_s


    MODULE PROCEDURE set_bc_wall_map_v
        USE tools_bc
        USE class_bc_math
        IMPLICIT NONE

        SELECT CASE(bc%id(2))
        CASE(bc_vel_free_sliding_)
            CALL bc%vel(1)%set_bc_math_map(i,a,b,c%x_())
            CALL bc%vel(2)%set_bc_math_map(i,a,b,c%y_())
            CALL bc%vel(3)%set_bc_math_map(i,a,b,c%x_())
        CASE DEFAULT
            WRITE(*,100)
            CALL abort_psblas
        END SELECT

100     FORMAT(' ERROR! Unsupported BC type in SET_BC_WALL_MAP')

    END PROCEDURE set_bc_wall_map_v

    ! ----- Boundary Values Updater -----

    MODULE PROCEDURE update_boundary_wall_s
        USE class_dimensions
        USE class_material
        USE class_mesh
        USE tools_bc
        IMPLICIT NONE
        !
        INTEGER :: i, id, IF, info, ib_offset, n
        REAL(psb_dpk_), ALLOCATABLE :: a(:), b(:), c(:)

        ! WARNING!
        ! - Use intent(inout) for BX(:).

        ! Number of boundary faces with flag < IB
        ib_offset = COUNT(msh%faces%flag_() > 0 .AND. msh%faces%flag_() < ib)

        n = COUNT(msh%faces%flag_() == ib)

        ALLOCATE(a(n),b(n),c(n),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        IF(dim == temperature_) THEN
            CALL get_abc_math(bc%temp,id,a,b,c)
            IF(  id == bc_neumann_flux_ .OR. &
                id == bc_robin_convection_) THEN
                DO i = 1, n
                    IF = ib_offset + i
                    CALL matlaw(mats,im(i),bx(IF),conductivity_,b(i))
                END DO
            END IF
        ELSEIF(dim == density_) THEN
            CALL get_abc_math(bc%temp,id,a,b,c)
            IF(  id == bc_neumann_flux_ .OR. &
                id == bc_robin_convection_) THEN
                DO i = 1, n
                    IF = ib_offset + i
                    CALL matlaw(mats,im(i),bx(IF),conductivity_,b(i))
                END DO
            END IF
        ELSE

            WRITE(*,200)
            CALL abort_psblas
        END IF

        CALL apply_abc_to_boundary(id,a,b,c,ib,msh,x,bx)

        DEALLOCATE(a,b,c,stat=info)
        IF(info /= 0) THEN
            WRITE(*,300)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Memory allocation failure in UPDATE_BOUNDARY_WALL_S')
200     FORMAT(' ERROR! Unsupported field dimensions in UPDATE_BOUNDARY_WALL_S')
300     FORMAT(' ERROR! Memory deallocation failure in UPDATE_BOUNDARY_WALL_S')

    END PROCEDURE update_boundary_wall_s


    MODULE PROCEDURE update_boundary_wall_v
        USE class_dimensions
        USE class_mesh
        USE tools_bc
        IMPLICIT NONE
        !
        INTEGER :: id, info, ib_offset, n
        REAL(psb_dpk_), ALLOCATABLE :: a(:), b(:)
        TYPE(vector),     ALLOCATABLE :: c(:)

        ! WARNING!
        ! - Use intent(inout) for BX(:).

        ! Number of boundary faces with flag < IB
        ib_offset = COUNT(msh%faces%flag_() > 0 .AND. msh%faces%flag_() < ib)

        n = COUNT(msh%faces%flag_() == ib)

        ALLOCATE(a(n),b(n),c(n),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        IF(dim == velocity_ .OR. dim == length_) THEN
            CALL get_abc_math(bc%vel,id,a,b,c)
            CALL apply_abc_to_boundary(id,a,b,c,ib,msh,x,bx)
        ELSE IF(dim == pressure_) THEN
            CALL get_abc_math(bc%stress,id,a,b,c)
            CALL apply_abc_to_boundary(id,a,b,c,ib,msh,x,bx)
        ELSE
            WRITE(*,200)
            CALL abort_psblas
        END IF

        DEALLOCATE(a,b,c,stat=info)
        IF(info /= 0) THEN
            WRITE(*,300)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Memory allocation failure in UPDATE_BOUNDARY_WALL_V')
200     FORMAT(' ERROR! Unsupported field dimensions in UPDATE_BOUNDARY_WALL_V')
300     FORMAT(' ERROR! Memory deallocation failure in UPDATE_BOUNDARY_WALL_V')

    END PROCEDURE update_boundary_wall_v

SUBROUTINE rd_inp_bc_wall(input_file,sec,nbf,id,bc_temp,bc_conc,bc_vel,bc_stress)
    USE class_psblas
    USE json_module, ONLY : json_file
    USE class_bc_math
    USE tools_bc
    USE tools_input

    IMPLICIT NONE
    !
    TYPE(json_file) :: nemo_json
    CHARACTER(len=*), INTENT(IN) :: input_file
    CHARACTER(len=*), INTENT(IN) :: sec
    INTEGER, INTENT(IN) :: nbf
    INTEGER, INTENT(OUT) :: id(:)
    TYPE(bc_math), INTENT(OUT) :: bc_temp
    TYPE(bc_math), INTENT(OUT) :: bc_conc
    TYPE(bc_math), INTENT(OUT) :: bc_vel(3)
    TYPE(bc_math), INTENT(OUT) :: bc_stress(3)
    !
    LOGICAL, PARAMETER :: debug = .FALSE.
    LOGICAL :: found
    !
    INTEGER :: mypnum, icontxt
    INTEGER :: id_sec, inp
    REAL(psb_dpk_) :: work(3,5), wtemp
    CHARACTER(len=15) :: par

    icontxt = icontxt_()
    mypnum  = mypnum_()

    id   = -1 ! DEFAULT value
    work = 0.d0

    IF(mypnum == 0) THEN

        CALL open_file(input_file,nemo_json)
        WRITE(*,*) '- Reading ', TRIM(sec), ' section: type WALL'

        CALL nemo_json%get(trim(sec)//'.temperature.id', id_sec, found)

        IF (.NOT.found) THEN
            WRITE(*,*) 'Temperature BC not found in RD_INP_BC_WALL'
            CALL abort_psblas
        ELSE
            id(bc_temp_) = id_sec
            SELECT CASE(id_sec)
            CASE(bc_temp_fixed_)      ! Fixed temperature
                CALL nemo_json%get(trim(sec)//'.temperature.value', work(1,bc_temp_), found)
            CASE(bc_temp_adiabatic_)  ! Adiabatic wall
                !READ(inp,'()')

            CASE(bc_temp_flux_)       ! Fixed heat flux
                CALL nemo_json%get(trim(sec)//'.temperature.value', work(1,bc_temp_), found)

            CASE(bc_temp_convection_) ! Convection
                CALL nemo_json%get(trim(sec)//'.temperature.value', work(1,bc_temp_), found)
                CALL nemo_json%get(trim(sec)//'.temperature.value2', work(2,bc_temp_), found)

            CASE(bc_temp_convection_map_) ! Convection map
                !READ(inp,'()')         ! To be actually set at the first mapping

            CASE DEFAULT
                WRITE(*,210)
                CALL abort_psblas
            END SELECT
        END IF

        CALL nemo_json%get(trim(sec)//'.concentration.id', id_sec, found)

        IF (.NOT.found) THEN
            WRITE(*,*) 'Concentration BC not found in RD_INP_BC_WALL'
            ! CALL abort_psblas
        ELSE
            id(bc_conc_) = id_sec
            SELECT CASE(id_sec)
            CASE(bc_conc_fixed_)      ! Fixed temperature
                CALL nemo_json%get(trim(sec)//'.concetration.value', work(1,bc_conc_), found)

            CASE(bc_conc_adiabatic_)  ! Adiabatic wall
                !READ(inp,'()')

            CASE DEFAULT
                WRITE(*,210)
                CALL abort_psblas
            END SELECT
        END IF

        CALL nemo_json%get(trim(sec)//'.velocity.id', id_sec, found)

        IF (.NOT.found) THEN
            WRITE(*,*) 'Velocity BC not found in RD_INP_BC_WALL'
        ELSE
            id(bc_vel_) = id_sec

            SELECT CASE(id_sec)
            CASE(bc_vel_no_slip_, bc_vel_free_slip_, bc_vel_free_sliding_)
                !READ(inp,'()')
            CASE(bc_vel_sliding_)
                CALL nemo_json%get(trim(sec)//'.velocity.v1', work(1,bc_vel_), found)
                CALL nemo_json%get(trim(sec)//'.velocity.v2', work(2,bc_vel_), found)
                CALL nemo_json%get(trim(sec)//'.velocity.v3', work(3,bc_vel_), found)
                WRITE(0,*) 'rd_inp_bc_wall: Debug: work', work(1:3,bc_vel_)
            CASE(bc_vel_moving_)
                !READ(inp,'()')
            CASE DEFAULT
                WRITE(*,220)
                CALL abort_psblas
            END SELECT
        END IF


        CALL nemo_json%get(trim(sec)//'.stress.id', id_sec, found)

        IF (.NOT.found) THEN
            WRITE(*,*) 'Stress BC not found in RD_INP_BC_WALL'
        ELSE
            id(bc_stress_) = id_sec

            SELECT CASE(id_sec)
            CASE(bc_stress_free_)
                !READ(inp,'()')
            CASE(bc_stress_prescribed_)
                CALL nemo_json%get(trim(sec)//'.stress.v1', work(1,bc_stress_), found)
                CALL nemo_json%get(trim(sec)//'.stress.v2', work(2,bc_stress_), found)
                CALL nemo_json%get(trim(sec)//'.stress.v3', work(3,bc_stress_), found)
                WRITE(0,*) 'rd_inp_bc_wall: Debug: work', work(1:3,bc_stress_)
            CASE DEFAULT
                WRITE(*,230)
                CALL abort_psblas
            END SELECT
        END IF
    END IF


    ! Broadcast
    CALL psb_bcast(icontxt,id)
    CALL psb_bcast(icontxt,work)


    ! TEMPERATURE section
    SELECT CASE(id(bc_temp_))
    CASE(bc_temp_fixed_)      ! Fixed temperature
        CALL bc_temp%alloc_bc_math(bc_dirichlet_,nbf,&
            & a=(/1.d0/),b=(/0.d0/),c=(/work(1,bc_temp_)/))

    CASE(bc_temp_adiabatic_)  ! Adiabatic wall
        CALL bc_temp%alloc_bc_math(bc_neumann_,nbf,&
            & a=(/0.d0/),b=(/1.d0/),c=(/0.d0/))

    CASE(bc_temp_flux_)       ! Fixed heat flux
        CALL bc_temp%alloc_bc_math(bc_neumann_flux_,nbf,&
            & a=(/0.d0/),b=(/1.d0/),c=(/work(1,bc_temp_)/))

    CASE(bc_temp_convection_) ! Convection
        CALL bc_temp%alloc_bc_math(bc_robin_convection_,nbf,&
            & a=(/work(1,bc_temp_)/),b=(/1.d0/),&
            & c=(/work(1,bc_temp_)*work(2,bc_temp_)/))

    CASE(bc_temp_convection_map_) ! Convection map
        CALL bc_temp%alloc_bc_math(bc_robin_map_,nbf,&
            & a=(/(0.d0, inp=1,nbf)/),b=(/(1.d0, inp=1,nbf)/),&
            & c=(/(0.d0, inp=1,nbf)/))

    END SELECT


    ! VELOCITY section
    SELECT CASE(id(bc_vel_))
    CASE(bc_vel_no_slip_)

        CALL bc_vel(1)%alloc_bc_math(bc_dirichlet_,nbf,&
            & a=(/1.d0/),b=(/0.d0/),c=(/0.d0/))
        CALL bc_vel(2)%alloc_bc_math(bc_dirichlet_,nbf,&
            & a=(/1.d0/),b=(/0.d0/),c=(/0.d0/))
        CALL bc_vel(3)%alloc_bc_math(bc_dirichlet_,nbf,&
            & a=(/1.d0/),b=(/0.d0/),c=(/0.d0/))

    CASE(bc_vel_free_slip_)

        CALL bc_vel(1)%alloc_bc_math(bc_neumann_,nbf,&
            & a=(/0.d0/),b=(/1.d0/),c=(/0.d0/))
        CALL bc_vel(2)%alloc_bc_math(bc_neumann_,nbf,&
            & a=(/0.d0/),b=(/1.d0/),c=(/0.d0/))
        CALL bc_vel(3)%alloc_bc_math(bc_neumann_,nbf,&
            & a=(/0.d0/),b=(/1.d0/),c=(/0.d0/))

    CASE(bc_vel_sliding_)
        WRITE(0,*) 'rd_inp_bc_wall: Debug: setting', work(1:3,bc_vel_)
        CALL bc_vel(1)%alloc_bc_math(bc_dirichlet_,nbf,&
            & a=(/1.d0/),b=(/0.d0/),c=(/work(1,bc_vel_)/))
        CALL bc_vel(2)%alloc_bc_math(bc_dirichlet_,nbf,&
            & a=(/1.d0/),b=(/0.d0/),c=(/work(2,bc_vel_)/))
        CALL bc_vel(3)%alloc_bc_math(bc_dirichlet_,nbf,&
            & a=(/1.d0/),b=(/0.d0/),c=(/work(3,bc_vel_)/))
        WRITE(0,*) 'rd_inp_bc_wall: Debug: set'

    CASE(bc_vel_free_sliding_)

        CALL bc_vel(1)%alloc_bc_math(bc_dirichlet_map_,nbf,&
            & a=(/1.d0/),b=(/0.d0/),c=(/(0.d0, inp=1,nbf)/))
        CALL bc_vel(2)%alloc_bc_math(bc_dirichlet_map_,nbf,&
            & a=(/1.d0/),b=(/0.d0/),c=(/(0.d0, inp=1,nbf)/))
        CALL bc_vel(3)%alloc_bc_math(bc_dirichlet_map_,nbf,&
            & a=(/1.d0/),b=(/0.d0/),c=(/(0.d0, inp=1,nbf)/))
    CASE(-1)
    CASE DEFAULT
        WRITE(0,*) 'Fix here unimplemented BC setup !!'

    END SELECT

    ! Stress section
    SELECT CASE(id(bc_stress_))
    CASE(bc_stress_free_)
        CALL bc_stress(1)%alloc_bc_math(bc_neumann_,nbf,&
            & a=(/0.d0/),b=(/1.d0/),c=(/0.d0/))
        CALL bc_stress(2)%alloc_bc_math(bc_neumann_,nbf,&
            & a=(/0.d0/),b=(/1.d0/),c=(/0.d0/))
        CALL bc_stress(3)%alloc_bc_math(bc_neumann_,nbf,&
            & a=(/0.d0/),b=(/1.d0/),c=(/0.d0/))

    CASE(bc_stress_prescribed_)
        CALL bc_stress(1)%alloc_bc_math(bc_dirichlet_,nbf,&
            & a=(/1.d0/),b=(/0.d0/),c=(/work(1,bc_stress_)/))
        CALL bc_stress(2)%alloc_bc_math(bc_dirichlet_,nbf,&
            & a=(/1.d0/),b=(/0.d0/),c=(/work(2,bc_stress_)/))
        CALL bc_stress(3)%alloc_bc_math(bc_dirichlet_,nbf,&
            & a=(/1.d0/),b=(/0.d0/),c=(/work(3,bc_stress_)/))

    CASE DEFAULT
        WRITE(0,*) 'Fix here unimplemented BC setup !!', id(bc_stress_)

    END SELECT


    ! ***
    ! If ID = -1 (Default)  the corresponding BC will not be allocated
    ! ***

    IF(debug) THEN
        WRITE(*,*)
        WRITE(*,400) mypnum
        WRITE(*,500) TRIM(sec),' - Type: Wall'

        WRITE(*,600) ' * TEMPERATURE Section *'
        WRITE(*,700) '  BC%id(bc_temp_) = ', id(bc_temp_)
        IF (id(bc_temp_) > 0) CALL bc_temp%debug_bc_math()

        WRITE(*,600) ' * VELOCITY Section *'
        WRITE(*,700) '  BC%id(bc_vel_) = ', id(bc_vel_)
        IF (id(bc_vel_) > 0) THEN
            CALL bc_vel(1)%debug_bc_math()
            CALL bc_vel(2)%debug_bc_math()
            CALL bc_vel(3)%debug_bc_math()
            WRITE(*,*)
        END IF
    END IF


100 FORMAT(a,i1)
210 FORMAT(' ERROR! Unsupported ID(BC_TEMP_) in RD_INP_BC_WALL')
220 FORMAT(' ERROR! Unsupported ID(BC_VEL_) in RD_INP_BC_WALL')
230 FORMAT(' ERROR! Unsupported ID(BC_STRESS_) in RD_INP_BC_WALL')
999 FORMAT(' ERROR! Missing boundary condition in RC_INP_BC_WALL')

400 FORMAT(' ----- Process ID = ',i2,' -----')
500 FORMAT(1x,a,a)
600 FORMAT(1x,a)
700 FORMAT(1x,a,i2)

END SUBROUTINE rd_inp_bc_wall

END SUBMODULE class_bc_wall_procedures
