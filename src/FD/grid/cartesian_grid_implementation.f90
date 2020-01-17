!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
submodule(cartesian_grid_interface) cartesian_grid_implementation
  !! author: Damian Rouson and Karla Morris
  use kind_parameters, only : i4k, r8k
  use assertions_interface,only : assert, max_errmsg_len, assertions
  use plate_3D_interface, only : plate_3D
  use surfaces_interface, only : backward, forward
  use package_interface, only : package, null_sender_id
  implicit none

contains

  module procedure build_surfaces
    integer, parameter :: first=1, last=2, success=0, x_dir=1, y_dir=2, z_dir=3, vec_components=3, space_dimensions=3, num_faces=2
    integer, parameter :: displacement(x_dir:z_dir, backward:forward, x_dir:z_dir) = &
      reshape( [ [-1,0,0], [1,0,0], [0,-1,0], [0,1,0], [0,0,-1], [0,0,1] ], [space_dimensions, num_faces, vec_components ] )
    type(package), allocatable, dimension(:,:,:) :: bare
    character(len=max_errmsg_len) error_message
    integer alloc_stat, b, coord_dir, face_dir
    class(package), allocatable, dimension(:,:,:) :: surface_packages

    select type(problem_geometry)
      type is(plate_3D)
        !! this verifies correct geometry for the assumptions in building the surfaces below
      class default
        error stop "cartesian_grid%build_surfaces: unsupported problem_geometry type"
    end select

    call assert(size(my_blocks)==2, "cartesian_grid%build_surfaces: size(my_blocks)==2")

    allocate( bare(my_blocks(first):my_blocks(last), space_dimension, backward:forward), &
      stat=alloc_stat, errmsg=error_message)
    call assert(alloc_stat==success, "cartesian_grid%build_surfaces: allocate(bare)", error_message)

    call bare%set_sender_block_id(null_sender_id)
    call bare%set_step(0)

     loop_over_blocks: &
     do b=my_blocks(first), my_blocks(last)
       loop_over_coordinate_directions: &
       do coord_dir = x_dir, z_dir
         loop_over_face_directions: &
         do face_dir = backward, forward
           associate( ijk_displaced => this%block_indicial_coordinates(b) + displacement(coord_dir, face_dir, :) )
             if (this%block_in_bounds(ijk_displaced)) then
               call bare(b, coord_dir, face_dir)%set_sender_block_id( this%block_identifier(ijk_displaced) )
             end if
           end associate
         end do loop_over_face_directions
       end do loop_over_coordinate_directions
     end do loop_over_blocks

    call block_faces%set_halo_data(bare, my_blocks)

  end procedure

  module procedure assign_structured_grid
    call assert( same_type_as(this, rhs), "cartesian_grid%assign_structured_grid: consistent types" )
    call this%clone(rhs)
  end procedure

  module procedure set_up_div_scalar_flux
    !! Evaluate 2nd partial derivatives in each coordinate direction according to equations (2.2) & (2.4) in
    !! Sundqvist & Veronis (1969) "A simple finite-difference grid with non-constant intervals", Tellus 22:1

    integer(i4k) i, j, k, alloc_stat
    integer(i4k), parameter :: success=0
    real(r8k), parameter :: half=0.5_r8k
    character(len=max_errmsg_len) :: alloc_error
    real(r8k), allocatable, dimension(:,:,:) :: div_flux_x, div_flux_y, div_flux_z

    call assert( same_type_as(this, vertices), "cartesian_grid%set_up_div_scalar_flux: same_type_as(this, vertices)" )

    associate( positions => vertices%vectors(), s=>this%get_scalar() )
      associate( npoints => shape(positions(:,:,:,1)) )

        allocate(div_flux_x, div_flux_y, div_flux_z, mold=positions(:,:,:,1), stat=alloc_stat, errmsg=alloc_error )
        call assert( alloc_stat==success, "cartesian_grid%set_up_div_scalar_flux: allocate(div_flux_{x,y,z})", alloc_error )

        div_flux_x = 0._r8k
        div_flux_y = 0._r8k
        div_flux_z = 0._r8k

        associate( x=>positions(:,:,:,1) )
          do concurrent(k=1:npoints(3), j=1:npoints(2), i=2:npoints(1)-1)

            associate( &
              dx_m => half*(x(i+1,j,k) - x(i-1,j,k)), &!! half*(x(i+1,j,k) - x(i,j,k) + x(i,j,k) - x(i-1,j,k))
              dx_f =>       x(i+1,j,k) - x(i,j,k), &
              dx_b =>       x(i,j,k)   - x(i-1,j,k), &
              s_f => half*( s(i+1,j,k) + s(i,j,k)  ), &
              s_b => half*( s(i,j,k)   + s(i-1,j,k) ))
              associate( &
                D_f => this%diffusion_coefficient( s_f ), &
                D_b => this%diffusion_coefficient( s_b) )

                div_flux_x(i,j,k) = &
                  D_f*(s(i+1,j,k) - s(i,j,k)  )/(dx_f*dx_m) - &
                  D_b*(s(i,j,k)   - s(i-1,j,k))/(dx_b*dx_m)
              end associate
            end associate
          end do
        end associate

        associate( y=>positions(:,:,:,2))
          do concurrent(k=1:npoints(3), j=2:npoints(2)-1, i=1:npoints(1))
            associate( &
              dy_m => (y(i,j+1,k) - y(i,j-1,k))*half, &
              dy_f =>  y(i,j+1,k) - y(i,j,k ), &
              dy_b =>  y(i,j,k)   - y(i,j-1,k), &
              s_f => half*(s(i,j+1,k) + s(i,j,k)), &
              s_b =>  half*(s(i,j,k) + s(i,j-1,k)) )

              associate( &
                D_f => this%diffusion_coefficient( s_f ), &
                D_b => this%diffusion_coefficient( s_b) )

                div_flux_y(i,j,k) = &
                  D_f*(s(i,j+1,k) - s(i,j,k))/(dy_f*dy_m) - &
                  D_b*(s(i,j,k) - s(i,j-1,k))/(dy_b*dy_m)
              end associate
            end associate
          end do
        end associate

        associate( z=>positions(:,:,:,3))
          do concurrent(k=2:npoints(3)-1, j=1:npoints(2), i=1:npoints(1))
            associate( &
              dz_m => (z(i,j,k+1) - z(i,j,k-1))*half, &
              dz_f =>  z(i,j,k+1) - z(i,j,k), &
              dz_b =>  z(i,j,k)   - z(i,j,k-1), &
              s_f => half*(s(i,j,k+1) + s(i,j,k)), &
              s_b => half*(s(i,j,k) + s(i,j,k-1)) )
              associate( &
                D_f => this%diffusion_coefficient( s_f ), &
                D_b => this%diffusion_coefficient( s_b) )

                div_flux_z(i,j,k) = &
                  D_f*(s(i,j,k+1) - s(i,j,k))/(dz_f*dz_m) - &
                  D_b*(s(i,j,k) - s(i,j,k-1))/(dz_b*dz_m)
              end associate
            end associate
          end do
        end associate

        ! TODO
        ! 1. Each block sets scalar_flux packages on halo blocks

        call div_flux_internal_points%set_scalar( div_flux_x + div_flux_y + div_flux_z )

      end associate
    end associate
  end procedure set_up_div_scalar_flux

  module procedure div_scalar_flux
    integer alloc_stat
    integer, parameter :: success=0
    character(len=max_errmsg_len) :: alloc_error
    real(r8k), allocatable, dimension(:,:,:) :: div_flux_x, div_flux_y, div_flux_z

    call assert( same_type_as(this, vertices), "cartesian_grid%div_scalar_flux: same_type_as(this, vertices)" )

    associate( positions => vertices%vectors(), s=>this%get_scalar() )
      associate( npoints => shape(positions(:,:,:,1)) )

        allocate(div_flux_x, div_flux_y, div_flux_z, mold=positions(:,:,:,1), stat=alloc_stat, errmsg=alloc_error )
        call assert( alloc_stat==success, "cartesian_grid%div_scalar_flux: allocate(div_flux_{x,y,z})", alloc_error )

        div_flux_x = 0._r8k
        div_flux_y = 0._r8k
        div_flux_z = 0._r8k

        ! 2. Each block gets scalar_flux packages from its halo
        ! 3. Each block uses its halo data to compute surface fluxes

        hardwire_known_boundary_values: &
        block
          real(r8k), parameter ::  x_center = 3*(0.25E-01) + 1.E-01/2., x_max = 2*x_center
          real(r8k), parameter ::  y_center = 0.5E-01 + 2*(0.25E-01) + 3.E-01/2., y_max = 2*y_center
          real(r8k), parameter ::  z_center = 20.E-01/2., z_max=2*z_center

          ! r_sq = ((x-x_center)/(x_max-x_center))**2 + ((y-y_center)/(y_max-y_center))**2 + ((z-z_center)/(z_max-z_center))**2)
          ! scalar = 1. - r_sq

          real(r8k), parameter ::  div_x_expected = -2./(x_max-x_center)
          real(r8k), parameter ::  div_y_expected = -2./(y_max-y_center)
          real(r8k), parameter ::  div_z_expected = -2./(z_max-z_center)
          real(r8k), parameter ::  tolerance=1.E-06

          div_flux_x(1,:,:) = div_x_expected
          div_flux_x(npoints(1),:,:) = div_x_expected

          div_flux_y(:,1,:) = div_y_expected
          div_flux_y(:,npoints(2),:) = div_y_expected

          div_flux_z(:,:,1) = div_z_expected
          div_flux_z(:,:,npoints(3)) = div_z_expected

          associate( div_flux_total => div_flux%get_scalar() + div_flux_x + div_flux_y + div_flux_z )
            associate( div_expected => div_x_expected + div_y_expected + div_z_expected )
              call assert( all( abs((div_flux_total - div_expected)/div_expected) < tolerance ), "div_scalar_flux: div_expected")
            end associate
          end associate
        end block hardwire_known_boundary_values
      end associate
    end associate

    call div_flux%increment_scalar( div_flux_x + div_flux_y + div_flux_z )

  end procedure

  module procedure block_indicial_coordinates

    call assert(n>0 .and. n<=product(this%get_global_block_shape()), "cartesian_grid%block_indicial_coordinates: identifier bounds")

    associate( extents=>this%get_global_block_shape() )
      associate( nx=>extents(1), ny=>extents(2), nz=>extents(3) )
        ijk = [ mod(n-1,nx)+1, mod( (n-1)/nx, ny ) + 1, (n-1)/(nx*ny) + 1]
      end associate
    end associate

    call assert(all(ijk>[0,0,0]) .and. all(ijk<=this%get_global_block_shape()), "block_indicial_coordinates: coordinates bounds")

  end procedure

  module procedure block_identifier

    character(len=256) :: diagnostic_string

    ! Requires
    if (assertions) then
      associate(assertion => all(ijk>[0,0,0]).and.all(ijk<=this%get_global_block_shape()))
        if (assertions) then
          write(diagnostic_string,*) "all(",ijk,">[0,0,0]), all(",ijk,"<=",this%get_global_block_shape(),")"
          call assert( assertion , "block_identifier: indicial coordinates in bounds", diagnostic_data=diagnostic_string)
        end if
      end associate
    end if

    associate( extents=>this%get_global_block_shape() )
      associate( i=>ijk(1), j=>ijk(2), k=>ijk(3),nx=>extents(1), ny=>extents(2) )
        n = (k-1)*(ny*nx) + (j-1)*nx + i
      end associate

      ! Assures
      if (assertions) then
        write(diagnostic_string,*) "ijk->n:",ijk,"->",n," extents=",extents
        call assert( n>0 .and. n<=product(extents), "structured_grid%block_identifier: identifier in bounds", diagnostic_string )
      end if
    end associate

  end procedure

  module procedure block_identifier_in_bounds
    associate( extents=>this%get_global_block_shape() )
      in_bounds = id>0 .and. id<=product(extents)
    end associate
  end procedure

  module procedure block_coordinates_in_bounds
    in_bounds = all(ijk>[0,0,0]) .and. all(ijk<=this%get_global_block_shape())
  end procedure

end submodule
