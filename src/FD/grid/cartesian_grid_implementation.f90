!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
submodule(cartesian_grid_interface) cartesian_grid_implementation
  !! author: Damian Rouson and Karla Morris
  use kind_parameters, only : i4k, r8k
  use assertions_interface,only : assert, max_errmsg_len, assertions
  use plate_3D_interface, only : plate_3D
  use surfaces_interface, only : surfaces, backward, forward
  use package_interface, only : package
  implicit none

  type(surfaces) block_faces

contains

  module procedure build_surfaces
    integer, parameter :: first=1, last=2, success=0, x_dir=1, y_dir=2, z_dir=3
    type(package), allocatable, dimension(:,:,:) :: bare
    character(len=max_errmsg_len) error_message
    integer alloc_stat, b

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

    loop_over_blocks: &
    do b=my_blocks(first), my_blocks(last)
      associate(ijk => this%block_indicial_coordinates(b))
        associate( &
          ijk_i_backward => ijk + [-1, 0, 0], &
          ijk_i_forward  => ijk + [ 1, 0, 0], &
          ijk_j_backward => ijk + [ 0,-1, 0], &
          ijk_j_forward  => ijk + [ 0, 1, 0], &
          ijk_k_backward => ijk + [ 0, 0,-1], &
          ijk_k_forward  => ijk + [ 0, 0, 1])

          if (this%block_in_bounds(ijk_i_backward)) &
             call bare(b, x_dir, backward)%set_sender_block_id( this%block_identifier(ijk_i_backward) )

          if (this%block_in_bounds(ijk_i_forward)) &
             call bare(b, x_dir, forward)%set_sender_block_id( this%block_identifier(ijk_i_forward) )

          if (this%block_in_bounds(ijk_j_backward)) &
             call bare(b, y_dir, backward)%set_sender_block_id( this%block_identifier(ijk_j_backward) )

          if (this%block_in_bounds(ijk_j_forward)) &
             call bare(b, y_dir, forward)%set_sender_block_id( this%block_identifier(ijk_j_forward) )

          if (this%block_in_bounds(ijk_k_backward)) &
             call bare(b, z_dir, backward)%set_sender_block_id( this%block_identifier(ijk_k_backward) )

          if (this%block_in_bounds(ijk_k_forward)) &
             call bare(b, z_dir, forward)%set_sender_block_id( this%block_identifier(ijk_k_forward) )

        end associate
      end associate
    end do loop_over_blocks

    call block_faces%set_halo_data(bare)
  end procedure

  module procedure assign_structured_grid
    call assert( same_type_as(this, rhs), "cartesian_grid%assign_structured_grid: consistent types" )
    call this%clone(rhs)
  end procedure

  module procedure div_scalar_flux
    !! Evaluate 2nd partial derivatives in each coordinate direction according to equations (2.2) & (2.4) in
    !! Sundqvist & Veronis (1969) "A simple finite-difference grid with non-constant intervals", Tellus 22:1

    integer(i4k) i, j, k, alloc_stat
    integer(i4k), parameter :: success=0
    real(r8k), parameter :: half=0.5_r8k
    character(len=max_errmsg_len) :: alloc_error
    real(r8k), allocatable, dimension(:,:,:) :: div_flux_x, div_flux_y, div_flux_z

    call assert( same_type_as(this, vertices), "div_scalar_flux: consistent types" )

    allocate(div_flux, mold=this, stat=alloc_stat, errmsg=alloc_error )
    call assert( alloc_stat==success, "div_scalar_flux (cartesian): result allocation fails with message '"//alloc_error//"'" )

    associate( positions => vertices%vectors(), s=>this%get_scalar() )
      associate( npoints => shape(positions(:,:,:,1)) )

        allocate(div_flux_x, div_flux_y, div_flux_z, mold=positions(:,:,:,1), stat=alloc_stat, errmsg=alloc_error )
        call assert( alloc_stat==success, "div_scalar_flux (cartesian): allocate(div_flux_{x,y,z}) (error: "//alloc_error//")" )


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
              s_f => half*(y(i,j+1,k) + y(i,j,k)), &
              s_b =>  half*(y(i,j,k) + y(i,j-1,k)) )

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
             s_f => half*(z(i,j,k+1)+z(i,j,k)), &
             s_b => half*(z(i,j,k)+z(i,j,k-1)) )
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

          call assert( all( abs((div_flux_x - div_x_expected)/div_x_expected) < tolerance ), "div_scalar_flux: div_x_expected")
          call assert( all( abs((div_flux_y - div_y_expected)/div_y_expected) < tolerance ), "div_scalar_flux: div_y_expected")
          call assert( all( abs((div_flux_z - div_z_expected)/div_z_expected) < tolerance ), "div_scalar_flux: div_z_expected")
        end block hardwire_known_boundary_values

        call div_flux%set_scalar( div_flux_x + div_flux_y + div_flux_z )

      end associate
    end associate
  end procedure div_scalar_flux

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
