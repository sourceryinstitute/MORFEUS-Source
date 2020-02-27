!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
submodule(spherical_grid_interface) spherical_grid_implementation
  !! author: Damian Rouson and Karla Morris
  use kind_parameters, only : i4k, r8k
  use assertions_interface,only : assert, max_errmsg_len, assertions
  use plate_3D_interface, only : plate_3D
  use surfaces_interface, only : backward, forward, x_dir, y_dir, z_dir
  use package_interface, only : package, null_neighbor_id
  implicit none

  integer, parameter :: success=0
    !! allocation stat value indicating success
  integer, parameter :: max_coord_dirs = size([x_dir,y_dir,z_dir])
    !! maximum number of coordinate directions
  integer, parameter :: max_vec_components = max_coord_dirs
    !! maximum number of vector components
  integer, parameter :: num_faces = size([backward, forward])
    !! number of faces in each coordinate direction of a hexahedral volume with faces orthogonal to coordinate directions
  integer, parameter, dimension(max_coord_dirs, num_faces , max_vec_components) :: displacement = &
    reshape([ [-1,0,0], [1,0,0], [0,-1,0], [0,1,0], [0,0,-1], [0,0,1] ], [max_coord_dirs, num_faces, max_vec_components])
    !! displacement vectors in indicial coordinates for structured_grid blocks

contains

  module procedure build_surfaces
    integer, parameter :: first=1, last=2
    type(package), allocatable, dimension(:,:,:) :: bare
    character(len=max_errmsg_len) error_message
    integer alloc_stat, b, coord_dir, face_dir

    select type(problem_geometry)
      type is(plate_3D)
        !! this verifies correct geometry for the assumptions in building the surfaces below
      class default
        error stop "spherical_grid%build_surfaces: unsupported problem_geometry type"
    end select

    define_bare_package: & !! initialize packages with only the time step & neighbor block_id (no halo-exchange data yet)
    associate( me => this_image() )
      associate( my_blocks => [block_partitions(me), block_partitions(me+1)-1] )

        allocate( bare(my_blocks(first):my_blocks(last), vertices(my_blocks(1))%space_dimension(), backward:forward), &
          stat=alloc_stat, errmsg=error_message)
        call assert(alloc_stat==success, "spherical_grid%build_surfaces: allocate(bare)", error_message)

        call bare%set_neighbor_block_id(null_neighbor_id)
        call bare%set_step(0)

        loop_over_blocks: &
        do b=my_blocks(first), my_blocks(last)

          loop_over_coordinate_directions: &
          do coord_dir = x_dir, z_dir
            loop_over_face_directions: &
            do face_dir = backward, forward
              associate( ijk_displaced => this%block_indicial_coordinates(b) + displacement(coord_dir, face_dir, :) )
                if (this%block_in_bounds(ijk_displaced)) then
                  call bare(b, coord_dir, face_dir)%set_neighbor_block_id( this%block_identifier(ijk_displaced) )
                end if
              end associate
            end do loop_over_face_directions
          end do loop_over_coordinate_directions

          set_surface_vertices: &
          associate(positions => vertices(b)%vectors())
            associate( nx => size(positions, x_dir), ny => size(positions, y_dir), nz => size(positions, z_dir))
              call bare(b, x_dir, backward)%set_surface_positions(positions(      2:2,:,:,:))
              call bare(b, x_dir, forward )%set_surface_positions(positions(nx-1:nx-1,:,:,:))

              call bare(b, y_dir, backward)%set_surface_positions(positions(:,      2:2,:,:))
              call bare(b, y_dir, forward )%set_surface_positions(positions(:,ny-1:ny-1,:,:))

              call bare(b, z_dir, backward)%set_surface_positions(positions(:,:,      2:2,:))
              call bare(b, z_dir, forward )%set_surface_positions(positions(:,:,nz-1:nz-1,:))
            end associate
          end associate set_surface_vertices

        end do loop_over_blocks
      end associate
    end associate define_bare_package

    call block_faces%set_halo_outbox(bare, block_partitions)
    call block_faces%set_num_scalars(num_scalars)

  end procedure

  module procedure set_up_div_scalar_flux
    !! Evaluate 2nd partial derivatives in each coordinate direction according to equations (2.2) & (2.4) in
    !! Sundqvist & Veronis (1969) "A simple finite-difference grid with non-constant intervals", Tellus 22:1

    integer(i4k) i, j, k, alloc_stat
    real(r8k), parameter :: half=0.5_r8k
    character(len=max_errmsg_len) :: alloc_error
    real(r8k), allocatable, dimension(:,:,:,:) :: div_flux
    real(r8k), allocatable, dimension(:,:,:) :: surface_fluxes

    call assert( same_type_as(this, vertices), "spherical_grid%set_up_div_scalar_flux: same_type_as(this, vertices)" )

    associate( positions => vertices%vectors(), s=>this%get_scalar() )
      associate( npoints => shape(positions(:,:,:,1)) )

        allocate(div_flux, mold=positions(:,:,:,:), stat=alloc_stat, errmsg=alloc_error )
        call assert( alloc_stat==success, "spherical_grid%set_up_div_scalar_flux: allocate(div_flux)", alloc_error )

        div_flux = 0._r8k

        x_direction_fluxes: &
        associate( x=>positions(:,:,:,1) )

          all_but_x_dir_boundaries: &
          do concurrent(k=1:npoints(3), j=1:npoints(2), i=2:npoints(1)-1)
            !! compute scalar flux divergence at internal points and boundary points excluding x-direction boundaries

            associate( &
              dx_m => half*(x(i+1,j,k) - x(i-1,j,k)), & !! (dx_b + dx_f)/2
              dx_f =>       x(i+1,j,k) - x(i,j,k), &
              dx_b =>       x(i,j,k)   - x(i-1,j,k), &
              s_f => half*( s(i+1,j,k) + s(i,j,k)  ), &
              s_b => half*( s(i,j,k)   + s(i-1,j,k) ))
              associate( &
                D_f => this%diffusion_coefficient( s_f ), &
                D_b => this%diffusion_coefficient( s_b) )

                div_flux(i,j,k,x_dir) = ( &
                  D_f*(s(i+1,j,k) - s(i,j,k)  )/dx_f - & !! forward flux in x direction
                  D_b*(s(i,j,k)   - s(i-1,j,k))/dx_b &   !! backward flux in x direction
                  ) / dx_m
              end associate
            end associate
          end do all_but_x_dir_boundaries

          allocate( surface_fluxes(npoints(y_dir), npoints(z_dir), backward:forward), stat = alloc_stat, errmsg = alloc_error)
          call assert(alloc_stat==success, "spherical_grid%set_up_div_scalar_flux: allocate(surface_fluxes) x_dir", alloc_error)

          x_normal_surface_fluxes: &
          do k=1, npoints(z_dir)
            do j=1, npoints(y_dir)

              i=1
              forward_difference_at_backward_face: &
              associate( &
                dx_f =>      x(i+1,j,k) - x(i,j,k), &
                s_f => half*(s(i+1,j,k) + s(i,j,k)) )
                associate( D_f => this%diffusion_coefficient( s_f ) )
                  surface_fluxes(j,k,backward) = D_f*(s(i+1,j,k) - s(i,j,k))/dx_f
                end associate
              end associate forward_difference_at_backward_face

              i=npoints(x_dir)
              backward_difference_at_forward_face: &
               associate( &
                 dx_b =>       x(i,j,k)   - x(i-1,j,k), &
                 s_b => half*( s(i,j,k)   + s(i-1,j,k) ))
                 associate( D_b => this%diffusion_coefficient( s_b) )
                   surface_fluxes(j,k,forward) = D_b*(s(i,j,k) - s(i-1,j,k))/dx_b
                 end associate
               end associate backward_difference_at_forward_face
            end do
          end do x_normal_surface_fluxes

          call block_surfaces%set_normal_scalar_fluxes( &
            this%get_block_identifier(), x_dir, backward, surface_fluxes(:,:,backward), this%get_scalar_identifier())
          call block_surfaces%set_normal_scalar_fluxes( &
            this%get_block_identifier(), x_dir,  forward, surface_fluxes(:,:, forward), this%get_scalar_identifier())

        end associate x_direction_fluxes

        y_direction_fluxes: &
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

                div_flux(i,j,k,y_dir) = ( &
                  D_f*(s(i,j+1,k) - s(i,j,k))/dy_f - & !! forward flux in y direction
                  D_b*(s(i,j,k) - s(i,j-1,k))/dy_b &   !! backward flux in y direction
                  ) / dy_m
              end associate
            end associate
          end do

          if ( any(shape(surface_fluxes) /= [npoints(x_dir), npoints(z_dir), forward]) ) then
            deallocate(surface_fluxes)
            allocate( surface_fluxes(npoints(x_dir), npoints(z_dir), backward:forward), stat = alloc_stat, errmsg = alloc_error)
            call assert( alloc_stat==success, "spherical_grid%set_up_div_scalar_flux: allocate(surface_fluxes) y_dir", alloc_error )
          end if

          y_normal_surface_fluxes: &
          do k=1, npoints(z_dir)
            do i=1, npoints(x_dir)

              j=1
              forward_difference_at_backward_face: &
              associate( &
                dy_f =>      y(i,j+1,k) - y(i,j,k), &
                s_f => half*(s(i,j+1,k) + s(i,j,k)) )
                associate( D_f => this%diffusion_coefficient( s_f ) )
                  surface_fluxes(i,k,backward) = D_f*(s(i,j+1,k) - s(i,j,k))/dy_f
                end associate
              end associate forward_difference_at_backward_face

              j=npoints(y_dir)
              backward_difference_at_forward_face: &
              associate( &
                dy_b =>       y(i,j,k)   - y(i,j-1,k), &
                s_b => half*( s(i,j,k)   + s(i,j-1,k) ))
                associate( D_b => this%diffusion_coefficient( s_b) )
                  surface_fluxes(i,k,forward) = D_b*(s(i,j,k) - s(i,j-1,k))/dy_b
                end associate
              end associate backward_difference_at_forward_face
            end do
          end do y_normal_surface_fluxes

          call block_surfaces%set_normal_scalar_fluxes( &
            this%get_block_identifier(), y_dir, backward, surface_fluxes(:,:,backward), this%get_scalar_identifier())
          call block_surfaces%set_normal_scalar_fluxes( &
            this%get_block_identifier(), y_dir,  forward, surface_fluxes(:,:, forward), this%get_scalar_identifier())

        end associate y_direction_fluxes

        z_direction_fluxes: &
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

                div_flux(i,j,k,z_dir) = ( &
                  D_f*(s(i,j,k+1) - s(i,j,k))/dz_f - &  !! forward flux in z direction
                  D_b*(s(i,j,k) - s(i,j,k-1))/dz_b &    !! backward flux in z direction
                  ) / dz_m
              end associate
            end associate
          end do

          if ( any(shape(surface_fluxes) /= [npoints(x_dir), npoints(y_dir), forward]) ) then
            deallocate(surface_fluxes)
            allocate( surface_fluxes(npoints(x_dir), npoints(y_dir), backward:forward), stat = alloc_stat, errmsg = alloc_error)
            call assert( alloc_stat==success, "spherical_grid%set_up_div_scalar_flux: allocate(surface_fluxes) z_dir", alloc_error )
          end if

          z_normal_surface_fluxes: &
          do j=1, npoints(y_dir)
            do i=1, npoints(x_dir)

              k=1
              forward_difference_at_backward_face: &
              associate( &
                dz_f =>      z(i,j,k+1) - z(i,j,k), &
                s_f => half*(s(i,j,k+1) + s(i,j,k)) )
                associate( D_f => this%diffusion_coefficient( s_f ) )
                  surface_fluxes(i,j,backward) = D_f*(s(i,j,k+1) - s(i,j,k))/dz_f
                end associate
              end associate forward_difference_at_backward_face

              k=npoints(z_dir)
              backward_difference_at_forward_face: &
              associate( &
                dz_b =>       z(i,j,k)   - z(i,j,k-1), &
                s_b => half*( s(i,j,k)   + s(i,j,k-1) ))
                associate( D_b => this%diffusion_coefficient( s_b) )
                  surface_fluxes(i,j,forward) = D_b*(s(i,j,k) - s(i,j,k-1))/dz_b
                end associate
              end associate backward_difference_at_forward_face
            end do
          end do z_normal_surface_fluxes

          call block_surfaces%set_normal_scalar_fluxes( &
            this%get_block_identifier(), z_dir, backward, surface_fluxes(:,:,backward), this%get_scalar_identifier())
          call block_surfaces%set_normal_scalar_fluxes( &
            this%get_block_identifier(), z_dir,  forward, surface_fluxes(:,:, forward), this%get_scalar_identifier())

        end associate z_direction_fluxes

        call div_flux_internal_points%set_scalar( div_flux(:,:,:,x_dir) + div_flux(:,:,:,y_dir) + div_flux(:,:,:,z_dir) )

      end associate
    end associate
  end procedure set_up_div_scalar_flux

  module procedure div_scalar_flux
    real(r8k), parameter :: half=0.5_r8k
    integer alloc_stat, i, j, k
    character(len=max_errmsg_len) :: alloc_error
    real(r8k), allocatable, dimension(:,:,:,:) :: div_flux_increment

    call assert( same_type_as(this, vertices), "spherical_grid%div_scalar_flux: same_type_as(this, vertices)" )

    associate( &
      positions => vertices%vectors(), &
      s => this%get_scalar(), &
      b => this%get_block_identifier() )

      allocate(div_flux_increment, mold=positions, stat=alloc_stat, errmsg=alloc_error )
      call assert( alloc_stat==success, "spherical_grid%div_scalar_flux: allocate(div_flux_increment)", alloc_error )

      div_flux_increment = 0._r8k

      associate( npoints => shape(positions(:,:,:,1)) )

        x_direction_fluxes: &
        associate( x=>positions(:,:,:,1) )

          i=1
          !backward_face: &
          !associate(neighbor_image => block_surfaces%get_block_image( block_surfaces%get_neighbor_block_id(b, x_dir, backward) ) )
          !  associate( &
          !    dx_b => block_surfaces%get_surface_normal_spacing(neighbor_image, b, x_dir, forward), &
          !    surface_fluxes => &
          !    block_surfaces%get_normal_scalar_fluxes(neighbor_image, b, x_dir, forward, this%get_scalar_identifier()) )

          !    do concurrent(k=1:npoints(z_dir), j=1:npoints(y_dir))
          !      associate( dx_f => x(i+1,j,k) - x(i,j,k) )
          !        associate( &
          !          dx_m => half*(dx_f + dx_b), &
          !          s_f => half*( s(i+1,j,k) + s(i,j,k)  ) )
          !          associate( D_f => this%diffusion_coefficient( s_f ) )

          !             div_flux_increment(i,j,k,x_dir) = &
          !               ( D_f*(s(i+1,j,k) - s(i,j,k)  )/dx_f - surface_fluxes(j,k) ) / dx_m

          !          end associate
          !        end associate
          !      end associate
          !    end do
          !  end associate
          !end associate backward_face

        end associate x_direction_fluxes
      end associate
    end associate

    ! TODO
    ! Write y & z directions for the following & test with Intel compiler on Linux:
    ! 1. Each block gets block_surfaces packages from its halo
    ! 2. Each block uses its halo data to compute surface fluxes

!      call div_flux%increment_scalar( &
!        div_flux_increment(:,:,:,x_dir) + div_flux_increment(:,:,:,y_dir) + div_flux_increment(:,:,:,z_dir))

    hardwire_known_boundary_values: &
    block
      real(r8k), parameter ::  x_center = 3*(0.25E-01) + 1.E-01/2., x_max = 2*x_center
      real(r8k), parameter ::  y_center = 0.5E-01 + 2*(0.25E-01) + 3.E-01/2., y_max = 2*y_center
      real(r8k), parameter ::  z_center = 20.E-01/2., z_max=2*z_center
      real(r8k), allocatable, dimension(:,:,:) :: div_flux_x, div_flux_y, div_flux_z

      ! r_sq = ((x-x_center)/(x_max-x_center))**2 + ((y-y_center)/(y_max-y_center))**2 + ((z-z_center)/(z_max-z_center))**2)
      ! scalar = 1. - r_sq

      real(r8k), parameter ::  div_x_expected = -2./(x_max-x_center)**2
      real(r8k), parameter ::  div_y_expected = -2./(y_max-y_center)**2
      real(r8k), parameter ::  div_z_expected = -2./(z_max-z_center)**2
      real(r8k), parameter ::  tolerance=1.E-06

      associate(positions => vertices%vectors())
        associate( npoints => shape(positions(:,:,:,1)) )

          allocate(div_flux_x, div_flux_y, div_flux_z, mold=positions(:,:,:,1), stat=alloc_stat, errmsg=alloc_error )
          call assert( alloc_stat==success, "spherical_grid%div_scalar_flux: allocate(div_flux_x/y/z)", alloc_error )

          div_flux_x = 0._r8k
          div_flux_y = 0._r8k
          div_flux_z = 0._r8k

          div_flux_x(1,:,:) = div_x_expected
          div_flux_x(npoints(1),:,:) = div_x_expected

          div_flux_y(:,1,:) = div_y_expected
          div_flux_y(:,npoints(2),:) = div_y_expected

          div_flux_z(:,:,1) = div_z_expected
          div_flux_z(:,:,npoints(3)) = div_z_expected

        end associate

        associate( div_flux_total => div_flux%get_scalar() + div_flux_x + div_flux_y + div_flux_z )
          associate( div_expected => div_x_expected + div_y_expected + div_z_expected )
            call assert( all( abs((div_flux_total - div_expected)/div_expected) < tolerance ), "div_scalar_flux: div_expected")
          end associate
        end associate

        call div_flux%increment_scalar( div_flux_x + div_flux_y + div_flux_z )

      end associate

    end block hardwire_known_boundary_values

  end procedure

  module procedure block_indicial_coordinates

    call assert(n>0 .and. n<=product(this%get_global_block_shape()), "spherical_grid%block_indicial_coordinates: identifier bounds")

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
