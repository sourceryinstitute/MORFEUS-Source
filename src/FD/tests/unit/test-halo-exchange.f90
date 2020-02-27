!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
program main
  !! author: Damian Rouson and Karla Morris
  !! date: 12/26/2019
  !!
  !! verify the setting of halo data on structured_grid block surfaces
  use assertions_interface, only : assert
  use problem_discretization_interface, only :  problem_discretization
  use surfaces_interface, only : surfaces, forward, backward, face_name, coordinate_name, x_dir, z_dir
  use package_interface, only : package
  use plate_3D_interface, only : plate_3D
  use ellipsoidal_field_interface, only : ellipsoidal_field
  implicit none

  type(plate_3D) plate_geometry
  type(problem_discretization) global_grid
  integer, parameter :: max_digits=9
  character(len=*), parameter :: input = "3Dplate-low-resolution-halo.json"
  character(len=*), parameter:: base_name = "3Dplate-low-resolution-halo"
  character(len=max_digits) image_number
  type(ellipsoidal_field) ellipsoidal_function
  integer b, coord_dir, face_dir
  type(package), allocatable, dimension(:,:,:) ::  surface_packages
  type(surfaces) block_surfaces
  character(len=2) id_string

  call plate_geometry%build(input)

  call global_grid%initialize_from_geometry(plate_geometry)
  call global_grid%set_scalars( [ellipsoidal_function] ) ! includes build_surfaces call
  call global_grid%set_scalar_flux_divergence( exact_result=[ellipsoidal_function] )
  write(image_number,'(i4)') this_image()
  call global_grid%write_output(base_name //"-image-"// trim(adjustl(image_number)) // ".vtu")
  surface_packages = global_grid%get_surface_packages()

  associate( &
    my_blocks => global_grid%my_blocks(), &
    num_blocks => size(surface_packages,1) )

    call assert( &
       assertion = size(surface_packages,2) == z_dir .and. size(surface_packages,3)==forward, &
       description = "test-halo-exchange: surface_packages shape" )

    loop_over_blocks: &
    do b = 1, num_blocks
      loop_over_coordinate_directions: &
      do coord_dir = x_dir, z_dir
        loop_over_face_directions: &
        do face_dir = backward, forward
          associate(block_id => my_blocks(1) + b - 1)
            write(id_string,'(i2)') block_id
            if (block_surfaces%is_external_boundary(block_id, coord_dir, face_dir)) then
              call assert( &
                assertion = surface_packages(b, coord_dir, face_dir)%neighbor_block_id_null(), &
                description = "test-halo-exchange: surface_packages(b, coord_dir, face_dir)%neighbor_block_id_null()", &
                diagnostic_data = face_name(face_dir) // "-" // coordinate_name(coord_dir) // " face on block " // id_string )
            else
              call assert( &
                assertion = .not. surface_packages(b, coord_dir, face_dir)%neighbor_block_id_null(), &
                description = "test-halo-exchange: .not. surface_packages(b, coord_dir, face_dir)%neighbor_block_id_null()", &
                diagnostic_data = face_name(face_dir) // "-" // coordinate_name(coord_dir) // " face on block " // id_string )
            end if
          end associate
        end do loop_over_face_directions
      end do loop_over_coordinate_directions
    end do loop_over_blocks
  end associate

  sync all
  if (this_image()==1) print *,"Test passed"
end program
