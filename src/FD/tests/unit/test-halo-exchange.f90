!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program main
  !! author: Damian Rouson and Karla Morris
  !! date: 12/26/2019
  !!
  !! verify the setting of halo data on structured_grid block surfaces
  use assertions_interface, only : assert
  use problem_discretization_interface, only :  problem_discretization
  use surfaces_interface, only : forward, backward
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
  class(package), allocatable, dimension(:,:,:) ::  surface_packages

  call plate_geometry%build(input)

  call global_grid%initialize_from_geometry(plate_geometry)
  call global_grid%set_scalars( [ellipsoidal_function] )
  call global_grid%set_scalar_flux_divergence( exact_result=[ellipsoidal_function] )
  write(image_number,'(i4)') this_image()
  call global_grid%write_output(base_name //"-image-"// trim(adjustl(image_number)) // ".vtu")
  call global_grid%get_surface_packages(surface_packages)

 associate( block_surfaces => global_grid%get_block_surfaces(), my_blocks => global_grid%my_blocks() )
    loop_over_blocks: &
    do b=my_blocks(1), my_blocks(2)
      loop_over_coordinate_directions: &
      do coord_dir = 1, 3
        loop_over_face_directions: &
        do face_dir = backward, forward
          if (block_surfaces%is_external_boundary(b, coord_dir, face_dir)) then
            call assert( surface_packages(b, coord_dir, face_dir)%sender_block_id_null(), &
              "test-halo-exchange: surface_packages(b, coord_dir, face_dir)%sender_block_id_null()")
          else
            call assert( .not. surface_packages(b, coord_dir, face_dir)%sender_block_id_null(), &
              "test-halo-exchange: .not. surface_packages(b, coord_dir, face_dir)%sender_block_id_null()")
          end if
        end do loop_over_face_directions
      end do loop_over_coordinate_directions
    end do loop_over_blocks
  end associate

  sync all
  print *,"Test passed"
end program
