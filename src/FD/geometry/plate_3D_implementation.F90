!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
submodule(plate_3D_interface) plate_3D_implementation
  !! author: Damian Rouson
  !! date: 8/16/2019
  use assertions_interface, only : assert
  implicit none

  character(len=*), parameter :: base_object = "MORFEUS_FD.GEOMETRY"
  character(len=*), parameter :: csv_format = '(*(G0,:,","))'
  integer, parameter :: success=0

contains

  module procedure set_grid_specification
    use string_functions_interface, only : file_extension
    use units_interface, only : units_system_names
    character(len=*), parameter :: expected_geometry_type="3D_plate"
    character(len=:), allocatable :: geometry_type, units_system
    logical found

    call assert( file_extension(grid_description_file)=="json", 'file_extension(grid_description_file)=="json"' )

    call this%grid_specification%load_file( grid_description_file )
    call assert( .not. this%grid_specification%failed(), "set_grid_specification: .not. this%grid_specification%failed()" )

    call this%grid_specification%get( base_object//".type",  geometry_type, found )
    call assert( found , base_object//".type found" )
    call assert(geometry_type==expected_geometry_type, "geometry_type==expected_geometry_type" )

    call this%grid_specification%get( base_object//".units_system",  units_system, found )
    call assert( found , base_object//".units_system found" )
    call assert(any(units_system==units_system_names), "any(units_system==units_system_names)" )

  end procedure

  module procedure get_block_metadata_shape
    shape_ = shape(this%metadata)
  end procedure

  module procedure get_block_domain
    associate( ix=>indicial_coordinates(1), iy=>indicial_coordinates(2), iz=>indicial_coordinates(3))
      call assert( all( [ [ix,iy,iz] >= lbound(this%metadata), [ix,iy,iz] <= ubound(this%metadata) ] ), &
      "get_block_domain: indicial_coordinates in bounds" )
      this_domain = this%metadata(ix,iy,iz)%get_subdomain()
    end associate
  end procedure

  module procedure get_block_metadata
    associate( ix=>indicial_coordinates(1), iy=>indicial_coordinates(2), iz=>indicial_coordinates(3))
      call assert( all( [ [ix,iy,iz] >= lbound(this%metadata), [ix,iy,iz] <= ubound(this%metadata) ] ), &
      "get_block_metadata: indicial_coordinates in bounds" )
      this_metadata = this%metadata(ix, iy, iz)
    end associate
  end procedure

  module procedure set_block_metadata
    integer ix, iy, iz, tag, alloc_stat
    integer, parameter :: max_stringlen=32 !! 3*10 digits + 2 commas
    character(len=max_stringlen) element
    character(len=:), allocatable :: key, label
    logical found
    integer, dimension(:), allocatable :: block_metadata_shape
    real, dimension(:), allocatable :: subdomain
    integer, parameter :: num_end_points=2

    call this%grid_specification%get( base_object // ".global_shape", block_metadata_shape, found)
    call assert( found, "set_block_metadata: found" )
    call assert( size(block_metadata_shape)==space_dimension, "size(block_metadata_shape)==space_dimension")

    associate(nx=>block_metadata_shape(1), ny=>block_metadata_shape(2), nz=>block_metadata_shape(3))
      allocate(this%metadata(nx, ny, nz), stat=alloc_stat )
      call assert( alloc_stat==success, "set_block_metadata: allocate(this%metadata(nx,ny,nz),...)" )

      do ix=1, nx
        do iy=1, ny
          do iz=1, nz
            write(element, csv_format) ix, iy, iz

            key = trim( base_object // ".structured_grid_blocks." // trim(element) // ".subdomain" )
            call this%grid_specification%get( key, subdomain, found)
            call assert( found, "set_block_metadata: found key " // key )
            call this%metadata(ix,iy,iz)%set_subdomain( reshape(subdomain, [space_dimension, num_end_points]) )

            key = trim( base_object // ".structured_grid_blocks." // trim(element) // ".tag" )
            call this%grid_specification%get( key, tag, found)
            call assert( found, "set_block_metadata: found key " // key )
            call this%metadata(ix,iy,iz)%set_tag(tag)

            key = trim( base_object // ".structured_grid_blocks." // trim(element) // ".label" )
            call this%grid_specification%get( key, label, found)
            call assert( found, "set_block_metadata: found key " // key )
            call this%metadata(ix,iy,iz)%set_label(label)
          end do
        end do
      end do
    end associate
  end procedure

end submodule plate_3D_implementation
