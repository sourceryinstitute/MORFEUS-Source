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
  use Kinds, only : r8k
  implicit none

  character(len=*), parameter :: base_object = "MORFEUS_FD.layers"
  character(len=*), parameter :: csv_format = '(*(G0,:,","))'
  integer, parameter :: success=0, max_stringlen=32

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

    !! Encapsulate json key/value pair hierarchy to be read from a 3Dplate*.json file

    type material_name
      !! Facilitate ragged-edged array of strings
      character(len=:), allocatable :: string
    end type

    type thickness
      real, allocatable, dimension(:) :: x, y, z
      character(len=:), allocatable :: dimensions
    end type

    type num_grid_blocks
      integer, allocatable, dimension(:) :: x, y, z
    end type

    type core
      type(material_name) material_name_
      type(thickness) thickness_
      type(num_grid_blocks) num_grid_blocks_
    end type

    type wrappers
      type(material_name), allocatable, dimension(:) ::  material_name_
      type(thickness) thickness_
      type(num_grid_blocks) num_grid_blocks_
    end type

    type layers
      character(len=:), allocatable :: type_, units_system_, max_spacing_
      type(core) core_
      type(wrappers) wrappers_
    end type

    type(layers) layers_

    integer ix, iy, iz, tag, alloc_stat
    character(len=max_stringlen) element
    character(len=:), allocatable :: key, label
    logical found
    integer, dimension(:), allocatable :: block_metadata_shape
    real(r8k), dimension(:), allocatable :: subdomain
    integer, parameter :: num_end_points=2
    real(r8k) max_spacing

    call this%grid_specification%get( base_object // ".core.num_grid_blocks.x", layers_%core_%num_grid_blocks_%x, found)
    call assert( found, "set_block_metadata: core_%num_grid_blocks_%x found" )

    call this%grid_specification%get( base_object // ".core.num_grid_blocks.y", layers_%core_%num_grid_blocks_%y, found)
    call assert( found, "set_block_metadata: core_%num_grid_blocks_%y found" )

    call this%grid_specification%get( base_object // ".core.num_grid_blocks.z", layers_%core_%num_grid_blocks_%z, found)
    call assert( found, "set_block_metadata: core_%num_grid_blocks_%z found" )

    call this%grid_specification%get( base_object // ".wrappers.num_grid_blocks.x", layers_%wrappers_%num_grid_blocks_%x, found)
    call assert( found, "set_block_metadata: wrappers_%num_grid_blocks_%x found" )

    call this%grid_specification%get( base_object // ".wrappers.num_grid_blocks.y", layers_%wrappers_%num_grid_blocks_%y, found)
    call assert( found, "set_block_metadata: wrappers_%num_grid_blocks_%y found" )

    call this%grid_specification%get( base_object // ".wrappers.num_grid_blocks.z", layers_%wrappers_%num_grid_blocks_%z, found)
    call assert( found, "set_block_metadata: wrappers_%num_grid_blocks_%z found" )

    call this%grid_specification%get( base_object//".max_spacing",  max_spacing, found )
    call assert( found , base_object//".max_spacing found" )
    call assert(max_spacing>0., "max_spacing>0." )

    associate( &
      num_core_blocks => layers_%core_%num_grid_blocks_, &
      num_wrapper_blocks => layers_%wrappers_%num_grid_blocks_ )
      associate(  &
        nx => sum(num_core_blocks%x) + sum(num_wrapper_blocks%x) , &
        ny => sum(num_core_blocks%y) + sum(num_wrapper_blocks%y) , &
        nz => sum(num_core_blocks%z) + sum(num_wrapper_blocks%z) )

        block_metadata_shape = [nx, ny, nz]

        allocate(this%metadata(nx, ny, nz), stat=alloc_stat )
        call assert( alloc_stat==success, "set_block_metadata: allocate(this%metadata(nx,ny,nz),...)" )

        do concurrent( ix=1:nx, iy=1:ny, iz=1:nz )

          write(element, csv_format) ix, iy, iz
          call this%metadata(ix,iy,iz)%set_label(element)
          call this%metadata(ix,iy,iz)%set_max_spacing(max_spacing)
         !call this%metadata(ix,iy,iz)%set_tag()
         !call this%metadata(ix,iy,iz)%set_subdomain( reshape(subdomain, [space_dimension, num_end_points]) )
        end do
      end associate
    end associate
  end procedure

end submodule plate_3D_implementation
