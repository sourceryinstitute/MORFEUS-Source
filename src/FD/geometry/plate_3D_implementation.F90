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

  character(len=*), parameter :: base_object = "MORFEUS_FD.layers"
  character(len=*), parameter :: csv_format = '(*(G0,:,","))'
  integer, parameter :: success=0

  !! Encapsulate json key/value pair hierarchy to be read from a 3Dplate*.json file

  type thickness_t
    real, allocatable, dimension(:) :: x, y, z
    character(len=:), allocatable :: dimensions
  end type

  type num_grid_blocks_t
    integer, allocatable, dimension(:) :: x, y, z
  end type

  type material_t
    !character(len=:), allocatable, dimension(:) :: material_name ! gfortran 8.3 bug causes subsequent memory corruption
    character(len=max_name_length), allocatable, dimension(:) :: material_name
    type(thickness_t) thickness
    type(num_grid_blocks_t) num_grid_blocks
  end type

  type layers_t
    character(len=:), allocatable :: type_, units_system
    real max_spacing
    type(material_t) core, wrappers
  end type

  type(layers_t) layers
contains

  module procedure set_grid_specification
    use string_functions_interface, only : file_extension
    use units_interface, only : units_system_names
    character(len=*), parameter :: expected_geometry_type="3D_plate"
    logical found

    call assert( file_extension(grid_description_file)=="json", 'file_extension(grid_description_file)=="json"' )

    call this%grid_specification%load_file( grid_description_file )
    call assert( .not. this%grid_specification%failed(), "set_grid_specification: .not. this%grid_specification%failed()" )

    call this%grid_specification%get( base_object//".type",  layers%type_, found )
    call assert( found , base_object//".type found" )
    call assert(layers%type_==expected_geometry_type, "layers%type_==expected_geometry_type" )

    call this%grid_specification%get( base_object//".units_system",  layers%units_system, found )
    call assert( found , base_object//".units_system found" )
    call assert(any(layers%units_system==units_system_names), "any(layers%units_system==units_system_names)" )

    call this%grid_specification%get( base_object//".max_spacing",  layers%max_spacing, found )
    call assert( found , base_object//".max_spacing found" )

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
#ifndef HAVE_FINDLOC
    use emulated_intrinsics_interface, only : findloc
#endif

    integer i, j, ix, iy, iz, alloc_stat, supremum
    logical found
    integer, parameter :: num_end_points=2, symmetry=2
    character(len=max_name_length), allocatable, dimension(:) :: names

    call read_core_components
    call verify_core_components

    call read_wrappers_components
    call verify_wrappers_components

    call verify_layers

    call set_metadata

  contains

    subroutine read_core_components

      call this%grid_specification%get( base_object // ".core.num_grid_blocks.x", layers%core%num_grid_blocks%x, found)
      call assert( found, "set_block_metadata: core_%num_grid_blocks_%x found" )

      call this%grid_specification%get( base_object // ".core.num_grid_blocks.y", layers%core%num_grid_blocks%y, found)
      call assert( found, "set_block_metadata: core_%num_grid_blocks_%y found" )

      call this%grid_specification%get( base_object // ".core.num_grid_blocks.z", layers%core%num_grid_blocks%z, found)
      call assert( found, "set_block_metadata: core_%num_grid_blocks_%z found" )

      call this%grid_specification%get( base_object//".core.material_name", names, found )
      call assert( found , base_object//".core.material_name found" )

      supremum = maxval( [( len( trim(names(i)) ), i=1,size(names) )] )
     !allocate( character(len=supremum) :: layers%core%material_name(size(names)) ) !! precluded by gfortran 8.3 bug
      allocate(  layers%core%material_name(size(names)) )
      layers%core%material_name = names
      call assert(size(layers%core%material_name)==1, "size(layers%core%material_name)==1" )

      call this%grid_specification%get( base_object//".core.thickness.x", layers%core%thickness%x, found )
      call assert( found , base_object//".core.thickness.x found" )

      call this%grid_specification%get( base_object//".core.thickness.y", layers%core%thickness%y, found )
      call assert( found , base_object//".core.thickness.y found" )

      call this%grid_specification%get( base_object//".core.thickness.z", layers%core%thickness%z, found )
      call assert( found , base_object//".core.thickness.z found" )
    end subroutine

    subroutine verify_core_components
      !! apply constraints specified in the json file

      associate( &
        num_grid_blocks=>layers%core%num_grid_blocks, thickness=>layers%core%thickness, material_name=>layers%core%material_name )

        call assert( all( [num_grid_blocks%x,num_grid_blocks%y,num_grid_blocks%z] > 0 ), &
                    "all( [num_grid_blocks%x,num_grid_blocks%y,num_grid_blocks%z] > 0" )
        call assert( all( [thickness%x,thickness%y,thickness%z] > 0. ), &
                    "all( [thickness%x,thickness%y,thickness%z] > 0." )
        call assert( size(material_name)==1, "size(material_name)==1" )
        call assert( all( [size(thickness%x), size(thickness%y), size(thickness%z)] ==1 ), &
                    "all( [size(thickness%x), size(thickness%y), size(thickness%z)] ==1 )" )
        call assert( all( [size(num_grid_blocks%x), size(num_grid_blocks%y), size(num_grid_blocks%z)]==1 ), &
                    "all( [size(num_grid_blocks%x), size(num_grid_blocks%y), size(num_grid_blocks%z)]==1 )" )

      end associate
    end subroutine

    subroutine read_wrappers_components
      call this%grid_specification%get( base_object // ".wrappers.num_grid_blocks.x", layers%wrappers%num_grid_blocks%x, found)
      call assert( found, "set_block_metadata: wrappers_%num_grid_blocks_%x found" )

      call this%grid_specification%get( base_object // ".wrappers.num_grid_blocks.y", layers%wrappers%num_grid_blocks%y, found)
      call assert( found, "set_block_metadata: wrappers_%num_grid_blocks_%y found" )

      call this%grid_specification%get( base_object // ".wrappers.num_grid_blocks.z", layers%wrappers%num_grid_blocks%z, found)
      call assert( found, "set_block_metadata: wrappers_%num_grid_blocks_%z found" )

      call this%grid_specification%get( base_object//".wrappers.material_name", names, found )
      call assert( found , base_object//".wrappers.material_name found" )

      supremum = maxval( [( len( trim(names(i)) ), i=1,size(names) )] )
     !allocate( character(len=supremum) :: layers%wrappers%material_name(size(names)) ) !! precluded by gfortran 8.3 bug
      allocate(  layers%wrappers%material_name(size(names)) )
      layers%wrappers%material_name = names

      call this%grid_specification%get( base_object//".wrappers.thickness.x", layers%wrappers%thickness%x, found )
      call assert( found , base_object//".wrappers.thickness.x found" )

      call this%grid_specification%get( base_object//".wrappers.thickness.y", layers%wrappers%thickness%y, found )
      call assert( found , base_object//".wrappers.thickness.y found" )

      call this%grid_specification%get( base_object//".wrappers.thickness.z", layers%wrappers%thickness%z, found )
      call assert( found , base_object//".wrappers.thickness.z found" )
    end subroutine

    subroutine verify_wrappers_components
      !! apply constraints specified in the json file

      associate( &
        num_grid_blocks=>layers%wrappers%num_grid_blocks, thickness=>layers%wrappers%thickness, &
        material_name=>layers%wrappers%material_name, max_spacing=>layers%max_spacing )

        call assert( max_spacing > 0., "max_spacing > 0." )
        call assert( all([num_grid_blocks%x,num_grid_blocks%y,num_grid_blocks%z]>0), &
                    "all([num_grid_blocks%x,num_grid_blocks%y,num_grid_blocks%z]>0)" )
        call assert( all([thickness%x,thickness%y,thickness%z]>0), &
                    "all([thickness%x,thickness%y,thickness%z]>0)" )

        call assert( all( [size(thickness%x), size(thickness%y), size(thickness%z)] &
          == [size(num_grid_blocks%x), size(num_grid_blocks%y), size(num_grid_blocks%z)] ), &
          "all( [size(thickness%x), size(thickness%y), size(thickness%z)] &
          & == [size(num_grid_blocks%x), size(num_grid_blocks%y), size(num_grid_blocks%z)] )" )

        call assert( all( [size(thickness%x), size(thickness%y)] == size(material_name) ), &
                    "all( [size(thickness%x), size(thickness%y)] == size(material_name) )" )

      end associate
    end subroutine

    subroutine verify_layers

       associate( wrappers => layers%wrappers, core=> layers%core )
         call assert( all( wrappers%thickness%z >= core%thickness%z ), "all( wrappers%thickness%z >= core%thickness%z )" )
         call assert( &
           all( wrappers%num_grid_blocks%z >= core%num_grid_blocks%z ), &
           "all( wrappers%num_grid_blocks%z >= core%num_grid_blocks%z )")
       end associate

    end subroutine

    subroutine set_metadata

      character(len=max_name_length), parameter :: cavity="cavity"
      character(len=max_name_length) block_material

      associate( &
        nx_wrappers => layers%wrappers%num_grid_blocks%x, &
        ny_wrappers => layers%wrappers%num_grid_blocks%y, &
        nz_wrappers => layers%wrappers%num_grid_blocks%z, &
        nx_core => layers%core%num_grid_blocks%x(1), &
        ny_core => layers%core%num_grid_blocks%y(1), &
        nz_core => layers%core%num_grid_blocks%z(1), &
        wrappers_thickness_x => layers%wrappers%thickness%x, &
        wrappers_thickness_y => layers%wrappers%thickness%y, &
        wrappers_thickness_z => layers%wrappers%thickness%z, &
        core_thickness_x => layers%core%thickness%x, &
        core_thickness_y => layers%core%thickness%y, &
        core_thickness_z => layers%core%thickness%z )

        associate(  &
          nx => nx_core + symmetry*sum(nx_wrappers), &
          ny => ny_core + symmetry*sum(ny_wrappers), &
          nz => sum(nz_wrappers) )

          allocate(this%metadata(nx, ny, nz), stat=alloc_stat )
          call assert( alloc_stat==success, "set_block_metadata: allocate(this%metadata(nx,ny,nz),...)" )

          associate( &
            nx_layers => [nx_wrappers, nx_core, nx_wrappers(size(nx_wrappers):1:-1) ], &
            ny_layers => [ny_wrappers, ny_core, ny_wrappers(size(ny_wrappers):1:-1) ], &
            nz_layers => [nz_core, nz_wrappers - nz_core], &
            wrappers_material => layers%wrappers%material_name, &
            core_material => layers%core%material_name(1) )

            associate( &
              material => [wrappers_material, core_material, wrappers_material(size(wrappers_material):1:-1)], &
              thickness_x => [wrappers_thickness_x, core_thickness_x, wrappers_thickness_x(size(wrappers_thickness_x):1:-1)], &
              thickness_y => [wrappers_thickness_y, core_thickness_y, wrappers_thickness_y(size(wrappers_thickness_y):1:-1)], &
              thickness_z => [ core_thickness_z, wrappers_thickness_z - core_thickness_z ] )

              associate( &
                block_material_x => [( [(material(i), j=1,nx_layers(i))], i=1,size(nx_layers) )], &
                block_material_y => [( [(material(i), j=1,ny_layers(i))], i=1,size(ny_layers) )], &
                block_thickness_x => [( [( thickness_x(i)/nx_layers(i), j=1,nx_layers(i))], i=1,size(nx_layers) )], &
                block_thickness_y => [( [( thickness_y(i)/ny_layers(i), j=1,ny_layers(i))], i=1,size(ny_layers) )], &
                block_thickness_z => [( [( thickness_z(i)/nz_layers(i), j=1,nz_layers(i))], i=1,size(nz_layers) )], &
                max_spacing => layers%max_spacing, &
                first_core_y_block => sum(ny_wrappers) + 1, &
                last_core_y_block => sum(ny_wrappers) + ny_core )

                call assert( all( lbound(this%metadata)==[1,1,1] .and. ubound(this%metadata)==[nx,ny,nz]), &
                "all( lbound(this%metadata)==[1,1,1] .and. ubound(this%metadata)==[nx,ny,nz] ) ")

               !set_spacing_and_subdomain: &
               !do concurrent( iz=1:nz, iy=1:ny, ix=1:nx )
                do iz=1,nz
                do iy=1,ny
                do ix=1,nx

                   !block_material => merge( block_material_x(ix), &
                   !  merge(block_material_x(ix),block_material_y(iy), block_material_x(ix)/=core_material ), &
                   !  block_material_x(ix)==block_material_y(iy) ), &

                  call this%metadata(ix,iy,iz)%set_max_spacing(real(max_spacing,r8k))
                  associate( &
                    x_domain =>  [ sum( block_thickness_x(1:ix-1) ), sum( block_thickness_x(1:ix) ) ], &
                    y_domain =>  [ sum( block_thickness_y(1:iy-1) ), sum( block_thickness_y(1:iy) ) ], &
                    z_domain =>  [ sum( block_thickness_z(1:iz-1) ), sum( block_thickness_z(1:iz) ) ], &
                    tag => [wrappers_material, core_material, cavity] )

                    if ( block_material_y(iy) == core_material ) then
                      block_material = block_material_x(ix)
                    else if ( block_material_x(ix) == core_material )
                      block_material = block_material_y(iy)
                    else
                      block_material = block_material_y(iy)
                    end if


                   !associate( block_material_ =>merge(block_material, cavity, iz<=nz_core .and. block_material==core_material) )

                    call this%metadata(ix,iy,iz)%set_label( block_material )

                    call this%metadata(ix,iy,iz)%set_tag( findloc(tag, block_material, 1, back=.true.) )

                    print *,"ix, iy, iz, tag:",ix, iy, iz, block_material, findloc(tag, block_material, 1, back=.true.)

                   !end associate

                    call this%metadata(ix,iy,iz)%set_subdomain( &
                     subdomain_t( reshape([x_domain(1), y_domain(1), z_domain(1), x_domain(2), y_domain(2), z_domain(2)], &
                     [space_dimension,num_end_points]) ) )

                  end associate
               !end do set_spacing_and_subdomain
                end do; end do; end do

              end associate
            end associate
          end associate
        end associate
      end associate
    end subroutine

  end procedure

end submodule plate_3D_implementation
