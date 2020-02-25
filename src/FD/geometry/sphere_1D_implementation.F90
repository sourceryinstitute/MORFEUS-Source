!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
submodule(sphere_1D_interface) sphere_1D_implementation
  !! author: Damian Rouson and Karla Morris
  !! date: 2/24/2020
  use assertions_interface, only : assert
  use string_functions_interface, only : csv_format
#ifndef HAVE_FINDLOC
    use emulated_intrinsics_interface, only : findloc
#endif
  implicit none

  character(len=*), parameter :: base_object = "MORFEUS_FD"
  character(len=*), parameter :: layers_object = "MORFEUS_FD.layers"
  integer, parameter :: success=0

  !! Encapsulate json key/value pair hierarchy to be read from a 3Dplate*.json file

  type thickness_t
    real, allocatable, dimension(:) :: r, theta, phi
    character(len=:), allocatable :: dimensions
  end type

  type num_grid_blocks_t
    integer, allocatable, dimension(:) :: r, theta, phi
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
    character(len=*), parameter :: expected_geometry_type="1D_sphere"
    logical found

    call assert( file_extension(grid_description_file)=="json", 'file_extension(grid_description_file)=="json"' )

    call this%grid_specification%load_file( grid_description_file )
    call assert( .not. this%grid_specification%failed(), "set_grid_specification: .not. this%grid_specification%failed()" )

    call this%grid_specification%get( layers_object//".type",  layers%type_, found )
    call assert( found , layers_object//".type found" )
    call assert(layers%type_==expected_geometry_type, "layers%type_==expected_geometry_type" )

    call this%grid_specification%get( base_object//".units_system",  layers%units_system, found )
    call assert( found , base_object//".units_system found" )
    call assert(any(layers%units_system==units_system_names), "any(layers%units_system==units_system_names)" )

    call this%grid_specification%get( layers_object//".max_spacing",  layers%max_spacing, found )
    call assert( found , layers_object//".max_spacing found" )

  end procedure

  module procedure get_block_metadata_shape
    shape_ = shape(this%metadata)
  end procedure

  module procedure get_block_domain
    associate( ir=>indicial_coordinates(1), itheta=>indicial_coordinates(2), iphi=>indicial_coordinates(3))
      associate(lower => lbound(this%metadata), upper => ubound(this%metadata))
        call assert( all( [ [ir] >= lower(1), [ir] <= upper(1) ] ), "get_block_domain: indicial_coordinates in bounds" )
      end associate
      this_domain = this%metadata(ir,itheta,iphi)%get_subdomain()
    end associate
  end procedure

  module procedure get_block_metadatum
    associate( ir=>indicial_coordinates(1), itheta=>indicial_coordinates(2), iphi=>indicial_coordinates(3))
      associate(lower => lbound(this%metadata), upper => ubound(this%metadata))
        call assert( all( [ [ir] >= lower(1), [ir] <= upper(1) ] ), "get_block_metadatum: indicial_coordinates in bounds" )
      end associate
      this_metadata_xyz = this%metadata(ir, itheta, iphi)
    end associate
  end procedure

  module procedure get_block_metadata
    this_metadata = this%metadata
  end procedure

  module procedure set_block_metadata
#ifndef HAVE_FINDLOC
    use emulated_intrinsics_interface, only : findloc
#endif

    integer i, j, ir, itheta, iphi, alloc_stat, supremum
    logical found
    integer, parameter :: num_end_points=2, symmetry=2
    character(len=max_name_length), allocatable, dimension(:) :: names

    call read_core_components
    call verify_core_components

    call read_wrappers_components
    call verify_wrappers_components

    call verify_layers

    call set_metadata

#ifdef FORD
  end procedure ! compilers never see this line; when generating documentation, run "ford -m FORD..."  to circumvent a ford bug
#else
  contains ! compilers see this line
#endif

    subroutine read_core_components

      call this%grid_specification%get( layers_object // ".core.num_grid_blocks.r", layers%core%num_grid_blocks%r, found)
      call assert( found, "set_block_metadata: core_%num_grid_blocks_%r found" )

      call this%grid_specification%get( layers_object // ".core.num_grid_blocks.theta", layers%core%num_grid_blocks%theta, found)
      call assert( .not. found, "set_block_metadata: core_%num_grid_blocks_%theta found" )
      layers%core%num_grid_blocks%theta = [0]

      call this%grid_specification%get( layers_object // ".core.num_grid_blocks.phi", layers%core%num_grid_blocks%phi, found)
      call assert( .not. found, "set_block_metadata: core_%num_grid_blocks_%phi found" )
      layers%core%num_grid_blocks%phi = [0]

      call this%grid_specification%get( layers_object//".core.material_name", names, found )
      call assert( found , layers_object//".core.material_name found" )

      supremum = maxval( [( len( trim(names(i)) ), i=1,size(names) )] )
      !allocate( character(len=supremum) :: layers%core%material_name(size(names)) ) !! precluded by gfortran 8.3 bug
      if (allocated(layers%core%material_name)) deallocate(layers%core%material_name)
      allocate(  layers%core%material_name(size(names)) )
      layers%core%material_name = names
      call assert(size(layers%core%material_name)==1, "size(layers%core%material_name)==1" )

      call this%grid_specification%get( layers_object//".core.thickness.r", layers%core%thickness%r, found )
      call assert( found , layers_object//".core.thickness.r found" )

      call this%grid_specification%get( layers_object//".core.thickness.theta", layers%core%thickness%theta, found )
      call assert( .not. found , layers_object//".core.thickness.theta found" )
      layers%core%thickness%theta = [0.0]

      call this%grid_specification%get( layers_object//".core.thickness.phi", layers%core%thickness%phi, found )
      call assert( .not. found , layers_object//".core.thickness.phi found" )
      layers%core%thickness%phi = [0.0]
    end subroutine

    subroutine verify_core_components
      !! apply constraints specified in the json file

      associate( &
        num_grid_blocks=>layers%core%num_grid_blocks, thickness=>layers%core%thickness, material_name=>layers%core%material_name )

        call assert( all( [num_grid_blocks%r] > 0 ), "all( [num_grid_blocks%r] > 0" )
        call assert( all( [num_grid_blocks%theta,num_grid_blocks%phi] == 0 ), &
                   " all( [num_grid_blocks%theta,num_grid_blocks%phi] == 0")
        call assert( all( [thickness%r] > 0. ), "all( [thickness%r] > 0." )
        call assert( all( [thickness%theta,thickness%phi] == 0. ), "all( [thickness%theta,thickness%phi] == 0." )
        call assert( size(material_name)==1, "size(material_name)==1" )
        call assert( all( [size(thickness%r), size(thickness%theta), size(thickness%phi)] ==1 ), &
                    "all( [size(thickness%r), size(thickness%theta), size(thickness%phi)] ==1 )" )
        call assert( all( [size(num_grid_blocks%r), size(num_grid_blocks%theta), size(num_grid_blocks%phi)]==1 ), &
                    "all( [size(num_grid_blocks%r), size(num_grid_blocks%theta), size(num_grid_blocks%phi)]==1 )" )

      end associate
    end subroutine

    subroutine read_wrappers_components
      call this%grid_specification%get( layers_object // ".wrappers.num_grid_blocks.r", layers%wrappers%num_grid_blocks%r, found)
      call assert( found, "set_block_metadata: wrappers_%num_grid_blocks_%r found" )

      call this%grid_specification%get( layers_object // ".wrappers.num_grid_blocks.theta", &
        layers%wrappers%num_grid_blocks%theta, found)
      call assert( .not. found, "set_block_metadata: wrappers_%num_grid_blocks_%theta found" )
      layers%wrappers%num_grid_blocks%theta = [0]

      call this%grid_specification%get( layers_object // ".wrappers.num_grid_blocks.phi", layers%wrappers%num_grid_blocks%phi,found)
      call assert( .not. found, "set_block_metadata: wrappers_%num_grid_blocks_%phi found" )
      layers%wrappers%num_grid_blocks%phi = [0]

      call this%grid_specification%get( layers_object//".wrappers.material_name", names, found )
      call assert( found , layers_object//".wrappers.material_name found" )

      supremum = maxval( [( len( trim(names(i)) ), i=1,size(names) )] )
      !allocate( character(len=supremum) :: layers%wrappers%material_name(size(names)) ) !! precluded by gfortran 8.3 bug
      if (allocated(layers%wrappers%material_name)) deallocate(layers%wrappers%material_name)
      allocate(  layers%wrappers%material_name(size(names)) )
      layers%wrappers%material_name = names

      call this%grid_specification%get( layers_object//".wrappers.thickness.r", layers%wrappers%thickness%r, found )
      call assert( found , layers_object//".wrappers.thickness.r found" )

      call this%grid_specification%get( layers_object//".wrappers.thickness.theta", layers%wrappers%thickness%theta, found )
      call assert( .not. found , layers_object//".wrappers.thickness.theta found" )
      layers%wrappers%thickness%theta = [0.0]

      call this%grid_specification%get( layers_object//".wrappers.thickness.phi", layers%wrappers%thickness%phi, found )
      call assert( .not. found , layers_object//".wrappers.thickness.phi found" )
      layers%wrappers%thickness%phi = [0.0]
    end subroutine

    subroutine verify_wrappers_components
      !! apply constraints specified in the json file

      associate( &
        num_grid_blocks=>layers%wrappers%num_grid_blocks, thickness=>layers%wrappers%thickness, &
        material_name=>layers%wrappers%material_name, max_spacing=>layers%max_spacing )

        call assert( max_spacing > 0., "max_spacing > 0." )
        call assert( all([num_grid_blocks%r]>0), "all([num_grid_blocks%r]>0)" )
        call assert( all([num_grid_blocks%theta,num_grid_blocks%phi]==0), &
                    "all([num_grid_blocks%theta,num_grid_blocks%phi]==0)" )
        call assert( all([thickness%r]>0), "all([thickness%r]>0)" )
        call assert( all([thickness%theta,thickness%phi]==0), "all([thickness%theta,thickness%phi]==0)" )

        call assert( all( [size(thickness%r), size(thickness%theta), size(thickness%phi)] &
          == [size(num_grid_blocks%r), size(num_grid_blocks%theta), size(num_grid_blocks%phi)] ), &
          "all( [size(thickness%r), size(thickness%theta), size(thickness%phi)] " // &
          "== [size(num_grid_blocks%r), size(num_grid_blocks%theta), size(num_grid_blocks%phi)] )" )

        call assert( all( [size(thickness%r)] == size(material_name) ), &
                    "all( [size(thickness%r)] == size(material_name) )" )

      end associate
    end subroutine

    subroutine verify_layers

      associate( wrappers => layers%wrappers, core=> layers%core )
        call assert( all( wrappers%thickness%phi >= core%thickness%phi ), "all( wrappers%thickness%phi >= core%thickness%phi )" )
        call assert( &
          all( wrappers%num_grid_blocks%phi >= core%num_grid_blocks%phi ), &
          "all( wrappers%num_grid_blocks%phi >= core%num_grid_blocks%phi )")
      end associate

    end subroutine

    subroutine set_metadata

      character(len=max_name_length), parameter :: cavity="cavity"

      associate( &
        nr_wrappers => layers%wrappers%num_grid_blocks%r, &
        ntheta_wrappers => layers%wrappers%num_grid_blocks%theta, &
        nphi_wrappers => layers%wrappers%num_grid_blocks%phi, &
        nr_core => layers%core%num_grid_blocks%r(1), &
        ntheta_core => layers%core%num_grid_blocks%theta(1), &
        nphi_core => layers%core%num_grid_blocks%phi(1), &
        wrappers_thickness_r => layers%wrappers%thickness%r, &
        wrappers_thickness_theta => layers%wrappers%thickness%theta, &
        wrappers_thickness_phi => layers%wrappers%thickness%phi, &
        core_thickness_r => layers%core%thickness%r, &
        core_thickness_theta => layers%core%thickness%theta, &
        core_thickness_phi => layers%core%thickness%phi )

        associate(  &
          nr => nr_core + symmetry*sum(nr_wrappers), &
          ntheta => 1, &
          nphi => 1 )

          allocate(this%metadata(nr, ntheta, nphi), stat=alloc_stat )
          call assert( alloc_stat==success, "set_block_metadata: allocate(this%metadata(nr,ntheta,nphi),...)" )

          associate( &
            nr_layers => [nr_wrappers, nr_core, nr_wrappers(size(nr_wrappers):1:-1) ], &
            ntheta_layers => [ntheta_wrappers, ntheta_core, ntheta_wrappers(size(ntheta_wrappers):1:-1) ], &
            nphi_layers => [nphi_core, nphi_wrappers - nphi_core], &
            wrappers_material => layers%wrappers%material_name, &
            core_material => layers%core%material_name(1) )

            associate( &
              thickness_r => [wrappers_thickness_r, core_thickness_r, wrappers_thickness_r(size(wrappers_thickness_r):1:-1)], &
              thickness_theta => [wrappers_thickness_theta, core_thickness_theta, &
                wrappers_thickness_theta(size(wrappers_thickness_theta):1:-1)], &
              thickness_phi => [ core_thickness_phi, wrappers_thickness_phi - core_thickness_phi ] )

              associate( &
                block_thickness_r => [( [( thickness_r(i)/nr_layers(i), j=1,nr_layers(i))], i=1,size(nr_layers) )], &
                block_thickness_theta => [0.0], &
                block_thickness_phi => [0.0] )

                call assert( all( lbound(this%metadata)==[1,1,1] .and. ubound(this%metadata)==[nr,ntheta,nphi]), &
                "all( lbound(this%metadata)==[1,1,1] .and. ubound(this%metadata)==[nr,ntheta,nphi] ) ")

                do iphi=1,nphi
                do itheta=1,ntheta
                do ir=1,nr

                  call this%metadata(ir,itheta,iphi)%set_max_spacing(real(layers%max_spacing, r8k))
                  associate( &
                    r_domain =>  [ sum( block_thickness_r(1:ir-1) ), sum( block_thickness_r(1:ir) ) ], &
                    theta_domain =>  [0.0, 0.0] , &
                    phi_domain =>  [0.0, 0.0], &
                    tag => [wrappers_material, core_material, cavity], &
                    block_material => material(ir,itheta,iphi,layers%core,layers%wrappers) )

                    call this%metadata(ir,itheta,iphi)%set_label( block_material )
                    call this%metadata(ir,itheta,iphi)%set_tag( findloc(tag, block_material, 1, back=.true.) )
                    call this%metadata(ir,itheta,iphi)%set_subdomain( &
                      subdomain_t( reshape([r_domain(1), theta_domain(1), phi_domain(1), &
                        r_domain(2), theta_domain(2), phi_domain(2)], [space_dimension,num_end_points]) ) )

                  end associate
                end do; end do; end do

              end associate
            end associate
          end associate
        end associate
      end associate
    end subroutine

#ifndef FORD
  end procedure ! compilers never see this line; when generating documentation, run "ford -m FORD ..." to circumvent a ford bug
#endif

  function material(ir, itheta, iphi, core_, wrapper_) result(material_ir_itheta)
    integer, intent(in) :: ir, itheta, iphi
    type(material_t), intent(in) :: core_, wrapper_
    integer i, j
    character(len=max_name_length), allocatable :: material_ir_itheta
    character(len=max_name_length), parameter :: void_name="cavity"

    associate( core_material_name => merge( core_%material_name, void_name, iphi <= core_%num_grid_blocks%phi ) )
      associate( &
        nr_layers => [wrapper_%num_grid_blocks%r, core_%num_grid_blocks%r, &
          wrapper_%num_grid_blocks%r(size(wrapper_%num_grid_blocks%r):1:-1) ], &
        ntheta_layers => [wrapper_%num_grid_blocks%theta, core_%num_grid_blocks%theta, &
          wrapper_%num_grid_blocks%theta(size(wrapper_%num_grid_blocks%theta):1:-1) ], &
        material => [wrapper_%material_name, core_material_name, wrapper_%material_name(size(wrapper_%material_name):1:-1)] )

        associate( &
          block_material_r => [( [(material(i), j=1,nr_layers(i))], i=1,size(nr_layers) )], &
          block_material_theta => [( [(material(i), j=1,ntheta_layers(i))], i=1,size(ntheta_layers) )] )

          associate( wrapper_material_r => replace_layers(block_material_r, block_material_theta, itheta) )
            material_ir_itheta = wrapper_material_r(ir)
          end associate
        end associate
      end associate
    end associate
  end function

  function replace_layers(full_delineation_r, full_delineation_theta, itheta) result(wrapper_layer_r)
    character(len=*), intent(in) :: full_delineation_r(:), full_delineation_theta(:)
    character(len=len(full_delineation_r)) :: wrapper_layer_r( size(full_delineation_r) )
    integer itheta

    associate( layer=>full_delineation_theta(itheta) )
      associate( &
        first => findloc(full_delineation_r, layer, dim=1, back=.false.), &
        last => findloc(full_delineation_r, layer, dim=1, back=.true.) )

        wrapper_layer_r(1:first-1) = full_delineation_r(1:first-1)
        wrapper_layer_r(first:last) = layer
        wrapper_layer_r(last+1:) = full_delineation_r(last+1:)
      end associate
    end associate

  end function
end submodule sphere_1D_implementation
