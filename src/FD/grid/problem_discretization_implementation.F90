!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
submodule(problem_discretization_interface) define_problem_discretization
  !! author: Damian Rouson and Karla Morris
  !! date: 9/9/2019
  use assertions_interface, only : assert,assertions
  use iso_fortran_env, only : error_unit
  use kind_parameters, only : i4k, r8k
  implicit none

  integer, parameter :: space_dimensions=3

contains

  module procedure write_formatted
    use string_functions_interface, only : file_extension
    integer subdomain

    character(len=:), allocatable :: file_name, suffix
    integer, parameter :: max_name_length=128

    iostat=0

    allocate(character(len=max_name_length)::file_name)
    inquire(unit=unit,name=file_name)

    suffix = file_extension( file_name )
    select case(suffix)
      case('csv')
        call csv_output
      case('vtk')
        CALL VTK_output
      case('json')
        CALL json_output
      case default
        error stop "problem_discretization%write_formatted: unsupported file name extension" // &
#ifdef HAVE_NON_CONSTANT_ERROR_STOP
        & file_name ! Fortran 2018
#else
        & "."       ! Fortran 2008
#endif
      end select

    contains

      subroutine json_output
        error stop "json_output: not yet implemented"
      end subroutine

      subroutine csv_output
        integer ix,iy,iz
        associate( nx=>this%global_block_shape_(1), ny=>this%global_block_shape_(2), nz=>this%global_block_shape_(3) )
          write(unit,'(4("      ",a,:,",",5x))') "x","y","z","layer (phony)"
          write(unit,*) new_line('a')
          do iz=1,nz
            do iy=1,ny
              do ix=1,nx
#ifdef HAVE_UDDTIO
                write(unit,*) this%vertices( this%block_identifier([ix,iy,iz])  )
#else
                CALL this%vertices( this%block_identifier([ix,iy,iz])  )%write_formatted(unit,iotype, v_list, iostat, iomsg)
#endif
              end do
            end do
          end do
        end associate
      end subroutine

      SUBROUTINE define_scalar( s, vals, dataname )
          USE vtk_attributes, ONLY : scalar, attributes
          CLASS(attributes), INTENT(INOUT) :: s
          REAL(r8k),         INTENT(IN)    :: vals(:)
          CHARACTER(LEN=*),  INTENT(IN)    :: dataname

          IF (.NOT. ALLOCATED(s%attribute)) ALLOCATE(scalar::s%attribute)
          CALL s%attribute%init (dataname, numcomp=1, real1d=vals)
      END SUBROUTINE

      SUBROUTINE  VTK_output
            USE vtk_datasets,   ONLY : unstruct_grid
            USE vtk,            ONLY : vtk_legacy_write
            USE vtk_cells, ONLY : voxel, vtkcell_list
            USE vtk_attributes, ONLY : attributes
            USE array_functions_interface, ONLY : OPERATOR(.catColumns.), OPERATOR(.columnVectors.)
            TYPE(voxel) :: voxel_cell  !! Voxel cell type
            TYPE(vtkcell_list), DIMENSION(:), ALLOCATABLE :: cell_list !! Full list of all cells
            TYPE (unstruct_grid) :: vtk_grid
            INTEGER(i4k) :: ip, jp, kp !! block-local point id
            INTEGER(i4k) :: ic, jc, kc !! block-local cell id
            REAL(r8k), DIMENSION(:,:), ALLOCATABLE :: points
            INTEGER(i4k), DIMENSION(:), ALLOCATABLE :: block_cell_material, point_block_id
            TYPE(attributes) :: point_values, cell_values
            INTEGER :: b, first_point_in_block, first_cell_in_block
            INTEGER, PARAMETER :: vector_indices=4, writer=1

            CALL assert(this_image()==writer, "problem_discretization%write_formatted: output from image 1")

            ALLOCATE( cell_list(SUM( this%vertices%num_cells())) )
            ALLOCATE( points(space_dimensions,0), block_cell_material(0), point_block_id(0) )

            first_point_in_block = 0
            first_cell_in_block = 1

            loop_over_grid_blocks: DO b=lbound(this%vertices,1), ubound(this%vertices,1)

              ASSOCIATE( vertex_positions => this%vertices(b)%vectors() ) !! 4D array of nodal position vectors
                ASSOCIATE( npoints => shape(vertex_positions(:,:,:,1)) )  !! shape of 1st 3 dims specifies # points in ea. direction
                  ASSOCIATE( ncells => npoints-1 )

                    points = points .catColumns. ( .columnVectors. vertex_positions )
                    block_cell_material = [ block_cell_material, (this%vertices(b)%get_tag(), ic=1,PRODUCT(ncells)) ]
                    point_block_id = [ point_block_id, [( b, ip=1, PRODUCT(npoints) )] ]

                    DO ic=1,ncells(1)
                      DO jc=1,ncells(2)
                        DO kc=1,ncells(3)
                          ASSOCIATE( block_local_point_id => & !! 8-element array of block-local ID's for voxel corners
                            [( [( [( kp*PRODUCT(npoints(1:2)) + jp*npoints(1) + ip, ip=ic,ic+1 )], jp=jc-1,jc )], kp=kc-1,kc )] &
                          )
                            CALL voxel_cell%setup ( first_point_in_block + block_local_point_id-1 )
                          END ASSOCIATE
                          ASSOCIATE( local_cell_id => (kc-1)*PRODUCT(ncells(1:2)) + (jc-1)*ncells(1) + ic )
                            ASSOCIATE(  cell_id => first_cell_in_block + local_cell_id - 1 )
                              ALLOCATE(cell_list(cell_id)%cell, source=voxel_cell)!GCC8 workaround for cell_list(i)%cell= voxel_cell
                            END ASSOCIATE
                          END ASSOCIATE
                        END DO
                      END DO
                    END DO

                    first_cell_in_block  = first_cell_in_block  + PRODUCT( ncells  )
                    first_point_in_block = first_point_in_block + PRODUCT( npoints )
                  END ASSOCIATE
                END ASSOCIATE
              END ASSOCIATE
            END DO loop_over_grid_blocks

            CALL assert( SIZE(point_block_id,1) == SIZE(points,2), "VTK block point data set & point set size match" )
            CALL assert( SIZE(block_cell_material,1) == SIZE(cell_list,1), "VTK cell data set & cell set size match" )

            CALL define_scalar(  cell_values, REAL( block_cell_material, r8k),  'material' )
            CALL define_scalar(  point_values, REAL( point_block_id, r8k),  'block' )

            CALL vtk_grid%init (points=points, cell_list=cell_list )
            CALL vtk_legacy_write( &
              unit=unit, geometry=vtk_grid, title='Morfeus-FD voxels', multiple_io=.TRUE., &
              celldatasets=[cell_values], pointdatasets=[point_values] )

        END SUBROUTINE VTK_output

  end procedure write_formatted

  pure function evenly_spaced_points( boundaries, resolution, direction ) result(grid_nodes)
    !! Define grid point coordinates with uniform spacing in the chosen subdomain
    real(r8k), intent(in) :: boundaries(:,:)
      !! subdomain boundaries of each coordinate direction
    integer, intent(in) :: resolution(:)
      !! number of grid points in each direction
    integer, intent(in) :: direction
      !! coordinate direction to define
    real, allocatable :: grid_nodes(:,:,:)
      !! grid node locations and spacing in each coordination direction
    real(r8k) dx(space_dimensions)
    integer alloc_status
    character(len=128) :: alloc_error

    integer, parameter :: lo_bound=1, up_bound=2, success=0, num_boundaries=2
    integer ix,iy,iz

    allocate(grid_nodes(resolution(1),resolution(2),resolution(3)), stat=alloc_status, errmsg=alloc_error )
    CALL assert( alloc_status==success, "evenly_spaced_points allocation ("//alloc_error//")", PRODUCT(resolution) )

    associate( num_intervals => resolution - 1 )
      dx = ( boundaries(:,up_bound) - boundaries(:,lo_bound) ) / num_intervals(:)
    end associate

    associate( nx=>resolution(1), ny=>resolution(2), nz=>resolution(3) )

    select case(direction)
      case(1)
        do concurrent(iy=1:ny,iz=1:nz)
          associate( internal_points => boundaries(direction,lo_bound) + [(ix*dx(direction),ix=1,nx-num_boundaries)] )
            grid_nodes(:,iy,iz) = [ boundaries(direction,lo_bound), internal_points , boundaries(direction,up_bound) ]
          end associate
        end do
      case(2)
        do concurrent(ix=1:nx,iz=1:nz)
          associate( internal_points => boundaries(direction,lo_bound) + [(iy*dx(direction),iy=1,ny-num_boundaries)] )
            grid_nodes(ix,:,iz) = [ boundaries(direction,lo_bound), internal_points , boundaries(direction,up_bound) ]
          end associate
        end do
      case(3)
        do concurrent(ix=1:nx,iy=1:ny)
          associate( internal_points => boundaries(direction,lo_bound) + [(iz*dx(direction),iz=1,nz-num_boundaries)] )
            grid_nodes(ix,iy,:) = [ boundaries(direction,lo_bound), internal_points, boundaries(direction,up_bound) ]
          end associate
        end do
      case default
        error stop "evenly_spaced_points: invalid direction"
    end select

    end associate

  end function

  module procedure initialize_from_plate_3D
    integer, parameter :: lo_bound=1, up_bound=2 !! array indices corresponding to end points on 1D spatial interval
    integer, parameter :: nx_min=2, ny_min=2, nz_min=2
    integer n

    call this%partition( plate_3D_geometry%get_block_metadata_shape() )
      !! partition a block-structured grid into subdomains with connectivity implied by the indexing of the 3D array of blocks

      associate( my_subdomains => this%my_subdomains() )
        do n = my_subdomains(lo_bound) , my_subdomains(up_bound) ! TODO: make concurrent after Intel supports co_sum

          associate( ijk => this%block_indicial_coordinates(n) )

            associate( metadata => plate_3D_geometry%get_block_metadatum(ijk))

              call this%vertices(n)%set_metadata( metadata  )

              associate( &
                subdomain => plate_3D_geometry%get_block_domain(ijk), &
                max_spacing => metadata%get_max_spacing() &
              )
                associate( &
                  nx => max( nx_min, floor( abs(subdomain(1,up_bound) - subdomain(1,lo_bound))/max_spacing ) + 1 ), &
                  ny => max( ny_min, floor( abs(subdomain(2,up_bound) - subdomain(2,lo_bound))/max_spacing ) + 1 ), &
                  nz => max( nz_min, floor( abs(subdomain(3,up_bound) - subdomain(3,lo_bound))/max_spacing ) + 1 ) &
                )
                  associate( &
                    x => evenly_spaced_points(  subdomain, [nx,ny,nz], direction=1 ), &
                    y => evenly_spaced_points(  subdomain, [nx,ny,nz], direction=2 ), &
                    z => evenly_spaced_points(  subdomain, [nx,ny,nz], direction=3 ) )
                    call this%set_vertices(x,y,z,block_identifier=n)
                  end associate
                end associate
              end associate
            end associate
          end associate
        end do
      end associate

  end procedure

  module procedure partition

    integer :: alloc_status, image
    !! error checking code, image number

    ! Requires
    if (assertions) call assert(size(global_block_shape)==3,"partition: 3D structured_grid blocks")

    this%global_block_shape_ =  global_block_shape

    associate( num_blocks => product(this%global_block_shape_) )
    associate( me => this_image() )
    associate( ni => num_images() )
    associate( remainder => mod(num_blocks,ni) )
    associate( quotient => num_blocks/ni )
    associate( my_first => 1 + sum([(quotient+overflow(image,remainder),image=1,me-1)]) )
    associate( my_last => my_first + quotient + overflow(me,remainder) - 1 )

    allocate( this%vertices(my_first:my_last), stat=alloc_status )
    !! allocate this image's subset of the vertices

    if (assertions) then
      call assert(alloc_status==0,"partition: data distribution established")
      block
#ifndef HAVE_COLLECTIVE_SUBROUTINES
        use emulated_intrinsics_interface, only : co_sum
#endif
        integer total_blocks
        total_blocks = size(this%vertices)
        call co_sum(total_blocks,result_image=1)
        if (me==1) call assert(total_blocks==num_blocks,"all blocks have been distributed amongst the images")
        sync all
      end block
    end if

    end associate; end associate; end associate; end associate; end associate; end associate; end associate


    ! Assures
    call this%mark_as_defined

  contains

    pure function overflow(image,remainder) result(filler)
      integer, intent(in) :: image,remainder
      integer :: filler
      filler = merge(1,0,image<=remainder)
    end function

  end procedure

  module procedure my_subdomains

    ! Requires
    if (assertions) call assert(size(this%global_block_shape_)==3,"partition: 3D structured_grid blocks")

    block_identifier_range = [ lbound(this%vertices), ubound(this%vertices) ]

  end procedure

  module procedure block_identifier

    character(len=256) :: diagnostic_string

    ! Requires
    if (assertions) then
      associate(assertion => all(ijk>[0,0,0]).and.all(ijk<=this%global_block_shape_))
        if (assertions) write(diagnostic_string,*) "all(",ijk,">[0,0,0]), all(",ijk,"<=",this%global_block_shape_,")"
        call assert( assertion , "block_identifier: indicial coordinates in bounds", diagnostic_data=diagnostic_string)
      end associate
    end if

    associate( i=>ijk(1), nx=>this%global_block_shape_(1) )
    associate( j=>ijk(2), ny=>this%global_block_shape_(2) )
    associate( k=>ijk(3))
      n = (k-1)*(ny*nx) + (j-1)*ny + i
    end associate; end associate; end associate

    ! Assures
    if (assertions) &
      call assert( n>0 .and. n<=product(this%global_block_shape_), "block_identifier: block identifier in bounds" )

  end procedure

  module procedure block_indicial_coordinates

    ! Requires
    if (assertions) &
      call assert( n>0 .and. n<=product(this%global_block_shape_), "block_indicial_coordinates: identifier in bounds" )

    associate( nx=>this%global_block_shape_(1) )
    associate( ny=>this%global_block_shape_(2) )
    associate( nz=>this%global_block_shape_(3) )
      ijk = [ mod(n-1,nx)+1, mod( (n-1)/nx, ny ) + 1, (n-1)/(nx*ny) + 1]
    end associate; end associate; end associate

    ! Assures
    if (assertions) &
      call assert( all(ijk>[0,0,0]).and.all(ijk<=this%global_block_shape_),"block_indicial_coordinates: coordinates in bounds")

  end procedure

  module procedure block_load

    ! Requires
    call assert(this%user_defined(),"block_load: this object defined")

    if (allocated(this%vertices)) then
      num_blocks = size(this%vertices)
    else
      num_blocks = 0
    end if

  end procedure

  module procedure user_defined_vertices

    ! Requires
    call assert( &
      block_identifier>=lbound(this%vertices,1) .and. block_identifier<=ubound(this%vertices,1), &
      "set_vertices: block identifier in bounds" &
    )

    associate(n=>block_identifier)
      call this%vertices(n)%set_vector_components(x_nodes,y_nodes,z_nodes)
    end associate

  end procedure

end submodule define_problem_discretization
