!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program main
  !!author:  Damian Rouson
  !!
  !! At inception, this is a standalone test of the logic for iterating through the
  !! structured-grid blocks defined in 3Dplate.json.  This might be later refactored
  !! to test the corresponding logic in the plate_3D_implmentation module procedure
  !! define_global_blocks, at which point it will be set up for automated execution
  !! via ctest.
  implicit none

  type material_t
    character(len=len('burrito')), allocatable :: material_name(:)
    integer, dimension(:), allocatable :: nx_blocks, ny_blocks, nz_blocks
  end type

  type(material_t) core, wrapper
  integer, parameter :: symmetry=2
  character(len=len(core%material_name)) , allocatable :: mat(:,:,:)

  core%material_name = ['burrito']
  wrapper%material_name = [character(len=len(core%material_name)) :: 'bag', 'air', 'foil']
  core%nx_blocks = [1]
  core%ny_blocks = [2]
  core%nz_blocks = [1]
  wrapper%nx_blocks = [3,2,1]
  wrapper%ny_blocks = [1,2,3]
  wrapper%nz_blocks = [2]

  if (wrapper%nz_blocks(1) < core%nz_blocks(1)) error stop "burrito not fully covered"

  associate( &
    nx => core%nx_blocks(1) + symmetry*sum(wrapper%nx_blocks), &
    ny => core%ny_blocks(1) + symmetry*sum(wrapper%ny_blocks), &
    nz => wrapper%nz_blocks(1) )

    allocate( mat(nx, ny, nz) )

    block
      character(len=*), parameter :: csv_format = '(*(G0,:,","))'
      integer i, j, k
      do k=1,nz
        do j=1,ny
          do i=1,nx
            mat(i,j,k) = material(i, j, k, core, wrapper)
          end do
          write(*,csv_format) mat(:,j,k)
        end do
        write(*,*)
      end do
    end block

  end associate

contains

  function material(ix, iy, iz, core_, wrapper_) result(material_ix_iy)
    integer, intent(in) :: ix, iy, iz
    type(material_t), intent(in) :: core_, wrapper_
    integer i, j, k
    character(len=:), allocatable :: material_ix_iy
    character(len=len(core_%material_name)), parameter :: void_name="void"

    associate( core_material_name => merge( core_%material_name, void_name, iz <= core_%nz_blocks ) )
      associate( &
        nx_layers => [wrapper_%nx_blocks, core_%nx_blocks, wrapper_%nx_blocks(size(wrapper_%nx_blocks):1:-1) ], &
        ny_layers => [wrapper_%ny_blocks, core_%ny_blocks, wrapper_%ny_blocks(size(wrapper_%ny_blocks):1:-1) ], &
        material => [wrapper_%material_name, core_material_name, wrapper_%material_name(size(wrapper_%material_name):1:-1)] )

        associate( &
          block_material_x => [( [(material(i), j=1,nx_layers(i))], i=1,size(nx_layers) )], &
          block_material_y => [( [(material(i), j=1,ny_layers(i))], i=1,size(ny_layers) )] )

          associate( wrapper_material_x => replace_layers(block_material_x, block_material_y, iy) )
            material_ix_iy = wrapper_material_x(ix)
          end associate
        end associate
      end associate
    end associate
  end function

  function replace_layers(full_delineation_x, full_delineation_y, iy) result(wrapper_layer_x)
    character(len=*), intent(in) :: full_delineation_x(:), full_delineation_y(:)
    character(len=len(full_delineation_x)) :: wrapper_layer_x( size(full_delineation_x) )
    integer iy

    associate( layer=>full_delineation_y(iy) )
      associate( &
        first => findloc(full_delineation_x, layer, dim=1, back=.false.), &
        last => findloc(full_delineation_x, layer, dim=1, back=.true.) )

        wrapper_layer_x(1:first-1) = full_delineation_x(1:first-1)
        wrapper_layer_x(first:last) = layer
        wrapper_layer_x(last+1:) = full_delineation_x(last+1:)
      end associate
    end associate

  end function

#ifdef __GNUC__
#if __GNUC__ >= 9
#define HAVE_FINDLOC
#endif
#endif

#ifndef HAVE_FINDLOC

  pure function findloc(array, value, dim, back) result(location)
    implicit none
    integer, intent(in) :: array(:), value, dim
    logical, intent(in) :: back
    integer location

    integer, parameter :: loop_increment=-1, base_index=1
    integer index_

    if ( (.not. back) .or. dim/=1) error stop "findloc: unsupported use case"

    associate( lower_bound=>lbound(array,dim) )
      do index_=ubound(array,dim), lower_bound, loop_increment
        if (array(index_)==value) then
          location = index_ - lower_bound + base_index
          return
        end if
      end do
    end associate
    location=0

  end function
#endif

end program main
