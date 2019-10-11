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

  integer i,j

  type material_t
    character(len=len('burrito')), dimension(:), allocatable :: material_name
    integer, dimension(:), allocatable :: nx_blocks, ny_blocks
  end type

  type(material_t) core, wrapper

  integer, parameter :: symmetry=2

  core%material_name = ['burrito']
  wrapper%material_name = [character(len=len('burrito')) :: 'bag', 'foil']
  core%nx_blocks = [2]
  core%ny_blocks = [1]
  wrapper%nx_blocks = [2,1]
  wrapper%ny_blocks = [1,1]

  associate( &
    nx => core%nx_blocks(1) + symmetry*sum(wrapper%nx_blocks), &
    ny => core%ny_blocks(1) + symmetry*sum(wrapper%ny_blocks) )

    do i=1,nx
      do j=1,ny
        call material(i, j, core, wrapper)
      end do
    end do

  end associate

contains

  subroutine material(ix, iy, core_, wrapper_)
    integer, intent(in) :: ix, iy
    type(material_t), intent(in) :: core_, wrapper_
    integer m

    associate( &
      nx_layers => [wrapper_%nx_blocks, core_%nx_blocks, wrapper_%nx_blocks(size(wrapper_%nx_blocks):1:-1) ], &
      ny_layers => [wrapper_%ny_blocks, core_%ny_blocks, wrapper_%ny_blocks(size(wrapper_%ny_blocks):1:-1) ], &
      material => [wrapper_%material_name, core_%material_name, wrapper_%material_name(size(wrapper_%material_name):1:-1)] )

      associate( &
        block_material_x => [( [(material(i), j=1,nx_layers(i))], i=1,size(nx_layers) )], &
        block_material_y => [( [(material(i), j=1,ny_layers(i))], i=1,size(ny_layers) )] )

       if ( (block_material_x(ix) == core_%material_name(1)) .and. (block_material_y(iy) == core_%material_name(1)) ) then
         print *,ix,",", iy, ":", block_material_x(ix), " | ", block_material_y(iy), " => ",core_%material_name(1)
       else if ( &
         (ix<=sum(wrapper_%nx_blocks)) .and. (iy>sum(wrapper_%ny_blocks)) .and. (iy<=sum(wrapper_%ny_blocks)+core%ny_blocks(1)))then
         print *,ix,",", iy, ":", block_material_x(ix), " | ", block_material_y(iy), " => ",block_material_x(ix)
       else if ( &
         (ix>sum(wrapper_%nx_blocks)+core_%nx_blocks(1)) .and. &
         (iy>sum(wrapper_%ny_blocks)) .and. (iy<=sum(wrapper_%ny_blocks)+core%ny_blocks(1))) then
         print *,ix,",", iy, ":", block_material_x(ix), " | ", block_material_y(iy), " => ",block_material_x(ix)
       else if ( &
         (iy<=sum(wrapper_%ny_blocks)) .and. (ix>sum(wrapper_%nx_blocks)) .and. (ix<=sum(wrapper_%nx_blocks)+core%nx_blocks(1)))then
         print *,ix,",", iy, ":", block_material_x(ix), " | ", block_material_y(iy), " => ",block_material_y(iy)
       else if ( &
         (iy>sum(wrapper_%ny_blocks)+core_%ny_blocks(1)) .and. &
         (ix>sum(wrapper_%nx_blocks)) .and. (ix<=sum(wrapper_%nx_blocks)+core%nx_blocks(1)))then
         print *,ix,",", iy, ":", block_material_x(ix), " | ", block_material_y(iy), " => ",block_material_y(iy)
       else if ( &
         (ix<=sum(wrapper_%nx_blocks)) .and. (iy<=sum(wrapper_%ny_blocks)) ) then
         print *,ix,",", iy, ":", block_material_x(ix), " | ", block_material_y(iy), " => ------- ",block_material_y(iy)
       else
         print *,ix,",", iy, ":", block_material_x(ix), " | ", block_material_y(iy)
       end if

      end associate
    end associate

  end subroutine


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
