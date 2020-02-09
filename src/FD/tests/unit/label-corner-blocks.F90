!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
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

  enum, bind(C)
    enumerator &
      absent, core, dressing, gap, wrap, ring, foil, skin, fabric, chamber
  end enum

  character(len=len('wrap')), parameter, dimension(absent:chamber) :: layer_name = &
    [character(len=len('wrap')) :: 'absent','core','dressing','gap','wrap','ring',&
    'foil','skin','fabric','chamber']

  associate( x_blocks => [5,0,1,2,1,0,1,2] )
  associate( y_blocks => x_blocks )
  associate( ny_core => y_blocks(1))
  associate(all_layers=>[core, dressing, gap,wrap, ring, foil, skin, fabric])

  if (size(x_blocks)/=size(all_layers)) error stop "size mismatch"

  associate( full_delineation => [( [(all_layers(i),j=1,x_blocks(i))],i=1,size(x_blocks) )] )

    print *,"x_blocks:        ",x_blocks
    print *,"all_layers:      ",all_layers
    print *

    print *,"full_delineation:       ",full_delineation

    associate( chamber_delineation =>merge(chamber,full_delineation,full_delineation==core) )
      print *,"chamber_delineation:     ",chamber_delineation
    end associate

    print *

    block
      use iso_c_binding, only : c_int
      integer(c_int), allocatable :: progressive_delineation(:)
      integer(c_int) layer, layer_block
      integer iy, iz
      integer, parameter :: nz_all=0, nz_chamber=1

      progressive_delineation = pack(full_delineation,full_delineation/=core)

      associate( ny_half=>size(full_delineation) )
      do concurrent(iz=1:nz_all+nz_chamber)
        iy=ny_core
        do layer=dressing,fabric-1
          associate(innermost_layer=>progressive_delineation(1))
            if (innermost_layer==layer) then
              associate( end_of_innermost => findloc(progressive_delineation, innermost_layer, dim=1, back=.true.) )
              associate( next_layer => progressive_delineation(end_of_innermost+1) )
                do layer_block=1,y_blocks(layer)
                  iy=iy+1
                  print *,"iy | iz | progressive_delineation (",layer_name(layer),"):",iy," | ",iz," | ",progressive_delineation
                  progressive_delineation = merge(progressive_delineation,next_layer,progressive_delineation/=layer)
                end do
              end associate; end associate
            else
              print *,"===> No ",layer_name(layer)," <==="
            end if
          end associate
        end do
        do layer_block=1,y_blocks(fabric)
          iy=iy+1
          print *,"iy | iz | progressive_delineation (",layer_name(fabric),"):",iy," | ",iz," | ",progressive_delineation
        end do
      end do

      end associate

    end block

  end associate; end associate; end associate; end associate; end associate


#ifdef __GNUC__
#if __GNUC__ >= 9
#define HAVE_FINDLOC
#endif
#endif

#ifndef HAVE_FINDLOC
contains

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
