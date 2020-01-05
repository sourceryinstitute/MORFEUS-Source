!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program main
  !! author: Damian Rouson
  !!
  !! Verify the relevant use cases for Fortran 2008/2018 intrinsic procedures
  !! that require emulation in the absence of compiler support.

#if ! (defined HAVE_FINDLOC) && ! (defined HAVE_COLLECTIVE_SUBROUTINES)
  use emulated_intrinsics_interface, only : co_sum, co_broadcast, findloc
#else
#ifndef HAVE_FINDLOC
  use emulated_intrinsics_interface, only : findloc
#endif
#ifndef HAVE_COLLECTIVE_SUBROUTINES
  use emulated_intrinsics_interface, only : co_sum, co_broadcast
#endif
#endif

  use assertions_interface, only : assert
  implicit none

  associate( me=>this_image(), ni=>num_images() )


    test_collective_broadcast: block
      integer ::i, messenger, message

      call assert(ni>1,"test-emulated-intrinsics: at least 2 images required")

      message=333
      messenger=1
      if (me==messenger) i=message
      call co_broadcast(i,source_image=messenger)
      call assert(i==message,"integer message broadcast from default binary tree root")

      message=666
      messenger=2
      if (me==messenger) i=message
      call co_broadcast(i,source_image=messenger)
      call assert(i==message,"integer message broadcast non-default binary tree root")

    end block test_collective_broadcast

    test_collective_sum: block
      integer ::i, image
      integer, parameter :: recipient=2

      image = me
      call co_sum(image,result_image=recipient)
      if (me==recipient) call assert(image==sum([(i,i=1,ni)]),"collective integer sum reduction accumulated on recipient image ")

      image = me
      call co_sum(image)
      call assert(image==sum([(i,i=1,ni)]),"collective integer sum accumulated on all images")

    end block test_collective_sum

    if (me==1) then
      test_findloc: block
        integer i,j
        character(len=*), parameter, dimension(*) :: names = &
          [character(len=len("dressing"))::"absent", "core", "dressing", "gap", "wrap", "ring", "foil", "skin", "fabric" ]

        character(len=*), parameter, dimension(*) :: empty = [character(len=len("dressing"))::]

        integer, parameter, dimension(*) :: zero_sized = [ integer :: ]

        enum, bind(C)
          enumerator &
          absent, core, dressing, gap, wrap, ring, foil, skin, fabric
        end enum

        associate( unused=>absent) !! eliminate unused variable warning
        end associate

        associate( nothing_to_search => findloc(zero_sized, value=99, dim=1, back=.true.) )
          call assert(nothing_to_search==0, "findloc: handles zero-sized array")
        end associate

        associate( non_existent_value => findloc( [1,2,3,4], value=99, dim=1, back=.true. ) )
          call assert(non_existent_value==0, "findloc: handles non-existent value")
        end associate

        associate( x_blocks => [5,0,1,2,1,0,1,2] )
        associate( all_layers => [core,dressing,gap,wrap,ring,foil,skin,fabric])
        associate( full_delineation => [( [(all_layers(i),j=1,x_blocks(i))],i=1,size(x_blocks) )] )
        associate( outermost_core => findloc(full_delineation, core, dim=1, back=.true.) )
          call assert( outermost_core==x_blocks(core), "findloc: right outermost core layer" )
        end associate; end associate; end associate; end associate

        associate( non_existent_value=> findloc(names, value="dreffing", dim=1, back=.true.) )
          call assert(non_existent_value==0, "findloc: handles non-existent string")
        end associate

        associate( dressing_location => findloc(names, value="dressing", dim=1, back=.true.) )
          call assert(dressing_location==3, "findloc: finds string location", dressing_location)
        end associate

        associate( empty_array => findloc(empty, value="dressing", dim=1, back=.true.) )
          call assert(empty_array==0, "findloc: handles empty character array")
        end associate

        associate( dressing_location => findloc(names, value="dressing", dim=1, back=.false.) )
          call assert(dressing_location==3, "findloc: finds string location from front", dressing_location)
        end associate

        associate( empty_array => findloc(empty, value="dressing", dim=1, back=.false.) )
          call assert(empty_array==0, "findloc: handles empty character array from front")
        end associate

      end block test_findloc
    end if

    sync all
    if (me==1) print *,"Test passed."

    end associate

end program main
