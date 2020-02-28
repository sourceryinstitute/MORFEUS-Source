!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
submodule(emulated_intrinsics_interface) emulated_intrinsics_implementation
  use assertions_interface, only : assert
  implicit none

contains

#ifndef HAVE_COLLECTIVE_SUBROUTINES
  module procedure co_sum_integer
    !! Binary tree collective sum reduction

    integer, save :: total[*]
    integer, parameter :: root_node=1, default_root_image=1
    integer even_child_image, odd_child_image, parent_image, my_image
    integer even_child_node , odd_child_node , parent_node , my_node, default_parent_image
    integer, allocatable :: relatives(:)

    my_image=this_image()

    if (.not. present(result_image)) then
      my_node = my_image
    else
      if (my_image==result_image) then
        my_node = root_node
      else
        my_node = merge(result_image, my_image, my_image==default_root_image)
      end if
    end if

    parent_node = my_node/2
    default_parent_image = my_image/2

    if (.not. present(result_image)) then
      parent_image = parent_node
    else
      if (parent_node==root_node) then
        parent_image = result_image
      else if (default_parent_image==result_image) then
        parent_image = default_root_image
      else
        parent_image = parent_node
      end if
    end if

    total = a

    even_child_node = 2*my_node
    if (.not. present(result_image)) then
      even_child_image = even_child_node
    else
      even_child_image = merge(even_child_node, default_root_image, even_child_node/=result_image)
    end if
    if ( exists(even_child_node) ) then
      sync images(even_child_image)
      total = total + total[even_child_image]
    end if

    odd_child_node = 2*my_node+1
    if (.not. present(result_image)) then
      odd_child_image = odd_child_node
    else
      odd_child_image = merge(odd_child_node, default_root_image, odd_child_node/=result_image)
    end if
    if ( exists(odd_child_node) )  then
      sync images(odd_child_image)
      total = total + total[odd_child_image]
    end if

    if (exists(parent_image)) sync images(parent_image)

    if (present(result_image)) then
      if (my_image==result_image) a = total
    else
      call co_broadcast(total,source_image=default_root_image)
      a = total
    end if

    if (present(errmsg)) errmsg=""
    if (present(stat)) stat=0

    relatives = [parent_image,even_child_image,odd_child_image]
    relatives = pack(relatives,exists(relatives))
    sync images(relatives)

  end procedure

  module procedure co_broadcast_integer
    !! Binary tree collective broadcast

    integer, save :: message[*]
    integer, parameter :: root_node=1, default_root_image=1
    integer even_child_image, odd_child_image, parent_image, my_image
    integer even_child_node , odd_child_node , parent_node , my_node, default_parent_image
    integer, allocatable :: relatives(:)

    my_image=this_image()

    if (my_image==source_image) then
      my_node = root_node
    else
      my_node = merge(source_image, my_image, my_image==default_root_image)
    end if

    parent_node = my_node/2
    default_parent_image = my_image/2

    if (parent_node==root_node) then
      parent_image = source_image
    else if (default_parent_image==source_image) then
      parent_image = default_root_image
    else
      parent_image = parent_node
    end if

    if (exists(parent_image)) then
      sync images(parent_image)
      message = message[parent_image]
    else
      message = a
        !! define the message for broadcasting down the tree when the orphaned root image falls through this branch
    end if

    even_child_node = 2*my_node
    even_child_image = merge(even_child_node, default_root_image, even_child_node/=source_image)
    if ( exists(even_child_node) ) sync images(even_child_image)

    odd_child_node = 2*my_node+1
    odd_child_image = merge(odd_child_node, default_root_image, odd_child_node/=source_image)
    if ( exists(odd_child_node) )  sync images(odd_child_image)

    a = message

    relatives = [parent_image,even_child_image,odd_child_image]
    relatives = pack(relatives,exists(relatives))
    sync images(relatives)

    if (present(errmsg)) errmsg=""
    if (present(stat)) stat=0

  end procedure

  elemental function exists(image) result(image_exists)
    !! Result true if image number is within the closed range [1,num_images()]
    integer, intent(in) :: image
    logical image_exists
    image_exists = (image>0 .and. image<=num_images())
  end function
#endif

#ifndef HAVE_FINDLOC
  module procedure findloc_integer_dim1

    if ( .not. present(back)) then

      location = minloc(array, dim, array == value)

    else if (back .eqv. .false.) then

      location = minloc(array, dim, array == value)

    else ! back is present and .true. so work around GCC 8 lack of support for the "back" argument
      block
        integer, parameter :: loop_increment=-1, base_index=1
        integer index_

        call assert(dim==1,"findloc_integer_dim1: unsupported use case")

        associate( lower_bound=>lbound(array,dim) )
          do index_=ubound(array,dim), lower_bound, loop_increment
            if (array(index_)==value) then
              location = index_ - lower_bound + base_index
              return
            end if
          end do
        end associate
        location=0
      end block
    end if

  end procedure

  module procedure findloc_logical_dim1
    integer, parameter :: loop_increment=-1, base_index=1
    integer index_

    call assert(back .and. dim==1,"findloc_logical_dim1_backtrue: unsupported use case")

    associate( lower_bound=>lbound(array,dim) )
      do index_=ubound(array,dim), lower_bound, loop_increment
        if (array(index_) .eqv. value) then
          location = index_ - lower_bound + base_index
          return
        end if
      end do
    end associate
    location=0

  end procedure

  module procedure findloc_character_dim1
    integer, parameter :: base_index=1
    integer index_, loop_increment, start, finish

    call assert(dim==1,"findloc_character_dim1: unsupported use case")

    associate( lower_bound=>lbound(array,dim) )
      select case(back)
        case(.true.)
          start = ubound(array,dim)
          finish = lower_bound
          loop_increment=-1
        case(.false.)
          start = lower_bound
          finish = ubound(array,dim)
          loop_increment=1
      end select

      do index_=start, finish, loop_increment
        if (array(index_)==value) then
          location = index_ - lower_bound + base_index
          return
        end if
      end do
    end associate
    location=0

  end procedure
#endif

end submodule emulated_intrinsics_implementation
