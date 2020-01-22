!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
submodule(package_interface) package_implementation
  !! author: Damian Rouson and Karla Morris
  !! date: 1/2/2019
  !! Encapsulate information and procedures for structured_grid block halo exchanges
  use assertions_interface, only: assert, assertions, max_errmsg_len
#ifndef HAVE_FINDLOC
  use emulated_intrinsics_interface, only : findloc
#endif
  implicit none

  integer, parameter :: success = 0

contains

  module procedure set_neighbor_block_id
    this%neighbor_block_id = neighbor_block_id
  end procedure

  module procedure set_step
    this%step = step
  end procedure

  module procedure set_surface_flux_positions
    if (assertions) call assert( count(shape(positions)==1) == 1 , "package%set_positions: count(shape(positions)==1) == 1")
    this%positions = positions
  end procedure

  module procedure set_num_scalars
    integer alloc_stat
    character(len=max_errmsg_len) error_message
    if (assertions) &
      call assert(.not.allocated(this%surface_normal_fluxes), "package%set_num_scalars: .not.allocated(this%surface_normal_fluxes)")
    allocate( this%surface_normal_fluxes(num_scalars), stat = alloc_stat, errmsg = error_message)
    if (assertions) call assert(alloc_stat==success, "package%set_num_scalars: allocate(this%surface_normal_fluxes)", error_message)
  end procedure

  module procedure set_normal_scalar_fluxes
    if (assertions) then
      call assert( allocated(this%positions), "package%normal_scalar_fluxes: allocated(positions)" )
      call assert(&
        findloc(shape(fluxes), dim=1, value=1, back=.true.) == findloc(shape(this%positions), dim=1, value=1, back=.true.), &
        "surfaces%set_normal_scalar_fluxes: matching surface flux/position normal direction")
      call assert( lbound(this%surface_normal_fluxes,1) <= scalar_id .and. scalar_id <= ubound(this%surface_normal_fluxes,1), &
        "package%set_normal_scalar_fluxes: lbound(surface_normal_fluxes) <= scalar_id <= ubound(surface_normal_fluxes)" )
    end if
    this%surface_normal_fluxes(scalar_id)%fluxes = fluxes
  end procedure

  module procedure get_neighbor_block_id
    this_neighbor_block_id = this%neighbor_block_id
  end procedure

  module procedure neighbor_block_id_null
    is_null = (this%neighbor_block_id == null_neighbor_id)
  end procedure

  module procedure copy

    if (assertions) then
      block
        integer i

        call assert( allocated(rhs%positions) .and. allocated(rhs%surface_normal_fluxes), &
          "package%copy: all([allocated(rhs%fluxes) .and. allocated(surface_normal_fluxes)])" )
        do concurrent (i = lbound(rhs%surface_normal_fluxes,1): ubound(rhs%surface_normal_fluxes,1) )
          call assert( allocated(rhs%surface_normal_fluxes(i)%fluxes), "package%copy: allocated(rhs%surface_normal_fluxes%fluxes)" )
        end do
      end block
    end if

    this%neighbor_block_id = rhs%neighbor_block_id
    this%step = rhs%step
    this%surface_normal_fluxes = rhs%surface_normal_fluxes
    this%positions = rhs%positions
  end procedure

end submodule package_implementation
