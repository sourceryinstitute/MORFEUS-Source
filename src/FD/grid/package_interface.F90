!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
module package_interface
  !! author: Damian Rouson and Karla Morris
  !! date: 1/2/2019
  !! Encapsulate halo data to be communicated across structured_grid block boundaries
  use kind_parameters, only : r8k
  implicit none

  private
  public :: package, null_sender_id

  integer, parameter :: null_sender_id=-1

  type package
    !! basic transmission data. extend this type to add coordinate-specific data
    private
    integer :: sender_block_id
    integer :: step
    real(r8k), allocatable, dimension(:,:,:,:) :: s_flux_normal
      !! surface-normal scalar flux components: shape = [Nx, Ny, Nz, scalar_field ID], where one of Nx|Ny|Nz is 1
    real(r8k), allocatable, dimension(:,:,:) :: positions
      !! flux locations: shape = [Nx, Ny, Nz], where one of Nx|Ny|Nz is 1
  contains
    procedure get_sender_block_id
    procedure set_sender_block_id
    procedure set_step
    procedure set_normal_scalar_fluxes
    procedure sender_block_id_null
    procedure copy
    generic :: assignment(=) => copy
  end type

  interface

    module function get_sender_block_id(this) result(this_sender_block_id)
      !! result is sender_block_id for this package
      implicit none
      class(package), intent(in) :: this
      integer :: this_sender_block_id
    end function

    elemental module subroutine set_sender_block_id(this, sender_block_id)
      !! set sender_bock_id in order to verify sender identity
      implicit none
      class(package), intent(inout) :: this
      integer, intent(in) :: sender_block_id
    end subroutine

    elemental module subroutine set_step(this, step)
      !! set time step in order to verify time step on package data
      implicit none
      class(package), intent(inout) :: this
      integer, intent(in) :: step
    end subroutine

    module subroutine set_normal_scalar_fluxes(this, s_flux_normal, positions)
      !! set datum to be communicated across structured_grid block internal surfaces
      implicit none
      class(package), intent(inout) :: this
      real(r8k), intent(in), dimension(:,:,:,:) :: s_flux_normal
        !! surface-normal scalar flux components: shape = [Nx, Ny, Nz, scalar_field ID], where one of Nx|Ny|Nz is 1
      real(r8k), intent(in), dimension(:,:,:) :: positions
        !! flux locations: shape = [Nx, Ny, Nz], where one of Nx|Ny|Nz is 1
    end subroutine

    elemental module function sender_block_id_null(this) result(is_null)
      !! set sender_block_id to invalid value for external boundaries (no block sends halo data to a boundary)
      implicit none
      class(package), intent(in) :: this
      logical is_null
    end function

    module subroutine copy(this, rhs)
      !! copy rhs package components into this package
      class(package), intent(inout) :: this
      type(package), intent(in) :: rhs
    end subroutine

  end interface

end module package_interface
