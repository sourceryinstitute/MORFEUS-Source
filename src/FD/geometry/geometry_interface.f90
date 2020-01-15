!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
module geometry_interface
  !! author: Damian Rouson
  !! date: 8/16/2019
  !! summary: abstract representation of problem geometry
  implicit none

  private
  public  :: geometry

  type, abstract :: geometry
    !! Abstract representation of the problem geometry in a problem-independent fashion;
    !! defer case-specific steps to child classes
    private
  contains
    procedure build
    procedure(set_json_file), deferred :: set_grid_specification
    procedure(set_metadata), deferred :: set_block_metadata
  end type


  interface

    module subroutine build(this, grid_description_file)
      !! template method pattern for building a geometry description from a json_file
      use json_module, only : json_file
      implicit none
      class(geometry), intent(out) :: this
      character(len=*), intent(in) :: grid_description_file
    end subroutine

  end interface

  abstract interface

    subroutine set_json_file(this, grid_description_file)
      !! define json_file problem description
      import geometry
      implicit none
      class(geometry), intent(out) :: this
      character(len=*), intent(in) :: grid_description_file
    end subroutine

    subroutine set_metadata(this)
      !! read grid metadata from json_file stored in child class
      import geometry
      implicit none
      class(geometry), intent(inout) :: this
    end subroutine

  end interface

end module geometry_interface
