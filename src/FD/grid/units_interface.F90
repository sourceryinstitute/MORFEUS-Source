!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
module units_interface
  !! author: Damian Rouson
  !! date: 9/9/2019
  !!
  !! Define SI and British units of measurement and associated arithmetic operators
  use object_interface, only : object
  implicit none

  private
  public :: units
  public :: K, m, kg, sec
  public :: R, ft, lbm
  public :: dimensionless, SI, British
  public :: num_fundamental
  public :: units_system_names
  public :: SI_units_names, British_units_names

  enum, bind(C)
    !! Enumerate the fundamental units for dimensional units quantities
    !! (meters, kilograms, seconds, and degrees Kelvin)
    enumerator :: K=1, m, kg, sec
    enumerator :: R=1, ft, lbm
    enumerator :: dimensionless=0, SI, British
  end enum

  integer, parameter :: num_fundamental=4, num_systems=2
  character(len=*), parameter :: units_system_names(num_systems)=[character(len=len("British")) :: "SI","British"]
  character(len=*), parameter :: SI_units_names(num_fundamental)=[character(len=len("sec")) :: "K", "m", "kg", "sec"]
  character(len=*), parameter :: British_units_names(num_fundamental)=[character(len=len("sec")) :: "R", "ft", "lbm", "sec"]

  type, extends(object) :: units
    !! Morfeus universal base type for all unitss
    private
    integer :: exponents_(num_fundamental)=dimensionless  !! Store the exponents for fundamental units
    integer :: system=dimensionless  !! Default to SI units
    character(len=:), allocatable :: description
  contains
    procedure :: set_units
    procedure :: get_units
    procedure :: get_system
    procedure :: is_dimensionless
    procedure :: has_length_units
    procedure :: has_mass_units
    procedure :: has_time_units
    procedure :: has_temperature_units
    procedure :: has_velocity_units
    procedure :: has_energy_units
    procedure :: has_density_units
    procedure :: has_specific_energy_units
    procedure :: has_stress_units
    procedure :: has_power_units
    procedure :: add
    procedure :: multiply
    procedure :: divide
    procedure :: subtract,negate
    procedure :: integer_power
    procedure :: real_power
    procedure :: assign_units
#ifndef FORD
    generic :: operator(+)=>add
    generic :: operator(*)=>multiply
    generic :: operator(/)=>divide
    generic :: operator(-)=>subtract,negate
    generic :: operator(**)=>integer_power,real_power
    generic :: assignment(=)=>assign_units
#endif
  end type

  interface

    pure module subroutine set_units(this,exponents,system)
      !! define units
        implicit none
        class(units), intent(inout) :: this
        integer, intent(in) :: exponents(num_fundamental)
        integer, intent(in) :: system
    end subroutine

    pure module subroutine assign_units(lhs,rhs)
      !! copy units information
        implicit none
        class(units), intent(inout) :: lhs
        class(units), intent(in) :: rhs
    end subroutine

#ifndef HAVE_ERROR_STOP_IN_PURE
    impure &
#endif
    elemental module function integer_power(this,exponent_) result(this_raised)
      !! result has units of the opearand raised to the power "exponent_"
        implicit none
        class(units), intent(in) :: this
        integer, intent(in) :: exponent_
        type(units) :: this_raised
    end function

    module function get_units(this) result(exponents)
      !! result holds the exponents of each unit in the argument (most useful when the actual argument is an expression)
        implicit none
        class(units), intent(in) :: this
        integer :: exponents(num_fundamental)
    end function

#ifndef HAVE_ERROR_STOP_IN_PURE
        impure &
#endif
    elemental module function get_system(this) result(system_of_units)
      !! result is enumerated value designating units system
        implicit none
        class(units), intent(in) :: this
        integer :: system_of_units
    end function


#ifndef HAVE_ERROR_STOP_IN_PURE
    impure &
#endif
    elemental module function real_power(this,exponent_) result(this_raised)
      !! result is the units of the operand raised to the power "exponent_"; includes check that operand is dimensionless
        implicit none
        class(units), intent(in) :: this
        real, intent(in) :: exponent_
        type(units) :: this_raised
    end function

#ifndef HAVE_ERROR_STOP_IN_PURE
    impure &
#endif
    elemental module function add(lhs,rhs) result(total)
      !! result is the units of the sum of two dimensional quantities; includes operand consistency check
        implicit none
        class(units), intent(in) :: lhs,rhs
        type(units) :: total
    end function

#ifndef HAVE_ERROR_STOP_IN_PURE
    impure &
#endif
    elemental module function subtract(lhs,rhs) result(difference)
      !! result is the units of the difference of two dimensional quantities; includes operand consistency check
        implicit none
        class(units), intent(in) :: lhs,rhs
        type(units) :: difference
    end function

    elemental module function negate(this) result(negative_this)
      !! result is the units of the negative of a dimensional quantities
        implicit none
        class(units), intent(in) :: this
        type(units) :: negative_this
    end function

    elemental module function multiply(lhs,rhs) result(product_)
      !! result is the units of the product of two dimensional quantities; includes units-system consistency check
        implicit none
        class(units), intent(in) :: lhs,rhs
        type(units) :: product_
    end function

    elemental module function divide(numerator,denominator) result(ratio)
      !! result is the units of the ratio of two dimensional quantities; includes units-sysetm consistency check
      implicit none
      class(units), intent(in) :: numerator,denominator
      type(units) :: ratio
    end function

    elemental module function is_dimensionless(this) result(nondimensional)
      !! Return true if all units exponents are zero; false otherwise.
      implicit none
      class(units), intent(in) :: this
      logical :: nondimensional
    end function

    elemental module function has_length_units(this) result(length_units)
      !! Return true if units match meters (m)
      implicit none
      class(units), intent(in) :: this
      logical :: length_units
    end function

    elemental module function has_mass_units(this) result(mass_units)
      !! Return true if units match kilograms (kg)
      implicit none
      class(units), intent(in) :: this
      logical :: mass_units
    end function

    elemental module function has_time_units(this) result(time_units)
      !! Return true if units match seconds (s)
      implicit none
      class(units), intent(in) :: this
      logical :: time_units
    end function

    elemental module function has_temperature_units(this) result(temperature_units)
      !! Return true if units match degrees Kelvin (K)
      implicit none
      class(units), intent(in) :: this
      logical :: temperature_units
    end function

    elemental module function has_velocity_units(this) result(velocity_units)
      !! Return true if units match meters/second^2 (m/s^2)
      implicit none
      class(units), intent(in) :: this
      logical :: velocity_units
    end function

    elemental module function has_energy_units(this) result(energy_units)
      !! Return true if units match joules (J)
      implicit none
      class(units), intent(in) :: this
      logical :: energy_units
    end function

    elemental module function has_density_units(this) result(density_units)
      !! Return true if units match kilograms (kg/m^3)
      implicit none
      class(units), intent(in) :: this
      logical :: density_units
    end function

    elemental module function has_specific_energy_units(this) result(specific_energy_units)
      !! Return true if units match Joules per kilogram (J/kg)
      implicit none
      class(units), intent(in) :: this
      logical :: specific_energy_units
    end function

    elemental module function has_stress_units(this) result(stress_units)
      !! Return true if units match Newtons per square meter (N/m^2)
      implicit none
      class(units), intent(in) :: this
      logical :: stress_units
    end function

    elemental module function has_power_units(this) result(power_units)
      !! Return true if units match Watts (W)
      implicit none
      class(units), intent(in) :: this
      logical :: power_units
    end function

  end interface

end module
