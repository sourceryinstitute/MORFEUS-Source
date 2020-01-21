!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
program test_units
  !! author:  Damian Rouson
  !!
  !! Verify the units-handling capabilities of the units class

  use units_interface, only : units, num_fundamental, K, m, kg, sec, SI
  use assertions_interface, only : assert
  implicit none

  call dimensionless_set_units
  call mass_units_setter_and_getter
  call velocity_energy_operations
  call length_operations
  call time_units_checks
  call temperature_units_checks
  call density_units_checks
  call stress_units_checks
  call power_units_checks

  sync all
  if (this_image()==1) print *,"Test passed."

contains

  pure subroutine dimensionless_set_units
    !! Verify default initialization and user-defined arithmetic for dimensionless quantities

    type(units) :: reynolds, mach

    call assert( reynolds%is_dimensionless(), "test-units: default initialization is dimensionless" )
      !! Verify dimensionless default initialization

    mach = reynolds**2.5
    call assert( mach%is_dimensionless(), "test-units: real power of dimensionless quantity" )
      !! Verify real power of dimensionless quantity

  end subroutine

  pure subroutine mass_units_setter_and_getter
    !! Verify setting and getting of mass units
    type(units) mass

    integer mass_units(num_fundamental)
    mass_units = 0      !! initially dimensionless
    mass_units(kg) = 1  !! set dimensions of mass

    call mass%set_units(exponents=mass_units,system=SI)
    call assert( mass%has_mass_units() .and. mass%get_system()==SI, "test-units: mass%has_mass_units() in SI" )

  end subroutine

  pure subroutine velocity_energy_operations
    !! Verify velocity, energy, and specific energy units, defined assignment, and  defined product and exponentiation operators

    type(units) mass, velocity, specific_energy, energy
    integer, dimension(num_fundamental) :: velocity_units,  mass_units
    velocity_units = 0
    velocity_units([m,sec])=[1,-1]
    mass_units = 0      !! initially dimensionless
    mass_units(kg) = 1  !! set dimensions of mass

    call velocity%set_units(exponents=velocity_units,system=SI)

    call assert( velocity%user_defined(), "test-units: velocity%user_defined()")
    call assert( velocity%has_velocity_units() .and. velocity%get_system()==SI,"test-units: velocity%has_velocity_units()")
      !! Verify setting and getting of velocity units

    specific_energy = velocity**2
    call assert(specific_energy%has_specific_energy_units() .and. specific_energy%get_system()==SI, &
      "test-units: specific_energy%has_specific_energy_units()")
      !! Verify integer exponentation of dimensional quantity

    call mass%set_units(exponents=mass_units,system=SI)
    energy = mass*specific_energy
    call assert( energy%has_energy_units() .and. energy%get_system()==SI,"test-units: energy%has_energy_units()")
      !! Verify multiplication of dimensional quantities

  end subroutine

  pure subroutine length_operations
    !! Verify setting and getting of length dimensions and addition, subtraction and negation operators
    type(units) height, width, depth, difference, negative, ratio
    integer length_units(num_fundamental)
    length_units=0
    length_units(m)=1

    call height%set_units(exponents=length_units,system=SI)
    call assert( height%has_length_units() .and. height%get_system()==SI,"test-units: height%has_length_units()")
      !! Verify setting and getting of length units

    width = height
    call assert( width%has_length_units() .and. width%get_system()==SI,"test-units: width%has_length_units()")

    depth = width + height
    call assert( depth%has_length_units() .and. depth%get_system()==SI,"test-units: depth%has_length_units()")
      !! Verify units of additive length

    difference = width - depth
    call assert( difference%has_length_units() .and. difference%get_system()==SI,"test-units: difference%has_length_units()")
      !! Verify units of differential length

    negative = -depth
    call assert( negative%has_length_units() .and. difference%get_system()==SI,"test-units: difference%has_length_units()")
      !! Verify units of negated length

    ratio = width/depth
    call assert( ratio%is_dimensionless(),"test-units: ratio%is_dimensionless()")
      !! Verify dimensionless ratio of like quantities

  end subroutine

  pure subroutine time_units_checks
    type(units) waiting
    integer time_units(num_fundamental)
    time_units=0
    time_units(sec)=1

    call waiting%set_units(exponents=time_units,system=SI)
    call assert( waiting%has_time_units() .and. waiting%get_system()==SI,"test-units: waiting%has_time_units()")
      !! Verify units of time

  end subroutine

  pure subroutine temperature_units_checks
    type(units) preheat
    integer temperature_units(num_fundamental)
    temperature_units=0
    temperature_units(K)=1

    call preheat%set_units(exponents=temperature_units,system=SI)
    call assert( preheat%has_temperature_units() .and. preheat%get_system()==SI,"test-units: preheat%has_temperature_units()")
      !! Verify units of temperature

  end subroutine

  pure subroutine density_units_checks
    type(units) mud
    integer density_units(num_fundamental)
    density_units=0
    density_units([kg,m])=[1,-3]

    call mud%set_units(exponents=density_units,system=SI)
    call assert( mud%has_density_units() .and. mud%get_system()==SI,"test-units: mud%has_density_units()")
      !! Verify units of density

  end subroutine

  pure subroutine stress_units_checks
    type(units) :: anxiety
    integer stress_units(num_fundamental)
    stress_units=0
    stress_units([kg,m,sec])=[1,-1,-2]

    call anxiety%set_units(exponents=stress_units,system=SI)
    call assert( anxiety%has_stress_units() .and. anxiety%get_system()==SI,"test-units: mud%has_stress_units()")
      !! Verify units of stress

  end subroutine

  pure subroutine power_units_checks
    type(units) :: solar
    integer power_units(num_fundamental)
    power_units=0
    power_units([kg,m,sec])=[1,2,-3]

    call solar%set_units(exponents=power_units,system=SI)
    call assert( solar%has_power_units() .and. solar%get_system()==SI,"test-units: mud%has_power_units()")
      !! Verify units of power

  end subroutine

end program test_units
