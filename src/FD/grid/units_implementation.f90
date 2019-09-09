!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
submodule(units_interface) units_implementation
    use assertions_interface, only : assertions,assert
    implicit none

    contains

        module procedure set_units
            !! define units exponents
            this%exponents_ = exponents
            this%system= system
            call this%mark_as_defined
        end procedure

        module procedure get_system
            system_of_units = this%system
        end procedure

        module procedure get_units
            exponents = this%exponents_
        end procedure

        module procedure assign_units
            lhs%exponents_ = rhs%exponents_
            lhs%system = rhs%system
        end procedure

        module procedure is_dimensionless
            nondimensional = all([this%exponents_==0]) .and. this%system==dimensionless
        end procedure

        module procedure has_length_units
            length_units = this%exponents_(m)==1 .and. all(this%exponents_([kg,sec,K])==0)
        end procedure

        module procedure has_mass_units
            mass_units = this%exponents_(kg)==1 .and. all(this%exponents_([m,sec,K])==0)
        end procedure

        module procedure has_time_units
            time_units = this%exponents_(sec)==1 .and. all(this%exponents_([m,kg,K])==0)
        end procedure

        module procedure has_temperature_units
            temperature_units = this%exponents_(K)==1 .and. all(this%exponents_([m,kg,sec])==0)
        end procedure

        module procedure has_velocity_units
            velocity_units = all(this%exponents_([m,sec])==[1,-1]) .and. all(this%exponents_([kg,K])==0)
        end procedure

        module procedure has_energy_units
            energy_units = all(this%exponents_([kg,m,sec])==[1,2,-2]) .and. this%exponents_(K)==0
        end procedure

        module procedure has_density_units
            density_units = all(this%exponents_([kg,m])==[1,-3]) .and. all(this%exponents_([sec,K])==0)
        end procedure

        module procedure has_specific_energy_units
            specific_energy_units = all(this%exponents_([m,sec])==[2,-2]) .and. all(this%exponents_([kg,K])==0)
        end procedure

        module procedure has_stress_units
            stress_units = all(this%exponents_([kg,m,sec])==[1,-1,-2]) .and. this%exponents_(K)==0
        end procedure

        module procedure has_power_units
            power_units = all(this%exponents_([kg,m,sec])==[1,2,-3]) .and. this%exponents_(K)==0
        end procedure

        module procedure integer_power
            this_raised%system = this%system
            this_raised%exponents_ = exponent_*this%exponents_
        end procedure

        module procedure real_power
            if (assertions) call assert(this%is_dimensionless(), &
                &   "units%real_power: an entity raised to a real power must be dimensionless")
            !! Require dimensionless operand
        end procedure

        module procedure negate
            negative_this%exponents_ = this%exponents_
            negative_this%system= this%system
        end procedure

        module procedure add
            if (assertions) then
                !! Require consistent operand units
                associate(preconditions => [lhs%system==rhs%system, lhs%exponents_==rhs%exponents_] )
                    call assert( all(preconditions), "units%add: consistent operands units")
                end associate
            end if
            total%exponents_  = lhs%exponents_
            total%system = lhs%system
        end procedure

        module procedure subtract
            if (assertions) then
                !! Require consistent operand units
                associate(preconditions => [lhs%system==rhs%system, lhs%exponents_==rhs%exponents_] )
                    call assert( all(preconditions), "units%subtract: consistent operand units")
                end associate
            end if
            difference%exponents_  = lhs%exponents_
            difference%system = lhs%system
        end procedure

        module procedure multiply

            if (assertions) call assert( lhs%system==rhs%system, "units%multiply: consistent operand units" )

            product_%exponents_  = lhs%exponents_ + rhs%exponents_
            product_%system = lhs%system
        end procedure

        module procedure divide

            if (assertions) call assert( numerator%system==denominator%system, "units%divide: consistent operand units" )

            ratio%exponents_  = numerator%exponents_ - denominator%exponents_
            ratio%system = merge(numerator%system,dimensionless,any(ratio%exponents_/=0))
        end procedure

end submodule units_implementation
