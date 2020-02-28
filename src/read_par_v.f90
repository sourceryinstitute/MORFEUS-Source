!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
!
!    NEMO - Numerical Engine (for) Multiphysics Operators
! Copyright (c) 2007, Stefano Toninel
!                     Gian Marco Bianchi  University of Bologna
!              David P. Schmidt    University of Massachusetts - Amherst
!              Salvatore Filippone University of Rome Tor Vergata
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification,
! are permitted provided that the following conditions are met:
!
!     1. Redistributions of source code must retain the above copyright notice,
!        this list of conditions and the following disclaimer.
!     2. Redistributions in binary form must reproduce the above copyright notice,
!        this list of conditions and the following disclaimer in the documentation
!        and/or other materials provided with the distribution.
!     3. Neither the name of the NEMO project nor the names of its contributors
!        may be used to endorse or promote products derived from this software
!        without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!---------------------------------------------------------------------------------
!
! $Id: $
!
! Description:
!    To be added...
!
SUBMODULE (tools_input) read_par_v_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE read_par_v
            USE class_psblas
            USE class_vector
            USE tools_input, ONLY : get_par, open_file, find_section
            IMPLICIT NONE
            !
            INTEGER :: icontxt, mypnum
            INTEGER :: inp
            REAL(psb_dpk_) :: x(3)

            icontxt = icontxt_()
            mypnum  = mypnum_()


            IF(mypnum == 0) THEN

                CALL open_file(input_file,inp)
                CALL find_section(sec,inp)

                read_par_v = get_par(inp,sec,par,default)
                x(1) = read_par_v%x_()
                x(2) = read_par_v%y_()
                x(3) = read_par_v%z_()

                CALL psb_bcast(icontxt,x)

                CLOSE(inp)
            ELSE
                CALL psb_bcast(icontxt,x)

                read_par_v = vector_(x(1),x(2),x(3))
            END IF

        END PROCEDURE read_par_v

END SUBMODULE read_par_v_implementation
