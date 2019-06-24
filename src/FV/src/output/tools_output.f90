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
! $Id: tools_output.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    To be added...
!
MODULE tools_output
    USE class_mesh
    USE class_iterating
    USE class_output
    USE class_cell, ONLY : cell
    USE class_psblas, ONLY : psb_dpk_
    USE class_scalar_field
    USE class_vector_field
    USE class_connectivity
    USE class_face
    USE class_vertex
    USE class_vector
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: write_field, write_mesh
    ! ----- Mesh & Field -----

    INTERFACE write_field
        MODULE PROCEDURE :: write_scalar_field
        MODULE PROCEDURE :: write_vector_field
    END INTERFACE

    INTERFACE

        MODULE SUBROUTINE write_scalar_field(fld,field,out,iter)
            IMPLICIT NONE
            TYPE(scalar_field), INTENT(IN) :: fld
            CHARACTER(len=*),   INTENT(IN) :: field
            TYPE(output),       INTENT(INOUT) :: out
            TYPE(iterating),    INTENT(IN), OPTIONAL :: iter
        END SUBROUTINE write_scalar_field

        MODULE SUBROUTINE write_vector_field(fld,field,out,iter)
            IMPLICIT NONE
            TYPE(vector_field), INTENT(IN) :: fld
            CHARACTER(len=*),   INTENT(IN) :: field
            TYPE(output),       INTENT(INOUT) :: out
            TYPE(iterating),    INTENT(IN), OPTIONAL :: iter
        END SUBROUTINE write_vector_field

        MODULE SUBROUTINE write_mesh(msh,out,iter)
            IMPLICIT NONE
            TYPE(mesh),   INTENT(IN) :: msh
            TYPE(output), INTENT(INOUT) :: out
            TYPE(iterating), INTENT(IN), OPTIONAL :: iter
        END SUBROUTINE write_mesh

    END INTERFACE

END MODULE tools_output
