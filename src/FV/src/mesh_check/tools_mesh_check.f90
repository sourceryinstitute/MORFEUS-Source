!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
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
! $Id: tools_mesh_check.f90 2469 2007-10-08 10:34:43Z sfilippo $
!
! Description:
!    To be added...
!
MODULE tools_mesh_check
    USE class_psblas, ONLY : psb_dpk_
    USE class_mesh, ONLY : mesh
      !! An Intel 18.0.5 bug precludes putting this in the interface bodies
    USE class_scalar_field, ONLY : scalar_field
      !! An Intel 18.0.5 bug precludes putting this in the interface bodies
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: check_tet_quality
    PUBLIC :: check_mesh_quality
    PUBLIC :: check_vertex_quality

    INTERFACE

        MODULE SUBROUTINE check_tet_quality(msh,ic,fmt)
            IMPLICIT NONE
            TYPE(mesh), INTENT(IN) :: msh
            INTEGER, INTENT(IN)    :: ic
            INTEGER, INTENT(IN), OPTIONAL :: fmt
        END SUBROUTINE check_tet_quality

       MODULE SUBROUTINE check_mesh_quality(msh,fquality,tol,bad_cells,quiet)
            IMPLICIT NONE
            TYPE(mesh), INTENT(IN) :: msh
            TYPE(scalar_field), INTENT(INOUT) :: fquality
              !! Note: fquality was changed from INTENT(IN) to INTENT(INOUT) b/c the value can be changed
              !!       in the implementation. IP 6/5/2019
            REAL(psb_dpk_), INTENT(IN),OPTIONAL :: tol
            INTEGER, ALLOCATABLE, INTENT(OUT), OPTIONAL :: bad_cells(:)
            LOGICAL,OPTIONAL :: quiet
        END SUBROUTINE check_mesh_quality

        MODULE FUNCTION check_vertex_quality(msh,ic2v) RESULT (quality)
            IMPLICIT NONE
            TYPE(mesh), INTENT(IN) :: msh
            INTEGER :: ic2v(:)
            REAL(psb_dpk_) :: quality
        END FUNCTION check_vertex_quality

    END INTERFACE

END MODULE tools_mesh_check
