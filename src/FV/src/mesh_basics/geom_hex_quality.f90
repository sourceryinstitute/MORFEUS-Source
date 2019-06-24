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
! $Id: geom_hex_quality.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description: Returns the quality of a hex, assuming CGNS order, ranging from 0 to 1
!
MODULE PROCEDURE geom_hex_quality

    USE class_psblas
    USE class_vector
    USE class_vertex

    IMPLICIT NONE
    !
    INTEGER :: iv
    REAL(psb_dpk_) :: sum_edge_sq

    sum_edge_sq = 0.d0

    DO iv = 2,4
        sum_edge_sq = sum_edge_sq + mag(verts(iv) - verts(iv-1) )**2
    ENDDO
    sum_edge_sq = sum_edge_sq + mag(verts(1) - verts(4) )**2

    DO iv = 6,8
        sum_edge_sq = sum_edge_sq + mag(verts(iv) - verts(iv-1) )**2
    ENDDO
    sum_edge_sq = sum_edge_sq + mag(verts(5) - verts(8) )**2

    DO iv = 5,8
        sum_edge_sq = sum_edge_sq + mag(verts(iv) - verts(iv-4) )**2
    ENDDO


    geom_hex_quality = SIGN((vol*vol)**(1.0d0/3.d0) * 12.0d0  / &
        &               sum_edge_sq, vol) ! signed quality, must handle neg. volume

END PROCEDURE geom_hex_quality
