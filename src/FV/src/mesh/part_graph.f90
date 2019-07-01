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
! $Id: part_graph.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    METIS partitioning
!
MODULE part_graph

    USE class_psblas

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: bld_part_graph, get_part_graph


    ! METIS
    INTEGER, ALLOCATABLE, SAVE :: part_glob(:)
    INTEGER, ALLOCATABLE, SAVE :: xadj(:), adjncy(:)
    INTEGER, ALLOCATABLE, SAVE :: perm(:), pinv(:)


    ! ----- Explicit Interfaces -----

    INTERFACE

        MODULE SUBROUTINE bld_part_graph(xadj_glob,adjncy_glob,ipart)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: xadj_glob(:), adjncy_glob(:)
            INTEGER, INTENT(IN) :: ipart
        END SUBROUTINE

        MODULE SUBROUTINE get_part_graph(part)
            IMPLICIT NONE
            INTEGER, ALLOCATABLE, INTENT(OUT) :: part(:)
        END SUBROUTINE

        SUBROUTINE METIS_PartGraphKway(n,xadj_glob,adjncy_glob,vwgt,adjwgt,&
            & wgtflag,numflag,nparts,options,edgecut,part)
            INTEGER, INTENT(IN):: n
            INTEGER, INTENT(IN), DIMENSION(*) :: xadj_glob, adjncy_glob
            INTEGER, INTENT(IN), DIMENSION(*) :: vwgt, adjwgt
            INTEGER, INTENT(IN) :: wgtflag, numflag,nparts
            INTEGER, INTENT(IN) :: options(*)
            INTEGER, INTENT(OUT) :: edgecut
            INTEGER, INTENT(OUT), DIMENSION(*) :: part
        END SUBROUTINE METIS_PartGraphKway

    END INTERFACE

END MODULE part_graph
