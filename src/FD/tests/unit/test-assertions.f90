!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program main
  use assertions_interface, only : assert
  implicit none
  call assert(assertion=(1==1),description="scalar tautology")
  call assert([(1==1),(2>0)] ,["vector tautology 1","vector tautology 2"])
  call assert([.true.,.true.,.true.],"group tautology")
  print *,"Test passed."
end program
