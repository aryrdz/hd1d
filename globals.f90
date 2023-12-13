!=======================================================================
!   This module contains global variables                                     
module globals
  use constants
  implicit none
  !========================================================================
  !                                                                           
  integer, parameter :: nx=3000
  !   This is the number of equation
  integer, parameter :: neq=3                             
  !   Here we set the extent of X and calculate $\Delta x$
  real, parameter :: xmax=1000.*AU
  real, parameter :: dx=xmax/real(nx)
  ! The simulation times
  real, parameter :: tmax= 1*YR             ! maximumn integration time
  real, parameter :: dtprint= 0.1*YR          ! interval between outputs
  ! simulation constants
  real, parameter :: gamma=5./3.
  real, parameter :: mu=1.4
  real, parameter :: eta = 0.0001
  ! Courant number                                                            
  real, parameter :: Co=0.5
  ! symmetry term
  real, parameter :: alpha=2.
  !                 
  ! Tolerance stuff
  real, parameter :: den_floor  = 1.e-15*(mu*mh) 
  real, parameter :: temp_floor = 100.
  !   This is a vector that contains u(x) 
  real,dimension(neq,0:nx+1) :: u,f
end module globals
