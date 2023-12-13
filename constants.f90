module constants
  implicit none
    !======================================================== 
  ! Astrophysical constants
  real, parameter :: boltz=1.38e-16
  real, parameter :: mh   = 1.67e-24
  real, parameter :: PC   = 3.085677588e+18   ! parsec(cm)
  real, parameter :: AU   = 1.495978707e+13     !unidad astronomica (cm)
  real, parameter :: MSUN = 1.988920000e+33   ! masa solar (g)
  real, parameter :: KYR  = 3.155673600e+10   ! Mil años en segundos (s)
  real, parameter :: AMU  = 1.660538782e-24   ! Unidad de masa atómica (g)
  real, parameter :: YR   = 3.155673600e+7      ! Año sideral terrestre (s)
  real, parameter :: PI = 4.D0*DATAN(1.D0)
  real, parameter :: KPS = 1.E5
end module constants
