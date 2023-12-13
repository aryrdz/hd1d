!=======================================================================
!   This program solves the hydrodynamics equation 1D
!=======================================================================
!   main program
program hd_1d
  use globals
  implicit none
  ! declaration of some variables needed by the main program
  real            :: time, dt             !  t, $\Delta t$
  real            :: tprint               ! time of next output
  integer         :: itprint              ! number of current output
  
  ! This subroutine generates the initial conditions
  call initflow(time, tprint, itprint)
  !   main loop, iterate until maximum time is reached
  do while (time.lt.tmax)
     !

     ! output at tprint intervals
     if(time.ge.tprint) then
        print*,'itprint=',itprint, time,tmax,dt
        call output(itprint)
        tprint=tprint+dtprint
        itprint=itprint+1
     end if
          ! Obtain the Delta t allowed by the CFL criterium
     call timestep(dt)
     !
     ! Integrate u fom t to t+dt
     call tstep(dt,time)
     !
     ! time counter increases
     print*,'time [yr]=',time/yr,'dt [yr]=',dt/yr
     time=time+dt
     !
  end do
  !
  stop
end program hd_1d

!=======================================================================
! generates initial condition
subroutine initflow(time, tprint, itprint)
  use globals
  use constants
  implicit none
  real, intent(out) :: time, tprint
  integer, intent (out) :: itprint
  ! The initial condition imposed here for the interestellar medium
  ! 
  real, parameter :: n_ism = 1.0           ! Numeric density (cm^-3)
  real, parameter :: mu0_ism = 1.3         ! Masa por partícula, gas neutro (mu)
  real, parameter :: mui_ism = 0.61        ! Masa por partícula, gas ionizado (mu)
  real, parameter :: T_ism = 100           ! Temperature (K)
  real, parameter :: u_ism = 0.0           ! x-velocity (cm/s)
  
  ! And the kinetical conditions for the SNR
  
!  real, parameter :: RSN = 20.*dx        !Initial radius of the explosion (cm)
!  real, parameter :: ESN = 1.0d+31       !Explosion energy (erg)
!  real, parameter :: MSN = 10*MSUN     ! Ejected mass (g)
!  real, parameter :: SN_x = 0.           ! Center of the explosion x-coordinate(cm)
  !In terms of the primitives, for the ISM 
  real, parameter :: rho_ism = n_ism * mu0_ism * mh
  real, parameter :: P_ism = boltz*n_ism*T_ism
  real, parameter :: E_ism = 0.5*rho_ism*(u_ism**2) + P_ism/(gamma-1)
  !
  !
!  real, parameter :: frac = 0.9     ! Fracción de energía cinética a total
!  real, parameter :: Ekin = frac*ESN         ! Energía cinética
!  real, parameter :: Eth = (1-frac)*ESN      ! Energía térmica
!  real, parameter :: rho_SN = MSN/(4.0*PI/3.0*RSN**3)   ! Densidad interior
!  real, parameter :: vmax = sqrt(10.0/3.0*Ekin/MSN)       ! Velocidad en el borde
!  real, parameter :: P_SN = (gamma-1)*Eth/(4.0*PI/3.0*RSN**3) ! Presión interior
!  real, parameter :: rho = rho_SN + rho_ism
!
  real :: x
  integer :: i
  real :: veloc
!
  !  fill the vector u
  do i=0, nx+1
    x=real(i)*dx   ! obtain the position $x_i$
!    if( x <= RSN )  then
!       veloc = (x/RSN)*vmax
!       u(1,i)=rho_SN
!       u(2,i)=rho_SN* veloc
!       u(neq,i)=0.5*rho_SN*veloc**2 + P_SN/(gamma-1)  
!    else
      u(1,i)=rho_ISM
      u(2,i)=rho_ISM*u_ism
      u(neq,i)=E_ism
!   end if
  end do
!
  !   reset the counters and time to 0
  time=0.
  tprint=0.
  itprint=0
  print*,'Condiciones iniciales ' 
  print*,'====================================='
  print*,'Medio interestelar ' 
  print*,'---------------------------------------'
  print*,'Densidad ISM = ', rho_ISM, 'g/cm^3'
  print*,'Velocidad ISM = ', u_ism, 'm/s'
  print*,'Presión ISM = ', P_ISM, 'dyn/cm^2'
  print*,'---------------------------------------'
!
! 
call jet(time)
  return
end subroutine initflow
!
!=======================================================================
! generates initial condition
subroutine jet(time)
  use globals
  use constants
  implicit none
  real, intent(in) :: time
  real :: n_j, rho_j, l_j, u_j, temp_j, mu_j, P_j, E_j
  real :: x
  integer :: i
!
  ! jet properies
  n_j = 1000.
  l_j = 3.*dx
  temp_j = 1000.   ! Kelvin
  u_j = 200.*kPS   !in km/s
  mu_j = 0.6
  !
  rho_j= n_j*mu_j*mh
  P_j = boltz*n_j*temp_j
  E_j = 0.5*rho_j*(u_j**2) + P_j/(gamma-1)
  !
  !
  !  fill the vector u
  do i=0, nx+1
     x=real(i)*dx   ! obtain the position $x_i$
     !
     if (x <= l_j) then
        u(1,i)=rho_j
        u(2,i)=rho_j*u_j
        u(neq,i)=E_j
     end if
  end do
  !
  return
end subroutine jet
!=======================================================================
! output to file
subroutine output(itprint)
  use globals
  implicit none
  integer, intent(in) :: itprint
  character (len=31) file1
  real                :: temp
  real,dimension(neq) :: prim
  integer :: i
  !
  ! open output file
  write(file1,'(a,i3.3,a)') './outputs/output-',itprint,'.dat'
  open(unit=10,file=file1,status='unknown')
  !
  ! writes x and u
  do i=1,nx
    call uprim(u(:,i),prim,temp)
    write(10,*) real(i)*dx,prim(1),prim(2),prim(neq)
  end do
  !
  ! closes output file
  close(10)
  !
  return
end subroutine output

!=======================================================================
! computes the timestep allowed by the CFL criterium
subroutine timestep(dt)
  use globals
  implicit none
  real, intent(out) ::dt
    ! Courant number 
  real :: temp,cs,csound,del
  real,dimension(neq) :: prim
  integer :: i
  !
  del=1.e+30
  do i=1,nx
    call uprim(u(:,i),prim,temp)
    cs=csound(prim(1),prim(neq))
    del=min(del,dx/(abs(prim(2))+cs))
!    print*,del,prim(2),cs
  enddo
  !
  dt=Co*del
  !
  return
end subroutine timestep
!
function csound(rho,press)
  use globals, only : gamma
  implicit none
  real :: rho,press,csound
  !
  csound=sqrt(gamma*press/rho)
  !
  return
end function
!
subroutine uprim(uu,prim,temp)
  use constants
  use globals
  implicit none
  real :: ek, et,temp
  real,dimension(neq) :: uu,prim
  real :: fl_r,fl_t
  !
  fl_r= den_floor
  fl_t= temp_floor
  !
  prim(1)=uu(1)
  prim(1)=max(fl_r,prim(1)) 
  prim(2)=uu(2)/prim(1)
  ek=0.5*prim(1)*prim(2)**2.
  et=uu(neq)-ek
  prim(neq)=et*(gamma-1.)
  !
  temp=max(fl_t,prim(neq)/(prim(1)*boltz/(mu*mh)))
  prim(neq)=prim(1)*boltz*temp/(mu*mh) 
  !
  return
end subroutine uprim
!=======================================================================
! integration from t to t+dt methods
subroutine tstep(dt,time)
  use globals
  implicit none
  real, dimension(neq,0:nx+1) :: up
  real, intent(in) :: dt, time
  real :: dtx
  integer :: i
  ! 
  call tstep_hll(u,f,dt)
  !call Mac(u,f,dt)
  ! call Lax(u,f,dt)
  !
  call jet(time)

  return
end subroutine tstep

!=======================================================================
! Obtain the fluxes F
subroutine fluxes(u,f)
  use globals, only :neq,nx,gamma
  implicit none
  real, dimension(neq) :: prim
  integer :: i
  real,dimension(neq,0:nx+1) :: u,f
  real :: temp,etot
!
   do i=0,nx+1
    call uprim(u(:,i),prim,temp)
    Etot=0.5*prim(1)*prim(2)**2.+prim(neq)/(gamma-1.)
    f(1,i)=prim(1)*prim(2)
    f(2,i)=prim(1)*prim(2)**2.+prim(neq)
    f(neq,i)=prim(2)*(etot+prim(neq))
  enddo
!
  return
end subroutine fluxes
!
  subroutine eulerfluxes(prim,ff)
  use globals, only : neq,nx,gamma
  implicit none
  real, dimension(neq) :: prim
  real,dimension(neq) :: ff
  real :: temp,etot
  integer :: i

    ! call uprim(u(:,i),prim,temp)
    Etot=0.5*prim(1)*prim(2)**2.+prim(neq)/(gamma-1.)
    ff(1)=prim(1)*prim(2)
    ff(2)=prim(1)*prim(2)**2.+prim(neq)
    ff(neq)=prim(2)*(etot+prim(neq))
    
  return
end subroutine eulerfluxes
!
!=======================================================================
! Set boundary conditions
subroutine boundaries(u)
  use globals, only : nx,neq
  implicit none
  real,dimension(neq,0:nx+1) :: u
  ! free outflow
  u(:,0)=u(:,1)
  u(:,nx+1)=u(:,nx)
  !  u(2,nx+1)=-u(2,nx)
  !
  return
end subroutine boundaries
!=======================================================================
subroutine tstep_hll(u,f,dt)
  use globals, only : nx,neq,dx,eta
  implicit none
  real,dimension(neq,0:nx+1) :: u,f
  real,dimension(neq,0:nx+1) :: up,upp
  real,dimension(neq) :: ss,prim
  real, intent(in)   :: dt
  real :: dtx, temp
  integer :: i
  dtx=dt/dx
!
  !PREDICTOR-----------------------
  !
!
  call hllfluxes(u,f)
  do i=1,nx
    !
     call uprim(u(:,i),prim,temp)
     call sources(i,prim,ss)
     !
     upp(:,i)=u(:,i)-(dtx/2.)*(f(:,i+1)-f(:,i))+eta*(u(:,i+1)+u(:,i-1)-2.*u(:,i))!-dt*ss(:)
!!
  end do
!!
  call boundaries(upp)
!!
  !corrector-----------------------
  !
!
call hllfluxes(upp,f)
  do i=1,nx
!   
     call uprim(upp(:,i),prim,temp)
     call sources(i,prim,ss)
    up(:,i)=(u(:,i)+upp(:,i))/2.-(dtx)*(f(:,i)-f(:,i-1))+eta*(upp(:,i+1)+upp(:,i-1)-2.*upp(:,i))!-dt*ss(:)
!
  end do
!
  call boundaries(up)

  u(:,:)=up(:,:)

  return
end subroutine tstep_hll

subroutine Lax(u,f,dt)
  use globals, only : nx,neq,dx,eta
  implicit none
  real,dimension(neq,0:nx+1) :: u,f
  real,dimension(neq,0:nx+1) :: up
  real,dimension(neq) :: ss,prim
  real, intent(in) :: dt
  real :: dtx, temp
  integer :: i
!
   dtx=dt/dx
!
   call fluxes(u,f)
   !
   do i=1,nx
      call uprim(u(:,i), prim, temp)
      ! flujos
      call sources(i,prim,ss)
      up(:,i)=0.5*(u(:,i-1)+u(:,i+1)-dtx*(f(:,i+1)-f(:,i-1)))+eta*(u(:,i+1)+u(:,i-1)-2.*u(:,i))-dt*ss(:)
   end do
   call boundaries(up)
   u(:,:)=up(:,:)
   return
 end subroutine Lax
!
subroutine sources(i,prim,ss)
use globals, only :neq,nx,gamma, dx 
  implicit none
  integer :: i
  real, dimension(neq,0:nx+1) :: u,f
  real, dimension(neq) :: prim, s,ss
  real ::  r, Etot

  ss(:) = 0.
  call symmetry(i,prim,s)
  ss(:) =+ s(:)
    !
    !   call cooling(i,prim,s)
    !   ss(:,i) = ss(:,i) + s(:,i)
    !
  return
end subroutine sources

subroutine symmetry(i,prim,ss)
  use globals, only :neq,nx,gamma, dx,u,alpha

  implicit none
  integer :: i
  real, dimension(neq) :: prim, ss
  real :: r, Etot, term,temp 
  !
  r = real(i)*dx
  Etot = 0.5*prim(1)*prim(2)**2 + prim(neq)/(gamma-1.)
  term = alpha/r 
  ss(1) = term*prim(1)*prim(2)
  ss(2) = term*prim(1)*prim(2)**2
  ss(neq) = term*prim(2)*(Etot+prim(neq))
end subroutine
!
 subroutine hllfluxes(u,f)
   use globals, only : neq,nx,gamma,dx,eta
   implicit none 
   real, dimension(neq,0:nx+1) :: u,f
   real,dimension(neq) :: uul, uur, ff, priml,primr
   real :: temp
   integer :: i
   !
  do i=1,nx
    call uprim(u(:,i),priml,temp)
    call uprim(u(:,i+1),primr,temp)
    uul=u(:,i)
    uur=u(:,i+1)
    call hll(priml,primr,uul,uur,ff)
    f(:,i)=ff(:)
  enddo
 return
end subroutine hllfluxes


 subroutine hll(priml,primr,uul,uur,ff)
  use globals, only : neq,nx,gamma,dx,eta
  implicit none
  real, dimension(neq) :: priml, primr
  !real,dimension(neq,0:nx+1) :: up
  real,dimension(neq) :: uul, uur, ffl, ffr,ff
  real :: csl, csr, sl, sr, T
  integer :: i,j
  real :: dtx

  ! Este método divide en dos, en izquiera (L) y derecha (R)
    !
    call eulerfluxes(priml,ffl) 
    call eulerfluxes(primr,ffr)
    !
    csl =sqrt(gamma*priml(neq)/priml(1))
    csr =sqrt(gamma*primr(neq)/primr(1))
!
!
    sl = dmin1(priml(2)-csl,primr(2)-csr)
    sr = dmax1(priml(2)+csl,primr(2)+csr)

    ! Ahora se evalua la dirección de los flujos

  if(sl .gt. 0.) then

    ff(:) = ffl(:)
    
  else if (sr .lt. 0.) then

    ff(:) = ffr(:)
    
  else
!
    ff(:) = (sr*ffl(:)-sl*ffr(:)+sl*sr*(uur(:)-uul(:)))/(sr-sl)

  end if

  end subroutine hll
!

!

!
function AVER(A,B)
  !
  !    Averaging function
  !
  AB=A*B
  IF(AB.LE.0) THEN
     AVER=0.
  ELSE
     AVER=(AB*B+AB*A)/(A*A+B*B)
  END IF
  !
  RETURN
END FUNCTION AVER
!---------------------------------

 subroutine Mac(u,f,dt)
   use globals, only : nx,neq,dx,eta
   implicit none
   real,dimension(neq,0:nx+1) :: u,f
   real,dimension(neq,0:nx+1) :: up, ut, ft
   real,dimension(neq) :: primr,ss
   real, intent(in) :: dt
   real :: dtx, temp
   integer :: i, j
   
   dtx=dt/dx
   !
   call fluxes(u,f)
   !
   do i=1,nx
      call uprim(u(:,i),primr, temp)
      call sources(i,primr,ss)
    !ss(:) =0.
      ut(:,i) = u(:,i)- dtx*(f(:,i+1)-f(:,i))-dt*ss(:)+eta*(u(:,i+1)+u(:,i-1)-2.*u(:,i))
   end do
   !
   call boundaries(ut)
   call fluxes(ut,ft)
   !
   do i=1, nx
      call uprim(u(:,i),primr, temp)
      call sources(i,primr,ss)
      !ss(:) =0.
      up(:,i) = (u(:,i)+ut(:,i))/2 - (dtx/2.)*(ft(:,i)-ft(:,i-1))+eta*(ut(:,i+1)+ut(:,i-1)-2.*ut(:,i))!-dt*ss(:)
   end do  
!   
 call boundaries(up)
 
    u(:,:)=up(:,:)
    !
  return
  !  
end subroutine Mac

