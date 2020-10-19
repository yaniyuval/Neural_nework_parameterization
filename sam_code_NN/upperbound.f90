    
subroutine upperbound

use vars
use params
use mse              ! peters
implicit none
real coef
integer i,j,k
real tau_nudging
parameter ( tau_nudging = 3600.)
real told(nx,ny,nzm), qold(nx,ny,nzm), tmp, coefp  ! peters store old t, q

if(doMSE) then
  told = t(1:nx,1:ny,:)
  qold = q(1:nx,1:ny,:)
end if

if(dolargescale) then

!
! if there is an "observed" sounding - nudge two highest levels to it
! to avoid problems with the upper boundary.
!

  coef = dtn / tau_nudging
  do k=nzm-1,nzm
    do j=1,ny
      do i=1,nx
         t(i,j,k)=t(i,j,k)-(t(i,j,k)-tg0(k)-gamaz(k))*coef
         q(i,j,k)=q(i,j,k)-(q(i,j,k)-qg0(k))*coef
      end do
    end do
  end do

else
!
!  otherwise, preserve the vertical gradients:
!
  coef = dz*adz(nzm)
  gamt0=(t0(nzm-1)-t0(nzm-2))/(z(nzm-1)-z(nzm-2))
  gamq0=(q0(nzm-1)-q0(nzm-2))/(z(nzm-1)-z(nzm-2))
  do j=1,ny
   do i=1,nx
     t(i,j,nzm)=t(i,j,nzm-1)+gamt0*coef
     q(i,j,nzm)=q(i,j,nzm-1)+gamq0*coef
   end do    
  end do 
           
end if

! peters integrate up t, q upperbound
if(doMSE) then
  coefp = 1./dtn*dtfactor/float(navgMSE)
  told = cp*(t(1:nx,1:ny,:)-told)*coefp                     ! J/kg/s
  qold = lcond*(q(1:nx,1:ny,:)-qold)*coefp                  ! J/kg/s
  do i=1,nx
    do j=1,ny
      call columnint(told(i,j,:),tmp)
      slidamp_mse(i,j) = slidamp_mse(i,j) + tmp                ! W/m^2
      call columnint(qold(i,j,:),tmp)
      hdamp_mse(i,j) = hdamp_mse(i,j) + tmp                    ! W/m^2
    end do
  end do
end if

end   
