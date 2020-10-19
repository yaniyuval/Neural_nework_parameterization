
subroutine kurant

use vars

implicit none

integer i, j, k, ncycle1(1),ncycle2(1)
real wm(nz)
real um(nz)
real vm(nz)
real tkhmax(nz)
real cfl
real coefx, coefy, coefz, dcoefx, dcoefy, dcoefz !cw changed dcoef to dcoefz
real, parameter :: cfl_limit = 0.68
!real, parameter :: diff_limit = 0.45 ! original value is 0.45
real, parameter :: diff_limit = 0.25
!zhiming
!real, parameter :: diff_limit = 0.45

real :: cflx, cfly, cflz, dcflz, dcflxy
real :: cflxmax, cflymax, cflzmax, dcflzmax, dcflxymax

!bloss ncycle = 1
	
!bloss  w_max=0.
!bloss  do k = 1,nzm
!bloss   wm(k) = 0.
!bloss   tkhmax(k) = 0.
!bloss   um(k) = 0.
!bloss   vm(k) = 0.
!bloss   do j = 1,ny
!bloss    do i = 1,nx
!bloss     wm(k) = max(wm(k),abs(w(i,j,k)))
!bloss     tkhmax(k)=max(tkhmax(k),tkh(i,j,k))
!bloss     w_max=max(w_max,w(i,j,k))
!bloss     um(k) = max(um(k),abs(u(i,j,k)))
!bloss    end do
!bloss   end do
!bloss  end do
!bloss  
!bloss  if(RUN3D) then
!bloss  
!bloss   do k=1,nzm
!bloss    do j = 1,ny
!bloss     do i = 1,nx
!bloss      vm(k) = max(vm(k),abs(v(i,j,k)))
!bloss     end do
!bloss    end do
!bloss   end do
!bloss  
!bloss  endif
!bloss  
!bloss  cfl = 0.
!bloss  do k=1,nzm
!bloss   cfl = max(cfl	&
!bloss       ,wm(k)*dt/(dz*adzw(k)), 1.5*tkhmax(k)*grdf_z(k)*dt/(dz*adzw(k))**2 &
!bloss       ,um(k)*dt/dx, 1.5*tkhmax(k)*grdf_x(k)*dt/dx**2 &
!bloss       ,vm(k)*dt/dy, 1.5*tkhmax(k)*grdf_y(k)*dt/dy**2)
!bloss  end do
!bloss  	
!bloss  10 continue
!bloss   
!bloss  !bloss   compute ncycle as next integer greater than cfl/0.7
!bloss  ncycle = ceiling(cfl/0.7)

!bloss compute cfl as sum of diffusion and advection contributions.
!       require: (u*dt/dx/0.5)^2 + (nu*dt/dx^2/0.5)^2 < 1.

cfl = 0.
cflx = 0.
cfly = 0.
cflz = 0.
dcflz = 0. 
dcflxy = 0.

coefx = dt/dx/cfl_limit 
coefy = dt/dy/cfl_limit
dcoefx = (dt/dx/dx/diff_limit)**2
dcoefy = (dt/dy/dy/diff_limit)**2

do k = 1,nzm
   coefz = dt/(dz*adzw(k))/cfl_limit
   if(RUN3D) then 
      dcoefz = (dt/(dz*adzw(k))**2/diff_limit)**2
      do j = 1,ny
         do i = 1,nx

            cfl = max(cfl, &
                 (coefx*u(i,j,k))**2 &
                 + (coefy*v(i,j,k))**2 &
                 + (coefz*w(i,j,k))**2 &
                 + (dcoefx+dcoefy)*(tk_xy(i,j,k)*ravefactor)**2 &
                 + dcoefz*(tk_z(i,j,k)/ravefactor)**2)
!                 + dcoefz*(tk_z(i,j,k)+tk_z_uw(i,j,k)*sqrt(betafactor))**2)

         end do
      end do
   else
      dcoefz = (dt/(dz*adzw(k))**2/diff_limit)**2
      do j = 1,ny
         do i = 1,nx

            cfl = max(cfl, &
                 (coefx*u(i,j,k))**2 &
                 + (coefz*w(i,j,k))**2 &
                 + dcoefx*(tk_xy(i,j,k)*ravefactor)**2 &
                 + dcoefz*(tk_z(i,j,k)/ravefactor)**2 )
         end do
      end do
   end if

end do

!fix ncycle based size of cfl (times 1.15 for margin of safety) and last ncycle
if(nstep.eq.1) then ! special condition for first step
   ncycle = max(ncycle0,ncyclemin)
else
      ncycle = max(ceiling(1.15*sqrt(cfl)),ncycle-1,ncyclemin)
endif

cfl = sqrt(cfl)/float(ncycle)
maxcfl = cfl

if(dompi) then
  ncycle1(1)=ncycle
  call task_max_integer(ncycle1,ncycle2,1)
  ncycle=ncycle2(1)

  call task_max_real(cfl,maxcfl,1)
end if

!!$if(ncycle.gt.10) then
!!$   save3Dbin = .true.
!!$   call write_fields3D() !bloss add output of 3D fields for diagnostics
!!$   call write_fields2D() !bloss add output of 3D fields for diagnostics
!!$end if

if(masterproc) then
   write(*,999) nstep, icycle, maxcfl,maxndiff, maxtkz
999 format('SUBNSTEP = ',i8,'   NCYCLE=',i4,'   MAX CFL=',f7.3,'   MAXNDIFF=',f4.0,'   MAXTKZ=',f8.2)
end if

 if(ncycle.gt.ncyclemax) then
    print *,'the number of cycles exceeded ', ncyclemax
  save3Dbin = .true.
  call write_fields3D() !bloss add output of 3D fields for diagnostics
 call task_abort()
 end if

end subroutine kurant	
