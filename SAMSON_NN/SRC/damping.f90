
subroutine damping

!  "Spange"-layer damping at the domain top region

use vars
use mse         ! peters
use params      ! peters
implicit none

real tau_min	! minimum damping time-scale (at the top)
real tau_max    ! maxim damping time-scale (base of damping layer)
real damp_depth ! damping depth as a fraction of the domain height
parameter(tau_min=1200., tau_max=72000., damp_depth=0.3)
real tau(nzm)   
integer i, j, k, n_damp
real    tau_ydamp, tau_xdamp   !kzm
integer it, jt, jlim, ilim,jj  !kzm
real told(nx,ny,nzm), qold(nx,ny,nzm), tmp, coef, t0x(nzm),q0x(nzm),u0x(nzm),w0x(nzm),v0x(nzm)  ! peters store old t, q; kzm variables for ydamp
double precision buffer(5*nzm),buffer1(5*nzm)

if(.not.douwpbl.and.tau_min*dxdzfactor*ravefactor.lt.2*dt) then ! changed by cw 3/9/06 -- dxdzfactor not used with uwpbl
   print*,'Error: in damping() tau_min is too small!'
   call task_abort()
end if

if(doMSE) then
  told = t(1:nx,1:ny,:)
  qold = q(1:nx,1:ny,:)
end if

do k=nzm,1,-1
 if(z(nzm)-z(k).lt.damp_depth*z(nzm)) then 
   n_damp=nzm-k+1
 endif
end do

do k=nzm,nzm-n_damp,-1
 tau(k) = tau_min *(tau_max/tau_min)**((z(nzm)-z(k))/(z(nzm)-z(nzm-n_damp)))*dxdzfactor
 tau(k)=1./tau(k)
end do

do k = nzm, nzm-n_damp, -1
   do j=1,ny
    do i=1,nx
      dudt(i,j,k,na)= dudt(i,j,k,na)-(u(i,j,k)-u0(k)) * tau(k)
      dvdt(i,j,k,na)= dvdt(i,j,k,na)-(v(i,j,k)-v0(k)) * tau(k)
      dwdt(i,j,k,na)= dwdt(i,j,k,na)-w(i,j,k) * tau(k)
      t(i,j,k)= t(i,j,k)-dtn*(t(i,j,k)-t0(k)) * tau(k)
      q(i,j,k)= q(i,j,k)-dtn*(q(i,j,k)-q0(k)) * tau(k)
      if(dokruegermicro) qv(i,j,k)= qv(i,j,k)-dtn*(qv(i,j,k)-qv0(k))*tau(k)
    end do! i 
   end do! j
end do ! k

!===========================================================================
! UW ADDITIONS

if (dobetafactor) then
   ! Nudge q back to sounding.
   do k = nzm, nzm-n_damp, -1
      do j=1,ny
         do i=1,nx
            q(i,j,k)= q(i,j,k)-dtn*(q0(k)-qg0(k)) * tau(k)
            if(dokruegermicro) qv(i,j,k)= qv(i,j,k)-dtn*(qv0(k)-qg0(k))*tau(k)
         end do
      end do
   end do
end if

if(dotro3) then
   do k = nzm, nzm-n_damp, -1
      do j=1,ny
         do i=1,nx
            tro3(i,j,k)=tro3(i,j,k)-dtn*(tro3(i,j,k)-o3g0(k))*tau(k)
         end do
      end do
   end do
endif

if (donudging_hadley) then
   coef=1./tauls*sqrt(betafactor)
   jlim = max(10,floor(0.3*ny_gl))     ! Size of damping zone in grid points
   call task_rank_to_index(rank,it,jt)  ! Get processor index in x,y
     do jj=1,nsubdomains_y 
      do j = 1,ny
            t0x=0.
            q0x=0.
            u0x=0.
            v0x=0.
            w0x=0.
            if((jj-1) * (ny_gl/nsubdomains_y).eq.jt) then
               do i = 1,nx
                  t0x=t0x+t(i,j,:)
                  q0x=q0x+q(i,j,:)
                  u0x=u0x+u(i,j,:)
                  v0x=v0x+v(i,j,:)
                  w0x=w0x+w(i,j,1:nzm)
               end do
               t0x=t0x/float(nx)
               q0x=q0x/float(nx)
               u0x=u0x/float(nx)
               v0x=v0x/float(nx)
               w0x=w0x/float(nx)
            endif
        if(dompi) then
            buffer(1:nzm) = t0x
            buffer(nzm+1:2*nzm) = q0x
            buffer(2*nzm+1:3*nzm) = u0x
            buffer(3*nzm+1:4*nzm) = w0x
            buffer(4*nzm+1:5*nzm) = v0x
            call task_sum_real8(buffer,buffer1,5*nzm)
	    t0x = buffer1(1:nzm) /float(nsubdomains_x)
	    q0x = buffer1(nzm+1:2*nzm) /float(nsubdomains_x)
	    u0x = buffer1(2*nzm+1:3*nzm) /float(nsubdomains_x)
	    w0x = buffer1(3*nzm+1:4*nzm) /float(nsubdomains_x)
	    v0x = buffer1(4*nzm+1:5*nzm) /float(nsubdomains_x)
        end if ! dompi 
        if((jj-1) * (ny_gl/nsubdomains_y).eq.jt) then
            do i = 1,nx
               dudt(i,j,:,na)= dudt(i,j,:,na)-(u0x-ug0x(j+jt,:)) * coef
               dvdt(i,j,:,na)= dvdt(i,j,:,na)-(v0x-vg0x(j+jt,:)) * coef
               dwdt(i,j,1:nzm,na)= dwdt(i,j,1:nzm,na)-(w0x-wg0x(j+jt,:)) * coef
               t(i,j,:)      = t(i,j,:)      -dtn*(t0x-tg0x(j+jt,:)-gamaz(:)) * coef
               q(i,j,:)      = q(i,j,:)      -dtn*(q0x-qg0x(j+jt,:)) * coef
            end do
            if(doydamping) then
            if(j+jt.le. jlim) then
               tau_ydamp=1./tau_min/2.*(1.-cos(3.1415926*float(jlim+1-j-jt)/float(jlim)))/2.
            do i = 1,nx
               dudt(i,j,:,na)= dudt(i,j,:,na)-(u(i,j,:)-u0x) * tau_ydamp
               dvdt(i,j,:,na)= dvdt(i,j,:,na)-(v(i,j,:)) * tau_ydamp
               dwdt(i,j,1:nzm,na)= dwdt(i,j,1:nzm,na)-(w(i,j,1:nzm)-w0x) * tau_ydamp
               t(i,j,:)      = t(i,j,:)      -dtn*(t(i,j,:)-t0x) * tau_ydamp
               q(i,j,:)      = q(i,j,:)      -dtn*(q(i,j,:)-q0x) * tau_ydamp
            end do
            endif
            if(j+jt.ge.ny_gl+1-jlim) then
               tau_ydamp=1./tau_min/2.*(1.-cos(3.1415926*float(ny_gl-jlim-j-jt)/float(jlim)))/2.
            do i = 1,nx
               dudt(i,j,:,na)= dudt(i,j,:,na)-(u(i,j,:)-u0x) * tau_ydamp
               dvdt(i,j,:,na)= dvdt(i,j,:,na)-(v(i,j,:)) * tau_ydamp
               dwdt(i,j,1:nzm,na)= dwdt(i,j,1:nzm,na)-(w(i,j,1:nzm)-w0x) * tau_ydamp
               t(i,j,:)      = t(i,j,:)      -dtn*(t(i,j,:)-t0x) * tau_ydamp
               q(i,j,:)      = q(i,j,:)      -dtn*(q(i,j,:)-q0x) * tau_ydamp
            end do
            endif
         endif
      endif
   end do
enddo
endif

if ((.not.donudging_hadley).and.doydamping) then
   call task_rank_to_index(rank,it,jt)  ! Get processor index in x,y
   jlim = max(10,floor(0.25*ny_gl))     ! Size of damping zone in grid points
!This requires the damping zone to be fitted in one processor (in y)
      do j = 1,jlim
!         tau_ydamp=tau_min*(6.*tau_max/tau_min)**(float(j-1)/float(jlim-1))
!         tau_ydamp=1./tau_ydamp
         tau_ydamp=1./tau_min/4.*(1.-cos(3.1415926*float(jlim+1-j)/float(jlim)))/2.
            t0x=0.
            q0x=0.
            u0x=0.
            v0x=0.
            w0x=0.
            if(jt.eq.0) then
               do i = 1,nx
                  t0x=t0x+t(i,j,:)
                  q0x=q0x+q(i,j,:)
                  u0x=u0x+u(i,j,:)
                  v0x=v0x+v(i,j,:)
                  w0x=w0x+w(i,j,1:nzm)
               end do
               t0x=t0x/float(nx)
               q0x=q0x/float(nx)
               u0x=u0x/float(nx)
               v0x=v0x/float(nx)
               w0x=w0x/float(nx)
            endif
        if(dompi) then
            buffer(1:nzm) = t0x
            buffer(nzm+1:2*nzm) = q0x
            buffer(2*nzm+1:3*nzm) = u0x
            buffer(3*nzm+1:4*nzm) = w0x
            buffer(4*nzm+1:5*nzm) = v0x
            call task_sum_real8(buffer,buffer1,5*nzm)
	    t0x = buffer1(1:nzm) /float(nsubdomains_x)
	    q0x = buffer1(nzm+1:2*nzm) /float(nsubdomains_x)
	    u0x = buffer1(2*nzm+1:3*nzm) /float(nsubdomains_x)
	    w0x = buffer1(3*nzm+1:4*nzm) /float(nsubdomains_x)
	    v0x = buffer1(4*nzm+1:5*nzm) /float(nsubdomains_x)
        end if ! dompi 
        if(jt.eq.0) then
            do i = 1,nx
               dudt(i,j,:,na)= dudt(i,j,:,na)-(u(i,j,:)-u0x) * tau_ydamp
               dvdt(i,j,:,na)= dvdt(i,j,:,na)-(v(i,j,:)) * tau_ydamp
               dwdt(i,j,1:nzm,na)= dwdt(i,j,1:nzm,na)-(w(i,j,1:nzm)-w0x) * tau_ydamp
               t(i,j,:)      = t(i,j,:)      -dtn*(t(i,j,:)-t0x) * tau_ydamp
               q(i,j,:)      = q(i,j,:)      -dtn*(q(i,j,:)-q0x) * tau_ydamp
            end do
         endif
      end do
!damping at high y
      do j = ny-jlim+1,ny
!            tau_ydamp = tau_min*(6.*tau_max/tau_min)**(float(ny-j)/float(jlim-1))
!            tau_ydamp=1./tau_ydamp
            tau_ydamp=1./tau_min/4.*(1.-cos(3.1415926*float(j-ny+jlim)/float(jlim)))/2.
            t0x=0.
            q0x=0.
            u0x=0.
            v0x=0.
            w0x=0.
            if(jt.eq.ny_gl-ny_gl/nsubdomains_y) then
               do i = 1,nx
                  t0x=t0x+t(i,j,:)
                  q0x=q0x+q(i,j,:)
                  u0x=u0x+u(i,j,:)
                  v0x=v0x+v(i,j,:)
                  w0x=w0x+w(i,j,1:nzm)
               end do
               t0x=t0x/float(nx)
               q0x=q0x/float(nx)
               u0x=u0x/float(nx)
               v0x=v0x/float(nx)
               w0x=w0x/float(nx)
            endif
        if(dompi) then
            buffer(1:nzm) = t0x
            buffer(nzm+1:2*nzm) = q0x
            buffer(2*nzm+1:3*nzm) = u0x
            buffer(3*nzm+1:4*nzm) = w0x
            buffer(4*nzm+1:5*nzm) = v0x
            call task_sum_real8(buffer,buffer1,5*nzm)
	    t0x = buffer1(1:nzm) /float(nsubdomains_x)
	    q0x = buffer1(nzm+1:2*nzm) /float(nsubdomains_x)
	    u0x = buffer1(2*nzm+1:3*nzm) /float(nsubdomains_x)
	    w0x = buffer1(3*nzm+1:4*nzm) /float(nsubdomains_x)
	    v0x = buffer1(4*nzm+1:5*nzm) /float(nsubdomains_x)
        end if ! dompi 
        if(jt.eq.ny_gl-ny_gl/nsubdomains_y) then
            do i = 1,nx
               dudt(i,j,:,na)= dudt(i,j,:,na)-(u(i,j,:)-u0x) * tau_ydamp
               dvdt(i,j,:,na)= dvdt(i,j,:,na)-(v(i,j,:)) * tau_ydamp
               dwdt(i,j,1:nzm,na)= dwdt(i,j,1:nzm,na)-(w(i,j,1:nzm)-w0x) * tau_ydamp
               t(i,j,:)      = t(i,j,:)      -dtn*(t(i,j,:)-t0x) * tau_ydamp
               q(i,j,:)      = q(i,j,:)      -dtn*(q(i,j,:)-q0x) * tau_ydamp
            end do
         endif
      end do
   end if

if (doxdamping) then
   call task_rank_to_index(rank,it,jt)  ! Get processor index in x,y
   ilim = max(10,floor(0.125*nx_gl))     ! Size of damping zone in grid points
   if (it.eq.1) then ! Damping zone at low x
      do i = 1,ilim
         tau_xdamp=tau_min*(tau_max/tau_min)**(float(ilim-i)/float(ilim-1))
         tau_xdamp=1./tau_xdamp
         do k = 1,nzm
            do j = 1,ny
               dudt(i,j,k,na)= dudt(i,j,k,na)-(u(i,j,k)-u0(k)) * tau_xdamp
               dvdt(i,j,k,na)= dvdt(i,j,k,na)-(v(i,j,k)-v0(k)) * tau_xdamp
               dwdt(i,j,k,na)= dwdt(i,j,k,na)-w(i,j,k) * tau_xdamp
               t(i,j,k)      = t(i,j,k)      -dtn*(t(i,j,k)-t0(k)) * tau_xdamp
               q(i,j,k)      = q(i,j,k)      -dtn*(q(i,j,k)-q0(k)) * tau_xdamp
               if(dokruegermicro) &
                    qv(i,j,k)= qv(i,j,k)     -dtn*(qv(i,j,k)-qv0(k))*tau_xdamp
            end do
         end do
      end do
   elseif (it.eq.nsubdomains_x) then ! Damping zone at high latitudes
      do i = nx-ilim+1,nx
         do k = 1,nzm
            tau_xdamp = tau_min &
                 *(tau_max/tau_min)**(float(i-(nx-ilim+1))/float(ilim-1))
            tau_xdamp=1./tau_xdamp
            do j = 1,ny
               dudt(i,j,k,na)= dudt(i,j,k,na)-(u(i,j,k)-u0(k)) * tau_xdamp
               dvdt(i,j,k,na)= dvdt(i,j,k,na)-(v(i,j,k)-v0(k)) * tau_xdamp
               dwdt(i,j,k,na)= dwdt(i,j,k,na)-w(i,j,k) * tau_ydamp
               t(i,j,k)      = t(i,j,k)      -dtn*(t(i,j,k)-t0(k)) * tau_xdamp
               q(i,j,k)      = q(i,j,k)      -dtn*(q(i,j,k)-q0(k)) * tau_xdamp
               if(dokruegermicro) &
                    qv(i,j,k)= qv(i,j,k)     -dtn*(qv(i,j,k)-qv0(k))*tau_xdamp
            end do
         end do
      end do
   end if
end if

! peters integrate up t, q damping
if(doMSE) then
  coef = 1./dtn*dtfactor/float(navgMSE)
  told = cp*(t(1:nx,1:ny,:)-told)*coef                     ! J/kg/s
  qold = lcond*(q(1:nx,1:ny,:)-qold)*coef                  ! J/kg/s
  do i=1,nx
    do j=1,ny
      call columnint(told(i,j,:),tmp)
      slidamp_mse(i,j) = slidamp_mse(i,j) + tmp                ! W/m^2
      call columnint(qold(i,j,:),tmp)
      hdamp_mse(i,j) = hdamp_mse(i,j) + tmp                    ! W/m^2
    end do
  end do
end if

! END UW ADDITIONS
!===========================================================================

end subroutine damping
