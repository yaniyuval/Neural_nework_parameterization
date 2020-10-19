subroutine surface()
	
use vars
use params
use simple_land, only : landmask_frac, xysurf_rough, nstart_surfrough, &
                        doxysurf_rough, xysoil_wet
implicit none
	
real t_s, q_s, u_h0
real taux0, tauy0, xlmo
real diag_ustar, coef, coef1
integer i,j,ic,jc,it,jt,jj,ii
real rrr,ranf_
double precision buffer(2), buffer1(2)
real delt,ssq,delq,cd,windspeed,wrk,wrk1,wrk2

real tland, lfluxbt(nx,ny), lfluxbq(nx,ny), lfluxbu(nx,ny), lfluxbv(nx,ny),z00
real zz00(0:nxp1,1-YES3D:nyp1)
! wrb:
real umin     ! Minimum wind speed at bottom level
real ubot    ! Bottom level u wind
real vbot    ! Bottom level v wind
parameter ( umin = 1.0 )
!
call task_rank_to_index(rank,it,jt) !Yani JY delete
! LES mode: 

if(.not.SFC_FLX_FXD) then

  if(.not.dosurfacefix) then

  lfluxbt = 0.0
  lfluxbq = 0.0
  lfluxbu = 0.0
  lfluxbv = 0.0
  fluxbt = 0.0
  fluxbq = 0.0
  fluxbu = 0.0
  fluxbv = 0.0

  do j=1,ny
    do i=1,nx

     tland = landmask_frac(i,j)  ! =0 if OCEAN, =1 if LAND

     if(tland.lt.0.999) then       ! THIS POINT HAS SOME OCEAN

       if(dobulksfc) then
          cd=1.1e-3
          !windspeed=5.
          ! wrb:
          jc=j+YES3D
          ic=i+1
          ubot=0.5*(u(ic,j,1)+u(i,j,1))+ug
!          ubot=0.5*(u(ic,j,1)+u(i,j,1)+u(ic,jc,1)+u(i,jc,1))+ug 
         vbot=0.5*(v(i,jc,1)+v(i,j,1))+vg
!          ubot=1.0*(u(i,j,1))+ug
!          vbot=1.0*(v(i,j,1))+vg
          windspeed=sqrt(ubot**2+vbot**2+umin**2)
          delt   = t(i,j,1)-gamaz(1) - sstxy(i,j)
          ssq = qsatw(sstxy(i,j),pres(1)) 
          delq   = q(i,j,1)  - ssq  
          wrk=(log(10/1.e-4)/log(z(1)/1.e-4))**2
          fluxbt(i,j) = -cd*windspeed*delt*wrk
          fluxbq(i,j) = -cd*windspeed*delq*wrk
          !ubot=1.0*(u(i,j,1))+ug
          !vbot=0.25*(v(i,j,1) + v(i,jc,1) + v(i-1,jc,1) + v(i-1,j,1))+vg
!          ubot=(0.2188*u(i+1,j-1,1)+0.2812*u(i+1,j,1)+0.2188*u(i,j-1,1)+0.2812*u(i,j,1))+ug
          ubot=(0.0625*u(i+1,j,1)+0.9375*u(i,j,1))+ug          
          vbot=(0.2812*v(i,j,1) + 0.2812*v(i,jc,1) + 0.2188*v(i-1,jc,1) + 0.2188*v(i-1,j,1))+vg
           windspeed=sqrt(ubot**2+vbot**2+umin**2)
    !      ubot=(0.0625*(u(i+1,j,1)**2)+0.9375*(u(i,j,1)**2))+ug**2
     !     vbot=(0.2812*(v(i,j,1)**2) + 0.2812*(v(i,jc,1)**2) + 0.2188*(v(i-1,jc,1)**2) + 0.2188*(v(i-1,j,1)**2))+vg**2
      !    windspeed=sqrt(ubot+vbot+umin**2)

          fluxbu(i,j) = -rho(1)*(u(i,j,1)+ug)*cd*windspeed*wrk
!          ubot=0.25*(u(i+1,j-1,1)+u(i+1,j,1)+u(i,j-1,1)+u(i,j,1))+ug
          !vbot=1.0*(v(i,j,1))+vg
          !windspeed=sqrt(ubot**2+vbot**2+umin**2)
          ubot=(0.2188*u(i+1,j-1,1)+0.2812*u(i+1,j,1)+0.2188*u(i,j-1,1)+0.2812*u(i,j,1))+ug
          vbot=(0.9375*v(i,j,1) + 0.0625*v(i,jc,1)) + vg
          windspeed=sqrt(ubot**2+vbot**2+umin**2)

       !   ubot=(0.2188*(u(i+1,j-1,1)**2)+0.2812*(u(i+1,j,1)**2)+0.2188*(u(i,j-1,1)**2)+0.2812*(u(i,j,1)**2))+ug**2
      !    vbot=(0.9375*(v(i,j,1)**2) + 0.0625*(v(i,jc,1)**2)) + vg**2
     !     windspeed=sqrt(ubot+vbot+umin**2)

          fluxbv(i,j) = -rho(1)*(v(i,j,1)+vg)*cd*windspeed*wrk
         ! if (i.eq.1) then
         ! print *, 'the SST:', sstxy(1,j),'The lat indices', j, jt
         ! end if                
       else

         jc=j+YES3D
         ic=i+1
         call oceflx(pres(1),0.5*(u(ic,j,1)+u(i,j,1))+ug, &
                       0.5*(v(i,jc,1)+v(i,j,1))+vg, &
                       t(i,j,1)-gamaz(1),q(i,j,1),t(i,j,1),z(1), &
                       sstxy(i,j), fluxt0, fluxq0, taux0, tauy0, q_s, sfcum)

         fluxbt(i,j) = fluxt0
         fluxbq(i,j) = fluxq0
            
         if (SFC_TAU_FXD) then
            u_h0 = max(1.,sqrt((u0(1)+ug)**2+(v0(1)+vg)**2))
            taux0 = -rho(1)*(u0(1)+ug)/u_h0*tau0
            tauy0 = -rho(1)*(v0(1)+vg)/u_h0*tau0
         endif
 
         fluxbu(i,j) = taux0
         fluxbv(i,j) = tauy0

       endif !bulksfc

     end if     ! if(tland.le.0.999)

     if(tland.gt.0.001) then               ! THIS POINT HAS SOME LAND

       coef = (1000./pres0)**(rgas/cp)
       jc=j+YES3D
       ic=i+1

       coef1 = (1000./pres(1))**(rgas/cp)
       t_s = sstxy(i,j)*coef
       q_s = xysoil_wet(i,j)*qsatw(sstxy(i,j),pres(1))

       ! limit surface roughness to some fraction of maxsurfrough if in spin up
       !  stage.
       if(nstep.le.nstart_surfrough .and. doxysurf_rough) then
         coef = float(nstep)/float(nstart_surfrough)
         z00 = max(xysurf_rough(i,j)*coef, 0.001) 
       else
         z00 = xysurf_rough(i,j)
       end if

       call landflx(pres(1),(t(i,j,1)-gamaz(1))*coef1, t_s,   &
                    q(i,j,1), q_s, 0.5*(u(ic,j,1)+u(i,j,1))+ug,     &
                    0.5*(v(i,jc,1)+v(i,j,1))+vg, z(1), z00,        &
                    fluxt0, fluxq0, taux0, tauy0, xlmo)

       lfluxbt(i,j) = fluxt0
       lfluxbq(i,j) = fluxq0
       lfluxbu(i,j) = taux0
       lfluxbv(i,j) = tauy0

     end if

     ! take weighted average of land/ocean fluxes
     fluxbt(i,j) = (1.-tland)*fluxbt(i,j) + tland*lfluxbt(i,j)
     fluxbq(i,j) = (1.-tland)*fluxbq(i,j) + tland*lfluxbq(i,j)
     fluxbu(i,j) = (1.-tland)*fluxbu(i,j) + tland*lfluxbu(i,j)
     fluxbv(i,j) = (1.-tland)*fluxbv(i,j) + tland*lfluxbv(i,j)

   end do
  end do

  else               ! dosurfacefix

    ! limit surface roughness to some fraction of maxsurfrough if in spin up
    !  stage.
    if(nstep.le.nstart_surfrough .and. doxysurf_rough) then
      coef = float(nstep)/float(nstart_surfrough)
      do j=1-YES3D,nyp1
        do i=0,nxp1
          zz00(i,j) = max(xysurf_rough(i,j)*coef, 0.001)
        end do
      end do
    else
      zz00 = xysurf_rough
    end if

    do j=1,ny
      do i=1,nx
        call surface_scheme(i,j,1,zz00)
        call surface_scheme(i,j,2,zz00)
        call surface_scheme(i,j,3,zz00)
      end do
    end do

  end if           ! if(dosurfacefix)

end if! .not.SFC_FLX_FXD


if(SFC_FLX_FXD) then

  if(doxy) then
    if(masterproc) print*,'doxy=',doxy,':SFC_FLX_FXD should be .F'
  end if

  u_h0 = max(1.,sqrt((u0(1)+ug)**2+(v0(1)+vg)**2))

  if(.not.SFC_TAU_FXD) then
    if(LAND) z0 = 0.03
    if(OCEAN) z0 = 0.0001

    tau0 = diag_ustar(z(1),  &
                bet(1)*(fluxt0+0.61*(t0(1)-gamaz(1))*fluxq0),u_h0,z0)**2  

  end if ! .not.SFC_TAU_FXD

  taux0 = -rho(1)*(u0(1)+ug)/u_h0*tau0
  tauy0 = -rho(1)*(v0(1)+vg)/u_h0*tau0

  if (ocean_type.eq.0) then 
     fluxbt(:,:) = fluxt0
     fluxbq(:,:) = fluxq0
  else
     fluxbt = fluxbt/sqrt(betafactor)
     fluxbq = fluxbq/sqrt(betafactor)
  endif
  fluxbu(:,:) = taux0
  fluxbv(:,:) = tauy0

end if ! SFC_FLX_FXD

!
! Homogenize the surface scalar fluxes if needed for sensitivity studies
!
   if(dosfchomo) then

	fluxt0 = 0.
	fluxq0 = 0.
	do j=1,ny
         do i=1,nx
	   fluxt0 = fluxt0 + fluxbt(i,j)
	   fluxq0 = fluxq0 + fluxbq(i,j)
         end do
        end do
	fluxt0 = fluxt0 / float(nx*ny)
	fluxq0 = fluxq0 / float(nx*ny)
        if(dompi) then
            buffer(1) = fluxt0
            buffer(2) = fluxq0
            call task_sum_real8(buffer,buffer1,2)
	    fluxt0 = buffer1(1) /float(nsubdomains)
	    fluxq0 = buffer1(2) /float(nsubdomains)
        end if ! dompi
	fluxbt(:,:) = fluxt0
	fluxbq(:,:) = fluxq0

   end if

!
! Homogenize the surface scalar fluxes in x if needed for sensitivity studies
!
   if(dosfchomox) then
     call task_rank_to_index(rank,it,jt)
     do jj=1,nsubdomains_y 
	do j=1,ny
           fluxt0 = 0.
           fluxq0 = 0.
        if((jj-1) * (ny_gl/nsubdomains_y).eq.jt) then
         do i=1,nx
	   fluxt0 = fluxt0 + fluxbt(i,j)
	   fluxq0 = fluxq0 + fluxbq(i,j)
         end do
	fluxt0 = fluxt0 / float(nx)
	fluxq0 = fluxq0 / float(nx)
        endif
        if(dompi) then
            buffer(1) = fluxt0
            buffer(2) = fluxq0
            call task_sum_real8(buffer,buffer1,2)
	    fluxt0 = buffer1(1) /float(nsubdomains_x)
	    fluxq0 = buffer1(2) /float(nsubdomains_x)
        end if ! dompi 
        if((jj-1) * (ny_gl/nsubdomains_y).eq.jt) then
           if(j+jt.le.ny_gl/2+hadley_lim.and.j+jt.ge.ny_gl/2-hadley_lim+1) then
              do i=1,nx/block_size
                 wrk1 = 0.
                 wrk2 = 0.
                 do ii=1,block_size
                    wrk1 = wrk1 + fluxbt((i-1)*block_size+ii,j)
                    wrk2 = wrk2 + fluxbq((i-1)*block_size+ii,j)
                 end do
                 fluxbt((i-1)*block_size+1:(i-1)*block_size+block_size,j) = fluxbt((i-1)*block_size+1:(i-1)*block_size+block_size,j)+(fluxt0-wrk1/float(block_size))
                 fluxbq((i-1)*block_size+1:(i-1)*block_size+block_size,j) = fluxbq((i-1)*block_size+1:(i-1)*block_size+block_size,j)+(fluxq0-wrk2/float(block_size))
              end do
           endif
        endif
     end do
  enddo
end if
!
! Homogenize the surface drag in x in certain latitudes if needed for sensitivity studies
!

   if(SFC_TAU_MASK) then
      call task_rank_to_index(rank,it,jt)
     do jj=1,nsubdomains_y 
	do j=1,ny
           taux0  = 0.
           tauy0  = 0.
        if((jj-1) * (ny_gl/nsubdomains_y).eq.jt) then
         do i=1,nx
	   taux0 = taux0 + fluxbu(i,j)
	   tauy0 = tauy0 + fluxbv(i,j)
         end do
	taux0 = taux0 / float(nx)
	tauy0 = tauy0 / float(nx)
        endif
        if(dompi) then
            buffer(1) = taux0
            buffer(2) = tauy0
            call task_sum_real8(buffer,buffer1,2)
	    taux0 = buffer1(1) /float(nsubdomains_x)
	    tauy0 = buffer1(2) /float(nsubdomains_x)
        end if ! dompi 
        if((jj-1) * (ny_gl/nsubdomains_y).eq.jt) then
         fluxbu(:,j)=fluxbu(:,j)*sfctaumask(j)+taux0*(1.-sfctaumask(j))
         fluxbv(:,j)=fluxbv(:,j)*sfctaumask(j)+tauy0*(1.-sfctaumask(j))
        endif
        end do
     enddo
  endif

   if (dobetafactor) then
      fluxbt = fluxbt * sqrt(betafactor)
      fluxbq = fluxbq * sqrt(betafactor)
      fluxbu = fluxbu * sqrt(betafactor)
      fluxbv = fluxbv * sqrt(betafactor)
   end if

if(doperturb_realtime) then
   do j=1,ny
      do i=1,nx
         rrr=1.-2.*ranf_()
         fluxbt(i,j)=fluxbt(i,j)*(1+0.1*rrr)
         rrr=1.-2.*ranf_()
         fluxbq(i,j)=fluxbq(i,j)*(1+0.1*rrr)
      end do
   end do
endif

shf_xy(:,:) = shf_xy(:,:) + fluxbt(:,:) * dtfactor
lhf_xy(:,:) = lhf_xy(:,:) + fluxbq(:,:) * dtfactor

end




! ----------------------------------------------------------------------
!
! DISCLAIMER : this code appears to be correct but has not been
!              very thouroughly tested. If you do notice any
!              anomalous behaviour then please contact Andy and/or
!              Bjorn
!
! Function diag_ustar:  returns value of ustar using the below 
! similarity functions and a specified buoyancy flux (bflx) given in
! kinematic units
!
! phi_m (zeta > 0) =  (1 + am * zeta)
! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
!
! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
!
! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface 
! Layer, in Workshop on Micormeteorology, pages 67-100.
!
! Code writen March, 1999 by Bjorn Stevens
!
! Code corrected 8th June 1999 (obukhov length was wrong way up,
! so now used as reciprocal of obukhov length)

      real function diag_ustar(z,bflx,wnd,z0)

      implicit none
      real, parameter      :: vonk =  0.4   ! von Karmans constant
      real, parameter      :: g    = 9.81   ! gravitational acceleration
      real, parameter      :: am   =  4.8   !   "          "         "
      real, parameter      :: bm   = 19.3   !   "          "         "
      real, parameter      :: eps  = 1.e-10 ! non-zero, small number

      real, intent (in)    :: z             ! height where u locates
      real, intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
      real, intent (in)    :: wnd           ! wind speed at z
      real, intent (in)    :: z0            ! momentum roughness height

      integer :: iterate
      real    :: lnz, klnz, c1, x, psi1, zeta, rlmo, ustar

      lnz   = log(z/z0) 
      klnz  = vonk/lnz              
      c1    = 3.14159/2. - 3.*log(2.)

      ustar =  wnd*klnz
      if (bflx /= 0.0) then 
        do iterate=1,4
          rlmo   = -bflx * vonk/(ustar**3 + eps)   !reciprocal of
                                                   !obukhov length
          zeta  = z*rlmo
          if (zeta > 0.) then
            ustar =  vonk*wnd  /(lnz + am*zeta)
          else
            x     = sqrt( sqrt( 1.0 - bm*zeta ) )
            psi1  = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + c1
            ustar = wnd*vonk/(lnz - psi1)
          end if
        end do
      end if

      diag_ustar = ustar

      return
      end function diag_ustar
! ----------------------------------------------------------------------

