subroutine precip_fall

!     positively definite monotonic advection with non-oscillatory option
!     and gravitational sedimentation (rain and snow advection)

use vars
use params
use mse, only : prec_inst_mse, doMSE, prec_inst_frozen_mse  ! peters
use microscaling
implicit none
 	
real mx(nzm),mn(nzm)
real www(nz),fz(nz)
real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real f0(nzm),df0(nzm)
real eps
integer i,j,k,kc,kb
logical nonos

real y 
real vrain, vsnow, vgrau, crain, csnow, cgrau
real lfac(nz), qrr, qss, qgg
real pp,pn, lat_heat
real wmax, omp, omg
real dtn_old

real wp(nzm), tmp_qp(nzm), irhoadz(nzm), iwmax(nzm), rhofac(nzm), prec_cfl
integer nprec, iprec
real fac_frozen, fac_total ! peters

pp(y)= max(0.,y)
pn(y)=-min(0.,y)


!--------------------------------------------------------
dtn_old = dtn

if(dokruegermicro) then
   call precip_fall_krueger()
   return
end if

if(dobetafactor.and.daremicrophysics.and..not.domicroscaling) dtn=dtn*mpbetafactor**0.5

eps = 1.e-10
nonos = .true.
  
crain = b_rain / 4.
csnow = b_snow / 4. 
cgrau = b_grau / 4. 
vrain = a_rain * gamr3 / 6. / (pi * rhor * nzeror) ** crain	  
vsnow = a_snow * gams3 / 6. / (pi * rhos * nzeros) ** csnow
vgrau = a_grau * gamg3 / 6. / (pi * rhog * nzerog) ** cgrau

! Initialize arrays that accumulate surface precipitation flux

 if(doMSE) then
    prec_inst_mse = 0.   !peters, need to accumulate over iprec steps
    prec_inst_frozen_mse = 0.
 end if

 if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
   do j=1,ny
    do i=1,nx
     precsfc(i,j)=0.
    end do
   end do
   do k=1,nzm
    precflux(k) = 0.
   end do
 end if

 do k = 1,nzm ! Initialize arrays which hold precipitation fluxes for stats.
    qpfall(k)=0. 
    tlat(k) = 0.
 end do

 do k = 1,nzm
    rhofac(k) = sqrt(1.29/rho(k)) ! Factor in precipitation velocity formula
    irhoadz(k) = 1./(rho(k)*adz(k)) ! Useful factor
    kb = max(1,k-1)
    wmax       = dz*adz(kb)/dtn   ! Velocity equivalent to a cfl of 1.0.
    iwmax(k)   = 1./wmax
 end do

! 	Add sedimentation of precipitation field to the vert. vel.

do j=1,ny
   do i=1,nx

      ! Compute precipitation velocity and flux column-by-column
      
      prec_cfl = 0.
      do k=1,nzm

	if(domicroscaling) dtn = dtn_scaled(i,j,k)	 

         wp(k)=0.
         omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
         lfac(k) = fac_cond+(1.-omp)*fac_fus
         if(k.eq.1) then      ! peters
           fac_frozen = (1.-omp)*lfus
           fac_total = omp*lcond + (1.-omp)*lfus
         end if
         if(qp(i,j,k).gt.qp_threshold) then
            if(omp.eq.1.) then
               wp(k)= rhofac(k)*vrain*(rho(k)*qp(i,j,k))**crain
            elseif(omp.eq.0.) then
               omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))
               qgg=omg*qp(i,j,k)
               qss=qp(i,j,k)-qgg
               wp(k)= rhofac(k)*(omg*vgrau*(rho(k)*qgg)**cgrau &
                                 +(1.-omg)*vsnow*(rho(k)*qss)**csnow)
            else
               omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))
               qrr=omp*qp(i,j,k)
               qss=qp(i,j,k)-qrr
               qgg=omg*qss
               qss=qss-qgg
               wp(k)=rhofac(k)*(omp*vrain*(rho(k)*qrr)**crain &
                     +(1.-omp)*(omg*vgrau*(rho(k)*qgg)**cgrau &
                          +(1.-omg)*vsnow*(rho(k)*qss)**csnow))
            endif
            prec_cfl = max(prec_cfl,wp(k)*iwmax(k)) ! Keep column maximum CFL
            wp(k) = -wp(k)*rhow(k)*dtn/dz 
         endif
      end do

      fz(nz)=0.
      www(nz)=0.
      lfac(nz)=0.

      ! If maximum CFL due to precipitation velocity is greater than 0.9,
      ! take more than one advection step to maintain stability.
!      if (prec_cfl.gt.0.9) then
!         nprec = CEILING(prec_cfl/0.9)
!zhiming
      if (prec_cfl.gt.0.3) then
         nprec = CEILING(prec_cfl/0.3)
         do k = 1,nzm
            ! wp already includes factor of dt, so reduce it by a
            ! factor equal to the number of precipitation steps.
            wp(k) = wp(k)/float(nprec) 
         end do
      else
         nprec = 1
      end if

      do iprec = 1,nprec

         do k = 1,nzm
            tmp_qp(k) = qp(i,j,k) ! Temporary array for qp in this column
         end do

         !-----------------------------------------

         if(nonos) then

            do k=1,nzm
               kc=min(nzm,k+1)
               kb=max(1,k-1)
               mx(k)=max(tmp_qp(kb),tmp_qp(kc),tmp_qp(k))
               mn(k)=min(tmp_qp(kb),tmp_qp(kc),tmp_qp(k))	  
            end do

         end if  ! nonos

         !  loop over iterations

         do k=1,nzm
            ! Define upwind precipitation flux
            fz(k)=tmp_qp(k)*wp(k)
         end do

         do k=1,nzm
            kc=k+1
            tmp_qp(k)=tmp_qp(k)-(fz(kc)-fz(k))*irhoadz(k) !Update temporary qp
         end do

         do k=1,nzm
            ! Also, compute anti-diffusive correction to previous
            ! (upwind) approximation to the flux
            kb=max(1,k-1)
            ! The precipitation velocity is a cell-centered quantity,
            ! since it is computed from the cell-centered
            ! precipitation mass fraction.  Therefore, a reformulated
            ! anti-diffusive flux is used here which accounts for
            ! this and results in reduced numerical diffusion.
            www(k) = 0.5*(1.+wp(k)*irhoadz(k)) &
                 *(tmp_qp(kb)*wp(kb) - tmp_qp(k)*wp(k)) ! works for wp(k)<0
         end do

         !---------- non-osscilatory option ---------------

         if(nonos) then

            do k=1,nzm
               kc=min(nzm,k+1)
               kb=max(1,k-1)
               mx(k)=max(tmp_qp(kb),tmp_qp(kc),tmp_qp(k),mx(k))
               mn(k)=min(tmp_qp(kb),tmp_qp(kc),tmp_qp(k),mn(k))	  
            end do

            do k=1,nzm
               kc=min(nzm,k+1)
               mx(k)=rho(k)*adz(k)*(mx(k)-tmp_qp(k)) &
                      /(pn(www(kc)) + pp(www(k))+eps)
               mn(k)=rho(k)*adz(k)*(tmp_qp(k)-mn(k)) &
                      /(pp(www(kc)) + pn(www(k))+eps)
            end do

            do k=1,nzm
               kb=max(1,k-1)
               ! Add limited flux correction to fz(k).
               fz(k) = fz(k) &                        ! Upwind flux
                    + pp(www(k))*min(1.,mx(k), mn(kb)) &
                    - pn(www(k))*min(1.,mx(kb),mn(k)) ! Anti-diffusive flux
            end do

         endif ! nonos

         ! Update precipitation mass fraction and liquid-ice static
         ! energy using precipitation fluxes computed in this column.
         do k=1,nzm
            
        if(domicroscaling) dtn = dtn_scaled(i,j,k)	 

            kc=k+1
            ! Update precipitation mass fraction.
            ! Note that fz is the total flux, including both the
            ! upwind flux and the anti-diffusive correction.
            qp(i,j,k)=qp(i,j,k)-(fz(kc)-fz(k))*irhoadz(k)
            qpfall(k)=qpfall(k)-(fz(kc)-fz(k))*irhoadz(k)  ! For qp budget
            lat_heat = -(lfac(kc)*fz(kc)-lfac(k)*fz(k))*irhoadz(k)
            t(i,j,k)=t(i,j,k)-lat_heat
            tlat(k)=tlat(k)-lat_heat            ! For energy budget
            precflux(k) = precflux(k) - fz(k)   ! For statistics
         end do
         precsfc(i,j) = precsfc(i,j) - fz(1) ! For statistics
         prec_xy(i,j) = prec_xy(i,j) - fz(1) ! For 2D output
         prec_inst(i,j) = -fz(1)*float(nprec)*dz/dtn
         if(doMSE) then ! peters
            prec_inst_mse(i,j) = prec_inst_mse(i,j) - fac_total*fz(1)
            prec_inst_frozen_mse(i,j) = prec_inst_frozen_mse(i,j) &
                                         - fac_frozen*fz(1)
         end if

         if (iprec.lt.nprec) then

            ! Re-compute precipitation velocity using new value of qp.
            do k=1,nzm

	if(domicroscaling) dtn = dtn_scaled(i,j,k)	 

               wp(k)=0.
               omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
               lfac(k) = fac_cond+(1.-omp)*fac_fus
               if(k.eq.1) then      ! peters
                 fac_frozen = (1.-omp)*lfus
                 fac_total = omp*lcond + (1.-omp)*lfus
               end if
               if(qp(i,j,k).gt.qp_threshold) then
                  if(omp.eq.1.) then
                     wp(k)= rhofac(k)*vrain*(rho(k)*qp(i,j,k))**crain
                  elseif(omp.eq.0.) then
                     omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))
                     qgg=omg*qp(i,j,k)
                     qss=qp(i,j,k)-qgg
                     wp(k)= rhofac(k)*(omg*vgrau*(rho(k)*qgg)**cgrau &
                                 +(1.-omg)*vsnow*(rho(k)*qss)**csnow)
                  else
                     omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))
                     qrr=omp*qp(i,j,k)
                     qss=qp(i,j,k)-qrr
                     qgg=omg*qss
                     qss=qss-qgg
                     wp(k)=rhofac(k)*(omp*vrain*(rho(k)*qrr)**crain &
                           +(1.-omp)*(omg*vgrau*(rho(k)*qgg)**cgrau &
                                +(1.-omg)*vsnow*(rho(k)*qss)**csnow))
                  endif
                  ! Decrease precipitation velocity by factor of nprec
                  wp(k) = -wp(k)*rhow(k)*dtn/dz/float(nprec)
                  ! Note: Don't bother checking CFL condition at each
                  ! substep since it's unlikely that the CFL will
                  ! increase very much between substeps when using
                  ! monotonic advection schemes.
               endif
            end do

            fz(nz)=0.
            www(nz)=0.
            lfac(nz)=0.

         end if

      end do !iprec

  end do
end do	

if(dostatis) then

 do k=1,nzm
   do j=dimy1_s,dimy2_s
     do i=dimx1_s,dimx2_s
        df(i,j,k) = t(i,j,k)
     end do
   end do
 end do
        
endif


 
do j=1,ny
  do i=1,nx
   if(qp(i,j,1).gt.1.e-6) s_ar=s_ar+dtfactor
  end do
end do


if(dostatis) then
          
 call stat_varscalar(t,df,f0,df0,t2leprec)
 call setvalue(twleprec,nzm,0.)
 call stat_sw2(t,df,twleprec)

endif


dtn = dtn_old
end subroutine precip_fall


