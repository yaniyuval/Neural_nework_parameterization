subroutine diagnose
	
! Diagnose some useful stuff

use vars
use params
use mse, only: doMSE, diagnoseMSE ! peters
implicit none
	
integer i,j,k,kb,kc,k200,k500
double precision coef, coef1, buffer(nzm,8), buffer1(nzm,8)
real omn, omp

! UW ADDITIONS
integer  k150,k850,k1000
real     wgt, pii, coef2, tmp(2), tmp2(2), tmp_lwp
real     press_total, press_top, tmp_wgt1, tmp_wgt2
real     p200, p200b, coef200, coef200b
real     p850, p850b, coef850, coef850b
real     p1000, p1000b, coef1000, coef1000b
real     h850, t850, th850, qv850
real     t1000, th1000, qv1000
double precision tmpw1_xy(nx,ny), tmpw2_xy(nx,ny)
double precision tmps1_xy(nx,ny), tmps2_xy(nx,ny)
double precision tmph1_xy(nx,ny), tmph2_xy(nx,ny)

coef = 1./float(nx*ny)

k200 = nzm

do k=1,nzm
  u0(k)=0.
  v0(k)=0.
  t01(k) = tabs0(k)
  q01(k) = q0(k)
  t0(k)=0.
  tabs0(k)=0.
  q0(k)=0.
  qv0(k)=0.
  p0(k)=0.
  tke0(k)=0.
  kc=min(nzm,k+1)
  kb=max(1,k-1)
  if(pres(kc).le.200..and.pres(kb).gt.200.) k200=k
  coef1 = rho(k)*dz*adz(k)*dtfactor
  do j=1,ny
    do i=1,nx
     if(.not.dokruegermicro) then ! bloss
     omn  = max(0.,min(1.,(tabs(i,j,k)-tbgmin)*a_bg))
     omp  = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
     tabs(i,j,k) = t(i,j,k)-gamaz(k)+ &
	      (fac_cond+(1.-omn)*fac_fus)*qn(i,j,k)+ &
	      (fac_cond+(1.-omp)*fac_fus)*qp(i,j,k)
     elseif(dokruegermicro) then
     tabs(i,j,k) = t(i,j,k)-gamaz(k)+ &
              fac_cond*(qc(i,j,k)+qr(i,j,k)) &
              + fac_sub*(qi(i,j,k)+qs(i,j,k)+qg(i,j,k))
     end if
     u(i,j,k) = dudt(i,j,k,nc)
     v(i,j,k) = dvdt(i,j,k,nc)
     w(i,j,k) = dwdt(i,j,k,nc)

     u0(k)=u0(k)+u(i,j,k)
     v0(k)=v0(k)+v(i,j,k)
     p0(k)=p0(k)+p(i,j,k)
     t0(k)=t0(k)+t(i,j,k)
     tabs0(k)=tabs0(k)+tabs(i,j,k)
     q0(k)=q0(k)+q(i,j,k)
     qv0(k)=qv0(k)+q(i,j,k)-qn(i,j,k)
     tke0(k)=tke0(k)+tke(i,j,k)

     pw_xy(i,j) = pw_xy(i,j)+q(i,j,k)*coef1
     if(.not.dokruegermicro) then ! bloss
     cw_xy(i,j) = cw_xy(i,j)+qn(i,j,k)*omn*coef1
     iw_xy(i,j) = iw_xy(i,j)+qn(i,j,k)*(1.-omn)*coef1
     elseif(dokruegermicro) then ! bloss
      cw_xy(i,j) = cw_xy(i,j)+qc(i,j,k)*coef1
      iw_xy(i,j) = iw_xy(i,j)+qi(i,j,k)*coef1
     end if
    end do
  end do
  u0(k)=u0(k)*coef
  v0(k)=v0(k)*coef
  t0(k)=t0(k)*coef
  tabs0(k)=tabs0(k)*coef
  q0(k)=q0(k)*coef
  qv0(k)=qv0(k)*coef
  p0(k)=p0(k)*coef
  tke0(k)=tke0(k)*coef
  rel0(k) = q0(k)/qsatw(tabs0(k),pres(k))

end do ! k

k500 = nzm
do k = 1,nzm
   kc = min(nzm,k+1)
   if((pres(kc).le.500.).and.(pres(k).gt.500.)) then
      if ((500.-pres(kc)).lt.(pres(k)-500.))then
         k500=kc
      else
         k500=k
      end if
   end if
end do

do j=1,ny
 do i=1,nx
  usfc_xy(i,j) = usfc_xy(i,j) + u(i,j,1)*dtfactor  
  vsfc_xy(i,j) = vsfc_xy(i,j) + v(i,j,1)*dtfactor  
  u200_xy(i,j) = u200_xy(i,j) + u(i,j,k200)*dtfactor  
  v200_xy(i,j) = v200_xy(i,j) + v(i,j,k200)*dtfactor  
  w500_xy(i,j) = w500_xy(i,j) -rho(k500)*ggr*w(i,j,k500)*dtfactor
 end do
end do

if(dompi) then

  coef1 = 1./float(nsubdomains)
  do k=1,nzm
    buffer(k,1) = u0(k)
    buffer(k,2) = v0(k)
    buffer(k,3) = t0(k)
    buffer(k,4) = q0(k)
    buffer(k,5) = p0(k)
    buffer(k,6) = tabs0(k)
    buffer(k,7) = tke0(k)
    buffer(k,8) = qv0(k)
  end do
  call task_sum_real8(buffer,buffer1,nzm*8)
  do k=1,nzm
    u0(k)=buffer1(k,1)*coef1
    v0(k)=buffer1(k,2)*coef1
    t0(k)=buffer1(k,3)*coef1
    q0(k)=buffer1(k,4)*coef1
    p0(k)=buffer1(k,5)*coef1
    tabs0(k)=buffer1(k,6)*coef1
    tke0(k)=buffer1(k,7)*coef1
    qv0(k)=buffer1(k,8)*coef1
  end do

end if ! dompi

!=====================================================
! UW ADDITIONS

pii = atan2(0.,-1.)

! peters update MSE fields, and output them if time
if(doMSE) call diagnoseMSE()

! FIND VERTICAL INDICES OF 150MB, 850MB AND 1000MB, COMPUTE SWVP
k150 = nzm
k850 = nzm
k1000 = 2
do k=1,nzm
  if((k.gt.2).and.(pres(k).le.1000..and.pres(kb).gt.1000.)) k1000=k
  if(pres(k).le.850..and.pres(kb).gt.850.) k850=k
  if(pres(k).le.150..and.pres(kb).gt.150.) k150=k
  coef1 = rho(k)*dz*adz(k)*dtfactor
  do j=1,ny
    do i=1,nx

     ! Saturated water vapor path with respect to water. Can be used
     ! with water vapor path (= pw-cw-iw) to compute column-average
     ! relative humidity.   
     swvp_xy(i,j) = swvp_xy(i,j)+qsatw(tabs(i,j,k),pres(k))*coef1
    end do
  end do
end do ! k

! ACCUMULATE AVERAGES OF TWO-DIMENSIONAL STATISTICS
do j=1,ny
 do i=1,nx
  psfc_xy(i,j) = psfc_xy(i,j) + p(i,j,1)*dtfactor  

  ! Compute coefficients to interpolate values at 200 millibar
  p200  = pres(k200)   + 0.01*p(i,j,k200)   ! Convert p(i,j,k) from Pa to mbar
  p200b = pres(k200-1) + 0.01*p(i,j,k200-1) 
  coef200  = (200.-p200b)/(p200-p200b)
  coef200b = (200.-p200 )/(p200b-p200)
  h200_xy(i,j) = h200_xy(i,j) + dtfactor*(z(k200-1)*coef200b + z(k200)*coef200)

  ! 850 mbar horizontal winds
  u850_xy(i,j) = u850_xy(i,j) + u(i,j,k850)*dtfactor  
  v850_xy(i,j) = v850_xy(i,j) + v(i,j,k850)*dtfactor  

  ! Compute coefficients to interpolate values at 850 millibar
  p850  = pres(k850)   + 0.01*p(i,j,k850)   ! Convert p(i,j,k) from Pa to mbar
  p850b = pres(k850-1) + 0.01*p(i,j,k850-1)
  coef850  = (850.-p850b)/(p850-p850b)
  coef850b = (850.-p850 )/(p850b-p850)

  ! Interpolate height, temperature, potential temp, water vapor at 850 mbar.
  h850  = z(k850-1)*coef850b + z(k850)*coef850
  t850  = tabs(i,j,k850-1)*coef850b + tabs(i,j,k850)*coef850
  th850 = tabs(i,j,k850-1)*prespot(k850-1)*coef850b &
         + tabs(i,j,k850)*prespot(k850)*coef850
  qv850 = (q(i,j,k850-1)-qn(i,j,k850-1))*coef850b &
         + (q(i,j,k850) -qn(i,j,k850))*coef850
  ! Keep running-average of 850 mbar height, theta_e and theta_es.
  h850_xy(i,j) = h850_xy(i,j) + dtfactor*h850
  the850_xy(i,j) = the850_xy(i,j) + dtfactor*th850*(1. + fac_cond*qv850/t850)
  thes850_xy(i,j) = the850_xy(i,j) &
       + dtfactor*th850*(1. + fac_cond*qsatw(t850,850.)/t850)
  
  ! Compute coefficients to interpolate values at 1000 millibar
  p1000  = pres(k1000)   + 0.01*p(i,j,k1000) ! Convert p(i,j,k) from Pa to mbar
  p1000b = pres(k1000-1) + 0.01*p(i,j,k1000-1)
  coef1000  = (1000.-p1000b)/(p1000-p1000b)
  coef1000b = (1000.-p1000 )/(p1000b-p1000)
  ! Interpolate height, temperature, potential temp, water vapor at 1000 mbar.
  t1000  = tabs(i,j,k1000-1)*coef1000b + tabs(i,j,k1000)*coef1000
  th1000 = tabs(i,j,k1000-1)*prespot(k1000-1)*coef1000b &
         + tabs(i,j,k1000)*prespot(k1000)*coef1000
  qv1000 = (q(i,j,k1000-1)-qn(i,j,k1000-1))*coef1000b &
         + (q(i,j,k1000) -qn(i,j,k1000))*coef1000
  ! Keep running-average of 1000 mbar height, theta_e and theta_es.
  the1000_xy(i,j) = the1000_xy(i,j) &
       + dtfactor*th1000*(1. + fac_cond*qv1000/t1000)
  thes1000_xy(i,j) = the1000_xy(i,j) &
       + dtfactor*th1000*(1. + fac_cond*qsatw(t1000,1000.)/t1000)
 end do
end do

! COMPUTE CLOUD/ECHO HEIGHTS AS WELL AS CLOUD TOP TEMPERATURE
! WHERE CLOUD TOP IS DEFINED AS THE HIGHEST MODEL LEVEL WITH A
! CONDENSATE PATH OF 0.01 kg/m2 ABOVE.  ECHO TOP IS THE HIGHEST LEVEL
! WHERE THE PRECIPITATE MIXING RATIO > 0.001 G/KG.
if(docloudechoheights) then
   ! initially, zero out heights and set cloudtoptemp to SST
   cloudtopheight = 0.
   cloudtoptemp = sstxy(1:nx,1:ny)
   echotopheight = 0.
   do j = 1,ny
      do i = 1,nx
         ! FIND CLOUD TOP HEIGHT
         tmp_lwp = 0.
         do k = nzm,1,-1
            tmp_lwp = tmp_lwp + qn(i,j,k)*rho(k)*dz*adz(k)
            if (tmp_lwp.gt.0.01) then
               cloudtopheight(i,j) = z(k)
               cloudtoptemp(i,j) = tabs(i,j,k)
               EXIT
            end if
         end do
         ! FIND ECHO TOP HEIGHT
         do k = nzm,1,-1
            if (qp(i,j,k).gt.1.e-6) then
               echotopheight(i,j) = z(k)
               EXIT
            end if
         end do
      end do
   end do
end if ! docloudechoheights


! COMPUTE STORAGE AND ADVECTION TERMS IN COLUMN-INTEGRATED BUDGETS
! OF FROZEN MOIST STATIC ENERGY AND LIQUID WATER ICE STATIC ENERGY.
if (docolumnbudgets) then 

   if (((nstep.eq.1).or.(nstep.eq.nsave2Dstart-nsave2D+1)) &
        .and.(icycle.eq.1)) then

      ! Compute the initial value of the frozen moist static energy
      ! and liquid-ice static energy.
      do j = 1,ny
         do i = 1,nx
            old_fmse(i,j) = 0.
            old_sli(i,j)  = 0.
         end do
      end do

      do k = 1,nzm
         coef2 = rho(k)*dz*adz(k)
         do j = 1,ny
            do i = 1,nx
               old_fmse(i,j) = old_fmse(i,j) &
                    + (cp*t(i,j,k) + lcond*(q(i,j,k) + qp(i,j,k)))*coef2
               old_sli(i,j)  = old_sli(i,j) + cp*t(i,j,k)*coef2
            end do
         end do
      end do
      
   elseif ((mod(nstep,nsave2D).eq.0).and.(icycle.eq.ncycle)) then

      ! Compute the final value of the frozen moist static energy
      ! and liquid-ice static energy.
      do j = 1,ny
         do i = 1,nx
            hstor_xy(i,j) = old_fmse(i,j)
            sstor_xy(i,j) = old_sli(i,j)

            old_fmse(i,j) = 0.
            old_sli(i,j)  = 0.
         end do
      end do

      do k = 1,nzm
         coef2 = rho(k)*dz*adz(k)
         do j = 1,ny
            do i = 1,nx
               old_fmse(i,j) = old_fmse(i,j) &
                    + (cp*t(i,j,k) + lcond*(q(i,j,k) + qp(i,j,k)))*coef2
               old_sli(i,j)  = old_sli(i,j) + cp*t(i,j,k)*coef2
            end do
         end do
      end do
      
      tmp(1) = 0.
      tmp(2) = 0.
      do j = 1,ny
         do i = 1,nx
            ! Compute storage term in budget
            hstor_xy(i,j) = &
                 (old_fmse(i,j) - hstor_xy(i,j))/(dt*float(nsave2D))
            sstor_xy(i,j) = &
                 (old_sli(i,j)  - sstor_xy(i,j))/(dt*float(nsave2D))

            ! Compute advection term in budget (scaled up by nsave2D).
            ! Note that the advection term is not computed direct,
            ! but only as a residual of the rest of the budget (which
            ! should ideally be in balance down to rounding error).
            hadv_xy(i,j) = hstor_xy(i,j) &
                 - (lwns_xy(i,j) + swntm_xy(i,j) &
                 - lwnt_xy(i,j) - swns_xy(i,j) &
                 + rhow(1)*cp*shf_xy(i,j) &
                 + rhow(1)*lcond*lhf_xy(i,j))/float(nsave2D)
            sadv_xy(i,j) = sstor_xy(i,j) &
                 - (lwns_xy(i,j) + swntm_xy(i,j) &
                 - lwnt_xy(i,j) - swns_xy(i,j) &
                 + rhow(1)*cp*shf_xy(i,j) &
                 + (lcond*dz/dt)*prec_xy(i,j))/float(nsave2D)

            tmp(1) = tmp(1) + hadv_xy(i,j)
            tmp(2) = tmp(2) + sadv_xy(i,j)
         end do
      end do

      if (dompi) then
         call task_sum_real(tmp/float(nx_gl*ny_gl),tmp2,2)
      else
         tmp2 = tmp/float(nx_gl*ny_gl)
      end if
      if ((masterproc).and.(.not.dowallx).and.(.not.dowally)) then
         ! Each advection term is computed as the residual of the
         ! rest of the budget terms.  If the domain is periodic, then
         ! the advection term should sum to zero (or be roughly the
         ! size of rounding error).  Any non-zero sum larger than
         ! rounding error indicates an imbalance in the budget.  Note
         ! that with individual terms in the budget on the order of
         ! 1000 W/m^2, the rounding error should be roughly 1e-4 to
         ! 1e-5 with single precision.
         write(*,*) 'FMSE Surface Budget Imbalance = ', tmp2(1), ' W/m^2'
         write(*,*) 'S_li Surface Budget Imbalance = ', tmp2(2), ' W/m^2'
      end if

   end if 
end if ! docolumnbudgets

! COMPUTE PROJECTION OF VARIOUS FIELDS ONTO FIRST/SECOND BAROCLINIC MODE
! SENSE OF PROJECTION IS SUCH THAT POSITIVE VALUES CORRESPOND TO
! POSITIVE ANOMALIES IN THE LOWER TROPOSPHERE.
if (dowaveoutput) then
   do j = 1,ny
      do i = 1,nx
         wmode1i_xy(i,j) = 0.
         wmode2i_xy(i,j) = 0.
      end do
   end do

   press_total = (presi(1)-presi(k150+1))
   press_top   = presi(k150+1)
   do k = 1,k150
      tmp_wgt1 = sin(pii*(pres(k)-press_top)/press_total) &
           *(presi(k)-presi(k+1))/press_total*dtfactor
      tmp_wgt2 = sin(-2.*pii*(pres(k)-press_top)/press_total) &
           *(presi(k)-presi(k+1))/press_total*dtfactor
      do j = 1,ny
         do i = 1,nx
            wmode1_xy(i,j) = wmode1_xy(i,j) + tmp_wgt1*w(i,j,k)
            wmode2_xy(i,j) = wmode2_xy(i,j) + tmp_wgt2*w(i,j,k)
            qmode1_xy(i,j) = qmode1_xy(i,j) + tmp_wgt1*(q(i,j,k)-q0(k))
            qmode2_xy(i,j) = qmode2_xy(i,j) + tmp_wgt2*(q(i,j,k)-q0(k))
            thmode1_xy(i,j) = thmode1_xy(i,j) + tmp_wgt1*(t(i,j,k)-t0(k))
            thmode2_xy(i,j) = thmode2_xy(i,j) + tmp_wgt2*(t(i,j,k)-t0(k))
            wmode1i_xy(i,j) = wmode1i_xy(i,j)  + tmp_wgt1*w(i,j,k)
            wmode2i_xy(i,j) = wmode2i_xy(i,j) + tmp_wgt2*w(i,j,k)
         end do
      end do
   end do
end if

! COMPUTE WAVE ENERGETICS

if (dowaveenergetics) then
   tmpw1_xy = 0.
   tmpw2_xy = 0.
   tmps1_xy = 0.
   tmps2_xy = 0.
   tmph1_xy = 0.
   tmph2_xy = 0.

   press_total = (presi(1)-presi(k150+1))
   do k = 1,k150
      tmp_wgt1 = 2.*sin(pii*(pres(k)-press_top)/press_total) &
           *(presi(k)-presi(k+1))/press_total
      tmp_wgt2 = 2.*sin(-2.*pii*(pres(k)-press_top)/press_total) &
           *(presi(k)-presi(k+1))/press_total
      do j = 1,ny
         do i = 1,nx
            tmpw1_xy(i,j) = tmpw1_xy(i,j) + tmp_wgt1*w(i,j,k)
            tmpw2_xy(i,j) = tmpw2_xy(i,j) + tmp_wgt2*w(i,j,k)
            tmps1_xy(i,j) = tmps1_xy(i,j) + tmp_wgt1*(t(i,j,k)-t0(k))
            tmps2_xy(i,j) = tmps2_xy(i,j) + tmp_wgt2*(t(i,j,k)-t0(k))
            tmph1_xy(i,j) = tmph1_xy(i,j) &
                 + tmp_wgt1*((t(i,j,k)-t0(k)) &
                             + fac_cond*(q(i,j,k)+qp(i,j,k)-q0(k)))
            tmph2_xy(i,j) = tmph2_xy(i,j) &
                 + tmp_wgt2*((t(i,j,k)-t0(k)) &
                             + fac_cond*(q(i,j,k)+qp(i,j,k)-q0(k)))
         end do
      end do
   end do

   wdwdt1_xy = wdwdt1_xy + dtfactor*tmpw1_xy*(tmpw1_xy - oldw1_xy)/dtn
   wdwdt2_xy = wdwdt2_xy + dtfactor*tmpw2_xy*(tmpw2_xy - oldw2_xy)/dtn

   sdsdt1_xy = sdsdt1_xy + dtfactor*tmps1_xy*(tmps1_xy - olds1_xy)/dtn
   sdsdt2_xy = sdsdt2_xy + dtfactor*tmps2_xy*(tmps2_xy - olds2_xy)/dtn

   hdhdt1_xy = hdhdt1_xy + dtfactor*tmph1_xy*(tmph1_xy - oldh1_xy)/dtn
   hdhdt2_xy = hdhdt2_xy + dtfactor*tmph2_xy*(tmph2_xy - oldh2_xy)/dtn

   oldw1_xy = tmpw1_xy
   oldw2_xy = tmpw2_xy
   olds1_xy = tmps1_xy
   olds2_xy = tmps2_xy
   oldh1_xy = tmph1_xy
   oldh2_xy = tmph2_xy

end if


if(douwpbl) then

   do i=1,nx
      do j=1,ny
         pblh_xy(i,j) = pblh_xy(i,j) + dtfactor*pblh(i,j)
      end do
   end do

end if

if(douwcu) then

   do i=1,nx
      do j=1,ny
         cin_xy(i,j) = cin_xy(i,j) + dtfactor*cin(i,j)
         cbmf_xy(i,j) = cbmf_xy(i,j) + dtfactor*cbmf(i,j)
      end do
   end do

end if

! END UW ADDITIONS
!=====================================================

!   recompute pressure levels:

call pressz()

total_water_after = 0.
do k=1,nzm
 total_water_after = total_water_after + &
            (sum(q(1:nx,1:ny,k))+sum(qp(1:nx,1:ny,k)))*adz(k)*dz *rho(k)
end do
	
end subroutine diagnose
