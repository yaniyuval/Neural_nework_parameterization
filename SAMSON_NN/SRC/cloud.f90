
subroutine cloud

!  Condensation of cloud water/cloud ice.

use vars
use params
use microscaling!, only: domicroscaling, dtn_scaled !Yani commented because it didn't finish the Makefile with this part

implicit none

integer i,j,k, kb, kc, kmax, kmin
real dtabs, tabs1, an, bn, ap, bp, om, ag, omp
real fac1,fac2  
real fff,dfff,qsat,dqsat
real lstarn,dlstarn,lstarp,dlstarp
real coef,dqi,lat_heat,vt_ice,qn0,coef_cl,vt_cl,fluxi,fluxc
real omnu, omnc, omnd, qiu, qic, qid, tmp_theta, tmp_phi
real fz(nx,ny,nz), fzt(nx,ny,nz)
real dtn_old
integer niter

if(dokruegermicro) then
   return ! saturation adjustment takes place in precip_proc_krueger()
end if

dtn_old = dtn
if(dobetafactor.and.daremicrophysics.and..not.domicroscaling) dtn=dtn*mpbetafactor**0.5

an = 1./(tbgmax-tbgmin)	
bn = tbgmin * an
ap = 1./(tprmax-tprmin)	
bp = tprmin * ap
fac1 = fac_cond+(1+bp)*fac_fus
fac2 = fac_fus*ap
ag = 1./(tgrmax-tgrmin)	

kmax=0
kmin=nzm+1

do k = 1,nzm
   qvcond(k) = 0.
   tlcond(k) = 0.
end do


do k = 1, nzm
 do j = 1, ny
  do i = 1, nx

     if(domicroscaling) dtn = dtn_scaled(i, j, k)

    qn0 = qn(i,j,k)

    q(i,j,k)=max(0.,q(i,j,k))


! Initail guess for temperature assuming no cloud water/ice:


    tabs(i,j,k) = t(i,j,k)-gamaz(k)
    tabs1=(tabs(i,j,k)+fac1*qp(i,j,k))/(1.+fac2*qp(i,j,k))

! Warm cloud:

    if(tabs1.ge.tbgmax) then

      tabs1=tabs(i,j,k)+fac_cond*qp(i,j,k)
      qsat = qsatw(tabs1,pres(k))

! Ice cloud:

    elseif(tabs1.le.tbgmin) then

      tabs1=tabs(i,j,k)+fac_sub*qp(i,j,k)
      qsat = qsati(tabs1,pres(k))

! Mixed-phase cloud:

    else

      om = an*tabs1-bn
      qsat = om*qsatw(tabs1,pres(k))+(1.-om)*qsati(tabs1,pres(k))

    endif


!  Test if condensation is possible:


    if(q(i,j,k) .gt. qsat) then

      niter=0
      dtabs = 100.
      do while(abs(dtabs).gt.0.01.and.niter.lt.10)
	if(tabs1.ge.tbgmax) then
	   om=1.
	   lstarn=fac_cond
	   dlstarn=0.
	   qsat=qsatw(tabs1,pres(k))
	   dqsat=dtqsatw(tabs1,pres(k))
        else if(tabs1.le.tbgmin) then
	   om=0.
	   lstarn=fac_sub
	   dlstarn=0.
	   qsat=qsati(tabs1,pres(k))
	   dqsat=dtqsati(tabs1,pres(k))
	else
	   om=an*tabs1-bn
	   lstarn=fac_cond+(1.-om)*fac_fus
	   dlstarn=an
	   qsat=om*qsatw(tabs1,pres(k))+(1.-om)*qsati(tabs1,pres(k))
	   dqsat=om*dtqsatw(tabs1,pres(k))+(1.-om)*dtqsati(tabs1,pres(k))
	endif
	if(tabs1.ge.tprmax) then
	   omp=1.
	   lstarp=fac_cond
	   dlstarp=0.
        else if(tabs1.le.tprmin) then
	   omp=0.
	   lstarp=fac_sub
	   dlstarp=0.
	else
	   omp=ap*tabs1-bp
	   lstarp=fac_cond+(1.-omp)*fac_fus
	   dlstarp=ap
	endif
	fff = tabs(i,j,k)-tabs1+lstarn*(q(i,j,k)-qsat)+lstarp*qp(i,j,k)
	dfff=dlstarn*(q(i,j,k)-qsat)+dlstarp*qp(i,j,k)-lstarn*dqsat-1.
	dtabs=-fff/dfff
	niter=niter+1
	tabs1=tabs1+dtabs
      end do   

      qsat = qsat + dqsat * dtabs
      qn(i,j,k) = max(0.,q(i,j,k)-qsat)
!bloss      if(tabs1.lt.tbgmax) then
      if(qn(i,j,k).gt.qp_threshold) then
        kmin = min(kmin,k)
        kmax = max(kmax,k)
      end if

    else

      qn(i,j,k) = 0.

    endif

    tabs(i,j,k) = tabs1
    qp(i,j,k) = max(0.,qp(i,j,k)) ! just in case

    qvcond(k) = qvcond(k) - qn(i,j,k) + qn0
    tlcond(k) = tlcond(k) + (fac_cond + fac_fus*(1.-omegan(tabs(i,j,k)))) &
         *(qn(i,j,k)-qn0)

  end do
 end do
end do


! Sedimentation of ice and water:

do k = 1,nzm
   qifall(k) = 0.
   tlatqi(k) = 0.
end do

fz = 0.
fzt = 0.

!
! Take into account sedimentation of cloud water which may be
! important for stratocumulus case. 
! Parameterization of sedimentation rate is taken from GCSS WG1
! DYCOMS2_RF2 case, and based on 
! Rogers and Yau, 1989

coef_cl = 1.19e8*(3./(4.*3.1415*1000.*Nc0*1.e6))**(2./3.)*exp(5.*log(1.5)**2)

! Compute cloud ice flux (using flux limited advection scheme, as in
! chapter 6 of Finite Volume Methods for Hyperbolic Problems by R.J.
! LeVeque, Cambridge University Press, 2002). 
do k = max(1,kmin-1),kmax
   ! Set up indices for x-y planes above and below current plane.
   kc = min(nzm,k+1)
   kb = max(1,k-1)
   ! CFL number based on grid spacing interpolated to interface i,j,k-1/2
!   coef = dtn/(0.5*(adz(kb)+adz(k))*dz)
   do j = 1,ny
      do i = 1,nx

         if(domicroscaling) dtn=dtn_scaled(i,j,k)
         coef = dtn/(0.5*(adz(kb)+adz(k))*dz)

         ! Compute cloud ice density in this cell and the ones above/below.
         ! Since cloud ice is falling, the above cell is u (upwind),
         ! this cell is c (center) and the one below is d (downwind). 
         omnu = max(0.,min(1.,(tabs(i,j,kc)-tbgmin)*a_bg))
         omnc = max(0.,min(1.,(tabs(i,j,k) -tbgmin)*a_bg))
         omnd = max(0.,min(1.,(tabs(i,j,kb)-tbgmin)*a_bg))         

         qiu = rho(kc)*qn(i,j,kc)*(1.-omnu)
         qic = rho(k) *qn(i,j,k) *(1.-omnc) 
         qid = rho(kb)*qn(i,j,kb)*(1.-omnd) 

         ! Ice sedimentation velocity depends on ice content. The fiting is
         ! based on the data by Heymsfield (JAS,2003). -Marat
         ! 0.1 m/s low bound was suggested by Chris Bretherton 
         vt_ice = max(0.1,0.5*log10(qic+1.e-12)+3.)

         ! Use MC flux limiter in computation of flux correction.
         ! (MC = monotonized centered difference).
         if (qic.eq.qid) then
            tmp_phi = 0.
         else
            tmp_theta = (qiu-qic)/(qic-qid)
            tmp_phi = max(0.,min(0.5*(1.+tmp_theta),2.,2.*tmp_theta))
         end if

         ! Compute limited flux.
         ! Since falling cloud ice is a 1D advection problem, this
         ! flux-limited advection scheme is monotonic.
         fluxi = -vt_ice*(qic - 0.5*(1.-coef*vt_ice)*tmp_phi*(qic-qid))

         if(doclouddropsed) then
            ! Compute cloud water density in this cell and the ones above/below
            ! Since cloud water is falling, the above cell is u (upwind),
            ! this cell is c (center) and the one below is d (downwind).
            qiu = rho(kc)*qn(i,j,kc)*omnu
            qic = rho(k) *qn(i,j,k) *omnc
            qid = rho(kb)*qn(i,j,kb)*omnd

            vt_cl = coef_cl*(qic+1.e-12)**(2./3.)

            ! Use MC flux limiter in computation of flux correction.
            ! (MC = monotonized centered difference).
            if (qic.eq.qid) then
               tmp_phi = 0.
            else
               tmp_theta = (qiu-qic)/(qic-qid)
               tmp_phi = max(0.,min(0.5*(1.+tmp_theta),2.,2.*tmp_theta))
            end if

            ! Compute limited flux.
            ! Since falling cloud water is a 1D advection problem, this
            ! flux-limited advection scheme is monotonic.
            fluxc = -vt_cl*(qic - 0.5*(1.-coef*vt_cl)*tmp_phi*(qic-qid))
         else
            fluxc = 0.
         end if

         fz(i,j,k) = fluxi + fluxc
         fzt(i,j,k) = -(fac_cond+fac_fus)*fluxi - fac_cond*fluxc
      end do
   end do
end do
fz(:,:,nz) = 0.

do k=max(1,kmin-2),kmax
!   coef=dtn/(dz*adz(k)*rho(k))
   do j=1,ny
      do i=1,nx

	if(domicroscaling) dtn=dtn_scaled(i,j,k)
	   coef=dtn/(dz*adz(k)*rho(k))

         ! The cloud ice increment is the difference of the fluxes.
         dqi=coef*(fz(i,j,k)-fz(i,j,k+1))
         ! Add this increment to both non-precipitating and total water.
         if (.not.do_sedimentation_rf) then    
          qn(i,j,k) = qn(i,j,k) + dqi
	  q(i,j,k)  = q(i,j,k)  + dqi
         end if
         ! Include this effect in the total moisture budget.
         qifall(k) = qifall(k) + dqi

         ! The latent heat flux induced by the falling cloud ice enters
         ! the liquid-ice static energy budget in the same way as the
         ! precipitation.  Note: use latent heat of sublimation. 
         lat_heat  = coef*(fzt(i,j,k)-fzt(i,j,k+1))
         ! Add divergence of latent heat flux to liquid-ice static energy.
         if (.not.do_sedimentation_rf) then    
            t(i,j,k)  = t(i,j,k)  + lat_heat
         end if
         ! Add divergence to liquid-ice static energy budget.
         tlatqi(k) = tlatqi(k) + lat_heat
      end do
   end do
end do

! peters.
! NOTE:  falling cloud water/ice through the bottom of the model is allowed
!  by advection scheme.  This is a source of surface precip but is not
!  included in MSE accumulated precip and budgets here.  It is hoped that the
!  error will be small since it requires a cloud a ground level, but
!  perhaps it will need to be taken into account

do j = 1,ny
   do i = 1,nx
	if(domicroscaling) dtn=dtn_scaled(i,j,1)
      if (.not.do_sedimentation_rf) then !JY - added not to count twice
       precsfc(i,j) = precsfc(i,j) - fz(i,j,1)*dtn/dz ! For statistics
       prec_xy(i,j) = prec_xy(i,j) - fz(i,j,1)*dtn/dz ! For 2D output
       prec_inst(i,j) = - fz(i,j,1)
      end if
   end do
end do

dtn = dtn_old

end subroutine cloud

