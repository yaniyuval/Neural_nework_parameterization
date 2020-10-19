subroutine precip_fall_krueger

  !     positively definite monotonic advection with non-oscillatory option
  !     and gravitational sedimentation (rain and snow advection)

  use vars
  use params
  use mse, only : prec_inst_mse, doMSE, prec_inst_frozen_mse  ! peters
  implicit none

  real, external :: vtermr, vtermg, vterms ! functions for computing fall speed

  real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
  real f0(nzm),df0(nzm)
  integer i,j,k,kc,kb

  real lat_heat, wmax, dq
  real, dimension(nzm) :: iwmax, rhofac, deltaz, dummy
  real, dimension(nzm) :: rqr, qr_initial               
  real, dimension(nzm) :: rqs, qs_initial
  real, dimension(nzm) :: rqg, qg_initial
  real, dimension(nz)  :: rflux, sflux, gflux
  logical :: dorainfall, dograufall, dosnowfall

! initialize constants for fall speed formulas -- taken from Kreuger's IMICRO

  real, parameter :: alin = 841.99667
  real, parameter :: clin = 4.836071224
  real, parameter :: rhoref = 1.226

  real, parameter :: zero = 0.
  real :: pie, gammafff

  !--------------------------------------------------------

  pie = atan2(0.,-1.)
  vconr = gamr3*alin /(6.*(pie*nzeror*rhor)**vexpr)
  vcons = gams3*clin /(6.*(pie*nzeros*rhos)**vexps)
  vcong = gamg3*40.74/(6.*(pie*nzerog*rhog)**vexpg)
  

  if(dobetafactor) dtn=dtn*betafactor**0.5

  ! Initialize arrays that accumulate surface precipitation flux

  if(doMSE) prec_inst_mse = 0.   !peters, need to accumulate over iprec steps

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
     qrfall(k)=0. 
     qgfall(k)=0. 
     qsfall(k)=0. 
     tlat(k) = 0.
  end do

  do k = 1,nzm
     deltaz(k) = adz(k)*dz
     ! Factor in formula for terminal velocity of precipitation
     rhofac(k) = sqrt(rhoref/rho(k)) ! Krueger uses rho_sfc = 1.226
     dummy(k) = 0. ! dummy argument for fall speed formulas
                   ! extra argument is included for T-dependent fall speeds.
  end do

  ! LOOP OVER VERTICAL COLUMNS
  do j=1,ny
     do i=1,nx

        ! =========== RAIN ===========
        dorainfall = .false.
        do k = 1,nzm
           qr_initial(k) = qr(i,j,k) ! SAVE INITIAL MIXTURE FRACTION PROFILE
           rqr(k) = rho(k)*qr(i,j,k) ! COMPUTE PROFILE OF RHO*QR
           if (qr(i,j,k).gt.zero) dorainfall = .true.
        end do
        rflux = 0.

        ! ONLY COMPUTE PRECIPITATION IF THERE IS RAIN IN THIS COLUMN
        if (dorainfall) then
           ! COMPUTE PRECIPITATION, TAKING MULTIPLE SUBSTEPS IF CFL > 0.9
           call precip_iter(vtermr,nzm,dtn,deltaz,rhofac,dummy,rqr,rflux)
           do k = 1,nzm
              qr(i,j,k) = max(rqr(k)/rho(k),0.) ! UPDATE MIXTURE FRACTION
              dq = qr(i,j,k) - qr_initial(k)
              t(i,j,k) = t(i,j,k) - fac_cond*dq !LIQUID WATER-ICE STATIC ENERGY
              qrfall(k) = qrfall(k) + dq        !moisture bugdet for statistics
              tlat(k)   = tlat(k)   - fac_cond*dq !energy budget for statistics
              precflux(k) = precflux(k) - rflux(k)*dtn/dz ! stats: precip flux
           end do
           precsfc(i,j) = precsfc(i,j) - rflux(1)*dtn/dz ! For statistics
           prec_xy(i,j) = prec_xy(i,j) - rflux(1)*dtn/dz ! For 2D output
           prec_inst(i,j) = - rflux(1)
           if(doMSE) prec_inst_mse(i,j) = &
                prec_inst_mse(i,j) - rflux(1)*dtn/dz ! peters
        end if

        ! =========== GRAUPEL ===========
        dograufall = .false.
        do k = 1,nzm
           qg_initial(k) = qg(i,j,k)
           rqg(k) = rho(k)*qg(i,j,k)
           if (qg(i,j,k).gt.zero) dograufall = .true.
        end do
        gflux = 0.

        ! Compute precipitation velocity and flux column-by-column
        if (dograufall) then
           ! compute precipitation, taking multiple substeps if cfl > 0.9
           call precip_iter(vtermg,nzm,dtn,deltaz,rhofac,dummy,rqg,gflux)
           do k = 1,nzm
              qg(i,j,k) = max(rqg(k)/rho(k),0.) ! update mixture fraction
              dq = qg(i,j,k) - qg_initial(k)
              t(i,j,k) = t(i,j,k) - fac_sub*dq  !LIQUID WATER-ICE STATIC ENERGY
              qgfall(k) = qgfall(k) + dq        !stats: moisture bugdet
              tlat(k)   = tlat(k)   - fac_sub*dq!stats: energy budget
              precflux(k) = precflux(k) - gflux(k)*dtn/dz !stats: precip flux
           end do
           precsfc(i,j) = precsfc(i,j) - gflux(1)*dtn/dz ! For statistics
           prec_xy(i,j) = prec_xy(i,j) - gflux(1)*dtn/dz ! For 2D output
           prec_inst(i,j) = prec_inst(i,j) - gflux(1)
           if(doMSE) prec_inst_frozen_mse(i,j) = &
                prec_inst_frozen_mse(i,j) - gflux(1)*dtn/dz ! peters
        end if


        ! =========== SNOW ===========
        dosnowfall = .false.
        do k = 1,nzm
           qs_initial(k) = qs(i,j,k)
           rqs(k) = rho(k)*qs(i,j,k)
           if (qs(i,j,k).gt.zero) dosnowfall = .true.
        end do
        sflux = 0.

        ! Compute precipitation velocity and flux column-by-column
        if (dosnowfall) then
           ! compute precipitation, taking multiple substeps if cfl > 0.9
           call precip_iter(vterms,nzm,dtn,deltaz,rhofac,dummy,rqs,sflux)
           do k = 1,nzm
              qs(i,j,k) = max(rqs(k)/rho(k),0.) ! update mixture fraction
              dq = qs(i,j,k) - qs_initial(k)
              t(i,j,k) = t(i,j,k) - fac_sub*dq  !LIQUID WATER-ICE STATIC ENERGY
              qsfall(k) = qsfall(k) + dq        !stats: moisture bugdet
              tlat(k)   = tlat(k)   - fac_sub*dq!stats: energy budget
              precflux(k) = precflux(k) - sflux(k)*dtn/dz !stats: precip flux
           end do
           precsfc(i,j) = precsfc(i,j) - sflux(1)*dtn/dz ! For statistics
           prec_xy(i,j) = prec_xy(i,j) - sflux(1)*dtn/dz ! For 2D output
           prec_inst(i,j) = prec_inst(i,j) - sflux(1)
           if(doMSE) prec_inst_frozen_mse(i,j) = &
                prec_inst_frozen_mse(i,j) - sflux(1)*dtn/dz ! peters
        end if

     end do
  end do

  if(dobetafactor) dtn=dtn/betafactor**0.5

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
        if(qr(i,j,1)+qg(i,j,1)+qs(i,j,1).gt.1.e-6) s_ar=s_ar+dtfactor
     end do
  end do


  if(dostatis) then

     call stat_varscalar(t,df,f0,df0,t2leprec)
     call setvalue(twleprec,nzm,0.)
     call stat_sw2(t,df,twleprec)

  endif

end subroutine precip_fall_krueger

real function vtermr(rq,rhofac,dummy)
  use params, only: vconr, vexpr, vmaxr
  implicit none

  real, intent(in) :: rq, rhofac, dummy

  vtermr = min(vmaxr,vconr*rhofac*max(rq,0.0)**vexpr)
end function vtermr

real function vtermg(rq,rhofac,dummy)
  use params, only: vcong, vexpg, vmaxg
  implicit none

  real, intent(in) :: rq, rhofac, dummy

  vtermg = min(vmaxg,vcong*rhofac*max(rq,0.0)**vexpg)
end function vtermg

real function vterms(rq,rhofac,dummy)
  use params, only: vcons, vexps
  implicit none

  real, intent(in) :: rq, rhofac, dummy

  vterms = vcons*rhofac*max(rq,0.0)**vexps
end function vterms



