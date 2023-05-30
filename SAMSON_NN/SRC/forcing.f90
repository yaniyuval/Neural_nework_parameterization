
subroutine forcing

  use vars
  use params
  use simple_ocean, only: sst_evolve, doevolveSSTforcing
  use simple_land                 ! peters


  implicit none

  integer i,j,k,k1,k2
  real coef, radtend, wrk1,wrk2,tmp1,tmp2,first_mode,second_mode,press_top,press_total
  real radfrac(nz), vadvtend !bloss

  do k=1,nzm
     total_water_ls = total_water_ls - &
          sum(q(1:nx,1:ny,k))*adz(k)*dz *rho(k)
  end do


  ! Large-scale sounding:

  k1=1
  k2=1
  do i=1,nsnd-1
     if(day.gt.daysnd(i)) then
        k2=i+1
        k1=i
     endif
  end do
  coef=(day-daysnd(k1))/(daysnd(k2)-daysnd(k1))
  do k=1,nzm
     tg0(k)=tsnd(k,k1)+(tsnd(k,k2)-tsnd(k,k1))*coef
     qg0(k)=qsnd(k,k1)+(qsnd(k,k2)-qsnd(k,k1))*coef
     qg0(k)=qg0(k)*1.e-3
     ! Note that ug0 and vg0 maybe reset if dolargescale is true)
     ug0(k)=usnd(k,k1)+(usnd(k,k2)-usnd(k,k1))*coef - ug
     vg0(k)=vsnd(k,k1)+(vsnd(k,k2)-vsnd(k,k1))*coef - vg
  end do

  ! Initialize tendencies:

  do k=1,nzm
     ttend(k)=0.
     qtend(k)=0.
  end do
  
  ! Large-Scale Advection Forcing:

  if(dolargescale.and.time.gt.timelargescale) then

     k1=1
     k2=1
     do i=1,nlsf-1
        if(day.gt.dayls(i)) then
           k2=i+1
           k1=i
        endif
     end do
     coef=(day-dayls(k1))/(dayls(k2)-dayls(k1))

     do k=1,nzm
        ttend(k)=dtls(k,k1)+(dtls(k,k2)-dtls(k,k1))*coef
        qtend(k)=(dqls(k,k1)+(dqls(k,k2)-dqls(k,k1))*coef)*qtendfactor
        ug0(k)=ugls(k,k1)+(ugls(k,k2)-ugls(k,k1))*coef - ug
        vg0(k)=vgls(k,k1)+(vgls(k,k2)-vgls(k,k1))*coef - vg
     end do

     if(.not.douwforcing) then !bloss: use an alternative method to
        !        specify large-scale forcing
        do k=1,nzm ! apply large-scale temperature/moisture tendencies
           do j=1,ny
              do i=1,nx
                 t(i,j,k)=t(i,j,k)+ttend(k) * dtn
                 q(i,j,k)=q(i,j,k)+qtend(k) * dtn 
              end do
           end do
        end do
           do k=1,nzm
              wsub(k)=wgls(k,k1)+(wgls(k,k2)-wgls(k,k1))*coef
           end do
        if(dosubsidence) call subsidence()

     else !bloss: use forcings appropriate to boundary layer simulations

        if(doblradforcing) then
           !bloss: apply fixed radiative cooling for q > blqmax
           !         and taper to zero where q < blqmin

           do k=1,nzm
              ttend(k) = 0.
              do j=1,ny
                 do i=1,nx
                    radtend = blqrad &
                         *max(0.,min(1.,(q(i,j,k)-blqmin)/(blqmax-blqmin)))
                    t(i,j,k)=t(i,j,k) + radtend * dtn
                    ttend(k) = ttend(k) + radtend
                 end do
              end do
           end do
           ttend = ttend/float(nx*ny) ! average ttend on this processor

           ! next, compute subsidence profile with fixed divergence rate.
           ! however, apply this subsidence in the same proportion as the
           ! radiative cooling above, so that they will approximately balance.

           ! compute fraction of maximum possible radiative cooling on each plane.
           radfrac = 0.
           radfrac(1:nzm) = ttend/blqrad/float(nsubdomains_x*nsubdomains_y)
           if (dompi) then
              call task_sum_real(radfrac,wsub,nz) ! average across all processors
           else
              wsub = radfrac
           end if
           wsub = -bldivg*zi*wsub ! multiply by divergence*height.

        elseif(doblforcing) then
           !bloss: boundary layer simulations with interactive radiation
           !       large-scale advective forcings only through a
           !       large-scale vertical velocity (based on user-specified
           !       divergence or subsidence) 

           !bloss: COMPUTE INVERSION HEIGHT as location of maximum 
           !        negative gradient in the total water profile. 
           coef = 0.
           do k=2,nzm
              if (-(q0(k)-q0(k-1))/adzw(k)/dz.gt.coef) then
                 k1 = k ! possible inversion height
                 coef = -(q0(k)-q0(k-1))/adzw(k)/dz ! new max neg gradient
              end if
           end do
           ! now, location of inversion is stored in k1. 

           ! SET UP INVERSE NUDGING TIMESCALE.
           donudging_tq = .true. ! make sure nudging is turned on.
           dotauz = .true. ! allow nudging timescale to be fn of height.
           itauz = 0.
           do k = k1+1,nzm
              ! nudging starts 150 meters above inversion.
              if (z(k).gt.z(k1)+150.) itauz(k) = 1./bltauz 
           end do

           ! SET UP SUBSIDENCE VELOCITY:
           !     subsidence velocity maximizes at inversion height,
           !     and nudging is applied starting 150m above the inversion.

           if(doblsubs) then ! user-specified subsidence
              wsub = -blsubs                       ! subsidence above inversion
              wsub(1:k1) = -blsubs*zi(1:k1)/zi(k1) ! assume divg rate const in bl
           elseif(doblsubs2) then 
              wsub = -blsubs2/rhow/ggr*(1.-exp(-1.e-2*(pres0-presi)))
           else ! user-specified divergence
              wsub = -bldivg*zi(k1)         ! subsidence above inversion
              wsub(1:k1) = -bldivg*zi(1:k1) ! assume divg rate const in bl
           end if

        elseif(doaddqvrestart) then
           if(firststep) then
              qg0=q0
              do i=1,nsnd
                 qsnd(:,i)=q0*1000.
              enddo
           endif

        elseif(dowtg) then

           ! apply large-scale temperature/moisture tendencies
           ! these should only include large-scale HORIZONTAL advection!!
           do k=1,nzm 
              do j=1,ny
                 do i=1,nx
                    t(i,j,k)=t(i,j,k)+ttend(k) * dtn
                    q(i,j,k)=q(i,j,k)+qtend(k) * dtn
                 end do
              end do
           end do

!!$         if(donudging_tq) then
!!$            dotauz = .true.
!!$            itauz = 0.
!!$         end if

           !use the current state as the reference profile so there is no initial bias
           if(firststep) then
              tg0=tabs0
              qg0=q0
              do i=1,nsnd
                 tsnd(:,i)=tabs0
                 qsnd(:,i)=q0*1000.
              enddo
           endif

           ! compute the wtg large-scale vertical velocity (kzm implementation)
           wrk1=0.
           wrk2=0.
           first_mode=0.
           second_mode=0.
           press_top=100.
           press_total=800.
           do k = 2,nzm-1 ! don't apply at surface
              !project the T anomaly onto the first mode
              if(pres(k).lt.(press_total+press_top).and.pres(k).gt.press_top) then
                 first_mode=first_mode+(tg0(k)-tabs0(k))*sin(3.1415926*(pres(k)-press_top)/press_total)*(presi(k)-presi(k+1))
                 wrk1=wrk1+(sin(3.1415926*(pres(k)-press_top)/press_total)**2)*(presi(k)-presi(k+1))
                 second_mode=second_mode+(tg0(k)-tabs0(k))*sin(2.*3.1415926*(pres(k)-press_top)/press_total)*(presi(k)-presi(k+1))
                 wrk2=wrk2+(sin(2.*3.1415926*(pres(k)-press_top)/press_total)**2)*(presi(k)-presi(k+1))
              end if
           enddo
           first_mode=first_mode/wrk1
           second_mode=second_mode/wrk2

           ! compute the wtg large-scale vertical velocity
           wsub = 0.
           do k = 2,nzm-1 ! don't apply at surface
              !            if(z(k).gt.1200.) then
              !               vadvtend = (tg0(k)-tabs0(k))/wtgtau - ttend(k)
              !kzm implementation 
              if(pres(k).lt.(press_total+press_top).and.pres(k).gt.press_top) then 
                 vadvtend = first_mode*sin(3.1415926*(pres(k)-press_top)/press_total)/wtgtau - ttend(k)
                 ! w_wtg = vadvtend / (dsli/dz)
                 ! NOTE: we impose a minimum static stability of 1 K/km
                 !        as in Raymond & Zeng, 2005. 
                 if(vadvtend.lt.0.) then ! wsub(k).gt.0.
                    wsub(k) = -vadvtend/max(1.e-3,(t0(k)-t0(k-1))/(dz*adzw(k)))
                 else ! wsub(k).lt.0.
                    wsub(k) = -vadvtend/max(1.e-3,(t0(k+1)-t0(k))/(dz*adzw(k+1)))
                 end if
                 !add w tendency due to the parameterization of the second mode
                 !so that second_mode_factor is 1/(wavelength in 1000km)
                 !note that press_total is in mb
                 wsub(k)=wsub(k)+dtn*1e-9*wavenumber_factor**2*(press_total/rho(k))**2*second_mode*sin(2*3.1415926*(pres(k)-press_top)/press_total)/tg0(k)
              else
!!$               if(dotauz) itauz(k) = 1./tauls
              end if
           end do

!!$         if(masterproc) write(*,*) &
!!$              'Stopping: WTG large-scale forcing not implemented yet'
!!$         call task_abort()

        elseif(doparameterizedwave) then

           ! apply large-scale temperature/moisture tendencies
           ! these should only include large-scale HORIZONTAL advection!!
           do k=1,nzm 
              do j=1,ny
                 do i=1,nx
                    t(i,j,k)=t(i,j,k)+ttend(k) * dtn
                    q(i,j,k)=q(i,j,k)+qtend(k) * dtn
                 end do
              end do
           end do

!!$         if(donudging_tq) then
!!$            dotauz = .true.
!!$            itauz = 0.
!!$         end if

           !use the current state as the reference profile so there is no initial bias
           if(firststep) then
              tg0=tabs0
              qg0=q0
              do i=1,nsnd
                 tsnd(:,i)=tabs0
                 qsnd(:,i)=q0*1000.
              enddo
           endif

           ! compute the wtg large-scale vertical velocity (kzm implementation)
           wrk1=0.
           wrk2=0.
           first_mode=0.
           second_mode=0.
           press_top=150.
           press_total=750.
           do k = 2,nzm-1 ! don't apply at surface
              !project the T anomaly onto the first mode
              tmp1=rho(k)/tg0(k)**3*(tg0(k)+gamaz(k)-tg0(k-1)-gamaz(k-1))/(dz*adzw(k))*(presi(k)-presi(k+1))
              first_mode=first_mode+(tg0(k)-tabs0(k))*tmp1*thmode1(k)
              wrk1=wrk1+thmode1(k)*thmode1(k)*tmp1
              second_mode=second_mode+(tg0(k)-tabs0(k))*tmp1*thmode2(k)
              wrk2=wrk2+thmode2(k)*thmode2(k)*tmp1
           enddo
           first_mode=first_mode/wrk1
           second_mode=second_mode/wrk2
           ! compute the wtg large-scale vertical velocity
           !         wsub = 0.
           do k = 2,nzm-1 ! don't apply at surface
              !            if(z(k).gt.1200.) then
              !               vadvtend = (tg0(k)-tabs0(k))/wtgtau - ttend(k)
              !kzm implementation 
              if(z(k).le.14000.) then 
                 !add w tendency due to the parameterization of the first/second mode
                 !so that second_mode_factor is 1/(wavelength in 1000km)
                 !the first mode factor is simply 2 times the second mode factor
                 !note that press_total is in mb
                 !2.e-3 is an estimate of the average stratification in K/km
                 vadvtend=wavenumber_factor**2*(2*3.1415926/1.e6)**2*(first_mode_speed**2*first_mode*thmode1(k)+second_mode_speed**2*second_mode*thmode2(k))

                 !let's not worry about upwind, so the computation of wsub is strictly linear
                 wsub(k) = wsub(k)-dtn*vadvtend/((tg0(k)+gamaz(k)-tg0(k-1)-gamaz(k-1))/(dz*adzw(k)))
!!$               if(vadvtend.lt.0.) then ! wsub(k).gt.0.
!!$                  wsub(k) = wsub(k)-dtn*vadvtend/max(1.e-3,(tg0(k)+gamaz(k)-tg0(k-1)-gamaz(k-1))/(dz*adzw(k)))
!!$               else ! wsub(k).lt.0.
!!$                  wsub(k) = wsub(k)-dtn*vadvtend/max(1.e-3,(tg0(k+1)+gamaz(k+1)-tg0(k)-gamaz(k))/(dz*adzw(k+1)))
!!$               end if
              elseif(z(k).gt.14000.) then !nudge higher altitudes back in 5 days
                 vadvtend = (tg0(k)-tabs0(k))/432000. - ttend(k)
                 ! w_wtg = vadvtend / (dsli/dz)
                 ! NOTE: we impose a minimum static stability of 1 K/km
                 !        as in Raymond & Zeng, 2005. 
                 if(vadvtend.lt.0.) then ! wsub(k).gt.0.
                    wsub(k) = -vadvtend/max(1.e-3,(tg0(k)+gamaz(k)-tg0(k-1)-gamaz(k-1))/(dz*adzw(k)))
                 else ! wsub(k).lt.0.
                    wsub(k) = -vadvtend/max(1.e-3,(tg0(k+1)+gamaz(k+1)-tg0(k)-gamaz(k))/(dz*adzw(k+1)))
                 end if

              else
!!$               if(dotauz) itauz(k) = 1./tauls
              end if
           end do

!!$         if(masterproc) write(*,*) &
!!$              'Stopping: WTG large-scale forcing not implemented yet'
!!$         call task_abort()

        end if ! if(doblradforcing)

        call subsidence() !bloss: all uw forcing options use subsidence.

     end if

  end if


  ! Prescribed Radiation Forcing:


  if(doradforcing.and.time.gt.timelargescale) then

     k1=1
     k2=1
     do i=1,nrfc-1
        if(day.gt.dayrfc(i)) then
           k2=i+1
           k1=i
        endif
     end do
     coef=(day-dayrfc(k1))/(dayrfc(k2)-dayrfc(k1))
     do k=1,nzm
        radtend=dtrfc(k,k1)+(dtrfc(k,k2)-dtrfc(k,k1))*coef
        radqrlw(k)=radtend*float(nx*ny)
        radqrsw(k)=0.
        do j=1,ny
           do i=1,nx
              t(i,j,k)=t(i,j,k)+radtend*dtn 
           end do
        end do
     end do

  endif



  ! Surface flux forcing:

  if(dosfcforcing.and.time.gt.timelargescale) then

     k1=1
     k2=1
     do i=1,nsfc-1
        if(day.gt.daysfc(i)) then
           k2=i+1
           k1=i
        endif
     end do
     if(doxy) then 
        if(masterproc)  print*,'doxy=',doxy,': no sfc forcing is allowed'
        call task_abort()
     end if
     coef=(day-daysfc(k1))/(daysfc(k2)-daysfc(k1))
     tabs_s=sstsfc(k1)+(sstsfc(k2)-sstsfc(k1))*coef
     fluxt0=(hsfc(k1)+(hsfc(k2)-hsfc(k1))*coef)/(rho(1)*cp)
     fluxq0=(lesfc(k1)+(lesfc(k2)-lesfc(k1))*coef)/(rho(1)*lcond)
     tau0=tausfc(k1)+(tausfc(k2)-tausfc(k1))*coef
     do j=1,ny
        do i=1,nx
           sstxy(i,j) = tabs_s
        end do
     end do

     if(dostatis) then
        sstobs = tabs_s  ! sst is not averaged over the sampling period
        lhobs = lhobs + fluxq0 * rho(1)*lcond
        shobs = shobs + fluxt0 * rho(1)*cp
     end if

  endif

  if(.not.dosfcforcing.and.dodynamicocean) then
     call sst_evolve()

  endif

  ! peters.  update the evolving SST if necessary
  if(.not.dosfcforcing.and.doevolveSSTforcing) then
    call set_sstevolveforcing()
  end if

  ! peters.  update the evolving soil wetness if necessary
  if(.not.dosfcforcing.and.doxysoil_wet.and.icycle.eq.1) then
    call set_evolveforcing(evolvesoilwet,xysoil_wet(1:nx,1:ny))
  end if

  do k=1,nzm
     total_water_ls = total_water_ls + &
          sum(q(1:nx,1:ny,k))*adz(k)*dz *rho(k)
  end do

end subroutine forcing
