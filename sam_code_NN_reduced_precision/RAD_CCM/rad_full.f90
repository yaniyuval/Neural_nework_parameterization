subroutine rad_full()

! Interface to the longwave and shortwave radiation code from the
! NCAR Community Climate Model (CCM3.5).
!
use rad
use vars
use params
use simple_land, only : landmask_frac, dolandseamask ! peters

implicit none

! Local space:

      	real pmid(nzm)		! Level pressure (Pa)
      	real pint(nz)		! Model interface pressure (Pa)
	real massl(nzm)		! Level mass (g/m2)
	real pmidrd(nzm)	! Level pressure (dynes/cm2)
	real pintrd(nz)		! Model interface pressure (dynes/cm2)
	real pmln(nzm)		! Natural Log of pmid
	real piln(nz)		! Natural Log of pint
	real tlayer(nzm)	! Temperature
	real qlayer(nzm)	! Specific humidity
	real plco2(nz)		! Prs weighted CO2 path
	real plh2o(nz)		! Prs weighted H2O path
	real tclrsf(nz)		! Total clear sky fraction
	real cld(nz)		! Fractional cloud cover
	real cldeff(nz)		! Effective cloud cover (for LW)
	real clwp(nzm)		! Cloud liquid water path
	real fice(nzm)		! Fractional ice content within cloud
	real rel(nzm)		! Liquid effective drop radius (micron)
	real rei(nzm)		! Ice effective drop size
	real irei(nzm)		! Inverse ice effective drop size
	real o3vmr(nzm)		! Ozone volume mixing ratio

	real qrl(nzm)		! Longwave heating rate (K/s)
	real qrs(nzm)		! Shortwave heating rate (K/s)
        real flu(nz)		! Longwave upward flux
	real fld(nz)		! Longwave downward flux
	real fsu(nz)		! Shortwave upward flux
	real fsd(nz)		! Shortwave downward flux

!	aerosol (mass mixing ratio):

	real aero(nzm)		! Aerosol
	real rh(nzm)		! relative humidity for aerorsol 

!       Diagnostics:

      	real flns          ! Surface cooling flux
      	real flnt          ! Net outgoing flux
      	real flnsc         ! Clear sky surface cooing
      	real flntc         ! Net clear sky outgoing flux
      	real flds         ! Down longwave flux at surface

      	real solin        ! Incident solar flux
      	real fsns         ! Surface absorbed solar flux
      	real fsnt         ! Total column absorbed solar flux
        real fsntm        ! Flux Shortwave Downwelling Top-of-model
      	real fsds         ! Flux Shortwave Downwelling Surface

      	real fsnsc        ! Clear sky surface absorbed solar flux
      	real fsntc        ! Clear sky total column absorbed solar flx
      	real sols         ! Direct solar rad incident on surface (< 0.7)
      	real soll         ! Direct solar rad incident on surface (>= 0.7)
      	real solsd        ! Diffuse solar rad incident on surface (< 0.7)
      	real solld        ! Diffuse solar rad incident on surface (>= 0.7)
      	real fsnirt       ! Near-IR flux absorbed at toa
      	real fsnrtc       ! Clear sky near-IR flux absorbed at toa
      	real fsnirtsq     ! Near-IR flux absorbed at toa >= 0.7 microns

		
	real qtot,lwupsfc,eccf,asdir,aldir,asdif,aldif
	real dayy
	integer i,j,k,m,ii,jj,i1,j1,count,nrad_call,it,jt
	integer iday, iday0
	real coef,factor,tmp(1)
	real(8) qradz(nzm),buffer(nzm)
	real perpetual_factor
        real perpetual_equinox,pii
        integer dayyear ! peters
        real tasdir,taldir,tasdif,taldif ! peters 

        if(icycle.ne.1) goto 999  ! ugly way to handle the subcycles. add rad heating.
	nrad_call = 3600./dt/sqrt(betafactor)
        pii = atan2(0.,-1.)

!-------------------------------------------------------
! Initialize some stuff
!


        if(initrad) then


           if(mod(nx,ndiv).ne.0.or.(RUN3D.and.mod(ny,ndiv).ne.0)) then
               if(masterproc) print*,'nx or ny is not divisible by ndiv'
               if(masterproc) print*,'set in RAD_CCM/rad.f90'
               if(masterproc) print*,'Stop.'
               call task_abort()
           end if


	   call tracesini()

           o3=o3*o3factor

           if(nrestart.eq.0) then

	     do k=1,nzm
	      do j=1,ny
	       do i=1,nx
	         tabs_rad(i,j,k)=0.
	         qv_rad(i,j,k)=0.
	         qc_rad(i,j,k)=0.
	         qi_rad(i,j,k)=0.
	         qs_rad(i,j,k)=0.
	         qrad(i,j,k)=0.
	       end do
	      end do
	     end do
	     nradsteps=0	  
             do k=1,nz
               radlwup(k) = 0.
               radlwdn(k) = 0.
               radswup(k) = 0.
               radswdn(k) = 0.
               radqrlw(k) = 0.
               radqrsw(k) = 0.
             end do
	 
             !bloss: interactive ozone initialization moved to settracer().

	   else

	       call read_rad()

	   endif

           ! peters.  updated to include doannual case
           if(doperpetual) then
              if(.not.doequinox .and. .not.doannual) then
            ! perpetual sun (no diurnal cycle)
              do j=1,ny
               do i=1,nx
                p_factor(i,j) = perpetual_factor(day0, latitude(i,j),longitude(i,j))
	       end do
	      end do
              elseif(doequinox) then
            ! perpetual sun (no diurnal cycle) symmetric about equator
              do j=1,ny
               do i=1,nx
                p_factor(i,j) = perpetual_equinox( &
                      latitude(i,j),longitude(i,j))
                if (i.eq.1) then !JY YAni
                 call task_rank_to_index(rank,it,jt)
                 print *, 'perpetual_equinox', j,jt , p_factor(1,j)
                end if 
               end do
              end do
              else
            ! peters perpetual sun (no diurnal cycle) averaged over ann. cycle
             if(masterproc) print*,'Initialize doannual'
              do dayyear = 1,365
               do j=1,ny
                do i=1,nx
                  p_factor(i,j) = p_factor(i,j) + &
                      perpetual_factor(float(dayyear), latitude(i,j),&
                                                       longitude(i,j))/365.
                end do
               end do
              end do
	     end if
           end if
	endif

        if (.not.dodoubleco2) co2factor = 1.
	call radini(ggr, cp, 0.622, 5.67e-8, co2factor)
	do k=1,nzm
	  m=nz-k
   	  pmid(m) = pres(k)*100.
	  pmln(m)=log(pmid(m))
	end do
	do m=2,nzm
          pint(m)=0.5*(pmid(m)+pmid(m-1))
	  piln(m) = log(pint(m)) 
	end do
        pint(1)=pmid(1)-0.5*(pmid(2)-pmid(1))
        pint(nz)=pres0*100.
	piln(1) = log(pint(1)) 
	piln(nz) = log(pint(nz)) 
        tmp(1) = pres0*100.

	call cldefr(tmp ,pmid, rel, rei)

	do k=1,nzm
          massl(k)=1000.*(pint(k+1)-pint(k))/ggr
	  aero(k)=0.
	  rh(k)=0.
	  o3vmr(k)=0.6034*o3(k)
	  irei(k)=1./rei(k)
	  qrl(k)=0.
	  qrs(k)=0.
	end do
	do k=1,nz
	  cld(k)=0.
	  cldeff(k)=0.
	  flu(k)=0.
	  fld(k)=0.
	  fsu(k)=0.
	  fsd(k)=0.
	end do

 
!-----------------------------------------------------------
! Check if it is time to compute gas absortion coefficients for
! longwave radiation. This coefficients are computed for
! horizontally average fields and storred internally, so
! they need not to be recalculated on every call of radclw() for
! efficiency reasons.

	if(initrad.or.mod(nstep,nrad_call).eq.0) then

          initrad=.false.

	  do jj=1,nydiv 
	  j1=(jj-1)*(ny/nydiv) 
	  do ii=1,nxdiv 
	   i1=(ii-1)*(nx/nxdiv) 
	   do k=1,nzm
	    tlayer(k)=0.
	    qlayer(k)=0.
	    cld(k) = 0.
	    cldeff(k) = 0.
	    clwp(k) = 0.
	    fice(k) = 0.
	    m=nz-k
	    
	    count = 0
	    do j=j1+1,j1+(ny/nydiv)
	     do i=i1+1,i1+(nx/nxdiv)	 
	       tlayer(k)=tlayer(k)+tabs(i,j,m)
	       qlayer(k)=qlayer(k)+q(i,j,m)-qn(i,j,m)
	       count = count+1
	     end do
	    end do	
	    tlayer(k)=tlayer(k)/float(count)
	    qlayer(k)=max(1.e-7,qlayer(k)/float(count))

              !--------------!interactive ozone
              if(dotro3.and.doactiveo3) then !interactive ozone
                 count = 0
                 o3(k) = 0.
                 do j=j1+1,j1+(ny/nydiv)
                    do i=i1+1,i1+(nx/nxdiv)	 
                       o3(k)=o3(k)+tro3(i,j,m)	!kzm interactive ozone
                       count = count+1
                    end do
                 end do
                 o3(k)=o3(k)/float(count)
                 o3vmr(k)=0.6034*o3(k)
              endif
              !--------------!interactive ozone

	   end do

	   call radinp(pmid, pint, qlayer, cldeff, pmidrd, pintrd, &
     			plco2, plh2o, tclrsf)

           if ((dodynamicocean).or.(ocean_type.ne.0)) then
              ! Compute local-averaged SST
              count = 0
              lwupsfc = 0.
              do j=j1+1,j1+(ny/nydiv)
                 do i=i1+1,i1+(nx/nxdiv)	 
                    lwupsfc = lwupsfc + 5.67e-8*sstxy(i,j)**4 * 1000. ! CGS units
                    count = count+1
                 end do
              end do
              lwupsfc=lwupsfc/float(count)
           else
	   lwupsfc = 5.67e-8*sstxy(1,1)**4 * 1000. ! CGS units
           end if
           tmp(1) = lwupsfc


	   call radclw(.true., tmp, tlayer, qlayer, o3vmr, &
     	          absnxt(1,1,ii,jj),abstot(1,1,ii,jj),emstot(1,ii,jj),	&
     	  	  pmidrd,pintrd, pmln, piln, plco2, plh2o, &
     		  n2o     ,ch4     ,cfc11       ,cfc12   ,  &
     		  cld, tclrsf, qrl , flu, fld, flns ,flnt , &
      		  flnsc ,flntc ,flds ) 
	  end do
	  end do

	endif
	

!------------------------------------------------------
!  Accumulate thermodynamical fields over nrad steps 
!
	
        if(.not.dokruegermicro) then
           do k=1,nzm
              do j=1,ny
                 do i=1,nx
                    tabs_rad(i,j,k)=tabs_rad(i,j,k)+tabs(i,j,k)
                    qv_rad(i,j,k)=qv_rad(i,j,k)+q(i,j,k)-qn(i,j,k)
                    qc_rad(i,j,k)=qc_rad(i,j,k)+qn(i,j,k)*omegan(tabs(i,j,k))
                    qi_rad(i,j,k)=qi_rad(i,j,k)+qn(i,j,k)* &
                         (1.-omegan(tabs(i,j,k)))
                 end do
              end do
           end do

           if(dokruegereffectiveradii) then ! snow is radiatively active
              do k=1,nzm
                 do j=1,ny
                    do i=1,nx
                       qs_rad(i,j,k)=qs_rad(i,j,k)+qp(i,j,k)* &
                            (1.-omegap(tabs(i,j,k)))*(1.-omegag(tabs(i,j,k)))
                    end do
                 end do
              end do
           end if

        else ! if using Lin/Lord microphysics (e.g. dokruegermicro)

           do k=1,nzm
              do j=1,ny
                 do i=1,nx
                    tabs_rad(i,j,k)=tabs_rad(i,j,k)+tabs(i,j,k)
                    qv_rad(i,j,k)=qv_rad(i,j,k)+qv(i,j,k)
                    qc_rad(i,j,k)=qc_rad(i,j,k)+qc(i,j,k)
                    qi_rad(i,j,k)=qi_rad(i,j,k)+qi(i,j,k)
                    qs_rad(i,j,k)=qs_rad(i,j,k)+qs(i,j,k)
                 end do
              end do
           end do

        end if
	nradsteps=nradsteps+1

!----------------------------------------------------
! Update radiation variables if the time is due
!

        !kzm Oct.14, 03 changed .eq.nrad to .ge.nrad to handle the
        ! case when a smaller nrad is used in restart  
	if(nstep.eq.1.or.nradsteps.ge.nrad) then 

! Compute radiation fields for averaged thermodynamic fields

	
	  coef=1./float(nradsteps)

          do k=1,nz
            radlwup(k) = 0.
            radlwdn(k) = 0.
            radswup(k) = 0.
            radswdn(k) = 0.
            radqrlw(k) = 0.
            radqrsw(k) = 0.
          end do


	  do k=1,nzm
	   do j=1,ny
	    do i=1,nx
	      tabs_rad(i,j,k)=tabs_rad(i,j,k)*coef
	      qv_rad(i,j,k)=qv_rad(i,j,k)*coef
	      qc_rad(i,j,k)=qc_rad(i,j,k)*coef
	      qi_rad(i,j,k)=qi_rad(i,j,k)*coef
	      qs_rad(i,j,k)=qs_rad(i,j,k)*coef
	    end do
	   end do
	  end do

	  do j=1,ny
	   jj = (j-1)/(ny/nydiv)+1
	   do i=1,nx
	   ii = (i-1)/(nx/nxdiv)+1

           if(.not.dokruegereffectiveradii) then

	     do k=1,nzm
	       m=nz-k
	       tlayer(k)=tabs_rad(i,j,m)
	       qtot = qc_rad(i,j,m)+qi_rad(i,j,m)
	       qlayer(k)=max(1.e-7,qv_rad(i,j,m))
	       if(qtot.gt.0.) then
	         clwp(k) = qtot*massl(k)
	         fice(k) = qi_rad(i,j,m)/qtot
	         cld(k) = 0.99
	         cldeff(k) = min(0.9,cld(k)*(1.-exp( &
		    -(0.15*(1.-fice(k))+ &
		       1.66 * (0.005+irei(k))*fice(k))*clwp(k))))
	       else
	         cld(k) = 0.
	         cldeff(k) = 0.
	         clwp(k) = 0.
	         fice(k) = 0.
	       endif
	     end do

            else ! dokruegereffectiveradii = .true.

              ! ============= DOKREUGEREFFECTIVERADII =================
              ! use effective radii for liquid/ice from Luo et al (2004)
              ! rel = 10 um, rei = 25 um, re_{snow} = 75 um
              do k=1,nzm
                 m=nz-k
                 tlayer(k)=tabs_rad(i,j,m)
                 qtot = qc_rad(i,j,m)+qi_rad(i,j,m)+qs_rad(i,j,m)
                 qlayer(k)=max(1.e-7,qv_rad(i,j,m))
                 rel(k) = 10.
                 rei(k) = 25.
                 if(qtot.gt.0.) then
                    clwp(k) = qtot*massl(k)
                    fice(k) = (qi_rad(i,j,m)+qs_rad(i,j,m))/qtot
                    cld(k) = 0.99
                    if (qs_rad(i,j,m).gt.0.) then
                       rei(k) = 25.*(qi_rad(i,j,m)+qs_rad(i,j,m)) &
                            /(qi_rad(i,j,m)+qs_rad(i,j,m)/3.)
                    end if
                    cldeff(k) = min(0.9,cld(k)*(1.-exp( &
                         -(0.15*(1.-fice(k))+ &
                         1.66 * (0.005+1./rei(k))*fice(k))*clwp(k))))
                 else
                    cld(k) = 0.
                    cldeff(k) = 0.
                    clwp(k) = 0.
                    fice(k) = 0.

                 endif
              end do

            end if

	     call radinp(pmid, pint, qlayer, cldeff, pmidrd, pintrd,      &
			plco2,plh2o,tclrsf)

	     lwupsfc = 5.67e-8*sstxy(i,j)**4 * 1000. ! CGS units
             tmp(1) = lwupsfc

	     if(dolongwave) then

	       call radclw(.false., tmp, tlayer, qlayer, o3vmr, &
	          absnxt(1,1,ii,jj),abstot(1,1,ii,jj),emstot(1,ii,jj),	&
	  	  pmidrd,pintrd, pmln, piln, plco2, plh2o, &
		  n2o     ,ch4     ,cfc11       ,cfc12   , &
 		  cldeff, tclrsf, qrl ,flu, fld, flns ,flnt , &
 		  flnsc ,flntc ,flds ) 
	     endif

	     if(doshortwave) then

               ! peters -- if necessary, scale time elapsed from day0 by
               ! darefactor.  Also use a 360 day year
               if(doDARErad) then
                 day = (day-day0)*sqrt(betafactor)*365./360. + day0
               end if

               if (doseasons) then
                  ! The diurnal cycle of insolation will vary
                  ! according to time of year of the current day.
                  dayy = day
               else
                  ! The diurnal cycle of insolation from the calendar
                  ! day on which the simulation starts (day0) will be
                  ! repeated throughout the simulation.
	          iday0 = day0
	          iday = day
	          dayy = day-iday
	          dayy = iday0 + dayy
	       end if
	       if(doperpetual) then
                 if (dosolarconstant) then
                    ! fix solar constant and zenith angle as specified
                    ! in prm file.
                    coszrs = cos(zenith_angle*pii/180.)
                    eccf = solar_constant/(1367.)
                 else
                  ! perpetual sun (no diurnal cycle) - Modeled after Tompkins
                  coszrs = 0.637 ! equivalent to zenith angle of 50.5 deg
                  eccf = p_factor(i,j)/coszrs ! Adjust solar constant
                 end if
	       else
	          call zenith(dayy,latitude(i,j),longitude(i,j),coszrs,eccf)
	       end if
               ! peters -- replace day with non-rescaled value
               if(doDARErad) then
                 day = (day-day0)/sqrt(betafactor)*360./365. + day0
               end if
	       tmp(1)=coszrs
               ! peters.  updated call to albedo for land/sea mask
               if(.not.dolandseamask) then
                 call albedo(OCEAN,tmp,asdir,aldir,asdif,aldif)
               else
                 ! albedo over ocean
                 call albedo(.true.,tmp,asdir,aldir,asdif,aldif)
                 ! albedo over land
                 call albedo(.false.,tmp,tasdir,taldir,tasdif,taldif)
                asdir = landmask_frac(i,j)*tasdir+(1.-landmask_frac(i,j))*asdir
                aldir = landmask_frac(i,j)*taldir+(1.-landmask_frac(i,j))*aldir
                asdif = landmask_frac(i,j)*tasdif+(1.-landmask_frac(i,j))*asdif
                aldif = landmask_frac(i,j)*taldif+(1.-landmask_frac(i,j))*aldif
               end if

	       call radcsw(pintrd, qlayer, o3, aero, rh, cld, clwp, &
	         rel, rei, fice, eccf, coszrs, &
                 asdir,aldir,asdif,aldif, solin ,qrs , fsu, fsd, &
	 	 fsns ,fsnt ,fsntm, fsds ,fsnsc ,fsntc ,sols ,soll , &
                 solsd ,solld ,fsnirt ,fsnrtc ,fsnirtsq)	 


	     endif

	     do k=1,nzm
	        m=nz-k
	        qrad(i,j,m)=qrl(k)+qrs(k)
	        radlwup(m)=radlwup(m)+flu(k)*1.e-3
	        radlwdn(m)=radlwdn(m)+fld(k)*1.e-3
	        radqrlw(m)=radqrlw(m)+qrl(k)
	        radswup(m)=radswup(m)+fsu(k)*1.e-3
	        radswdn(m)=radswdn(m)+fsd(k)*1.e-3
	        radqrsw(m)=radqrsw(m)+qrs(k)
	     enddo

	     lwnsxy(i,j) = flns
	     swnsxy(i,j) = fsns
	     lwntxy(i,j) = flnt
	     swntxy(i,j) = fsnt
             swntmxy(i,j)= fsntm
	     lwnscxy(i,j) = flnsc
	     swnscxy(i,j) = fsnsc
	     lwntcxy(i,j) = flntc
	     swntcxy(i,j) = fsntc
	     swdsxy(i,j) = fsds
	     lwdsxy(i,j) = flds
	     solinxy(i,j) = solin

	   end do
	  end do

	  do k=1,nzm
	   do j=1,ny
	    do i=1,nx
	     tabs_rad(i,j,k)=0.
	     qv_rad(i,j,k)=0.
	     qc_rad(i,j,k)=0.
	     qi_rad(i,j,k)=0.
	    end do
	   end do
	  end do
	  nradsteps=0

          !-----kzm scale the radiative fluxes up by sqrt(betafactor)-----
          factor=sqrt(betafactor)
          if (dobetafactor) then
             qrad=qrad*factor
             radqrlw=radqrlw*factor
             radqrsw=radqrsw*factor
             lwnsxy=lwnsxy*factor
             lwntxy=lwntxy*factor
             swnsxy=swnsxy*factor
             swntxy=swntxy*factor
             swntmxy=swntmxy*factor
             solinxy=solinxy*factor
          end if

	  if(masterproc.and.doshortwave.and..not.doperpetual) &
		print*,'radiation: coszrs=',coszrs,' solin=',solin
	  if(masterproc.and.doshortwave.and.doperpetual) &
		print*,'radiation: perpetual sun, solin=',solin
	  if(masterproc.and..not.doshortwave) &
		print*,'longwave radiation is called'

	endif ! (nradsteps.eq.nrad) 

!------------------------------------------------------
! Homogenize radiation:

       if(doradhomo) then    

	factor = 1./dble(nx*ny)
        do k=1,nzm
         qradz(k) = 0.
         do j=1,ny
          do i=1,nx
             qradz(k) = qradz(k) + qrad(i,j,k)
	  end do
	 end do
	 qradz(k) = qradz(k) * factor
	 buffer(k) = qradz(k)
	end do

        factor = 1./float(nsubdomains)
        if(dompi) call task_sum_real8(qradz,buffer,nzm)
	   
        do k=1,nzm
          qradz(k)=buffer(k)*factor
          do j=1,ny
           do i=1,nx
             qrad(i,j,k) = qradz(k) 
           end do
          end do
        end do

       end if
!------------------------------------------------------
! Homogenize radiation in x direction:

       if(doradhomox) then    

     call task_rank_to_index(rank,it,jt)

     do jj=1,nsubdomains_y 
         do j=1,ny
         qradz(:) = 0.
         if((jj-1) * (ny_gl/nsubdomains_y).eq.jt) then
            factor = 1./dble(nx)
            do k=1,nzm
               do i=1,nx
                  qradz(k) = qradz(k) + qrad(i,j,k)
               end do
               qradz(k) = qradz(k) * factor
               buffer(k) = qradz(k)
            end do
         endif

        factor = 1./float(nsubdomains_x)
        if(dompi) call task_sum_real8(qradz,buffer,nzm)
	   
        if((jj-1) * (ny_gl/nsubdomains_y).eq.jt) then
        do k=1,nzm
          qradz(k)=buffer(k)*factor
           do i=1,nx
             qrad(i,j,k) = qradz(k) 
           end do
        end do
        endif
        end do
       end do
       end if

	 do j=1,ny
	  do i=1,nx
	     lwns_xy(i,j) = lwns_xy(i,j) + lwnsxy(i,j) 
	     swns_xy(i,j) = swns_xy(i,j) + swnsxy(i,j)
	     lwnt_xy(i,j) = lwnt_xy(i,j) + lwntxy(i,j) 
	     swnt_xy(i,j) = swnt_xy(i,j) + swntxy(i,j)
	     swntm_xy(i,j)= swntm_xy(i,j)+ swntmxy(i,j)
	     lwnsc_xy(i,j) = lwnsc_xy(i,j) + lwnscxy(i,j) 
	     swnsc_xy(i,j) = swnsc_xy(i,j) + swnscxy(i,j)
	     lwntc_xy(i,j) = lwntc_xy(i,j) + lwntcxy(i,j) 
	     swntc_xy(i,j) = swntc_xy(i,j) + swntcxy(i,j)
	  !   solin_xy(i,j) = solin_xy(i,j) + solinxy(i,j)
             solin_xy(i,j) =  solinxy(i,j) ! JY Yani changed - wanted to get solin_xy to the RF. Should also do the same with SST? 
	  end do
	 end do
!----------------------------------------------------------------
	if(dostatis) then

	  do j=1,ny
	   do i=1,nx
	    s_flns = s_flns + lwnsxy(i,j) 
	    s_fsns = s_fsns + swnsxy(i,j) 
	    s_flnt = s_flnt + lwntxy(i,j) 
            ! Changed to reflect TOA net sw, rather than top of model
	    s_fsnt = s_fsnt + swntxy(i,j) 
	    s_flnsc = s_flnsc + lwnscxy(i,j) 
	    s_fsnsc = s_fsnsc + swnscxy(i,j) 
	    s_flntc = s_flntc + lwntcxy(i,j) 
	    s_fsntc = s_fsntc + swntcxy(i,j) 
	    s_fsds = s_fsds + swdsxy(i,j) 
	    s_flds = s_flds + lwdsxy(i,j) 
	    s_solin = s_solin + solinxy(i,j) 
	   end do
	  end do
	end if
!----------------------------------------------------------------
!  Write the radiation-restart file:

        ! peters
        if((mod(nstep,nstat).eq.0).and.(.not.dorestart_last) .and.    &
           (nwriterestart.eq.-1 .or. mod(nstep,nwriterestart).eq.0)) then

          call write_rad() ! write radiation restart file

        endif

!-------------------------------------------------------
! Update the temperature field:
999     continue
       if (rad_lev_pred<48) then
        do k=rad_lev_pred+1,nzm
         do j=1,ny
          do i=1,nx
           t(i,j,k)=t(i,j,k)+qrad(i,j,k)*dtn
           misc(i,j,k)=qrad(i,j,k)
          end do
         end do
        end do
       endif
!print *,'changed the radiation to work only in certain levels...'



end
