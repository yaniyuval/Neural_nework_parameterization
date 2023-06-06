	
subroutine setparm
	
!       initialize parameters:

use vars
use params
use isccp, only : isccp_zero
use isccpTables, only : isccp_tables_init
use mse   ! peters
use simple_ocean   ! peters
use simple_land    ! peters
use microscaling
implicit none
	
integer j,k, ios

NAMELIST /PARAMETERS/ dodamping, doupperbound, docloud, doprecip, &
                dorandomforest, do_rf_diffusion, do_rf_q_surf_flux, rf_filename, rf_uses_qp, rf_uses_rh, &
                do_wind_input_rf, do_z_diffusion_rf, do_q_T_surf_fluxes_rf, do_surf_wind_rf, do_q_surf_fluxes_out_rf, & !New parameters after matlab corrections.
                do_sedimentation_rf, do_fall_tend_rf, do_radiation_output_rf, do_flux_rf, do_hor_advection_rf, do_hor_diffusion_rf, &
                Tin_feature_rf, qin_feature_rf, Tin_z_grad_feature_rf, qin_z_grad_feature_rf, predict_tendencies_rf, do_tkh_z, do_tkh_xy, &
		do_vert_wind_input_rf, do_yin_input, diffusivity_levels_rf, &
                rf_filename_diffusion, rf_filename_tkh, &
                rad_lev_pred, input_ver_dim, int1_dum, int2_dum, &
                dolongwave, doshortwave, dosgs, &
                docoriolis, dosurface, dolargescale, doradforcing, &
		nadams,fluxt0,fluxq0,tau0,tabs_s,z0,tauls,nelapse, &
		dt, dx, dy, fcor, ug, vg, nstop, caseid, &
		nstat, nstatfrq, nprint, nrestart, doradsimple, &
		nsave3D, nsave3Dstart, nsave3Dend, dosfcforcing, &
		donudging_uv, donudging_tq, dosmagor, doscalar,  &
		timelargescale, longitude0, latitude0, day0, nrad, &
		CEM,LES,OCEAN,LAND,SFC_FLX_FXD,SFC_TAU_FXD, SFC_TAU_MASK,soil_wetness, &
                doensemble, nensemble, doxy, dowallx, dowally, &
                nsave2D, nsave2Dstart, nsave2Dend, qnsave3D, & 
                dt_cup, docup, docolumn, save2Dbin, save3Dbin, &
                save2Dsep, save3Dsep, dogzip2D, dogzip3D, restart_sep, &
	        doseasons, doperpetual, doradhomo, dosfchomo,doradhomox, dosfchomox, doisccp, &
	        dodynamicocean, ocean_type, &
		dosolarconstant, solar_constant, zenith_angle, raddir
	
NAMELIST /UWOPTIONS/  nelapsemin, dorestart_last, dooldrestart, &
                nwriterestart,                    & ! peters
                doequinox, dodoubleco2, co2factor, o3factor, iyear, &
                doreflectivity, doinstantoutput, dowaveenergetics, & 
                dosaturatedwvp, dorceintercomparison, dowaveoutput, &
                docolumnbudgets, do850mbarwinds, dooutputheights, &
                dolowlevelthetae, doclearsky, doinstantolr, & 
                docloudechoheights, docasefive, &
                hml, Szero, sst_target, & !uw ocean options
                deltaS, delta_sst, dofixedoceancooling, resetsst, &
                dokruegermicro, dokruegereffectiveradii, & !bloss
                tracer_option, dorestart_tracer, dotauz, doreset_tracer, & !kzm
		doactiveo3, dotro3, itauo3, itautrx,itautry,itautrz,itautrzz,&
                dosubsidence_tonly, dosubsidence, domoistadiabaticnudging, &
                dobetaplane, betafactor, dxdzfactor, ravefactor, dobetafactor, dowaterbubble,&
                doydamping,doxdamping, doperiodic_forcing, day0_forcing, period_forcing,&
                nsaveMSE,navgMSE,nsaveMSEstart,nsaveMSEend,nsaveMSE3D,&!peters
                domombudget, writeMSE3D, & !peters
                doblradforcing, blqrad, blqmin, blqmax, bldivg, &
                doblforcing, doblsubs, blsubs, bltauz, dobldivg, &
                doaddqvrestart, dowtg, wtgtau, doblsubs2, blsubs2, &
                ncycle0, ncyclemax, ncyclemin,do3Daverage, &
                doclouddropsed, donudging_q, nudging_q_factor,donudging_hadley, donoevap, donoqpload, doparameterizedwave,second_mode_factor, &
                douwpbl, douwcu, uwpbl_ndiv, uwcu_ndiv,mpbetafactor,daremicrophysics, fixedndiff, domicroscaling, cmin, cmax,wavenumber_factor, &
                 alfafactor,nstartlinearwave,nsteplinearwavebg,nsteplinearwave,&
                 first_mode_speed,second_mode_speed,wavedampingtime,doevolvewavebg,pbldampingfactor,&
                 qtendfactor,dosmoothwconv,nhwave_avg,dowavesubsidence_tonly,dosteadylinearwave,doperturb_realtime,&
                 dobulksfc,dolinearwavelid,dorave,dowavedampingt,doadvectbg,dozerowsub,evapfactor,tmin_evap, &
                dosubsidenceonbackground, & !JAA
                 doevolveSSTforcing,evolveSST,& ! peters
                dolandseamask, land_filename, land_ncfldname, & !peters
                dointeractiveland, hml_land, tau_landdamping, & !peters
                doxysurf_rough, surfrough_filename,surfrough_ncfldname,&!peters
                maxsurfrough, nstart_surfrough,doxysoil_wet, evolvesoilwet,&
                doannual, doDARErad, doDAREradlongitude, &
                dosurfacefix , block_size, hadley_lim , sfcum





!----------------------------
!  Set defaults:

	dodamping 	= .false.
	doupperbound   	= .false.
	docloud   	= .false.
	doprecip        = .false.
	dolongwave	= .false.
	doshortwave	= .false.
	doradsimple 	= .false.
	dosgs		= .false.
	dosmagor	= .false.
	doscalar	= .false.
	dosubsidence	= .false.
	docoriolis	= .false.
	dosurface	= .false.
	dolargescale    = .false.
	doradforcing    = .false.
	dosfcforcing    = .false.
	donudging_uv	= .false.
	donudging_tq	= .false.
	doensemble	= .false.
	doxy    	= .false.
	dowallx    	= .false.
	dowally    	= .false.
	docup		= .false.
	docolumn	= .false.
	doseasons	= .false.
	doperpetual	= .false.
	doradhomo	= .false.
	dosfchomo	= .false.
	doradhomox	= .false.
	dosfchomox	= .false.
	doisccp		= .false.
	dodynamicocean 	= .false.
        dosolarconstant = .false.
	CEM		= .false.
	LES		= .false.
	OCEAN		= .false.
	LAND		= .false.
	SFC_FLX_FXD	= .false.
	SFC_TAU_FXD	= .false.
	SFC_TAU_MASK	= .false.
        daremicrophysics= .true.
        block_size      = 1
        hadley_lim      = 100000
        ! JA modified by cw 5/7/07
        domicroscaling  = .false.
        cmin            = 0.2
        cmax            = 0.75
        
        doperturb_realtime=.false.

        ! pog yani
	dorandomforest  = .false.
    do_rf_diffusion= .true.
    do_rf_q_surf_flux = .false.
	rf_uses_qp  = .false.
	rf_uses_rh  = .false.
   !New parameters after matlab corrections. 
    do_wind_input_rf = .false.
    do_z_diffusion_rf = .false.
    do_q_T_surf_fluxes_rf = .false.
    do_surf_wind_rf = .false.
    do_q_surf_fluxes_out_rf = .false.
    do_sedimentation_rf = .false.
    do_fall_tend_rf = .false.
    do_radiation_output_rf = .false.
    rad_lev_pred = 48
    input_ver_dim = 48
    int1_dum = 0 
    int2_dum = 0
    do_flux_rf = .false.
    do_hor_advection_rf = .false.
    do_hor_diffusion_rf = .false.
    Tin_feature_rf = .true.
    qin_feature_rf = .true.
    Tin_z_grad_feature_rf = .false.
    qin_z_grad_feature_rf = .false.
    predict_tendencies_rf = .true.
    do_tkh_z = .false. 
    do_tkh_xy = .false.
    do_vert_wind_input_rf = .false. 
    do_yin_input = .false.
    diffusivity_levels_rf = 15


        rf_filename = 'INPUT/rf_trees.nc'
        rf_filename_diffusion = 'INPUT/rf_trees.nc'
        rf_filename_tkh = 'INPUT/rf_trees.nc'    

    
        !kzm linear wave stuff
         dowavesubsidence_tonly = .false.
         doevolvewavebg = .true.
         wavenumber_factor = 0. !1 is for 1000km, 2 is for 500km wavelength wave
         alfafactor=1.
         nstartlinearwave=999999999
         nsteplinearwavebg=100
         nsteplinearwave=100
         wavedampingtime=10000.
         dowavedampingt=.false.
         doadvectbg=.true.
         dozerowsub=.false.
         pbldampingfactor=1.
         qtendfactor=1.
         nhwave_avg=1
         dosmoothwconv=.true.
         dosteadylinearwave=.false.
         dolinearwavelid=.false.
         
         sfcum=0.
         !JAA
         dosubsidenceonbackground = .true.

         dobulksfc=.false. !use aerodynamic formula
         evapfactor  = 1. !enhance factor for evaporation of precip
         tmin_evap = 0. !only do evaporation of qp above this temperature

	nadams		= 3
	dt		= 0
	dt_cup		= 0
	dx		= 0
	dy		= 0
	longitude0	= 0.
	latitude0	= 0.
	fcor	        = -999.
	day0		= 0.
	nrad		= 1
	ug		= 0.
	vg		= 0.
	fluxt0		= 0.
	fluxq0		= 0.
	tau0		= 0.
	z0		= 0.035
        soil_wetness    = 1.
	timelargescale  = 0.
	tauls		= 7200.
	tabs_s 		= 0.
	nstop 		= 0
	nelapse		= 999999999
	caseid		= 'les00000'
	nstat		= 1000
	nstatfrq	= 50
	nprint		= 1000
	nrestart 	= 0
	restart_sep 	= .false.
	save3Dbin	= .false.
	save2Dsep	= .false.
	save3Dsep	= .false.
	nsave3D		= 1
	nsave3Dstart	= 99999999
	nsave3Dend	= 999999999
	dogzip2D	= .false.
	dogzip3D	= .false.  ! SAM6.2 default value is .false.
	save2Dbin	= .false.
	nsave2D		= 1
	nsave2Dstart	= 99999999
	nsave2Dend	= 999999999
	nensemble	= 0
 	qnsave3D	= 0.
	ocean_type	= 0 
        raddir          = '/home/kzm/Model/SAMson/RADDATA/' !SAM6.3 default:'./RADDATA'

        ! Specify solar constant and zenith angle for perpetual insolation.
        ! Note that if doperpetual=.true. and dosolarconstant=.false.
        ! the insolation will be set to the daily-averaged value on day0. 
        solar_constant = 685. ! Values from Tompkins & Craig, J. Climate (1998)
        zenith_angle = 51.7

        !==================================================================
        ! UW OPTIONS

        ! radiation options
        doequinox       = .false. ! Symmetric insolation about equator
        iyear           = 1999    ! Year for which insolation is calculated
        dodoubleco2  = .false.    ! Double CO2 levels in radiation calculation
        co2factor    = 2.         ! Scaling factor for CO2 (defalt = 2)
        o3factor     = 1.         ! Scaling factor for ozone (defalt = 1)
        doannual = .false.        ! peters.  annual mean solar rad
        doDARErad= .false.        ! peters
        doDAREradlongitude= .false. ! peters
        
        dokruegereffectiveradii = .false. ! Use Luo et al (2004) effective
                                          ! radii for cloud water/ice/snow
        !!!!!!!! NOTE THAT SNOW IS RADIATIVELY ACTIVE WHEN TRUE !!!!!!!

        ! nudging options
        dotauz     = .false.               !kzm: nudge T to mimic
                                           ! effect of radiation 
 	domoistadiabaticnudging = .false.  !bloss: Apply nudging so that it
                                           ! introduces a perturbation
                                           ! which is uniform in saturated  
                                           ! moist static energy and
                                           ! neutral in relative humidity. 

        ! forcing options
 	dosubsidence_tonly = .false.  !kzm: Only apply subsidence velocity 
                                      ! to temperature, so that subsidence 
                                      ! is only a dynamical heating term
	doperiodic_forcing= .false.   !kzm: use periodic forcing
	day0_forcing    = 0.
        period_forcing  = 3.

        ! MICROPHYSICS OPTIONS
        dokruegermicro = .false. ! use a port of Steve Krueger's
                                 ! version of the Lord-Lin scheme.
                                 ! Includes prognostic equations for
                                 ! vapor, cloud, ice, snow, graupel and rain.
        
        ! tracer options
 	tracer_option	= 0  
	itauo3		= 0.     !kzm tracer
	itautrz		= 0.     !kzm tracer
	itautrx		= 0.     !kzm tracer
	itautry		= 0.     !kzm tracer
	itautrzz		= 0.     !kzm tracer
	dotrz		= .false. !kzm tracer
	dotrx		= .false. !kzm tracer
	dotry		= .false. !kzm tracer
	dotrzz		= .false. !kzm tracer
	dotro3		= .false. !kzm tracer
	doactiveo3	= .false. !kzm tracer
        donudging_hadley= .false. !kzm nudget towards a given zonal mean circulation
        donudging_q     = .false.
        nudging_q_factor= 1.

       ! RESTART OPTIONS
        resetsst = .false. ! reset sst when restart
	nelapsemin = 999999999 ! stop run after nrestart wallclock minutes
        dorestart_last = .false. ! write restart file only after last timestep
        dooldrestart = .false. ! use restart file from old UW version of SAM
        nwriterestart = -1     ! peters

        ! SLAB OCEAN OPTIONS
        hml = 20. ! Ocean Mixed Layer (ML) depth in meters
        Szero = -70. ! Initial value of deep ocean cooling.
        sst_target = 300. ! Target horizontally-averaged SST (K)
        dofixedoceancooling = .true. ! Fix ocean cooling in time at Szero.
        deltaS = 50. ! Ocean ML heating change across domain (W/m^2)
        delta_sst = 2. ! Amplitude of SST variation (K) (static ocean)

        
        ! 2D OUTPUT OPTIONS       ! THIS OPTION OUTPUTS ...
        doreflectivity  = .false. ! simulated radar reflectivities
        doclearsky = .false. ! clear sky radiative fluxes/insolation
        doinstantolr = .false. ! instantaneous value of net OLR flux
        docolumnbudgets = .false. ! column-integrated mse/sli budgets
        dosaturatedwvp = .false. ! saturated water vapor path
        do850mbarwinds = .false. ! 850 mbar horizontal winds
        dooutputheights = .false. ! geopotential heights at 200/850 mbar
        dolowlevelthetae = .false. ! theta_e and theta_es at 1000/850 mbar
        dowaveoutput = .false. ! projection of w, q, theta onto 
                               ! first two gravity wave modes
        dowaveenergetics = .false. ! energetics of projections of w, q, theta 
                                   ! onto first two gravity wave modes
        doinstantoutput = .false. ! instantaneous slices: pcp,e,fsnt,pw,
                                  ! flnt,fsntc,flntc,cond,tsfc,qsfc,wssfc
        dorceintercomparison  = .false. ! additional for rce intercomparison
                                        ! also adds outputs to .stat file
        docloudechoheights = .false. ! output 2D cloud/echo hts & cldtoptemp
        docasefive = .false. ! GCSS case 5 intercomparison statistics
        do3Daverage= .false.

        ! BETA PLANE OPTIONS
        dobetaplane = .false. ! Turn on beta effect
        dobetafactor = .false. ! Turn on beta factor
        doydamping  = .false. ! Damp fluctuations at y boundaries
        doxdamping  = .false. ! Damp fluctuations at x boundaries
        betafactor  = 1. ! Multiply beta by this factor.
        mpbetafactor=10000.
        ravefactor  = 1. ! divide Dw/Dt by square of this factor
        dorave=.true.
        dxdzfactor  = 1. ! to reduce the anisotropy in the grid used in tke_full
        dowaterbubble = .false. !initiate with a humidity anomaly
        donoevap    = .false.   !disable evaporation of precip
        donoqpload=.false. !disable the density effect of condensates

        ! peters MSE options
        nsaveMSE = 0
        navgMSE = 0
        nsaveMSEstart = 1e9
        nsaveMSEend =  2e9
        nsaveMSE3D = 1
        domombudget=.true.
        writeMSE3D=.true.

        ! peters option for surface scheme
        dosurfacefix = .false.

        ! boundary layer specified radiative cooling
        doblradforcing = .false.
        blqrad = -2. ! K/day
        blqmin = 4. ! g/kg
        blqmax = 6. ! g/kg
        bldivg = 3.e-6 ! 1/s
        
        ! boundary layer forcing (subsidence + nudging of free troposphere)
        doblforcing = .false.
        doblsubs = .false. ! option for specified subsidence
        blsubs = 3.e-3 ! m/s
        bltauz = 10800. ! nudging timescale [s]
        ! default divergence rate [1/s], same as above

        ! option within doblforcing for the use of a specified divergence
        !  (rather than a specified subsidence above the bl)
        dobldivg = .false.
        bldivg = 3.e-6 ! 1/s

        doblsubs2 = .false.
        blsubs2 = 0.04 ! Pa/s

        ! use weak temperature gradient to specify large-scale
        ! vertical velocity
        dowtg = .false.
        wtgtau = 7200. ! seconds
        doaddqvrestart = .false.
        doparameterizedwave = .false.
        second_mode_factor = 0. !1 is for 1000km, 2 is for 500km wavelength wave

        ! allow user to specify initial/maximum value of ncycle.
        ! This should allow the user to run with a relatively large
        ! timestep which SAM will break up (into ncycle pieces) to
        ! keep the cfl number close to the maximum allowable.
        ncycle0 = 1
        ncyclemax = 4
        ncyclemin = 1

        ! cloud drop sedimentation
        doclouddropsed = .false.

        douwpbl = .false.     ! Bretherton McCaa PBL scheme
        uwpbl_ndiv = 1
        douwcu = .false. ! UW shallow convection scheme
        uwcu_ndiv = 1
        fixedndiff = 0

        !  peters.  options for evolving SST and LAND mask.  vars in either
        !   simple_ocean or simple_land
        doevolveSSTforcing = .false.
        call setdefault_evolveparms(evolveSST)
        doxysoil_wet = .false.
        call setdefault_evolveparms(evolvesoilwet)

        dolandseamask=.false.
        doxysurf_rough = .false.
        maxsurfrough = 9.e10        ! no maximum value
        nstart_surfrough = 0.       ! start at max at first time step
        dointeractiveland = .false.
        tau_landdamping = -1.        ! don't do any land damping
        hml_land = 0.1               ! mixed layer depth for land, 10 cm
                                     ! ref: Hartmann, p 85
        ! END UW OPTIONS
        !==================================================================


!----------------------------------
!  Read namelist variables from the standard input:
!------------

open(55,file='./'//trim(case)//'/prm', status='old',form='formatted') 
read (55,PARAMETERS)
close(55)

!----------------------------------
!  Read namelist for uw options from same prm file:
!------------

open(55,file='./'//trim(case)//'/prm', status='old',form='formatted') 
read (UNIT=55,NML=UWOPTIONS,IOSTAT=ios)
if (ios.ne.0) write(*,*) 'No UW options included in prm file'
close(55)

! Test for incompatible options
if ((dobetaplane).and.(doxy)) then
   print*,'Cannot specify both dobetaplane and doxy'
   call task_abort()
end if

if ((dodoubleco2).and.(doradforcing)) then
   print*,'dodoubleco2 only works with interactive radiation.'
   call task_abort()
end if

if ((domoistadiabaticnudging).and.(.not.donudging_tq)) then
   print*,'domoistadiabaticnudging requires donudging_tq=.true.'
   call task_abort()
end if

if ((doequinox).and.(.not.doperpetual)) then
   print*,'doequinox requires doperpetual=.true.'
   call task_abort()
end if

! peters  restart options
if(nwriterestart.ne.-1) then
  k = nstat/nstatfrq
  k = k * nstatfrq     ! k is updated value of nstat in main.f90, after
                       ! setparm is called

  if(mod(nwriterestart,k).ne.0 .or. mod(nstop,nwriterestart).ne.0) then
  print*,'Error: nwriterestart was set but does not divide both nstat and nstop.'
  print*,'  Restart file will not be written and run will be lost.  Quitting.'
  call task_abort()
  end if
end if

if (dorestart_last.and.(nelapse.eq.999999999) &
     .and.(nelapsemin.eq.999999999)) then
   print*,'Only use dorestart_last in combination with nelapse or nelapsemin.'
   print*,'Otherwise, you risk losing an entire run because of a crash.'
   call task_abort()
end if

if (douwpbl .and. (.not.dosmagor)) then
   print*,'UWPBL requires dosmagor=.true.'
   call task_abort()
endif
   
   
if (douwcu .and. .not.douwpbl) then
   print*,'douwcu requires douwpbl=.true.'
   call task_abort()
endif
   
!if (dorceintercomparison) doinstantoutput=.true.

!------------------------------------
!  Set parameters 

if (do3Daverage) then
   doMSE=.true.
   nsaveMSEstart = nsave3Dstart
   nsaveMSEend = nsave3Dend
   nsaveMSE = nsave3D	
   navgMSE = nsave3D	
endif
	if(RUN2D) dy=dx

	if(RUN2D.and.YES3D.eq.1) then
	  print*,'Error: 2D run and YES3D is set to 1. Exitting...'
	  call task_abort()
	endif
	if(RUN3D.and.YES3D.eq.0) then
	  print*,'Error: 3D run and YES3D is set to 0. Exitting...'
	  call task_abort()
	endif

        if(nstat.lt.nstatfrq) then
	  print*,'Error: nstat is less than nstatfrq. Exitting...'
	  call task_abort()
	endif

	pi = acos(-1.)
	if(fcor.eq.-999.) fcor= 4*pi/86400.*sin(latitude0*pi/180.)*sqrt(betafactor)
	fcorz = sqrt(4.*(2*pi/(3600.*24.))**2*betafactor-fcor**2)	  
	coszrs = 0.637 ! default mean solar zenith angle
	
	if(ny.eq.1) dy=dx

        !kzm give a default nstep, na,nb,nc only when it is a new run
        !otherwise the restart would be in error!!
	if(nrestart.eq.0) then
           na = 1
           nb = 2
           nc = 3
           nstep = 0
           time = 0.
           dtn = dt
        endif

	a_bg = 1./(tbgmax-tbgmin)
	a_pr = 1./(tprmax-tprmin)
	a_gr = 1./(tgrmax-tgrmin)

	notopened2D = .true.
	notopened3D = .true.

        call isccp_tables_init()   ! initialize isccp tables
	call isccp_zero()

        !===============================================================
        ! UW ADDITIONS

!bloss    INITIALIZE kzm TRACER FLAGS
!bloss      I removed this from settracer() and put it here so that the
!bloss    various dotr* flags would be set within setparm().  This
!bloss    allows the tracer arrays to be dynamically allocated,
!bloss    saving a big chunk of memory (and reducing restart file
!bloss    size) when they're not in use.

        dotrz=.false.
        dotrx=.false.
        dotry=.false.
        dotrzz=.false.

        select case (tracer_option)

        case(0) ! tracer_option=0: no tracers
        case(1) ! tracer_option=1: put tracer trz in the subcloud layer and
                !                 trzz in the bulk troposphere (z<=10.75km)
           dotrz=.true.
        case(2) ! tracer_option=2: put tracer trz in the subcloud layer 
                !                 trzz in the bulk troposphere (z<=10.75km)
           dotrz=.true.
           dotrzz=.true.
        case(3) ! tracer_option=3: put tracer trz and trzz in the subcloud
                !                    layer (z<600m)	  	 
           dotrz=.true.
           dotrzz=.true.
        case(4) ! tracer_option=4: put tracer trz, try and trzz in the  bulk
                !                    troposphere (z<=10km) 
           dotrz=.true.
           dotry=.true.
           dotrzz=.true.
        case(5) ! tracer_option=5: put tracer trx below 2km, try from 2 to
                !                    10km, and trz from 10 to 15 km.
           dotrz=.true.
           dotrx=.true.
           dotry=.true.
           dotrzz=.true.
        case(6) ! tracer_option=6: constant, horizontally-uniform
           ! flux of tracer into domain at surface.  Tracer decays
           ! with specified timescale, itautrx
           dotrx=.true.
        case(99)! tracer_option=99: use the full suit of tracers, 
                !                     each is set to the coordinate. 
                !                     trzz is set to k^2 to track mixing.
                !                     This is useful for
                !                     backtrajectory calculations  
           dotrz=.true.
           dotrx=.true.
           dotry=.true.
           dotrzz=.true.
        case default
           if(masterproc) write(*,*) 'Invalid value for tracer_option ', &
                'specified in prm file' 
           call task_abort()
        end select

        !radiatively-active ozone tracer
	if(doactiveo3) dotro3=.true. 
	
	!convert the inverse timescale for tracers from 1/day to 1/sec
	itauo3 = itauo3/86400.  ! kzm scalar
	itautrx = itautrx/86400.  ! kzm scalar
	itautry = itautry/86400.  ! kzm scalar
	itautrz = itautrz/86400.  ! kzm scalar
	itautrzz = itautrzz/86400.  ! kzm scalar
	
        !kzm Apr8, 04, read in the relaxation time scale profile
        if (dotauz) then
           open(8,file='./'//trim(case)//'/tau',status='old',form='formatted') 

           do k=1,nzm      
              read(8,*) itauz(k)
           end do
           close (8)

           itauz=itauz/3600./24.!change from 1/day to 1/seconds
        end if
        !kzm Aug17, 05, read in the zonal mean circulation
        if (donudging_hadley) then
           open(8,file='./'//trim(case)//'/hadley',status='old',form='formatted') 
           do j=1,ny_gl
           do k=1,nzm 
              read(8,*) ug0x(j,k),vg0x(j,k),wg0x(j,k),tg0x(j,k),qg0x(j,k)
           end do
           enddo
           close (8)
           qg0x=qg0x/1.e3 !convert to kg/kg
        end if
!!$        if (doparameterizedwave) then
!!$           open(8,file='./'//trim(case)//'/thmodes',status='old',form='formatted') 
!!$            do k=1,nzm 
!!$              read(8,*) thmode1(k),thmode2(k)
!!$            enddo
!!$           close (8)
!!$        end if

        !bloss: move re-definition of nsaveMSE to here from main
        if(nsaveMSE.eq.0) nsaveMSE=nstat  ! nsaveMSE=nstat if not initialized

        ! peters radiative options
        if ((doannual).and.(.not.doperpetual)) then
          print*,'doannual requires doperpetual=.true.'
          call task_abort()
        end if
        if ((doannual).and.(doequinox)) then
          print*,'doannual and doequinox both set!'
          call task_abort()
        end if

        ! peters allocate MSE variables if needed
        ! moved this stuff and other error checking associated with MSE
        ! input to checkMSEinput   May 2006
        call checkMSEinput()

        ! peters.  check simple land flags and initialize some stuff
        call checkSSTinput()
        call checkLANDinput()

        ! peters always initialize simple land, even after restart
        call initialize_simple_land()                 ! set defaults
        if(dolandseamask) call set_landmask()         !  initialize landmask
        if(doxysurf_rough) call set_surfrough()      ! initialize xy z0
        if(doxysoil_wet)                                                   &
               call set_evolveforcing(evolvesoilwet,xysoil_wet(1:nx,1:ny))

        ! peters check dosurfacefix
        if(dosurfacefix .and. (SFC_FLX_FXD .or. dobulksfc)) then
          print*,'dosurfacefix cannot be set with SFC_FLX_FXD or dobulksfc'
          call task_abort()
        end if

        !kzm check that number of xpoints per processor can be divided by block_size
        if(mod(nx_gl/nsubdomains_x,block_size).ne.0) then
          print*,'number of xpoints per processor cannot be divided by block_size'
          call task_abort()
        end if

        !bloss: if doblforcing=.true, set doblsubs=.true.
        if(doblforcing) doblsubs = .true.

        !bloss: if using an alternate scheme for large-scale forcing, 
        !         set douwforcing=.true., since the use of douwforcing as 
        !         a general flag simplifies things in forcing.f90 
        if(doblradforcing) douwforcing = .true.
        if(dobldivg) douwforcing = .true.
        if(doblsubs) douwforcing = .true.
        if(doblsubs2) douwforcing = .true.
        if(dowtg) douwforcing = .true.
!        if(doparameterizedwave) douwforcing = .true.

        !bloss: use doblforcing as a general flag for these two options.
        if(dobldivg) doblforcing = .true.
        if(doblsubs) doblforcing = .true.
        if(doblsubs2) doblforcing = .true.


        !bloss: boundary layer specified radiative cooling
        if (doblradforcing) then
           blqrad = blqrad/86400. ! convert to K/s
           blqmin = blqmin/1.e3 ! convert to kg/kg
           blqmax = blqmax/1.e3 ! convert to kg/kg
        end if
	mpbetafactor=min(betafactor,mpbetafactor);
        ! END UW ADDITIONS
        !===============================================================

        call allocate_uwvars() ! allocate variables for uw additions

end
