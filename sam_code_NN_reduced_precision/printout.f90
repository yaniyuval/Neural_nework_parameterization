	subroutine printout()

	use vars
	use params
        use simple_land
        use simple_ocean, only : doevolveSSTforcing
	implicit none
	
	write(6,*)
	write(6,*)'Case: ',case
	write(6,*)'Caseid: ',caseid
	if(nrestart.eq.0) then
	  write(6,*)'New Run.'
	else
	  write(6,*)'Restart. nrestart=',nrestart
	endif
	write(6,*)'Day:',day
	write(6,*)'Latitude:',latitude0
	write(6,*)'Longitude:',longitude0

        if((.not.(OCEAN.or.LAND)) .and. (.not.dolandseamask)) then
	 write(6,*) 'Neither OCEAN nor LAND are set. Exitting...'
	 call task_abort()
	endif
	if(OCEAN.and.LAND) then
	 write(6,*) 'Both OCEAN and LAND are set. Confused...'
	 call task_abort()
	endif
	if(.not.(CEM.or.LES)) then
	 write(6,*) 'Neither CEM nor LES are set. Exitting...'
	 call task_abort()
	endif
	if(CEM.and.LES) then
	 write(6,*) 'Both CEM and LES are set. Confused...'
	 call task_abort()
	endif

	if(LES) print*,'Model type: LES'
	if(CEM) print*,'Model type: CEM'
	
	write(6,*) 'Finish at timestep:',nstop
	write(6,*) 'Finish on day:',day+(nstop-nstep)*dt/3600./24.
	write(6,*)
	write(6,*) 'Statistics file ouput frequency: ',nstat,' steps'
	write(6,*) 'Statistics file sampling: every ',nstatfrq,' steps'
	write(6,*) 'printouts frequency:',nprint,' steps'
	if(nstop-nstep.lt.nstat) then
	  write(6,*) 'Error: job will finish before statistics is done'
	  call task_abort()
	endif
	
	if(nadams.eq.2.or.nadams.eq.3) then
	   write(6,*) 'Adams-Bashforth scheme order:',nadams
	else
	   write(6,*) 'Error: nadams =',nadams
	   call task_abort()
	endif 
	  	
	write(6,*)
	if(dx.gt.0.and.dy.gt.0.and.dz.gt.0. ) then
	   write(6,*) 'Global Grid:',nx_gl,ny_gl,nz_gl
	   write(6,*) 'Local Grid:',nx,ny,nzm
	   write(6,*) 'Grid spacing (m) dx, dy, dz:',dx, dy, dz
	   write(6,*) 'Domain dimensions (m):',	dx*nx_gl,dy*ny_gl,z(nzm)
	else
	   write(6,*) 'Error: grid spacings dx, dy, dz:',dx, dy, dz
	   call task_abort()
	endif
	write(6,*)
	if(dt.gt.0) then
	   write(6,*) 'Timestep (sec):',dt
	else
	   write(6,*) 'Error: dt =',dt
	   call task_abort()
	endif   	
	write(6,*)
	write(6,*) 'do convective parameterization', docup
	write(6,*) 'do column model', docolumn
	if(docup) then
	if(dt_cup.gt.0) then
	   write(6,*) 'Timestep for cconv parameterization (sec):',dt_cup
	else
	   write(6,*) 'Error: dt_cup =',dt_cup
	   call task_abort()
	endif   	
	end if

	write(6,*) 'do spange damping at the domain top: ', dodamping
	write(6,*) 'maintain grad. of scalars at the top:',doupperbound
	write(6,*) 'clouds are allowed:',docloud
	write(6,*) 'precipitation is allowed:',doprecip.and.docloud
	dosmagor=dosgs.and.dosmagor
	write(6,*) 'SGS scheme is on:',dosgs
	if(dosgs) then
	  write(6,*) 'Smagorinsky SGS closure is on:',dosmagor
	end if
	if(.not.dosmagor.and.doscalar) then
           print*,'Error: dosmagor should be TRUE for doscalar be TRUE.'
           print*,'Aborting...'
           call task_abort()       
        end if
	write(6,*) 'larger-scale subsidence is on:',dosubsidence
	write(6,*) 'larger-scale tendency is on:',dolargescale
	write(6,*) 'coriolis force is allowed:',docoriolis	
	if(docoriolis) then
		write(6,*) '   Coriolis parameter (1/s):',fcor
		write(6,*) '   Vertical Coriolis parameter (1/s):',fcorz
	endif	
        if(doradforcing.and.(dolongwave.or.doshortwave)) then
          write(6,*) 'prescribed rad. forcing and radiation '// &
          'calculations cannot be done at the same time.'
          call task_abort()
        endif
        if(dolongwave) write(6,*) 'longwave radiation:',dolongwave  
        if(doshortwave) then
            write(6,*) 'shortwave radiation:',doshortwave
            write(6,*) 'do seasonal solar cycle:',doseasons
            write(6,*) 'do perpetual sun:',doperpetual
	endif
        if(doradforcing) write(6,*) 'radiation forcing is prescribed'
	write(6,*) 'surface flux parameterization is on:',dosurface
	if(dosurface) then
	    if(LAND) then
               print*,'Surface type: LAND'
               print*,'soil_wetness=',soil_wetness
               print*,'z0=',z0
            end if
	    if(OCEAN) print*,'Surface type: OCEAN'
	    write(6,*) ' sensible heat flux prescribed:',SFC_FLX_FXD
	    write(6,*) ' latent heat flux prescribed:',SFC_FLX_FXD
	    write(6,*) ' surface stress prescribed:',SFC_TAU_FXD
	endif
	
        if(dolargescale.or.dosubsidence) then
          if(    day.lt.dayls(1) &
            .or.day+(nstop-nstep)*dt/86400..gt.dayls(nlsf)) then
             print*,'Error: simulation time (from start to stop)'// &
              'can be beyond the l.s. forcing intervals'
             print*,'current day=',day
             print*,'stop day=',day+(nstop-nstep)*dt/86400.
             print*,'ls forcing: start =',dayls(1)
             print*,'ls forcing:   end =',dayls(nlsf)
             call task_abort()
          endif
        endif
        if(dosurface.and.dosfcforcing) then
          print*,'surface temperature isss prescribed'
          if(dodynamicocean) then
	     print*, 'dodynamicocean cannot be set to T'// &
                     'when dosfcforcing is also T'
	     call task_abort()
	  end if
          if(    day.lt.daysfc(1) &
             .or.day+(nstop-nstep)*dt/86400..gt.daysfc(nsfc))then
             print*,'Error: simulation time (from start to stop)'// &
              'can be beyond the sfc forcing intervals'
             print*,'current day=',day
             print*,'stop day=',day+(nstop-nstep)*dt/86400.
             print*,'sfc forcing:start =',daysfc(1)
             print*,'sfc forcing:  end =',daysfc(nsfc)
             call task_abort()
          endif
        endif
        if(doradforcing) then
          if ( day.lt.dayrfc(1) &
            .or.day+(nstop-nstep)/86400.*dt.gt.dayrfc(nrfc))then
             print*,'Error: simulation time (from start to stop)'// &
              'can be beyond the rad. forcing intervals'
             print*,'current day=',day
             print*,'stop day=',day+(nstop-nstep-1)*dt/86400.
             print*,'rad forcing:start =',dayrfc(1)
             print*,'rad forcing:  end =',dayrfc(nrfc)
             call task_abort()
          endif
        endif

	if(donudging_uv) print*, 'Nudging of U and V:', donudging_uv
	if(donudging_uv) print*,'tauls = ',tauls
	if(donudging_tq) print*, 'Nudging of T and Q:', donudging_tq

        print*,'dodynamicocean =',dodynamicocean
	print*,'doxy = ', doxy

        print*,'douwpbl', douwpbl
        print*,'douwcu', douwcu
        print*,'uwpbl_ndiv', uwpbl_ndiv
        print*,'uwcu_ndiv', uwcu_ndiv

        print*,'doevolveSSTforcing=',doevolveSSTforcing
        print*,'dolandseamask=',dolandseamask
        print*,'dointeractiveland=',dointeractiveland
        print*,'doxysurf_rough=',doxysurf_rough
        print*,'doxysoil_wet=',doxysoil_wet
        print*,'dosurfacefix=',dosurfacefix
	return
	end

