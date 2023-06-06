	subroutine write_all()
	
	use vars
	use params
	implicit none
	character *4 rankchar
	character *256 filename
	integer irank


	if(restart_sep) then

          write(rankchar,'(i4)') rank

          filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart.bin'


          open(65,file=trim(filename), status='unknown',form='unformatted')
          write(65) nsubdomains

	  call write_statement
          if (.not.dooldrestart) call write_uwextras !bloss
          close(65)


	else
	  write(rankchar,'(i4)') nsubdomains
	  filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart.bin'

	  do irank=0,nsubdomains-1
	
	     call task_barrier()

	     if(irank.eq.rank) then

	       if(masterproc) then
	      
	        open(65,file=trim(filename), status='unknown',form='unformatted')
	        write(65) nsubdomains

	       else

                open(65,file=trim(filename), status='unknown',form='unformatted',&
                   position='append')

	       end if

               call write_statement

               close(65)

             end if
	  end do

          !===============================================================
          ! UW ADDITION

          if(.not.dooldrestart) then !bloss: append extra uw stuff to file
             do irank=0,nsubdomains-1
                call task_barrier()
                if(irank.eq.rank) then
                   open(65,file=trim(filename), status='unknown', &
                        form='unformatted', position='append')
                   call write_uwextras
                   close(65)
                end if
             end do
          end if
          !===============================================================

	end if ! restart_sep

	call task_barrier()

        return
        end
 
 
 
 
     
	subroutine read_all()
	
        use simple_ocean, only: readmlohflx
	use vars
	use params
	implicit none
	character *4 rankchar
	character *256 filename
	integer irank, ii, i, j, k

        integer k150, kb, kc !bloss for dowaveenergetics
        real press_total, tmp_wgt1, pii, press_top, tmp_wgt2


	if(restart_sep) then

           write(rankchar,'(i4)') rank

           filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart.bin'


           open(65,file=trim(filename), status='unknown',form='unformatted')
           read(65)

	   call read_statement
           if (.not.dooldrestart) call read_uwextras !bloss
           close(65)


	else

	  write(rankchar,'(i4)') nsubdomains

	  filename='./RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart.bin'
          open(65,file=trim(filename), status='unknown',form='unformatted')

	  do irank=0,nsubdomains-1
	
	     call task_barrier()

	     if(irank.eq.rank) then

	       read (65)
 
               do ii=0,irank-1 ! skip records
                 read(65)
	       end do

	       call read_statement

               !===============================================================
               ! UW ADDITION

               if(.not.dooldrestart) then 
                  do ii = 1,nsubdomains-1 ! bloss: skip additional records
                     read(65) 
                  end do
                   call read_uwextras
               end if
               !===============================================================

               close(65)

             end if

	  end do

	end if ! restart_sep

	call task_barrier()
	

        if(dooldrestart) then

	prec_xy(:,:) = 0.
	lhf_xy(:,:) = 0.
	shf_xy(:,:) = 0.
	lwns_xy(:,:) = 0.
	swns_xy(:,:) = 0.
	lwnsc_xy(:,:) = 0.
	swnsc_xy(:,:) = 0.
	lwnt_xy(:,:) = 0.
	swnt_xy(:,:) = 0.
	lwntc_xy(:,:) = 0.
	swntc_xy(:,:) = 0.
	pw_xy(:,:) = 0.
	cw_xy(:,:) = 0.
	iw_xy(:,:) = 0.
	u200_xy(:,:) = 0.
	v200_xy(:,:) = 0.
	usfc_xy(:,:) = 0.
	vsfc_xy(:,:) = 0.
        w500_xy(:,:) = 0.
        solin_xy(:,:) = 0.
        qocean_xy(:,:) = 0.

!==================================================================
! UW ADDITIONS

	psfc_xy = 0.
        h200_xy = 0.
        h850_xy = 0.
        u850_xy = 0.
        v850_xy = 0.
        the850_xy = 0.
        thes850_xy = 0.
        the1000_xy = 0.
        thes1000_xy = 0.

	swvp_xy(:,:) = 0.
	swntm_xy(:,:) = 0.
	wmode1_xy(:,:) = 0.
	wmode2_xy(:,:) = 0.
	qmode1_xy(:,:) = 0.
	qmode2_xy(:,:) = 0.
	thmode1_xy(:,:) = 0.
	thmode2_xy(:,:) = 0.
	wmode1i_xy(:,:) = 0.
	wmode2i_xy(:,:) = 0.
	hstor_xy(:,:) = 0.
	hadv_xy(:,:) = 0.
	sstor_xy(:,:) = 0.
	sadv_xy(:,:) = 0.

        cloudtopheight = 0.
        echotopheight = 0.
        cloudtoptemp = 0.

        if (docolumnbudgets) then 
           ! compute initial frozen moist static energy and liquid
           ! water ice static energy which will be used to compute
           ! the storage terms in the column-integrated hf and sli
           ! budgets 
           old_fmse = 0.
           old_sli  = 0.
           do k = 1,nzm
              do j = 1,ny
                 do i = 1,nx
                    old_fmse(i,j) = old_fmse(i,j) + rho(k)*dz*adz(k) &
                         *(cp*t(i,j,k) + lcond*(q(i,j,k) + qp(i,j,k)))
                    old_sli(i,j)  = old_sli(i,j) + rho(k)*dz*adz(k)*cp*t(i,j,k)
                 end do
              end do
           end do
        end if

        if (dowaveenergetics) then
           ! FIND 150 MBAR LEVEL -- ROUGH DEFINITION OF TOP OF TROPOSPHERE
           ! AND TOP OF BAROCLINIC MODES.
           k150 = nzm
           do k = 1,nzm
              kc=min(nzm,k+1)
              kb=max(1,k-1)
              if(pres(kc).le.150..and.pres(kb).gt.150.) k150=k
           end do

           ! INITIALIZE OLD MODE STRENGTHS TO ZERO
           oldw1_xy = 0.
           oldw2_xy = 0.
           olds1_xy = 0.
           olds2_xy = 0.
           oldh1_xy = 0.
           oldh2_xy = 0.

           pii = atan2(0.,-1.)
           press_total = (presi(1)-presi(k150+1))
           press_top   = presi(k150+1)

           do k = 1,k150
              ! COMPUTE PROJECTION OF W/S/H ONTO BAROCLINIC MODES
              tmp_wgt1 = 2.*sin(pii*(pres(k)-press_top)/press_total) &
                   *(presi(k)-presi(k+1))/press_total
              tmp_wgt2 = 2.*sin(-2.*pii*(pres(k)-press_top)/press_total) &
                   *(presi(k)-presi(k+1))/press_total
              do j = 1,ny
                 do i = 1,nx
                    oldw1_xy(i,j) = oldw1_xy(i,j) + tmp_wgt1*w(i,j,k)
                    oldw2_xy(i,j) = oldw2_xy(i,j) + tmp_wgt2*w(i,j,k)
                    olds1_xy(i,j) = olds1_xy(i,j) + tmp_wgt1*(t(i,j,k)-t0(k))
                    olds2_xy(i,j) = olds2_xy(i,j) + tmp_wgt2*(t(i,j,k)-t0(k))
                    oldh1_xy(i,j) = oldh1_xy(i,j) &
                         + tmp_wgt1*((t(i,j,k)-t0(k)) &
                         + fac_cond*(q(i,j,k)+qp(i,j,k)-q0(k)))
                    oldh2_xy(i,j) = oldh2_xy(i,j) &
                         + tmp_wgt2*((t(i,j,k)-t0(k)) &
                         + fac_cond*(q(i,j,k)+qp(i,j,k)-q0(k)))
                 end do
              end do
           end do

        end if

        if ((ocean_type.ge.3).and.(dodynamicocean)) call readmlohflx

        end if

! END UW ADDITIONS
!==================================================================

        return
        end
 
 
 
        subroutine write_statement()

        use vars
        use params
        implicit none

        if(dooldrestart) then

             write(65)  &
               rank, dx, dy, dz, adz, adzw, at, bt, ct, dt, dtn, dt3, time, &
               day, day0, nstep, na, nb, nc, nadams, rank, rank, caseid, case, &
!kzm add tracers
               u, v, w, t, q, qp, tro3, trx, try, trz, trzz, tke, tk_xy, tk_z, tkh_xy, tkh_z, p, tabs, qn, dudt, dvdt, dwdt,&
               fluxbu, fluxbv, fluxbt, fluxbq, fluxtu, fluxtv, fluxtt, fluxtq,&
               fzero, precsfc, t0, q0, qc0, qv0, qi0,tabs0, tv0, rel0, u0, v0,&
               p0, tg0, qg0, ug0, vg0, z, pres, rho, rhow, bet, gamaz, &
               t01, q01, wsub, qtend, ttend, prespot, tke0, betp,tauls, &
               dqls,dtls,ugls,vgls,wgls,pres0ls,dayls, zi, presi, &
               dtrfc,dayrfc,sstsfc,hsfc,lesfc,tausfc,daysfc, &
!               usnd,vsnd, & ! ADD WHEN CHANGING RESTART
               tsnd,qsnd,daysnd,nlsf,nrfc,nsfc,nsnd, &
               dodamping, doupperbound, docloud, doprecip, doradhomo, dosfchomo,&
               dolongwave, doshortwave, dosgs, dosubsidence, &
               docoriolis, dosurface, dolargescale,doradforcing, &
               dosfcforcing, doradsimple, donudging_uv, donudging_tq, &
               dosmagor, doscalar,doxy,dowallx,dowally, doperpetual, doseasons, &
               docup, docolumn, dt_cup, soil_wetness, &
               ttend_cup, qtend_cup, utend_cup, vtend_cup, &
               pres0, ug,vg,fcor,fcorz,tabs_s,z0,sstxy,fcory,fcorzy, &
               longitude,latitude,fluxt0,fluxq0,gamt0,gamq0, &
               tau0,timelargescale,gam3,gams1,gams2,gams3,gamr1,gamr2,gamr3,&
               gamg1, gamg2, gamg3,accrsc,accrsi,accrgc,accrgi,accrrc,coefice,&
               evaps1,evaps2,evapg1,evapg2,evapr1,evapr2, a_bg, a_pr, a_gr, &
               CEM, LES, OCEAN, LAND, SFC_FLX_FXD, SFC_FLX_FXD, SFC_TAU_FXD, &
               nx, ny, nz, grdf_x, grdf_y, grdf_z, &
               tkh_x_rf, tkh_y_rf, tkh_z_rf, u_i, v_i, w_i, t_i, q_i, qp_i !JY Yani janniy added
                

          else

             write(65)  &
               rank, dx, dy, dz, adz, adzw, at, bt, ct, dt, dtn, dt3, time, &
               day, day0, nstep, na, nb, nc, nadams, rank, rank, caseid, case, &
               u, v, w, t, q, qp, tke, tk_xy, tk_z, tkh_xy, tkh_z, p, tabs, qn, dudt, dvdt, dwdt,&
               fluxbu, fluxbv, fluxbt, fluxbq, fluxtu, fluxtv, fluxtt, fluxtq,&
               fzero, precsfc, t0, q0, qc0, qv0, qi0,tabs0, tv0, rel0, u0, v0,&
               p0, tg0, qg0, ug0, vg0, z, pres, rho, rhow, bet, gamaz, &
               t01, q01, wsub, qtend, ttend, prespot, tke0, betp,tauls, &
               dqls,dtls,ugls,vgls,wgls,pres0ls,dayls, zi, presi, &
               dtrfc,dayrfc,sstsfc,hsfc,lesfc,tausfc,daysfc, &
               usnd,vsnd,tsnd,qsnd,daysnd,nlsf,nrfc,nsfc,nsnd, &
               dodamping, doupperbound, docloud, doprecip, doradhomo, dosfchomo,&
               dolongwave, doshortwave, dosgs, dosubsidence, &
               docoriolis, dosurface, dolargescale,doradforcing, &
               dosfcforcing, doradsimple, donudging_uv, donudging_tq, &
               dosmagor, doscalar,doxy,dowallx,dowally, doperpetual, doseasons, &
               docup, docolumn, dt_cup, soil_wetness, dodynamicocean, ocean_type,&
               ttend_cup, qtend_cup, utend_cup, vtend_cup, &
               pres0, ug,vg,fcor,fcorz,tabs_s,z0,sstxy,fcory,fcorzy, &
               longitude,latitude,fluxt0,fluxq0,gamt0,gamq0, &
               tau0,timelargescale,gam3,gams1,gams2,gams3,gamr1,gamr2,gamr3,&
               gamg1, gamg2, gamg3,accrsc,accrsi,accrgc,accrgi,accrrc,coefice,&
               evaps1,evaps2,evapg1,evapg2,evapr1,evapr2, a_bg, a_pr, a_gr, &
               CEM, LES, OCEAN, LAND, SFC_FLX_FXD, SFC_FLX_FXD, SFC_TAU_FXD, &
               nx, ny, nz, grdf_x, grdf_y, grdf_z, &
               tkh_x_rf, tkh_y_rf, tkh_z_rf, u_i, v_i, w_i, t_i, q_i, qp_i !JY Yani janniy added


          end if

               if(rank.eq.nsubdomains-1) then
                  print *,'Restart file was saved. nstep=',nstep
               endif


        return
        end




        subroutine read_statement()

        use vars
        use params
        implicit none
        integer  nx1, ny1, nz1, rank1,k
        real addqv(nzm)
        character(40) case1,caseid1

        if(dooldrestart) then

             read(65) &
               rank1, dx, dy, dz, adz, adzw, at, bt, ct, dt, dtn, dt3, time, &
               day, day0, nstep, na, nb, nc, nadams, rank, rank, caseid1, case1, &
!kzm add tracers
               u, v, w, t, q, qp, tro3, trx, try, trz, trzz, tke, tk_xy, tk_z, tkh_xy, tkh_z, p, tabs, qn, dudt, dvdt, dwdt,&
               fluxbu, fluxbv, fluxbt, fluxbq, fluxtu, fluxtv, fluxtt, fluxtq,&
               fzero, precsfc, t0, q0, qc0, qv0, qi0,tabs0, tv0, rel0, u0, v0,&
               p0, tg0, qg0, ug0, vg0, z, pres, rho, rhow, bet, gamaz, &
               t01, q01, wsub, qtend, ttend, prespot, tke0, betp,tauls, &
               dqls,dtls,ugls,vgls,wgls,pres0ls,dayls, zi, presi, &
               dtrfc,dayrfc,sstsfc,hsfc,lesfc,tausfc,daysfc, &
!               usnd,vsnd, & ! ADD WHEN CHANGING RESTART
               tsnd,qsnd,daysnd,nlsf,nrfc,nsfc,nsnd, &
               dodamping, doupperbound, docloud, doprecip, doradhomo, dosfchomo, &
               dolongwave, doshortwave, dosgs, dosubsidence, &
               docoriolis, dosurface, dolargescale,doradforcing, &
               dosfcforcing, doradsimple, donudging_uv, donudging_tq, &
               dosmagor, doscalar,doxy,dowallx,dowally, doperpetual, doseasons, &
               docup, docolumn, dt_cup, soil_wetness, &
               ttend_cup, qtend_cup, utend_cup, vtend_cup, &
               pres0, ug,vg,fcor,fcorz,tabs_s,z0,sstxy,fcory,fcorzy, &
               longitude,latitude,fluxt0,fluxq0,gamt0,gamq0, &
               tau0,timelargescale,gam3,gams1,gams2,gams3,gamr1,gamr2,gamr3,&
               gamg1, gamg2, gamg3,accrsc,accrsi,accrgc,accrgi,accrrc,coefice,&
               evaps1,evaps2,evapg1,evapg2,evapr1,evapr2, a_bg, a_pr, a_gr, &
               CEM, LES, OCEAN, LAND, SFC_FLX_FXD, SFC_FLX_FXD, SFC_TAU_FXD, &
               nx1, ny1, nz1, grdf_x, grdf_y, grdf_z, &
               tkh_x_rf, tkh_y_rf, tkh_z_rf, u_i, v_i, w_i, t_i, q_i, qp_i !JY Yani janniy added


          else

             read(65) &
               rank1, dx, dy, dz, adz, adzw, at, bt, ct, dt, dtn, dt3, time, &
               day, day0, nstep, na, nb, nc, nadams, rank, rank, caseid1, case1, &
               u, v, w, t, q, qp, tke, tk_xy, tk_z, tkh_xy, tkh_z, p, tabs, qn, dudt, dvdt, dwdt,&
               fluxbu, fluxbv, fluxbt, fluxbq, fluxtu, fluxtv, fluxtt, fluxtq,&
               fzero, precsfc, t0, q0, qc0, qv0, qi0,tabs0, tv0, rel0, u0, v0,&
               p0, tg0, qg0, ug0, vg0, z, pres, rho, rhow, bet, gamaz, &
               t01, q01, wsub, qtend, ttend, prespot, tke0, betp,tauls, &
               dqls,dtls,ugls,vgls,wgls,pres0ls,dayls, zi, presi, &
               dtrfc,dayrfc,sstsfc,hsfc,lesfc,tausfc,daysfc, &
               usnd,vsnd,tsnd,qsnd,daysnd,nlsf,nrfc,nsfc,nsnd, &
               dodamping, doupperbound, docloud, doprecip, doradhomo, dosfchomo, &
               dolongwave, doshortwave, dosgs, dosubsidence, &
               docoriolis, dosurface, dolargescale,doradforcing, &
               dosfcforcing, doradsimple, donudging_uv, donudging_tq, &
               dosmagor, doscalar,doxy,dowallx,dowally, doperpetual, doseasons, &
               docup, docolumn, dt_cup, soil_wetness, dodynamicocean, ocean_type, &
               ttend_cup, qtend_cup, utend_cup, vtend_cup, &
               pres0, ug,vg,fcor,fcorz,tabs_s,z0,sstxy,fcory,fcorzy, &
               longitude,latitude,fluxt0,fluxq0,gamt0,gamq0, &
               tau0,timelargescale,gam3,gams1,gams2,gams3,gamr1,gamr2,gamr3,&
               gamg1, gamg2, gamg3,accrsc,accrsi,accrgc,accrgi,accrrc,coefice,&
               evaps1,evaps2,evapg1,evapg2,evapr1,evapr2, a_bg, a_pr, a_gr, &
               CEM, LES, OCEAN, LAND, SFC_FLX_FXD, SFC_FLX_FXD, SFC_TAU_FXD, &
               nx1, ny1, nz1, grdf_x, grdf_y, grdf_z, &
               tkh_x_rf, tkh_y_rf, tkh_z_rf, u_i, v_i, w_i, t_i, q_i, qp_i !JY Yani janniy added


          end if

             if(nx.ne.nx1.or.ny.ne.ny1.or.nz.ne.nz1) then
              print *,'Error: domain dims (nx,ny,nz) set by grid.f'
              print *,'do not correspond to ones in the restart file.'
              print *,'in executable:   nx, ny, nz:',nx,ny,nz
              print *,'in restart file: nx, ny, nz:',nx1,ny1,nz1
              print *,'Exiting...'
              call task_abort()
             endif
             if(doaddqvrestart) then 
            open(8,file='./'//trim(case)//'/addqv',status='old',form='formatted') 

!            if(rank.eq.nsubdomains/2-1) then
           do k=1,nzm      
              read(8,*) addqv(k)
              q(:,:,k)=q(:,:,k)+addqv(k)/1000.
           end do
           close (8)
                 print *, 'Add qv',addqv
!             endif
            endif
            if(rank.eq.nsubdomains-1) then
                 print *,'Case:',caseid
                 print *,'Restarting at step:',nstep
                 print *,'Time:',nstep*dt
              endif


        return
        end
        
        
        subroutine write_uwextras()

          use vars
          use params
          use mse    ! peters
          use simple_ocean ! peters
          use simple_land ! peters
          implicit none

          !bloss: If arrays are not allocated for optional uw features,
          ! allocate 1x1x1 arrays here as placeholders in the restart
          ! file.  This will make for a consistent restart file
          ! across versions of the model with different features, but
          ! will allow for much smaller restart files if not all
          ! features are enabled.
          call allocate_placeholder_arrays()

          write(65) qv,qc,qi,qr,qs,qg, &   ! Krueger microphysics
               iso_qv,iso_qc,iso_qi,iso_qr,iso_qs,iso_qg, &   ! isotopes
               iso2_qv,iso2_qc,iso2_qi,iso2_qr,iso2_qs,iso2_qg, &  
               trx, try, trz, trzz, tro3, & ! Tracers (kzm)
               prec_xy, lhf_xy, shf_xy, &
               lwns_xy, swns_xy, lwnsc_xy, swnsc_xy, &
               lwnt_xy, swnt_xy, lwntc_xy, swntc_xy, &
               pw_xy, cw_xy, iw_xy, &
               u200_xy, v200_xy, usfc_xy, vsfc_xy, &
               w500_xy, solin_xy, qocean_xy, psfc_xy, h200_xy, &
               h850_xy, u850_xy, v850_xy, the850_xy, thes850_xy, &
               the1000_xy, thes1000_xy, swvp_xy, swntm_xy, &
               wmode1_xy, wmode2_xy, qmode1_xy, qmode2_xy, &
               thmode1_xy, thmode2_xy, wmode1i_xy, wmode2i_xy, &
               hstor_xy, hadv_xy, sstor_xy, sadv_xy, old_fmse, old_sli, &
               oldw1_xy, oldw2_xy, wdwdt1_xy, wdwdt2_xy, &
               olds1_xy, olds2_xy, sdsdt1_xy, sdsdt2_xy, &
               oldh1_xy, oldh2_xy, hdhdt1_xy, hdhdt2_xy, &
               dodynamicocean, ocean_type, qoceanxy, &
               hml, Szero, deltaS, sst_target, delta_sst, &
               dofixedoceancooling, mlohflx, ssty, &
               dodoubleco2, co2factor, doequinox, iyear, &
               doclearsky, doinstantolr, doreflectivity, &      
               dosaturatedwvp, dowaveoutput, docolumnbudgets, & 
               do850mbarwinds, dooutputheights, dolowlevelthetae, &
               dorceintercomparison, doinstantoutput, &
               dotrx,dotry,dotrz,dotrzz, dotro3, doactiveo3, &
               doreset_tracer, dorestart_tracer,dotauz, tracer_option, &
               dokruegermicro,  dokruegereffectiveradii, &
               dosubsidence_tonly, domoistadiabaticnudging, &
               doperiodic_forcing, day0_forcing, period_forcing, &
               betafactor, dobetaplane, dowaterbubble,dobetafactor, &
               doydamping, dowaveenergetics, &
               douwpbl, douwcu, uwpbl_ndiv, uwcu_ndiv, pblh, tke_uw, &
               cush, tk_z_uw, tkh_z_uw, kvm_out_flip_3d, kvh_out_flip_3d, &
               domicroscaling, cmin, cmax,wavedampingtime,t_wavebg_hist, &
               q_wavebg_hist,heating_wavebg_hist,t_wave_local,q_wave_local, &
               h_wave_local,t_wave,q_wave,heating_wave,t_wavebg,q_wavebg, &
               heating_wavebg,nsteplinearwave,wavecounter,w_wave, &
               dtu_uu, dtu_uv, dtu_uw, dtv_vu, dtv_vv, dtv_vw, &  ! peters mse
               dtu_diff, dtv_diff, &   ! peters stuff for mse
               evolveSST%readflds, evolveSST%nctimestep, & ! peters
               evolvesoilwet%readflds, evolvesoilwet%nctimestep, & ! peters
               sstxy_dble   ! peters


          call deallocate_placeholder_arrays()

        end subroutine write_uwextras

        subroutine read_uwextras()

          use vars
          use params
          use mse
          use simple_ocean
          use simple_land
          implicit none

          !bloss: If arrays are not allocated for optional uw features,
          ! allocate 1x1x1 arrays here as placeholders in the restart
          ! file.  This will make for a consistent restart file
          ! across versions of the model with different features, but
          ! will allow for much smaller restart files if not all
          ! features are enabled.
          call allocate_placeholder_arrays()

          read(65) qv,qc,qi,qr,qs,qg, &   ! Krueger microphysics
               iso_qv,iso_qc,iso_qi,iso_qr,iso_qs,iso_qg, &   ! isotopes
               iso2_qv,iso2_qc,iso2_qi,iso2_qr,iso2_qs,iso2_qg, &  
               trx, try, trz, trzz, tro3, & ! Tracers (kzm)
               prec_xy, lhf_xy, shf_xy, &
               lwns_xy, swns_xy, lwnsc_xy, swnsc_xy, &
               lwnt_xy, swnt_xy, lwntc_xy, swntc_xy, &
               pw_xy, cw_xy, iw_xy, &
               u200_xy, v200_xy, usfc_xy, vsfc_xy, &
               w500_xy, solin_xy, qocean_xy, psfc_xy, h200_xy, &
               h850_xy, u850_xy, v850_xy, the850_xy, thes850_xy, &
               the1000_xy, thes1000_xy, swvp_xy, swntm_xy, &
               wmode1_xy, wmode2_xy, qmode1_xy, qmode2_xy, &
               thmode1_xy, thmode2_xy, wmode1i_xy, wmode2i_xy, &
               hstor_xy, hadv_xy, sstor_xy, sadv_xy, old_fmse, old_sli, &
               oldw1_xy, oldw2_xy, wdwdt1_xy, wdwdt2_xy, &
               olds1_xy, olds2_xy, sdsdt1_xy, sdsdt2_xy, &
               oldh1_xy, oldh2_xy, hdhdt1_xy, hdhdt2_xy, &
               dodynamicocean, ocean_type, qoceanxy, &
               hml, Szero, deltaS, sst_target, delta_sst, &
               dofixedoceancooling, mlohflx, ssty, &
               dodoubleco2, co2factor, doequinox, iyear, &
               doclearsky, doinstantolr, doreflectivity, &      
               dosaturatedwvp, dowaveoutput, docolumnbudgets, & 
               do850mbarwinds, dooutputheights, dolowlevelthetae, &
               dorceintercomparison, doinstantoutput, &
               dotrx,dotry,dotrz,dotrzz, dotro3, doactiveo3, &
               doreset_tracer, dorestart_tracer,dotauz, tracer_option, &
               dokruegermicro,  dokruegereffectiveradii, &
               dosubsidence_tonly, domoistadiabaticnudging, &
               doperiodic_forcing, day0_forcing, period_forcing, &
               betafactor, dobetaplane, dowaterbubble,dobetafactor, &
               doydamping, dowaveenergetics, &
               douwpbl, douwcu, uwpbl_ndiv, uwcu_ndiv, pblh, tke_uw, &
               cush, tk_z_uw, tkh_z_uw, kvm_out_flip_3d, kvh_out_flip_3d,&
               domicroscaling, cmin, cmax,wavedampingtime,t_wavebg_hist, &
               q_wavebg_hist,heating_wavebg_hist,t_wave_local,q_wave_local, &
               h_wave_local,t_wave,q_wave,heating_wave,t_wavebg,q_wavebg, &
               heating_wavebg,nsteplinearwave,wavecounter,w_wave,   &
               dtu_uu, dtu_uv, dtu_uw, dtv_vu, dtv_vv, dtv_vw, &  ! peters mse
               dtu_diff, dtv_diff, &   ! peters stuff for mse
               evolveSST%readflds, evolveSST%nctimestep, & ! peters
               evolvesoilwet%readflds, evolvesoilwet%nctimestep, & ! peters
               sstxy_dble   ! peters

          call deallocate_placeholder_arrays()

        end subroutine read_uwextras

        subroutine allocate_placeholder_arrays
          use vars
          use mse
          implicit none

          ! arrays for lin/lord-style microphysics
          if(.not.ALLOCATED(qv)) then
             ALLOCATE(qv(1,1,1)); qv=0.
          end if
          if(.not.ALLOCATED(qc)) then
             ALLOCATE(qc(1,1,1)); qc=0.
          end if
          if(.not.ALLOCATED(qi)) then
             ALLOCATE(qi(1,1,1)); qi=0.
          end if
          if(.not.ALLOCATED(qr)) then
             ALLOCATE(qr(1,1,1)); qr=0.
          end if
          if(.not.ALLOCATED(qs)) then
             ALLOCATE(qs(1,1,1)); qs=0.
          end if
          if(.not.ALLOCATED(qg)) then
             ALLOCATE(qg(1,1,1)); qg=0.
          ! arrays for isotopes
          end if
          if(.not.ALLOCATED(iso_qv)) then
             ALLOCATE(iso_qv(1,1,1)); iso_qv=0.
          end if
          if(.not.ALLOCATED(iso_qc)) then
             ALLOCATE(iso_qc(1,1,1)); iso_qc=0.
          end if
          if(.not.ALLOCATED(iso_qi)) then
             ALLOCATE(iso_qi(1,1,1)); iso_qi=0.
          end if
          if(.not.ALLOCATED(iso_qr)) then
             ALLOCATE(iso_qr(1,1,1)); iso_qr=0.
          end if
          if(.not.ALLOCATED(iso_qs)) then
             ALLOCATE(iso_qs(1,1,1)); iso_qs=0.
          end if
          if(.not.ALLOCATED(iso_qg)) then
             ALLOCATE(iso_qg(1,1,1)); iso_qg=0.
          end if
          if(.not.ALLOCATED(iso2_qv)) then
             ALLOCATE(iso2_qv(1,1,1)); iso2_qv=0.
          end if
          if(.not.ALLOCATED(iso2_qc)) then
             ALLOCATE(iso2_qc(1,1,1)); iso2_qc=0.
          end if
          if(.not.ALLOCATED(iso2_qi)) then
             ALLOCATE(iso2_qi(1,1,1)); iso2_qi=0.
          end if
          if(.not.ALLOCATED(iso2_qr)) then
             ALLOCATE(iso2_qr(1,1,1)); iso2_qr=0.
          end if
          if(.not.ALLOCATED(iso2_qs)) then
             ALLOCATE(iso2_qs(1,1,1)); iso2_qs=0.
          end if
          if(.not.ALLOCATED(iso2_qg)) then
             ALLOCATE(iso2_qg(1,1,1)); iso2_qg=0.
          ! arrays for kzm tracers
          end if
          if(.not.ALLOCATED(trx)) then
             ALLOCATE(trx(1,1,1)); trx=0.
          end if
          if(.not.ALLOCATED(try)) then
             ALLOCATE(try(1,1,1)); try=0.
          end if
          if(.not.ALLOCATED(trz)) then
             ALLOCATE(trz(1,1,1)); trz=0.
          end if
          if(.not.ALLOCATED(trzz)) then
             ALLOCATE(trzz(1,1,1)); trzz=0.
          end if
          if(.not.ALLOCATED(tro3)) then
             ALLOCATE(tro3(1,1,1)); tro3=0.
          end if
          ! arrays for peters mse
          if(.not.isallocateMSE) then
            allocate(dtu_uu(1,1,1,1), dtu_uv(1,1,1,1), dtu_uw(1,1,1,1), &
                     dtv_vu(1,1,1,1), dtv_vv(1,1,1,1), dtv_vw(1,1,1,1), &
                     dtu_diff(1,1,1,1), dtv_diff(1,1,1,1)             )
            dtu_uu=0; dtu_uv=0; dtu_uw=0
            dtv_vu=0; dtv_vv=0; dtv_vw=0
            dtu_diff=0; dtv_diff=0
          end if


        end subroutine allocate_placeholder_arrays

        subroutine deallocate_placeholder_arrays
          use vars
          use mse
          implicit none

          ! De-allocate any 1x1x1 arrays allocated above
          ! arrays for lin/lord-style microphysics
          if(SIZE(qv,1).eq.1) DEALLOCATE(qv)
          if(SIZE(qc,1).eq.1) DEALLOCATE(qc)
          if(SIZE(qi,1).eq.1) DEALLOCATE(qi)
          if(SIZE(qr,1).eq.1) DEALLOCATE(qr)
          if(SIZE(qs,1).eq.1) DEALLOCATE(qs)
          if(SIZE(qg,1).eq.1) DEALLOCATE(qg)
          ! arrays for isotopes
          if(SIZE(iso_qv,1).eq.1) DEALLOCATE(iso_qv)
          if(SIZE(iso_qc,1).eq.1) DEALLOCATE(iso_qc)
          if(SIZE(iso_qi,1).eq.1) DEALLOCATE(iso_qi)
          if(SIZE(iso_qr,1).eq.1) DEALLOCATE(iso_qr)
          if(SIZE(iso_qs,1).eq.1) DEALLOCATE(iso_qs)
          if(SIZE(iso_qg,1).eq.1) DEALLOCATE(iso_qg)
          if(SIZE(iso2_qv,1).eq.1) DEALLOCATE(iso2_qv)
          if(SIZE(iso2_qc,1).eq.1) DEALLOCATE(iso2_qc)
          if(SIZE(iso2_qi,1).eq.1) DEALLOCATE(iso2_qi)
          if(SIZE(iso2_qr,1).eq.1) DEALLOCATE(iso2_qr)
          if(SIZE(iso2_qs,1).eq.1) DEALLOCATE(iso2_qs)
          if(SIZE(iso2_qg,1).eq.1) DEALLOCATE(iso2_qg)
          ! arrays for kzm tracers
          if(SIZE(trx,1).eq.1) DEALLOCATE(trx)
          if(SIZE(try,1).eq.1) DEALLOCATE(try)
          if(SIZE(trz,1).eq.1) DEALLOCATE(trz)
          if(SIZE(trzz,1).eq.1) DEALLOCATE(trzz)
          if(SIZE(tro3,1).eq.1) DEALLOCATE(tro3)
          ! arrays for peters MSE
          if(.not.isallocateMSE) then
            deallocate(dtu_uu, dtu_uv, dtu_uw, &
                       dtv_vu, dtv_vv, dtv_vw, &
                       dtu_diff, dtv_diff      )
          end if

        end subroutine deallocate_placeholder_arrays

