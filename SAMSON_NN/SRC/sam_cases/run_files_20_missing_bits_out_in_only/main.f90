program crm

!       Main module.

use vars
use hbuffer
use mse       ! peters
use params, only : lcond   ! peters
use uwpbl, only : uwpbldrv
use microscaling ! JA
!use random_forest_param ! pog
use nn_convection_flux_mod
use nn_diffusion_mod
use random_forest_param_diffusion ! janniy
!use random_forest_param_tkh       ! janniy



!use random_forest_param_del ! pog

!use diffuse_scalar_vert_flag_mod !Yani

implicit none

!public :: diffuse_scalar_vert_flag !Yani

integer k, icyc, nn, nstatsteps, i,j
real rrr, ranf_
real dummy(nzm)
double precision dtime !bloss  to lose rounding error in long runs.
double precision cputime, oldtime, init_time, elapsed_time !bloss wallclocktime
real dtaccum

firststep = .true. !bloss: when (re)starting, initialize necessary stuff
!-------------------------------------------------------------------
! determine the rank of the current task and of the neighbour's ranks

call task_init() 
!------------------------------------------------------------------
! print time, version, etc


if(masterproc) call header()	
!------------------------------------------------------------------
! Get initial time of job

call secondsf(init_time)
!------------------------------------------------------------------

call init()     ! initialize some statistics arrays
call setparm()	! set all parameters and constants
if(doparameterizedwave) call init_linear_wave
!------------------------------------------------------------------
! Initialize or restart from the save-dataset:
if(nrestart.eq.0) then
   day=day0 
   call setgrid() ! initialize vertical grid structure
   call setdata() ! initialize all variables

   call setforcing()

elseif(nrestart.eq.1) then
   call read_all()
elseif(nrestart.eq.2) then
   call read_all()
   call setparm() ! overwrite the parameters
   call setforcing()
else
   print *,'Error: confused by value of NRESTART'
   call task_abort() 
endif

if(dorandomforest) then
! call random_forest_init()
 call nn_convection_flux_init()
 if(.not.rf_uses_qp) then
  qp = 0 ! qp is not predicted 
 endif
endif

if(do_z_diffusion_rf) then
 call random_forest_init_diffusion()
endif

if(do_tkh_z) then
 !call random_forest_init_tkh()
 call nn_diffusion_init()
endif


!===================================================================
!UW ADDITIONS

!kzm initialize the tracers..
if(tracer_option.ne.0) call settracer()
  
!bloss: initialize q, qn, qp if using krueger microphysics
if (dokruegermicro) then
   q(1:nx,1:ny,1:nzm) = qv(1:nx,1:ny,1:nzm) + qc(1:nx,1:ny,1:nzm) &
        + qi(1:nx,1:ny,1:nzm)
   qn(1:nx,1:ny,1:nzm) = qc(1:nx,1:ny,1:nzm) + qi(1:nx,1:ny,1:nzm)
   qp(1:nx,1:ny,1:nzm) = qr(1:nx,1:ny,1:nzm) + qs(1:nx,1:ny,1:nzm) &
        + qg(1:nx,1:ny,1:nzm)
end if
!===================================================================

if(masterproc) call printout()
!------------------------------------------------------------------
!  Initialize statistics buffer:


call hbuf_init()
	
!------------------------------------------------------------------
nstatis = nstat/nstatfrq
nstat = nstatis * nstatfrq
nstatsteps = 0
!------------------------------------------------------------------
!   Main time loop    
!------------------------------------------------------------------
call secondsf(oldtime)

dtime = dble(time) !bloss: model time in double precision
totalcycles = 0 !bloss: track total number of cycles as well as time steps

do while(nstep.lt.nstop.and.nelapse.gt.0) 
         
  nstep = nstep + 1
  dtime = dtime + dt  !bloss: add increment to double precision time counter
  time  = dtime
  day = day0 + dtime/86400.
!bloss  time = time + dt
!bloss  day = day0 + time/86400.
  nelapse = nelapse - 1
!!$  do i = 1, nx
!!$     do j = 1, ny
!!$        rrr=(1.-2.*ranf_())
!!$        uwcu_pert_cldsrc(i,j)=0.81*uwcu_pert_cldsrc(i,j)+0.19*rrr
!!$     enddo
!!$  enddo
!------------------------------------------------------------------
!  Check if the dynamical time step should be decreased 
!  to handle the cases when the flow being locally linearly unstable
!------------------------------------------------------------------

!kzm renew parts of the tracer fields
  if(tracer_option.ne.0) call renewtracer()
  
  if(doperturb_realtime) call setperturb_realtime()

  call kurant() !bloss: ncycle is set within kurant()

  do icyc=1,ncycle

      icycle = icyc
      dtn = dt/ncycle
      dt3(na) = dtn
      dtfactor = dtn/dt
      
      totalcycles = totalcycles + 1 !bloss

     if(mod(nstep,nstatis).eq.0.and.icycle.eq.ncycle) then
        nstatsteps = nstatsteps + 1
        dostatis = .true.
        if(masterproc) print *,'Collecting statistics...'
     else
        dostatis = .false.
     endif
     !bloss: make special statistics flag for radiation,
     !       since it's only updated at icycle==1.
     dostatisrad = .false.
     if(mod(nstep,nstatis).eq.0.and.icycle.eq.1) dostatisrad = .true.

!---------------------------------------------
! UW ADDITIONS

     ! peters decide if time to collect MSE variables
     if (nsaveMSE.ge.0) call docollectMSE()

!---------------------------------------------
!  	the Adams-Bashforth scheme in time

     call abcoefs()
 
!---------------------------------------------
!  	initialize stuff: 
	
     call zero()

!-----------------------------------------------------------
!       Buoyancy term:
	     
     call buoyancy()
!------------------------------
do k=1,nzm
     do j=1,ny
       do i=1,nx
         u_i(i,j,k) = u(i,j,k)
         v_i(i,j,k) = v(i,j,k)
         w_i(i,j,k) = w(i,j,k)
         t_i(i,j,k) = tabs(i,j,k) ! - commented becasue  I wanted to run rf with t and not tabs.... 
         q_i(i,j,k) = q(i,j,k)
	 qp_i(i,j,k) = qp(i,j,k)
       end do
     end do
    end do


!------------------------------------------------------------
!       Large-scale and surface forcing:

     call forcing()
     if(doparameterizedwave) call wavesubsidence()

!----------------------------------------------------------
!       Nadging:

     call nudging()

!----------------------------------------------------------
!   	suppress turbulence near the upper boundary (spange):

     if(dodamping) call damping()

!-----------------------------------------------------------
!       Handle upper boundary for scalars

     if(doupperbound) call upperbound()

!-----------------------------------------------------------
!bloss: forcing/nudging/damping/upperbound above applied to q.
!         when using krueger microphysics, apply these to qv
!         instead since qv (not q) has prognostic eqn.
     if (dokruegermicro) qv(1:nx,1:ny,1:nzm) = q(1:nx,1:ny,1:nzm) &
          - qc(1:nx,1:ny,1:nzm) - qi(1:nx,1:ny,1:nzm)

!----------------------------------------------------------
!      Update the subdomain's boundaries for velocity

     call boundaries(0)

!---------------------------------------------------------
!	SGS TKE equation:     	
      if(dosgs) call tke_full()


!---------------------------------------------------------
!        Update boundaries for the SGS exchange coefficients:

     call boundaries(4)

!-----------------------------------------------
!       advection of momentum:

     call advect_mom()

!-----------------------------------------------
!   	surface fluxes:

     if(dosurface) then

       if(dosurfacefix) call boundaries(11)

       call surface()
      
     end if

!----------------------------------------------------------
!	SGS diffusion of momentum:

     if(dosgs) call diffuse_mom()

!-----------------------------------------------------------
!       Coriolis force:
	     
     if(docoriolis) call coriolis()
	 
!----------------------------------------------------------

     if(douwpbl) then
!        call orig_uwpbldrv( )
        call uwpbldrv( )

!!$       if(douwpbl.and.(uwpbl_ndiv.gt.1)) then
!!$        tk_z(1:nx, 1:ny, :) = tk_z(1:nx, 1:ny, 1:nzm) + tk_z_uw(:, :, :)*sqrt(betafactor)
!!$        tkh_z(1:nx, 1:ny, :) = tkh_z(1:nx, 1:ny, 1:nzm) + tkh_z_uw(:, :, :)*sqrt(betafactor)
!!$       endif

! **** temp by cw --- needed if updating velocities in uwpbl
     call boundaries(0)
! ****

     endif

!---------------------------------------------------------
!       compute rhs of the Poisson equation and solve it for pressure. 

     call pressure()

!---------------------------------------------------------
!       find velocity field at n+1/2 timestep needed for advection of scalars:
	 
     call adams()


!----------------------------------------------------------
!     Update boundaries for velocity fields to use for advection of scalars:

     call boundaries(1)

!        Update boundaries for scalars:

     call boundaries(2)

!---------------------------------------------------------
!      advection of scalars :
     if(doMSE) call divergence()  ! peters compute divergence if needed
     mseflag=1 ! peters
     call advect_scalar(t,tadv,twle,t2leadv,t2legrad,twleadv,.true.)
     
     mseflag=2 ! peters
     if(.not.dokruegermicro) then !bloss
     call advect_scalar(q,qadv,qwle,q2leadv,q2legrad,qwleadv,.true.)
     !==============================================================
     elseif(dokruegermicro) then !bloss
      call advect_scalar(qv,qadv,qwle,q2leadv,q2legrad,qwleadv,.true.)
     !==============================================================
     end if

     mseflag=0 ! peters
     if(dosgs.and..not.dosmagor) then
      call advect_scalar(tke,dummy,tkewle,dummy,dummy,dummy,.false.)
     else if(doscalar) then
      call advect_scalar(tke,dummy,tkewle,s2leadv,s2legrad,swleadv,.true.)
     end if

     mseflag=2 ! peters
     if(docloud.and.doprecip) then

       if((.not.dorandomforest).or.(dorandomforest.and.rf_uses_qp)) then

        if(.not.dokruegermicro) then !bloss
        call advect_scalar(qp,qpadv,qpwle,dummy,dummy,dummy,.false.)
        !==============================================================
        elseif(dokruegermicro) then !bloss
         call advect_scalar(qc,qcadv,qcwle,dummy,dummy,dummy,.false.)
         call advect_scalar(qi,qiadv,qiwle,dummy,dummy,dummy,.false.)
         if(dostatis) qadv = qadv + qcadv + qiadv
         if(dostatis) qwle = qwle + qcwle + qiwle
         call advect_scalar(qr,qradv,qrwle,dummy,dummy,dummy,.false.)
         call advect_scalar(qg,qgadv,qgwle,dummy,dummy,dummy,.false.)
         call advect_scalar(qs,qsadv,qswle,dummy,dummy,dummy,.false.)

         ! update q, qp
         q(1:nx,1:ny,1:nzm) = qv(1:nx,1:ny,1:nzm) + qc(1:nx,1:ny,1:nzm) &
              + qi(1:nx,1:ny,1:nzm)
         qp(1:nx,1:ny,1:nzm) = qr(1:nx,1:ny,1:nzm) + qs(1:nx,1:ny,1:nzm) &
              + qg(1:nx,1:ny,1:nzm)
        !==============================================================
        end if


        do k=1,nzm
          total_water_prec = total_water_prec + &
                         sum(qp(1:nx,1:ny,k))*adz(k)*dz *rho(k)
        end do

       endif
 
       if((.not.dorandomforest).or.(dorandomforest.and.rf_uses_qp)) then
        call precip_fall()
       endif


       if((.not.dorandomforest).or.(dorandomforest.and.rf_uses_qp)) then
        ! update qp for computation of precip part of total water budget
        if (dokruegermicro) qp(1:nx,1:ny,1:nzm) = qr(1:nx,1:ny,1:nzm) &
             + qs(1:nx,1:ny,1:nzm) + qg(1:nx,1:ny,1:nzm)

        do k=1,nzm
          total_water_prec = total_water_prec - &
                         sum(qp(1:nx,1:ny,k))*adz(k)*dz *rho(k)
        end do
       endif

     endif	

!==============================================================
!    kzm advect some tracers

     mseflag=0 ! peters
     if(dotrz)  call advect_scalar(trz,trzadv,trzwle,trz2leadv&
          &,trz2legrad,trzwleadv,.true.) 
     if(dotrx) 	call advect_scalar(trx,trxadv,trxwle,trx2leadv&
          &,trx2legrad,trxwleadv,.true.) 
     if(dotry) 	call advect_scalar(try,tryadv,trywle,try2leadv&
          &,try2legrad,trywleadv,.true.) 
     if(dotrzz) 	call advect_scalar(trzz,trzzadv,trzzwle,trzz2leadv&
          &,trzz2legrad,trzzwleadv,.true.) 
     if(dotro3)  call advect_scalar(tro3,tro3adv,tro3wle,tro32leadv&
          &,tro32legrad,tro3wleadv,.true.) 
!==============================================================

!---------------------------------------------------------
!      diffusion of scalars :

!        Update boundaries for scalars:

      if(dosgs) call boundaries(3)


      if(do_tkh_z) then !Janniy yani 
       !call random_forest_tkh()
       call nn_diffusion()
      end if




      mseflag=1 ! peters
      call diffuse_scalar_vert_flag_surf_flag(t,fluxbt-uwpbl_fluxbt,fluxtt,tdiff,twsb, &
                           t2lediff,t2lediss,twlediff,.true.) !Yani change the functio call 

    !call diffuse_scalar(t,fluxbt-uwpbl_fluxbt,fluxtt,tdiff,twsb, &
    !                       t2lediff,t2lediss,twlediff,.true.) !Yani change the functio call 
      do k=1,nzm
         total_water_evap = total_water_evap - &
                         sum(q(1:nx,1:ny,k))*adz(k)*dz *rho(k)
      end do

      mseflag=2 ! peters
      if(.not.dokruegermicro) then !.and.((.not.do_rf_diffusion).or.(.not.dorandomforest))) then !bloss  Yani-I think the logic should be corrent - to consider to print something

        call diffuse_scalar_vert_flag_surf_flag(q,fluxbq-uwpbl_fluxbq,fluxtq,qdiff,qwsb, &
                         q2lediff,q2lediss,qwlediff,.true.) !Yani change the functio call 

      if(do_z_diffusion_rf) then !Janniy yani
        call random_forest_diffusion()
      end if 


       !call diffuse_scalar(q,fluxbq-uwpbl_fluxbq,fluxtq,qdiff,qwsb, &
       !                   q2lediff,q2lediss,qwlediff,.true.) !Yani change the functio call 
  !    elseif((do_rf_diffusion).and.(.not.do_rf_q_surf_flux))  !Yani 
  !     call diffuse_scalar_only_surf(q,fluxbq-uwpbl_fluxbq,fluxtq,qdiff,qwsb, &
  !                        q2lediff,q2lediss,qwlediff,.true.)
  !    elseif((do_rf_diffusion).and.(do_rf_q_surf_flux))  !Yani 
  !     call diffuse_scalar_no_vert_diff(q,fluxbq-uwpbl_fluxbq,fluxtq,qdiff,qwsb, &
  !                        q2lediff,q2lediss,qwlediff,.true.)

      !==============================================================
      elseif(dokruegermicro) then !bloss
       call diffuse_scalar(qv,fluxbq,fluxtq,qdiff,qwsb, &
                           q2lediff,q2lediss,qwlediff,.true.)

       ! update q for computation of evaporation part of total water budget
       q(1:nx,1:ny,1:nzm) = qv(1:nx,1:ny,1:nzm) + qc(1:nx,1:ny,1:nzm) &
            + qi(1:nx,1:ny,1:nzm)
      !==============================================================
      end if

      do k=1,nzm
         total_water_evap = total_water_evap + &
                         sum(q(1:nx,1:ny,k))*adz(k)*dz *rho(k)
      end do

      mseflag=0 ! peters
      if(.not.dosmagor) then
          call diffuse_scalar(tke,fzero,fzero,dummy,tkewsb, &
                                    dummy,dummy,dummy,.false.)
      else if(doscalar) then
          call diffuse_scalar(tke,fluxbq,fluxtq,dummy,tkewsb, &
                           s2lediff,s2lediss,swlediff,.true.)
      end if

      mseflag=2 ! peters
      if(docloud.and.doprecip) then
       if((.not.dorandomforest).or.(dorandomforest.and.rf_uses_qp)) then
        if(.not.dokruegermicro) then !bloss
           !call diffuse_scalar(qp,fzero,fzero,qpdiff,qpwsb, &   
           !                dummy,dummy,dummy,.false.)    !JY - changed to my diffusion scheme that can use the RF results
	   call diffuse_scalar_vert_flag_surf_flag(qp,fzero,fzero,qpdiff,qpwsb,dummy,dummy,dummy,.false.) 
        !==============================================================
        elseif(dokruegermicro) then !bloss
           call diffuse_scalar(qc,fzero,fzero,qcdiff,qcwsb, &
                dummy,dummy,dummy,.false.)
           call diffuse_scalar(qi,fzero,fzero,qidiff,qiwsb, &
                dummy,dummy,dummy,.false.)
           if(dostatis) qdiff = qdiff + qcdiff + qidiff
           if(dostatis) qwsb = qwsb + qcwsb + qiwsb
           call diffuse_scalar(qr,fzero,fzero,qrdiff,qrwsb, &
                dummy,dummy,dummy,.false.)
           call diffuse_scalar(qg,fzero,fzero,qgdiff,qgwsb, &
                dummy,dummy,dummy,.false.)
           call diffuse_scalar(qs,fzero,fzero,qsdiff,qswsb, &
                dummy,dummy,dummy,.false.)
        !==============================================================
        endif
       endif
      endif

!==============================================================
!kzm adding some tracers

     mseflag=0 ! peters
     if(dotrz) call diffuse_scalar(trz,fluxbtrz,fzero,trzdiff,trzwsb, &
          trz2lediff,trz2lediss,trzwlediff,.true.)
     if(dotrx) call diffuse_scalar(trx,fluxbtrx,fzero,trxdiff,trxwsb, &
          trx2lediff,trx2lediss,trxwlediff,.true.)
     if(dotry) call diffuse_scalar(try,fluxbtry,fzero,trydiff,trywsb, &
          try2lediff,try2lediss,trywlediff,.true.)
     if(dotrzz) call diffuse_scalar(trzz,fluxbtrzz,fzero,trzzdiff,trzzwsb, &
          trzz2lediff,trzz2lediss,trzzwlediff,.true.)
     if(dotro3) call diffuse_scalar(tro3,fluxbtro3,fzero,tro3diff,tro3wsb, &
          tro32lediff,tro32lediss,tro3wlediff,.true.)
!==============================================================

!-----------------------------------------------------------
!       Scaling microphysics rates based upon whether the motion is convective, 
!       stratiform or somewhere in between
!       JA 8/8/06

        if(domicroscaling) then
                call testforstratiform()
        end if


!-----------------------------------------------------------
!       Cloud condensation/evaporation and precipitation processes:


      if(docloud) then
        call cloud()
        if(doprecip) then
         if((.not.dorandomforest).or.(dorandomforest.and.rf_uses_qp)) then
          call precip_proc() !-TRIED TO DO ONLY CORRECTION... JY Yani at the moment I substitute it.... 
         end if
         if(dorandomforest) then 
         ! call random_forest()
          call nn_convection_flux()
         end if
     !    if(do_z_diffusion_rf) then !Janniy yani
     !     call random_forest_diffusion()
     !    end if 
        end if
      end if

      !bloss: initialize q, qn, qp if using krueger microphysics
      if (dokruegermicro) then
         q(1:nx,1:ny,1:nzm) = qv(1:nx,1:ny,1:nzm) + qc(1:nx,1:ny,1:nzm) &
              + qi(1:nx,1:ny,1:nzm)
         qn(1:nx,1:ny,1:nzm) = qc(1:nx,1:ny,1:nzm) + qi(1:nx,1:ny,1:nzm)
         qp(1:nx,1:ny,1:nzm) = qr(1:nx,1:ny,1:nzm) + qs(1:nx,1:ny,1:nzm) &
              + qg(1:nx,1:ny,1:nzm)
      end if

!-----------------------------------------------------------
!	Radiation

!      if(icycle.eq.ncycle.and.(dolongwave.or.doshortwave)) then 
      if((dolongwave.or.doshortwave)) then 
	call radiation()     
      end if

!-----------------------------------------------------------
!       Compute field diagnostics and update the velocity field:

      call diagnose()
!------compute parameterized linear large-scale wave response
      if(doparameterizedwave) call linear_wave()

!----------------------------------------------------------
! Rotate the dynamic tendency arrays for Adams-bashforth scheme:

      nn=na
      na=nc
      nc=nb
      nb=nn

      firststep = .false. ! no longer first step of run

   end do ! icycle	
          
!----------------------------------------------------------

! don't want to send both diffusivities through statistics

!!$   if(douwpbl.and.(uwpbl_ndiv.gt.1)) then
!!$      tk_z(1:nx, 1:ny, :) = tk_z(1:nx, 1:ny, 1:nzm) - tk_z_uw(:, :, :)*sqrt(betafactor)
!!$      tkh_z(1:nx, 1:ny, :) = tkh_z(1:nx, 1:ny, 1:nzm) - tkh_z_uw(:, :, :)*sqrt(betafactor)
!!$   endif
   if(mod(nstep,nstatis).eq.0) call statistics()
	
!----------------------------------------------------------
!  collect statistics, write save-file, etc.

   call stepout(nstatsteps)
  
! do want to send both diffusivities through kurant

!!$   if(douwpbl.and.(uwpbl_ndiv.gt.1)) then
!!$      tk_z(1:nx, 1:ny, :) = tk_z(1:nx, 1:ny, 1:nzm) + tk_z_uw(:, :, :)*sqrt(betafactor)
!!$      tkh_z(1:nx, 1:ny, :) = tkh_z(1:nx, 1:ny, 1:nzm) + tkh_z_uw(:, :, :)*sqrt(betafactor)
!!$   endif

!----------------------------------------------------------

   ! peters only end when MSE and stat is output
   !               or when stat is output and MSE is off (i.e. nsaveMSE==-1)
   if ((mod(nstep,nstat).eq.0).and. &
        ((nsaveMSE.eq.-1).or.(mod(nstep,nsaveMSE).eq.0))) then

      call secondsf(cputime)
      if (masterproc) write(*,999) cputime-oldtime, float(nstat)*dt/3600.
999   format('CPU TIME = ',f12.4, ' OVER ', f8.2,' MODEL HOURS')

      elapsed_time = cputime-init_time
      if ((elapsed_time+2.*(cputime-oldtime)).gt.60.*float(nelapsemin)) then
         nelapse=0 ! Job will stop when nelapse=0
         if (masterproc) write(*,*) 'Job stopping -- nelapsemin reached'
      end if

      oldtime = cputime
   end if

end do ! main loop

!----------------------------------------------------------
!----------------------------------------------------------

! ONLY WRITE RESTART FILES AT END OF LAST TIMESTEP
!    Useful to minimize disk usage during shorter runs over NFS.
!    Only use this option with nelapsemin, so that only short amounts 
!    of runtime will be lost if a run does crash. 	
if (dorestart_last) then
   call write_all() ! save restart file
   call write_rad() ! write radiation restart file
end if

call task_abort()

end program crm
