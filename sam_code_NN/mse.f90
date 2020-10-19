
! peters module mse in mse_budget.f90

module mse

implicit none

private ! make module variables/subroutines private by default
public :: divergence, diagnoseMSE, writeMSE, &             ! subroutines 
     columnint, allocateMSE, zeroMSE, deallocateMSE, &     ! subroutines
     checkMSEinput,           &                            ! subroutines
     oldsli_mse, oldhf_mse, div, prec_inst_mse, prec_inst_frozen_mse, & !output
     slidamp_mse, hdamp_mse, hdifh_mse, hdifs_mse, &       ! output
     hadvs_mse, vadvs_mse, hadvh_mse, vadvh_mse, &         ! output
     huwpbl_mse, suwpbl_mse,     &                          ! output
     dtu_uu, dtu_uv, dtu_uw, dtv_vu, dtv_vv, &             ! output
     dtv_vw, dtu_diff, dtv_diff,&                          ! output
     u_uwpbl, v_uwpbl,  &                                   ! output
     nsaveMSE, navgMSE, nsaveMSEstart, nsaveMSEend, nsaveMSE3D, & ! input 
     domombudget, writeMSE3D, &                                   ! input
     mseflag, doMSE, isallocateMSE, docollectMSE   ! internal switches

!  variables to output
real, dimension(:,:,:), allocatable :: u_mse, v_mse, w_mse   ! m/s
real, dimension(:,:,:), allocatable :: tabs_mse, qrad_mse    ! K, K/s
real, dimension(:,:,:), allocatable :: rh_mse                 ! [%]
real, dimension(:,:,:), allocatable :: qv_mse, qc_mse, qi_mse ! kg/kg
real, dimension(:,:,:), allocatable :: qr_mse, qsg_mse        ! kg/kg
real, dimension(:,:,:), allocatable :: pp_mse                 ! [Pa]

!  2D fields to output for budgets
!  All 2D output fields have units of W/m^2 except for sst=K
real, dimension(:,:), allocatable :: prec_mse, prec_frozen_mse
real, dimension(:,:), allocatable :: shf_mse, lhf_mse, sst_mse
real, dimension(:,:), allocatable :: hadvh_mse, vadvh_mse, hstor_mse
real, dimension(:,:), allocatable :: hadvs_mse, vadvs_mse, sstor_mse
real, dimension(:,:), allocatable :: hdifh_mse, hdifs_mse
real, dimension(:,:), allocatable :: rad_mse
real, dimension(:,:), allocatable :: lwns_mse, lwnt_mse, lwnsc_mse, lwntc_mse
real, dimension(:,:), allocatable :: swns_mse, swntm_mse, swnsc_mse, swntc_mse
real, dimension(:,:), allocatable :: slidamp_mse, hdamp_mse
real, dimension(:,:), allocatable :: huwpbl_mse, suwpbl_mse

! 1D field to output
real, dimension(:), allocatable   :: pres_mse

! momentum budget variables.  these are variables to output
real, dimension(:,:,:), allocatable :: u_uu, u_uv, u_uw     ! m/s^2
real, dimension(:,:,:), allocatable :: v_vu, v_vv, v_vw     ! m/s^2
real, dimension(:,:,:), allocatable :: u_diff, v_diff       ! m/s^2
real, dimension(:,:,:), allocatable :: ustor, vstor         ! m/s^2
real, dimension(:,:,:), allocatable :: u_uwpbl, v_uwpbl     ! m/s^2

! 3D variables accumulated over previous nsaveMSE3D time steps.  These are
!  output when the 3D fields are output
real, dimension(:,:,:), allocatable :: u_mse_out, v_mse_out, w_mse_out 
real, dimension(:,:,:), allocatable :: tabs_mse_out, qrad_mse_out, rh_mse_out
real, dimension(:,:,:), allocatable :: qv_mse_out, qc_mse_out, qi_mse_out
real, dimension(:,:,:), allocatable :: qr_mse_out, qsg_mse_out
real, dimension(:,:,:), allocatable :: pp_mse_out
real, dimension(:,:,:), allocatable :: u_uu_out, u_uv_out, u_uw_out
real, dimension(:,:,:), allocatable :: v_vu_out, v_vv_out, v_vw_out
real, dimension(:,:,:), allocatable :: u_diff_out, v_diff_out
real, dimension(:,:,:), allocatable :: ustor_out, vstor_out
real, dimension(:,:,:), allocatable :: u_uwpbl_out, v_uwpbl_out

! momentum budget variables;  these hold the previous steps for SAM's
!  adams bashforth scheme
real, dimension(:,:,:,:), allocatable :: dtu_uu, dtu_uv, dtu_uw     ! m/s^2
real, dimension(:,:,:,:), allocatable :: dtv_vu, dtv_vv, dtv_vw     ! m/s^2
real, dimension(:,:,:,:), allocatable :: dtu_diff, dtv_diff       ! m/s^2

! variables used by SAM
real, dimension(:,:,:), allocatable :: oldsli_mse, oldhf_mse   ! J/kg
real, dimension(:,:,:), allocatable :: oldu, oldv              ! m/s
real, dimension(:,:,:), allocatable :: div      ! divergence
real, dimension(:,:), allocatable   :: prec_inst_mse
real, dimension(:,:), allocatable   :: prec_inst_frozen_mse
integer nsaveMSE, navgMSE, nsaveMSEstart, nsaveMSEend, nsaveMSE3D
integer mseflag
logical doMSE, isallocateMSE
logical domombudget, writeMSE3D

!--------------------------------------------------------------------
!  end variable definition.  now start subroutine definition
CONTAINS

subroutine divergence()
 ! computes the divergence from the velocity field
 !  This is called after the velocity field is updated to the next half
 !  step, and u is replaced with u*rho(k)*dtn/dx, and v is replaced with
 !  v*rho(k)*dtn/dy
 !  Don't worry about the w term since d/dz(w*scalar) integrates to zero
 !
 ! on exit, div = du/dx + dv/dy
 
 use vars
 implicit none

 integer i, j, k, ic, jc
 real coef

 ! computes the divergence from the velocity field
  
 if(RUN3D) then
   do k=1,nzm
     coef = 1./dtn/rho(k)
     do j=1,ny
       jc=j+1 
       do i=1,nx
         ic=i+1
         div(i,j,k) = (u(ic,j,k)-u(i,j,k) + v(i,jc,k)-v(i,j,k))*coef
       end do
     end do
   end do

 else    ! 2D run
   j=1
   do k=1,nzm
     coef = dtn/rho(k)
     do i=1,nx
       ic=i+1
       div(i,j,k) = (u(ic,j,k)-u(i,j,k))*coef
     end do
   end do
 end if

end subroutine divergence

!-------------------------------------------------------------------
subroutine diagnoseMSE()
  ! updates the fields to output, and writes them if necessary

  use vars
  use params
  use rad
  implicit none

  integer i, j, k
  real coef, coef2, omn, omp, tmp
  logical local_writeMSE3D
  real relhum

  coef = dtfactor/float(navgMSE)
  coef2 = 1./float(navgMSE)

  do k=1,nzm
    ! dtfactor = dtn/dt = fraction of time step
    u_mse(:,:,k) = u_mse(:,:,k) + u(1:nx,1:ny,k)*coef
    v_mse(:,:,k) = v_mse(:,:,k) + v(1:nx,1:ny,k)*coef
    w_mse(:,:,k) = w_mse(:,:,k) + w(1:nx,1:ny,k)*coef

    tabs_mse(:,:,k) = tabs_mse(:,:,k) + tabs(1:nx,1:ny,k)*coef
    qrad_mse(:,:,k) = qrad_mse(:,:,k) + qrad(1:nx,1:ny,k)*coef
    pp_mse(:,:,k)   = pp_mse(:,:,k)   + p(1:nx,1:ny,k)*coef

    if(.not.dokruegermicro) then
       do j=1,ny !bloss: swapped loop order
          do i=1,nx
             omn = omegan(tabs(i,j,k))
             omp = omegap(tabs(i,j,k))
             qv_mse(i,j,k) = qv_mse(i,j,k) + (q(i,j,k)-qn(i,j,k))*coef
             qc_mse(i,j,k) = qc_mse(i,j,k) + omn*qn(i,j,k)*coef
             qi_mse(i,j,k) = qi_mse(i,j,k) + (1.-omn)*qn(i,j,k)*coef
             qr_mse(i,j,k) = qr_mse(i,j,k) + omp*qp(i,j,k)*coef
             qsg_mse(i,j,k) = qsg_mse(i,j,k) + (1.-omp)*qp(i,j,k)*coef

             ! compute rh
             relhum = omn*qsatw(tabs(i,j,k),pres(k))+ &
                          (1.-omn)*qsati(tabs(i,j,k),pres(k))
             relhum = (q(i,j,k)-qn(i,j,k))/relhum
             relhum = max(0.,min(1.,relhum))
             rh_mse(i,j,k) = rh_mse(i,j,k) + relhum*100.*coef
          end do
       end do
    else ! dokruegermicro = .true.
       qv_mse(:,:,k) = qv_mse(:,:,k) + qv(1:nx,1:ny,k)*coef
       qc_mse(:,:,k) = qc_mse(:,:,k) + qc(1:nx,1:ny,k)*coef
       qi_mse(:,:,k) = qi_mse(:,:,k) + qi(1:nx,1:ny,k)*coef
       qr_mse(:,:,k) = qr_mse(:,:,k) + qr(1:nx,1:ny,k)*coef
       qsg_mse(:,:,k) = qsg_mse(:,:,k) + qs(1:nx,1:ny,k)*coef &
            + qg(1:nx,1:ny,k)*coef

       ! compute rh
       do j=1,ny
          do i=1,nx
            omn = omegan(tabs(i,j,k))
            relhum = omn*qsatw(tabs(i,j,k),pres(k))+ &
                       (1.-omn)*qsati(tabs(i,j,k),pres(k))
            relhum = qv(i,j,k)/relhum
            relhum = max(0.,min(1.,relhum))
            rh_mse(i,j,k) = rh_mse(i,j,k) + relhum*100.*coef
          end do
       end do
    end if
  end do   ! end k=1,nzm

  ! don't scale prec by dtfactor since this is implicitly done in
  ! precip_fall.f90
  prec_mse = prec_mse + (dz/dt)*prec_inst_mse*coef2   ! W/m^2
  prec_frozen_mse = prec_frozen_mse + (dz/dt)*prec_inst_frozen_mse*coef2 !W/m^2
  shf_mse = shf_mse   + cp*rhow(1)*fluxbt*coef              ! W/m^2
  lhf_mse = lhf_mse   + lcond*rhow(1)*fluxbq*coef           ! W/m^2
  sst_mse = sst_mse   + sstxy(1:nx,1:ny)*coef                          ! K

  pres_mse = pres_mse + pres*coef

  ! accumulate radiation fluxes.  only need to do at end of time step,
  !  since this is the only time they change.  Don't need to worry about
  !  dtfactor either because of this
  if(icycle.eq.ncycle) then
    lwns_mse = lwns_mse + lwnsxy*coef2
    lwnt_mse = lwnt_mse + lwntxy*coef2
    lwnsc_mse = lwnsc_mse + lwnscxy*coef2
    lwntc_mse = lwntc_mse + lwntcxy*coef2
    swns_mse = swns_mse + swnsxy*coef2
    swntm_mse = swntm_mse + swntmxy*coef2
    swnsc_mse = swnsc_mse + swnscxy*coef2
    swntc_mse = swntc_mse + swntcxy*coef2
    do j=1,ny
      do i=1,nx
        call columnint(qrad(i,j,:),tmp)
        rad_mse(i,j) = rad_mse(i,j) + cp*tmp*coef2         ! W/m^2
      end do
    end do
  end if

  ! now compute the mom budget terms
  do k=1,nzm
    do j=1,ny
      do i=1,nx

        u_uu(i,j,k) = u_uu(i,j,k) + coef &
              *(at*dtu_uu(i,j,k,na)+bt*dtu_uu(i,j,k,nb)+ct*dtu_uu(i,j,k,nc))
        u_uv(i,j,k) = u_uv(i,j,k) + coef &
              *(at*dtu_uv(i,j,k,na)+bt*dtu_uv(i,j,k,nb)+ct*dtu_uv(i,j,k,nc))
        u_uw(i,j,k) = u_uw(i,j,k) + coef &
              *(at*dtu_uw(i,j,k,na)+bt*dtu_uw(i,j,k,nb)+ct*dtu_uw(i,j,k,nc))

        v_vu(i,j,k) = v_vu(i,j,k) + coef &
              *(at*dtv_vu(i,j,k,na)+bt*dtv_vu(i,j,k,nb)+ct*dtv_vu(i,j,k,nc))
        v_vv(i,j,k) = v_vv(i,j,k) + coef &
              *(at*dtv_vv(i,j,k,na)+bt*dtv_vv(i,j,k,nb)+ct*dtv_vv(i,j,k,nc))
        v_vw(i,j,k) = v_vw(i,j,k) + coef &
              *(at*dtv_vw(i,j,k,na)+bt*dtv_vw(i,j,k,nb)+ct*dtv_vw(i,j,k,nc))

        u_diff(i,j,k) = u_diff(i,j,k) + coef  &
           *(at*dtu_diff(i,j,k,na)+bt*dtu_diff(i,j,k,nb)+ct*dtu_diff(i,j,k,nc))
        v_diff(i,j,k) = v_diff(i,j,k) + coef &
           *(at*dtv_diff(i,j,k,na)+bt*dtv_diff(i,j,k,nb)+ct*dtv_diff(i,j,k,nc))

      end do
    end do
  end do


  !--------------------------------------------------
  ! all the fluxes are now accumulated.  if time to write the data, update
  ! the adv, dif fluxes, the storage terms, and write the fields
  if((icycle.eq.ncycle) .and. (mod(nstep,nsaveMSE).eq.0) ) then

    ! need to update the damping, adv and dif fluxes.  now, h is only
    ! q+qp, so need to add in sli terms
    hadvh_mse = hadvh_mse + hadvs_mse
    vadvh_mse = vadvh_mse + vadvs_mse
    hdifh_mse = hdifh_mse + hdifs_mse
    hdamp_mse = hdamp_mse + slidamp_mse
    huwpbl_mse = huwpbl_mse + suwpbl_mse

    ! update the storage term
    do i=1,nx
      do j=1,ny
        call columnint(cp*t(i,j,:)-oldsli_mse(i,j,:),sstor_mse(i,j))
        call columnint(cp*t(i,j,:)+lcond*(q(i,j,:)+qp(i,j,:))-       &
                         oldhf_mse(i,j,:), hstor_mse(i,j))
      end do
    end do
    sstor_mse = sstor_mse/dt/float(navgMSE)
    hstor_mse = hstor_mse/dt/float(navgMSE)

    ! update u,v storage
    ustor = (u(1:nx,1:ny,1:nzm) - oldu)/dt/float(navgMSE)
    vstor = (v(1:nx,1:ny,1:nzm) - oldv)/dt/float(navgMSE)

    ! accumulate the 3D out variables
    if(writeMSE3D) then
      call copyMSE3D_to_MSE3Dout(nsaveMSE)
    elseif(nsaveMSE3D.ne.1) then
      call copyMSE3D_to_MSE3Dout(nsaveMSE3D)
    end if

    ! write all the data to disk
    call writeMSE(local_writeMSE3D)

    ! reset all fields to zero.
    call zeroMSE()
    if(local_writeMSE3D) call zeroMSE_out()

  end if  ! if (icycle.eq.ncycle.and.mod(nstep,nsaveMSE).eq.0), ie, time2write


end subroutine diagnoseMSE


!-------------------------------------------------------------------
subroutine docollectMSE()
  ! decide if it is time to collect MSE statistics
  use vars
  use params
  implicit none

  if( (nstep.ge.(nsaveMSEstart-navgMSE+1)).and.(nstep.le.nsaveMSEend).and. &
    ((mod(nstep,nsaveMSE).ge.(nsaveMSE-navgMSE+1)).or. &
                                  (mod(nstep,nsaveMSE).eq.0)) ) then
     doMSE=.true.
  else
     doMSE=.false.
  end if

  ! save oldsli, oldhf if navgMSE steps before time to write
  if(isallocateMSE .and. mod(nstep,nsaveMSE).eq.(nsaveMSE-navgMSE+1) .and. &
                                  icycle.eq.1) then
    oldsli_mse = cp*t(1:nx,1:ny,1:nzm)
    oldhf_mse = oldsli_mse+lcond*(q(1:nx,1:ny,1:nzm)+qp(1:nx,1:ny,1:nzm))
    oldu = u(1:nx,1:ny,1:nzm)
    oldv = v(1:nx,1:ny,1:nzm)
    if(masterproc) print*,'Starting to collect MSE fluxes...'
  end if

end subroutine docollectMSE

!-------------------------------------------------------------------
subroutine writeMSE(local_writeMSE3D)

  use vars
  implicit none

  ! output
  logical local_writeMSE3D       ! write 3D file this time step??

  character *80 filename,long_name
  character *8 name
  character *10 timechar
  character *4 rankchar
  character *6 filetype
  character *10 units
  character *10 c_z(nzm),c_p(nzm),c_dx, c_dy, c_time
  integer i,j,k,nfields1,nfields2D,nfields3D
  real tmp(nx,ny,nzm)
  logical MSE3Dbin

  integer k200, k500, k850

  !--------------------
  MSE3Dbin = .true.

  ! decide if time to write 3D file
  if(writeMSE3D) then
    local_writeMSE3D = .true.
  elseif(nsaveMSE3D.ne.1 .and. mod(nstep,nsaveMSE3D).eq.0) then
    local_writeMSE3D = .true.
  else
    local_writeMSE3D = .false.
  end if

  nfields3D=24  ! yikes!
  nfields2D=48  ! yikes!
  if(.not.domombudget) nfields3D = 11
  if(.not.local_writeMSE3D) nfields3D = 0

  nfields1=0

  if(masterproc) then

    write(rankchar,'(i4)') nsubdomains
    write(timechar,'(i10)') nstep
    do k=1,11-lenstr(timechar)-1
      timechar(k:k)='0'
    end do

    if(RUN3D) then
      if(local_writeMSE3D) then
        filetype = '.MSE3D'
      else
        filetype = '.MSE2D'
      end if
      filename='./DATA3D/'//trim(case)//'_MSE_'//trim(caseid)//'_'// &
          rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype
      open(46,file=filename,status='unknown',form='unformatted')
    else
      filetype = '.bin2D'
      filename='./DATA3D/'//trim(case)//'_MSE_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//filetype
      if(nrestart.eq.0.and.notopened3D) then
         open(46,file=filename,status='unknown',form='unformatted')	
      else
         open(46,file=filename,status='unknown', &
                            form='unformatted', position='append')
      end if
      notopened3D=.false.
    end if  ! end if(run3D)

    ! write fields in 3Dbin format
    write(46) nx,ny,nzm,nsubdomains,nsubdomains_x,nsubdomains_y
    write(46) nfields3D,nfields2D
    do k=1,nzm
      write(46) z(k) 
    end do
    do k=1,nzm
      write(46) pres_mse(k)
    end do
    do k=1,nzm
      write(46) rho(k)
    end do
    do k=1,nz
      write(46) rhow(k)
    end do
    write(46) dx
    write(46) dy
    write(46) nstep*dt/(3600.*24.)+day0

  end if         ! if(masterproc)

  ! write 3D fields if necessary
 if(local_writeMSE3D) then

  nfields1=nfields1+1
   tmp=u_mse_out
  name='U'
  long_name='X Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=v_mse_out
  name='V'
  long_name='Y Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=w_mse_out
  name='W'
  long_name='Z Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=tabs_mse_out
  name='TABS'
  long_name='Absolute Temperature'
  units='K'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=rh_mse_out
  name='RH'
  long_name='Relative Humidity'
  units='%'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)


  nfields1=nfields1+1
   tmp=qrad_mse_out
  name='QRAD'
  long_name='Radiative heating rate'
  units='K/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=qv_mse_out
  name='QV'
  long_name='Water Vapor'
  units='kg/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=qc_mse_out
  name='QC'
  long_name='Cloud water'
  units='kg/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=qi_mse_out
  name='QI'
  long_name='Cloud Ice'
  units='kg/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=qr_mse_out
  name='QR'
  long_name='Rain'
  units='kg/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=qsg_mse_out
  name='QSG'
  long_name='Snow + Grapuel'
  units='kg/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=pp_mse_out
  name='PP'
  long_name='Pressure Perturbation'
  units='Pa'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  if(domombudget) then           ! only write mombudget if necessary

  nfields1=nfields1+1
   tmp=u_uu_out
  name='U_UU'
  long_name='Advective term -d(u*u)/dx in u-momentum equation'
  units='m/s^2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=u_uv_out
  name='U_UV'
  long_name='Advective term -d(u*v)/dy in u-momentum equation'
  units='m/s^2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=u_uw_out
  name='U_UW'
  long_name='Advective term -d(u*w)/dz in u-momentum equation'
  units='m/s^2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=v_vu_out
  name='V_VU'
  long_name='Advective term -d(v*u)/dx in v-momentum equation'
  units='m/s^2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=v_vv_out
  name='V_VV'
  long_name='Advective term -d(v*v)/dy in v-momentum equation'
  units='m/s^2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=v_vw_out
  name='V_VW'
  long_name='Advective term -d(v*w)/dz in v-momentum equation'
  units='m/s^2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=u_diff_out
  name='U_DIFF'
  long_name='Diffusive flux in u-momentum equation'
  units='m/s^2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=v_diff_out
  name='V_DIFF'
  long_name='Diffusive flux in v-momentum equation'
  units='m/s^2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=u_uwpbl_out
  name='U_UWPBL'
  long_name='Tendency from UWPBL in u-momentum equation'
  units='m/s^2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=v_uwpbl_out
  name='V_UWPBL'
  long_name='Tendency from UWPBL in v-momentum equation'
  units='m/s^2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=ustor_out
  name='USTOR'
  long_name='Storage of u-momentum du/dt'
  units='m/s^2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp=vstor_out
  name='VSTOR'
  long_name='Storage of v-momentum dv/dt'
  units='m/s^2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 MSE3Dbin,dompi,rank,nsubdomains)

  end if   ! if(domombudget)
 
  end if   ! if(local_writeMSE3D)


  !--------------------
  ! surface fluxes
  nfields1=nfields1+1
   tmp(:,:,1)=prec_mse
  name='PREC'
  long_name='Surface Precip. Rate'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=prec_frozen_mse
  name='PREC_FZN'
  long_name='Frozen Surface Precip. Rate'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=shf_mse
  name='SHF'
  long_name='Sensible Heat Flux'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=lhf_mse
  name='LHF'
  long_name='Latent Heat Flux'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=sst_mse
  name='SST'
  long_name='Sea Surface Temperature'
  units='K'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  !--------------------
  ! MSE budget stuff
  nfields1=nfields1+1
   tmp(:,:,1)=hadvh_mse
  name='HADVH'
  long_name='Column int. FMSE Horizontal Advection Source'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=vadvh_mse
  name='VADVH'
  long_name='Column int. FMSE Vertical Advection Source'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=hdifh_mse
  name='HDIFH'
  long_name='Column int. FMSE Horizontal Diffusive Source'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=hstor_mse
  name='HSTOR'
  long_name='Column int. FMSE Storage'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=hadvs_mse
  name='HADVS'
  long_name='Column int. S-LI Horizontal Advection Source'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)
                                                                                
  nfields1=nfields1+1
   tmp(:,:,1)=vadvs_mse
  name='VADVS'
  long_name='Column int. S-LI Vertical Advection Source'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)
                                                                                
  nfields1=nfields1+1
   tmp(:,:,1)=hdifs_mse
  name='HDIFS'
  long_name='Column int. S-LI Horizontal Diffusive Source'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)
                                                                                
  nfields1=nfields1+1
   tmp(:,:,1)=sstor_mse
  name='SSTOR'
  long_name='Column int. S-LI Storage'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  !------------------
  ! radiative fluxes
  nfields1=nfields1+1
   tmp(:,:,1)=rad_mse
  name='RAD'
  long_name='Net Column integrated Radiative Flux Divergence'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=lwns_mse
  name='LWNS'
  long_name='Net Surface Upward LW Flux'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=lwnt_mse
  name='LWNT'
  long_name='Net TOA Upward LW Flux (OLR)'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=lwnsc_mse
  name='LWNSC'
  long_name='Clear Sky Surface Upward LW Flux'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=lwntc_mse
  name='LWNTC'
  long_name='Clear Sky TOA Upward LW Flux'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=swns_mse
  name='SWNS'
  long_name='Net Surface Downward SW Flux'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)
                                                                                
  nfields1=nfields1+1
   tmp(:,:,1)=swntm_mse
  name='SWNTM'
  long_name='Net TOM Downward SW Flux'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)
                                                                                
  nfields1=nfields1+1
   tmp(:,:,1)=swnsc_mse
  name='SWNSC'
  long_name='Clear Sky Surface Downward SW Flux'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)
                                                                                
  nfields1=nfields1+1
   tmp(:,:,1)=swntc_mse
  name='SWNTC'
  long_name='Clear Sky TOA Downward SW Flux'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=slidamp_mse
  name='SLI_DAMP'
  long_name='Column int. S-LI Damping source'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=hdamp_mse
  name='HF_DAMP'
  long_name='Column int. FMSE Damping source'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=huwpbl_mse
  name='H_UWPBL'
  long_name='Column int. FMSE source from UWPBL scheme'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=suwpbl_mse
  name='S_UWPBL'
  long_name='Column int. S-LI source from UWPBL scheme'
  units='W/m^2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  !----------------------------------------------------------------------
  ! write some derived 2D fields
  ! find indices of 850, 500, 200 mb 
  call findindex(pres_mse,850.,k850)
  call findindex(pres_mse,500.,k500)
  call findindex(pres_mse,200.,k200)

  nfields1=nfields1+1
   tmp(:,:,1)=u_mse(:,:,1)
  name='USURF'
  long_name='U at surface'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=v_mse(:,:,1)
  name='VSURF'
  long_name='V at surface'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=u_mse(:,:,k850)
  name='U850'
  long_name='U at 850 hPa'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=v_mse(:,:,k850)
  name='V850'
  long_name='V at 850 hPa'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=u_mse(:,:,k500)
  name='U500'
  long_name='U at 500 hPa'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=v_mse(:,:,k500)
  name='V500'
  long_name='V at 500 hPa'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=u_mse(:,:,k200)
  name='U200'
  long_name='U at 200 hPa'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=v_mse(:,:,k200)
  name='V200'
  long_name='V at 200 hPa'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=w_mse(:,:,k850)
  name='W850'
  long_name='W at 850 hPa'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=w_mse(:,:,k500)
  name='W500'
  long_name='W at 500 hPa'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=w_mse(:,:,k200)
  name='W200'
  long_name='W at 200 hPa'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)


  nfields1=nfields1+1
    tmp(:,:,1)=tabs_mse(:,:,1)
  name='TSURF'
  long_name='TABS at surface'
  units='K'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=tabs_mse(:,:,k850)
  name='T850'
  long_name='T at 850 hPa'
  units='K'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=tabs_mse(:,:,k500)
  name='T500'
  long_name='TABS at 500 hPa'
  units='K'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=tabs_mse(:,:,k200)
  name='T200'
  long_name='TABS at 200 hPa'
  units='K'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)


  nfields1=nfields1+1
    tmp(:,:,1)=qv_mse(:,:,1) + &
               qc_mse(:,:,1) + &
               qi_mse(:,:,1)
  name='QSURF'
  long_name='Q at surface'
  units='kg/kg'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=qv_mse(:,:,k850)+qc_mse(:,:,k850)+qi_mse(:,:,k850)
  name='Q850'
  long_name='Q at 850 hPa'
  units='kg/kg'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=qv_mse(:,:,k500)+qc_mse(:,:,k500)+qi_mse(:,:,k500)
  name='Q500'
  long_name='Q at 500 hPa'
  units='kg/kg'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=qv_mse(:,:,k200)+qc_mse(:,:,k200)+qi_mse(:,:,k200)
  name='Q200'
  long_name='Q at 200 hPa'
  units='kg/kg'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=qc_mse(:,:,k850)+qi_mse(:,:,k850)
  name='QN850'
  long_name='QN at 850 hPa'
  units='kg/kg'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=qc_mse(:,:,k500)+qi_mse(:,:,k500)
  name='QN500'
  long_name='QN at 500 hPa'
  units='kg/kg'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   tmp(:,:,1)=qc_mse(:,:,k200)+qi_mse(:,:,k200)
  name='QN200'
  long_name='QN at 200 hPa'
  units='kg/kg'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               MSE3Dbin,dompi,rank,nsubdomains)
  !----------------------------------------------------------------------




  call task_barrier()

  if((nfields3D+nfields2D).ne.nfields1) then
    if(masterproc) print*,'write_fields3D error: nfields'
    call task_abort()
  end if
  if(masterproc) then
    close (46)
    if(RUN3D) then
       if(dogzip3D) call systemf('gzip -f '//filename)
       print*, 'Writing MSE data. file:'//filename
    else
       print*, 'Appending MSE data. file:'//filename
    end if
  endif

end subroutine writeMSE

!---------------------------------------------------------------------
subroutine findindex(f,val,indval)
  ! finds the index in the vector f closest to val and returns it in indval

use grid, only : nzm

implicit none

real, intent(in) :: f(1:nzm)
real, intent(in) :: val
integer, intent(out) :: indval

real fdiff(1:nzm)
integer k

fdiff = abs(f-val)
indval = 1
do k=2,nzm
  if(fdiff(k) .lt. fdiff(indval)) indval = k
end do

end subroutine findindex



!-------------------------------------------------------------------
subroutine columnint(f,value)
  !  mass weighted column integral of scalar field
  !        returns int_0^TOA (rho*f dz)


use grid, only : adz, nzm, dz
use vars, only : rho

implicit none

real, intent(in)  :: f(1:nzm)
real, intent(out) :: value

integer k

!-----------------------------
value = 0.
do k=1,nzm
   value = value + f(k)*(rho(k)*dz*adz(k))
end do

end subroutine columnint

!-------------------------------------------------------------------
subroutine checkMSEinput()
  ! checks the MSE input; this is called from setparm.  Also allocates
  !  the MSE vars if needed

  use vars
  implicit none

  if(nsaveMSE.eq.0) nsaveMSE=nstat  ! nsaveMSE=nstat if not initialized

  !  need to check if already allocated since setparm can get
  !  called twice if nrestart=2
  if((nsaveMSEstart.le.nstop) .and. (.not.isallocateMSE))  then
    call allocateMSE()
    isallocateMSE = .true.
  elseif (isallocateMSE) then
           !bloss: do nothing
  else
    isallocateMSE=.false.
    doMSE=.false. ! If MSE not needed
    nsaveMSE=-1 ! If MSE not needed
  end if

  !  check the value of nsaveMSE3D
  if(.not.writeMSE3D .and. nsaveMSE3D.ne.1 .and. &
    mod(nsaveMSE3D,nsaveMSE).ne.0) then
    print*,'Error in specification of nsaveMSE3D.  Quitting.'
    call task_abort()
  end if

  ! because of restarting, doMSE only works if dooldrestart=false
  if(nsaveMSEstart.le.nstop .and. dooldrestart) then
    print*,'ERROR: MSE budgets set to output, but dooldrestart.true.'
    print*,'       For MSE restart, dooldrestart must be false'
    call task_abort()
  end if

  ! print MSE stuff to screen
  if(masterproc) then
    print*,'nsaveMSE,navgMSE,nsaveMSEstart,nsaveMSEend,nsaveMSE3D='
    print*,nsaveMSE,navgMSE,nsaveMSEstart,nsaveMSEend,nsaveMSE3D
  end if

end subroutine checkMSEinput


!-------------------------------------------------------------------
subroutine allocateMSE()
  ! allocates all the MSE variables and sets their initial values to 0

  use grid, only : nx, ny, nzm, nxp1, nyp1
  implicit none

  allocate(u_mse(nx,ny,nzm), v_mse(nx,ny,nzm), w_mse(nx,ny,nzm))
  allocate(tabs_mse(nx,ny,nzm), qrad_mse(nx,ny,nzm), rh_mse(nx,ny,nzm))
  allocate(qv_mse(nx,ny,nzm), qc_mse(nx,ny,nzm), qi_mse(nx,ny,nzm))
  allocate(qr_mse(nx,ny,nzm), qsg_mse(nx,ny,nzm))
  allocate(pp_mse(nx,ny,nzm))

  allocate(prec_mse(nx,ny), prec_frozen_mse(nx,ny))
  allocate(shf_mse(nx,ny), lhf_mse(nx,ny), sst_mse(nx,ny))
  allocate(hadvh_mse(nx,ny), vadvh_mse(nx,ny), hstor_mse(nx,ny))
  allocate(hadvs_mse(nx,ny), vadvs_mse(nx,ny), sstor_mse(nx,ny))
  allocate(hdifh_mse(nx,ny), hdifs_mse(nx,ny))
  allocate(rad_mse(nx,ny))
  allocate(lwns_mse(nx,ny), lwnt_mse(nx,ny), lwnsc_mse(nx,ny), lwntc_mse(nx,ny))
  allocate(swns_mse(nx,ny), swntm_mse(nx,ny),swnsc_mse(nx,ny),swntc_mse(nx,ny))
  allocate(slidamp_mse(nx,ny), hdamp_mse(nx,ny))
  allocate(huwpbl_mse(nx,ny), suwpbl_mse(nx,ny))

  allocate(pres_mse(nzm))

  allocate(u_uu(nx,ny,nzm), u_uv(nx,ny,nzm), u_uw(nx,ny,nzm))
  allocate(v_vu(nx,ny,nzm), v_vv(nx,ny,nzm), v_vw(nx,ny,nzm))
  allocate(u_diff(nx,ny,nzm), v_diff(nx,ny,nzm))
  allocate(ustor(nx,ny,nzm), vstor(nx,ny,nzm))
  allocate(u_uwpbl(nx,ny,nzm), v_uwpbl(nx,ny,nzm))

  allocate(dtu_uu(nxp1,ny,nzm,3), dtu_uv(nxp1,ny,nzm,3), dtu_uw(nxp1,ny,nzm,3))
  allocate(dtv_vu(nxp1,ny,nzm,3), dtv_vv(nx,nyp1,nzm,3), dtv_vw(nx,nyp1,nzm,3))
  allocate(dtu_diff(nxp1,ny,nzm,3), dtv_diff(nx,nyp1,nzm,3))

  allocate(oldsli_mse(nx,ny,nzm), oldhf_mse(nx,ny,nzm))
  allocate(oldu(nx,ny,nzm), oldv(nx,ny,nzm))
  allocate(div(nx,ny,nzm))
  allocate(prec_inst_mse(nx,ny), prec_inst_frozen_mse(nx,ny))

  allocate(u_mse_out(nx,ny,nzm), v_mse_out(nx,ny,nzm), w_mse_out(nx,ny,nzm))
  allocate(tabs_mse_out(nx,ny,nzm), qrad_mse_out(nx,ny,nzm))
  allocate(rh_mse_out(nx,ny,nzm))
  allocate(qv_mse_out(nx,ny,nzm), qc_mse_out(nx,ny,nzm), qi_mse_out(nx,ny,nzm))
  allocate(qr_mse_out(nx,ny,nzm), qsg_mse_out(nx,ny,nzm))
  allocate(pp_mse_out(nx,ny,nzm))
  allocate(u_uu_out(nx,ny,nzm), u_uv_out(nx,ny,nzm), u_uw_out(nx,ny,nzm))
  allocate(v_vu_out(nx,ny,nzm), v_vv_out(nx,ny,nzm), v_vw_out(nx,ny,nzm))
  allocate(u_diff_out(nx,ny,nzm), v_diff_out(nx,ny,nzm))
  allocate(ustor_out(nx,ny,nzm), vstor_out(nx,ny,nzm))
  allocate(u_uwpbl_out(nx,ny,nzm), v_uwpbl_out(nx,ny,nzm))


  ! initialize all fields to zero
  call zeroMSE()
  call zeroMSE_out()

end subroutine allocateMSE

!--------------------------------------------------------------------------
subroutine zeroMSE()

  use grid, only : firststep, nrestart

  implicit none

  u_mse=0.
  v_mse=0.
  w_mse=0.
  tabs_mse=0.
  qrad_mse=0.
  rh_mse=0.
  qv_mse=0.
  qc_mse=0.
  qi_mse=0.
  qr_mse=0.
  qsg_mse=0.
  pp_mse=0.

  prec_mse=0.
  prec_frozen_mse=0.
  shf_mse=0.
  lhf_mse=0.
  sst_mse=0.
  hadvh_mse=0.
  vadvh_mse=0.
  hstor_mse=0.
  hadvs_mse=0.
  vadvs_mse=0.
  sstor_mse=0.
  hdifh_mse=0.
  hdifs_mse=0.
  rad_mse=0.
  lwns_mse=0.
  lwnt_mse=0.
  lwnsc_mse=0.
  lwntc_mse=0.
  swns_mse=0.
  swntm_mse=0.
  swnsc_mse=0.
  swntc_mse=0.
  slidamp_mse=0.
  hdamp_mse=0.
  huwpbl_mse=0.
  suwpbl_mse=0.

  pres_mse=0.

  u_uu=0.
  u_uv=0.
  u_uw=0.
  v_vu=0.
  v_vv=0.
  v_vw=0.
  u_diff=0.
  v_diff=0.
  ustor=0.
  vstor=0.
  u_uwpbl=0.
  v_uwpbl=0.

  if(nrestart.eq.0 .and. firststep) then
    dtu_uu=0.
    dtu_uv=0.
    dtu_uw=0.
    dtv_vu=0.
    dtv_vv=0.
    dtv_vw=0.
    dtu_diff=0.
    dtv_diff=0.
  end if

  oldsli_mse=0.
  oldhf_mse=0.
  oldu=0.
  oldv=0.
  div=0.
  prec_inst_mse=0.
  prec_inst_frozen_mse=0.

end subroutine zeroMSE

!----------------------------------------------------------------------
subroutine zeroMSE_out()

  implicit none

  u_mse_out=0.
  v_mse_out=0.
  w_mse_out=0.
  tabs_mse_out=0.
  qrad_mse_out=0.
  rh_mse_out=0.
  qv_mse_out=0.
  qc_mse_out=0.
  qi_mse_out=0.
  qr_mse_out=0.
  qsg_mse_out=0.
  pp_mse_out=0.

  u_uu_out=0.
  u_uv_out=0.
  u_uw_out=0.
  v_vu_out=0.
  v_vv_out=0.
  v_vw_out=0.
  u_diff_out=0.
  v_diff_out=0.
  ustor_out=0.
  vstor_out=0.

  u_uwpbl_out=0.
  v_uwpbl_out=0.

end subroutine zeroMSE_out

!----------------------------------------------------------------------
subroutine copyMSE3D_to_MSE3Dout(coeffac)

  implicit none
  integer coeffac   ! factor to scale by

  real coef

  coef = float(nsaveMSE)/float(coeffac)

  u_mse_out=u_mse_out + coef*u_mse
  v_mse_out=v_mse_out + coef*v_mse
  w_mse_out=w_mse_out + coef*w_mse
  tabs_mse_out=tabs_mse_out + coef*tabs_mse
  rh_mse_out=rh_mse_out + coef*rh_mse
  qrad_mse_out=qrad_mse_out + coef*qrad_mse
  qv_mse_out=qv_mse_out + coef*qv_mse
  qc_mse_out=qc_mse_out + coef*qc_mse
  qi_mse_out=qi_mse_out + coef*qi_mse
  qr_mse_out=qr_mse_out + coef*qr_mse
  qsg_mse_out=qsg_mse_out + coef*qsg_mse
  pp_mse_out=pp_mse_out + coef*pp_mse

  u_uu_out=u_uu_out + coef*u_uu
  u_uv_out=u_uv_out + coef*u_uv
  u_uw_out=u_uw_out + coef*u_uw
  v_vu_out=v_vu_out + coef*v_vu
  v_vv_out=v_vv_out + coef*v_vv
  v_vw_out=v_vw_out + coef*v_vw
  u_diff_out=u_diff_out + coef*u_diff
  v_diff_out=v_diff_out + coef*v_diff
  ustor_out=ustor_out + coef*ustor
  vstor_out=vstor_out + coef*vstor

  u_uwpbl_out = u_uwpbl_out + coef*u_uwpbl
  v_uwpbl_out = v_uwpbl_out + coef*v_uwpbl

end subroutine copyMSE3D_to_MSE3Dout



!----------------------------------------------------------------------
subroutine deallocateMSE()
  ! deallocates all the MSE variables
  implicit none

  deallocate(u_mse, v_mse, w_mse)
  deallocate(tabs_mse, qrad_mse, rh_mse)
  deallocate(qv_mse, qc_mse, qi_mse)
  deallocate(qr_mse, qsg_mse)
  deallocate(pp_mse)

  deallocate(prec_mse, prec_frozen_mse, shf_mse, lhf_mse, sst_mse)
  deallocate(hadvh_mse, vadvh_mse, hstor_mse)
  deallocate(hadvs_mse, vadvs_mse, sstor_mse)
  deallocate(hdifh_mse, hdifs_mse)
  deallocate(rad_mse, lwns_mse, lwnt_mse, lwnsc_mse, lwntc_mse)
  deallocate(swns_mse, swntm_mse,swnsc_mse,swntc_mse)
  deallocate(slidamp_mse, hdamp_mse)
  deallocate(huwpbl_mse, suwpbl_mse)

  deallocate(pres_mse)

  deallocate(u_uu, u_uv, u_uw)
  deallocate(v_vu, v_vv, v_vw)
  deallocate(u_diff, v_diff)
  deallocate(ustor, vstor)

  deallocate(oldsli_mse, oldhf_mse)
  deallocate(oldu, oldv)
  deallocate(div)
  deallocate(prec_inst_mse, prec_inst_frozen_mse)

  deallocate(u_mse_out, v_mse_out, w_mse_out)
  deallocate(tabs_mse_out, qrad_mse_out, rh_mse_out)
  deallocate(qv_mse_out, qc_mse_out, qi_mse_out)
  deallocate(qr_mse_out, qsg_mse_out)
  deallocate(pp_mse_out)
  deallocate(u_uu_out, u_uv_out, u_uw_out)
  deallocate(v_vu_out, v_vv_out, v_vw_out)
  deallocate(u_diff_out, v_diff_out)
  deallocate(ustor_out, vstor_out)
  deallocate(u_uwpbl, v_uwpbl, u_uwpbl_out, v_uwpbl_out)


end subroutine deallocateMSE


end module mse

