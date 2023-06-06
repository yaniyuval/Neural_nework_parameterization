     
subroutine write_fields3D
	
use vars
use microscaling

implicit none
character *80 filename,long_name
character *8 name
character *10 timechar
character *4 icyclechar
character *4 rankchar
character *6 filetype
character *10 units
character *10 c_z(nzm),c_p(nzm),c_dx, c_dy, c_time
integer i,j,k,nfields,nfields1
real tmp(nx,ny,nzm)

nfields=11 ! number of 3D fields to save
if(.not.(doshortwave.or.dolongwave).or.doradhomo) nfields=nfields-1
if(.not.docloud) nfields=nfields-1
if(.not.doprecip) nfields=nfields-1
!kzm added 4 tracers + tro3
if (dotrx) nfields=nfields+1
if (dotry) nfields=nfields+1
if (domicroscaling) nfields=nfields+2
if (dotrz) nfields=nfields+1
if (dotrzz) nfields=nfields+1
if (dotro3) nfields=nfields+1
!bloss added 5 fields for dokruegermicro
if (dokruegermicro) nfields = nfields + 3
!cw added temperature tendency from uwcu
if (douwcu) nfields = nfields + 1
if (douwpbl) nfields = nfields + 4
if (uwpbl_ndiv.gt.1) nfields = nfields + 2

nfields1=0

if(masterproc) then

  write(rankchar,'(i4)') nsubdomains
  write(timechar,'(i10)') nstep
  write(icyclechar,'(i4)') icycle
  do k=1,11-lenstr(timechar)-1
    timechar(k:k)='0'
  end do

  do k=1,5-lenstr(icyclechar)-1
    icyclechar(k:k)='0'
  end do

  if(RUN3D) then
    if(save3Dbin) then
      filetype = '.bin3D'
    else
      filetype = '.com3D'
    end if
    filename='./DATA3D/'//trim(case)//'_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)// &
        '_'//icyclechar(1:4)//filetype
    open(46,file=filename,status='unknown',form='unformatted')

  else
    if(save3Dbin) then
     if(save3Dsep) then
       filetype = '.bin3D'
     else
       filetype = '.bin2D'
     end if
    else
     if(save3Dsep) then
       filetype = '.com3D'
     else
       filetype = '.com2D'
     end if
    end if
    if(save3Dsep) then
      filename='./DATA3D/'//trim(case)//'_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype
      open(46,file=filename,status='unknown',form='unformatted')	
    else
      filename='./DATA3D/'//trim(case)//'_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//filetype
      if(nrestart.eq.0.and.notopened3D) then
         open(46,file=filename,status='unknown',form='unformatted')	
      else
         open(46,file=filename,status='unknown', &
                              form='unformatted', position='append')
      end if
      notopened3D=.false.
    end if  

  end if

  if(save3Dbin) then

    write(46) nx,ny,nzm,nsubdomains,nsubdomains_x,nsubdomains_y,nfields
    do k=1,nzm
      write(46) z(k) 
    end do
    do k=1,nzm
      write(46) pres(k)
    end do
    write(46) dx
    write(46) dy
    write(46) nstep*dt/(3600.*24.)+day0

  else

    write(long_name,'(8i4)') nx,ny,nzm,nsubdomains, &
                                   nsubdomains_x,nsubdomains_y,nfields
    do k=1,nzm
       write(c_z(k),'(f10.3)') z(k)
    end do
    do k=1,nzm
       write(c_p(k),'(f10.3)') pres(k)
    end do
    write(c_dx,'(f10.0)') dx
    write(c_dy,'(f10.0)') dy
    write(c_time,'(f10.2)') nstep*dt/(3600.*24.)+day0
	
    write(46) long_name(1:32)
    write(46) c_time,c_dx,c_dy, (c_z(k),k=1,nzm),(c_p(k),k=1,nzm)

  end if ! save3Dbin

end if ! masterproc

if(douwpbl) then
   nfields1=nfields1+1
   do k=1,nzm
      do j=1,ny
         do i=1,nx
            tmp(i,j,k)=ttend_uwpbl(i,j,k)*86400
         end do
      end do
   end do
   name='ttend_uwpbl'
   long_name='Temperature tendency from UW PBL scheme'
   units='K/day'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
        save3Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do k=1,nzm
      do j=1,ny
         do i=1,nx
            tmp(i,j,k)=qtend_uwpbl(i,j,k)*86400
         end do
      end do
   end do
   name='qtend_uwpbl'
   long_name='Total moisture tendency from UW PBL scheme'
   units='kg/kg/day'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
        save3Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do k=1,nzm
      do j=1,ny
         do i=1,nx
            tmp(i,j,k)=utend_uwpbl(i,j,k)
         end do
      end do
   end do
   name='uten_uwpbl'
   long_name='Zonal wind tendency from UW PBL scheme'
   units='m/s/s'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
        save3Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do k=1,nzm
      do j=1,ny
         do i=1,nx
            tmp(i,j,k)=vtend_uwpbl(i,j,k)
         end do
      end do
   end do
   name='vten_uwpbl'
   long_name='Meridional wind tendency from UW PBL scheme'
   units='m/s/s'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
        save3Dbin,dompi,rank,nsubdomains)

!!$   nfields1=nfields1+1
!!$   do k=1,nzm
!!$      do j=1,ny
!!$         do i=1,nx
!!$            tmp(i,j,k)=tke_uw(i,j,k)
!!$         end do
!!$      end do
!!$   end do
!!$   name='tke_uw'
!!$   long_name='TKE from UW PBL scheme'
!!$   units='m.m/s/s'
!!$   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
!!$        save3Dbin,dompi,rank,nsubdomains)
end if

if(uwpbl_ndiv.gt.1) then

   nfields1=nfields1+1
   do k=1,nzm
      do j=1,ny
         do i=1,nx
            tmp(i,j,k)=tk_z_uw(i,j,k)
         end do
      end do
   end do
   name='tkz_uw'
   long_name='TKZ_UW from patch-averaged UW PBL scheme'
   units='m.m/s'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
        save3Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do k=1,nzm
      do j=1,ny
         do i=1,nx
            tmp(i,j,k)=tkh_z_uw(i,j,k)
         end do
      end do
   end do
   name='tkhz_uw'
   long_name='TKHZ_UW from patch-averaged UW PBL scheme'
   units='m.m/s'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
        save3Dbin,dompi,rank,nsubdomains)

end if

if(douwcu) then
   nfields1=nfields1+1
   do k=1,nzm
      do j=1,ny
         do i=1,nx
            tmp(i,j,k)=ttend_uwcu(i,j,k)*86400.0
         end do
      end do
   end do
   name='ttend_uwcu'
   long_name='t tendency from uwcu'
   units='K/day'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
        save3Dbin,dompi,rank,nsubdomains)

end if


nfields1=nfields1+1
do k=1,nzm
   do j=1,ny
      do i=1,nx
         tmp(i,j,k)=tk_xy(i,j,k)
      end do
   end do
end do
name='tk_xy'
long_name='horizontal diffusivity'
units='m.m/s'
call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
     save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nzm
   do j=1,ny
      do i=1,nx
         tmp(i,j,k)=tk_z(i,j,k)
      end do
   end do
end do
name='tk_z'
long_name='vertical diffusivity'
units='m.m/s'
call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
        save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=u(i,j,k)
    end do
   end do
  end do
  name='U'
  long_name='X Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=v(i,j,k)
    end do
   end do
  end do
  name='V'
  long_name='Y Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=w(i,j,k)
    end do
   end do
  end do
  name='W'
  long_name='Z Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=p(i,j,k)
    end do
   end do
  end do
  name='PP'
  long_name='Pressure Perturbation'
  units='Pa'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)


if((dolongwave.or.doshortwave).and..not.doradhomo) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=misc(i,j,k)*86400.
    end do
   end do
  end do
  name='QRAD'
  long_name='Radiative heating rate'
  units='K/day'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=tabs(i,j,k)
    end do
   end do
  end do
  name='TABS'
  long_name='Absolute Temperature'
  units='K'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=(q(i,j,k)-qn(i,j,k))*1.e3
    end do
   end do
  end do
  name='Q'
  long_name='Water Vapor'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

if((docloud).and.(.not.dokruegermicro)) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qn(i,j,k)*1.e3
    end do
   end do
  end do
  name='QN'
  long_name='Non-precipitating Condensate (Water+Ice)'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if


if((doprecip).and.(.not.dokruegermicro)) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qp(i,j,k)*1.e3
    end do
   end do
  end do
  name='QP'
  long_name='Precipitating Water (Rain+Snow)'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if

!=================================================================
! UW ADDITIONS

!bloss krueger microphysics
if((docloud).and.(dokruegermicro)) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qc(i,j,k)*1.e3
    end do
   end do
  end do
  name='QC'
  long_name='Cloud Water'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qi(i,j,k)*1.e3
    end do
   end do
  end do
  name='QI'
  long_name='Cloud Ice'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if

if((doprecip).and.(dokruegermicro)) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qr(i,j,k)*1.e3
    end do
   end do
  end do
  name='QR'
  long_name='Rain'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qg(i,j,k)*1.e3
    end do
   end do
  end do
  name='QG'
  long_name='Graupel'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qs(i,j,k)*1.e3
    end do
   end do
  end do
  name='QS'
  long_name='Snow'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if

!kzm tracer
 if(dotrx) then 
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=trx(i,j,k)
    end do
   end do
  end do
  name='TRX'
  long_name='X coordinate tracer'
  units=''
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
 endif
 
 if(dotry) then 
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=try(i,j,k)
    end do
   end do
  end do
  name='TRY'
  long_name='Y coordinate tracer'
  units=''
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
 endif

 if(dotrz) then 
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=trz(i,j,k)
    end do
   end do
  end do
  name='TRZ'
  long_name='Z coordinate tracer'
  units=''
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
 endif

 if(dotrzz) then 
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=trzz(i,j,k)
    end do
   end do
  end do
  name='TRZ2'
  long_name='Z2 coordinate tracer'
  units=''
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
 endif

  if(dotro3) then 
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=tro3(i,j,k)*1.e6
    end do
   end do
  end do
  name='TRO3'
  long_name='Interactive O3'
  units='ppmm'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
 endif


! JA 8/14/06
 
if(domicroscaling) then
  nfields1=nfields1+2
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=real_stratiformnessc(i,j,k)
    end do
   end do
  end do
  name='STR'
  long_name='Stratiformness index (JA)'
  units=''
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=dtn_scaled(i,j,k)
    end do
   end do
  end do
  name='DTN'
  long_name='scaled time step(JA)'
  units=''
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

endif



! Finished outputing 3D fields.

 if (doreset_tracer) call settracer()  ! Reset tracer if necessary.

! END UW ADDITIONS
!=================================================================

  call task_barrier()

  if(nfields.ne.nfields1) then
    if(masterproc) print*,'write_fields3D error: nfields'
    call task_abort()
  end if
  if(masterproc) then
    close (46)
    if(RUN3D.or.save3Dsep) then
       if(dogzip3D) call systemf('gzip -f '//filename)
       print*, 'Writting 3D data. file:'//filename
    else
       print*, 'Appending 3D data. file:'//filename
    end if
  endif

end
