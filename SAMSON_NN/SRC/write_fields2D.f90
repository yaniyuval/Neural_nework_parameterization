     
subroutine write_fields2D
	
use vars
use params
use rad
use simple_ocean, only : doevolveSSTforcing       ! peters
use simple_land                                   ! peters
implicit none
character *80 filename,long_name
character *8 name
character *10 timechar
character *4 rankchar
character *6 filetype
character *10 units
character *10 c_dx, c_dy, c_time
integer i,j,k,nfields,nfields1
real tmp(nx,ny,nzm)
character*3 filestatus

nfields= 24 ! added PSFC SWNTM to standard output
if(douwpbl) nfields = nfields + 2 ! cw added pblh
if(douwcu) nfields = nfields + 3 ! cw added cin, cbmf, uwcu_err
if((.not.doxy).and.(.not.dodynamicocean).and.(ocean_type.eq.0)) nfields=nfields-1 ! no SST
if(.not.dolongwave) nfields = nfields-4  ! no LWNT LWNTC LWNS LWNSC
if(.not.doshortwave) nfields = nfields-6 ! no SOLIN LWNS LWNSC SWNS SWNSC SWNTM
if(.not.dodynamicocean) nfields=nfields-1! no QOCEAN
if(.not.doprecip) nfields = nfields-1    ! no PREC
if(.not.docloud) nfields = nfields-2     ! no CWP IWP
if(SFC_FLX_FXD) nfields = nfields-2      ! no SHF LHF
if(dosaturatedwvp) nfields = nfields + 1  ! Add SWVP
if(doreflectivity) nfields=nfields+3      ! Add dbZe1km 3km 6km
if(doinstantolr) nfields = nfields + 2    ! Add LWNTI LWNTCI
if(docolumnbudgets) nfields = nfields + 4 ! Add SADV SSTOR HADV HSTOR
if(do850mbarwinds) nfields = nfields + 2  ! Add U850 V850
if(dooutputheights) nfields = nfields + 2 ! Add H850 H200
if(dolowlevelthetae) nfields = nfields + 4! Add THE850 THES850 THE1000 THES1000
if(dowaveoutput) nfields = nfields + 8    ! proj. onto 1/2nd baroclinic modes
if(dowaveenergetics) nfields = nfields + 6 ! 1/2nd baroclinic modes energetics
if(doinstantoutput) nfields = nfields + 11! Add a bunch of inst. fields
if(docloudechoheights) nfields = nfields + 3! Add cloud/echo heights & cttemp
if(pres(nzm).gt.500.) nfields = nfields - 1 ! no W500
if(doxysoil_wet) nfields = nfields+1        ! peters output soilwetness
nfields1=0

if(masterproc) then

  write(rankchar,'(i4)') nsubdomains
  write(timechar,'(i10)') nstep
  do i=1,11-lenstr(timechar)-1
    timechar(i:i)='0'
  end do


! Make sure that the new run doesn't overwrite the file from the old run 

    filestatus='old'
    if(notopened2D.and.(nrestart.eq.0.or.nrestart.eq.2)) then
      filestatus='new'
    end if

    if(save2Dbin) then
      filetype = '.2Dbin'
    else
      filetype = '.2Dcom'
    end if
    if(save2Dsep) then
       filename='./DATA3D/'//trim(case)//'_'//trim(caseid)//'_'// &
          rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype 
          open(46,file=filename,status='unknown',form='unformatted')
    else
       filename='./DATA3D/'//trim(case)//'_'//trim(caseid)//'_'// &
           rankchar(5-lenstr(rankchar):4)//filetype
	      
       if(nrestart.eq.0.and.notopened2D) then
          open(46,file=filename,status=filestatus,form='unformatted')	
       else
          open(46,file=filename,status=filestatus,form='unformatted', &
                                                     position='append')
       end if
       notopened2D=.false.
    end if

     if(save2Dbin) then

       write(46) nx,ny,nzm,nsubdomains, nsubdomains_x,nsubdomains_y,nfields
       write(46) dx
       write(46) dy
       write(46) day  !nstep*dt/(3600.*24.)+day0

     else

       write(long_name,'(8i4)') nx,ny,nzm,nsubdomains,  &
                                        nsubdomains_x,nsubdomains_y,nfields
       write(c_dx,'(f10.0)') dx
       write(c_dy,'(f10.0)') dy
       write(c_time,'(f10.2)') day  !nstep*dt/(3600.*24.)+day0
       write(46) long_name(1:32)
       write(46) c_time,c_dx,c_dy

     end if ! save2Dbin

end if! masterproc


! 2D fields:

! PBLH added by cw 2/23/06:
if (douwpbl) then
   nfields1=nfields1+1
   do j=1,ny
      do i=1,nx
         tmp(i,j,1)=pblh_xy(i,j)/float(nsave2D)
         pblh_xy(i,j) = 0.0
      end do
   end do
   name='PBLH'
   long_name='PBL height from UW PBL scheme'
   units='m'
   call compress3D(tmp,nx,ny,1,name,long_name,units, &
        save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
      do i=1,nx
         tmp(i,j,1) = float(uwpbl_ndiff(i,j))
      end do
   end do
   name='uwpbl_ndiff'
   long_name='number of PBL iterations'
   units='none'
   call compress3D(tmp,nx,ny,1,name,long_name,units, &
        save2Dbin,dompi,rank,nsubdomains)
end if

if (douwcu) then
   nfields1=nfields1+1
   do j=1,ny
      do i=1,nx
         tmp(i,j,1)=cin_xy(i,j)/float(nsave2D)
         cin_xy(i,j) = 0.0
      end do
   end do
   name='cin'
   long_name='Convective inhibition'
   units='J/kg'
   call compress3D(tmp,nx,ny,1,name,long_name,units, &
        save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
      do i=1,nx
         tmp(i,j,1)=uwcu_err(i,j)
      end do
   end do
   name='uwcu_err'
   long_name='UWCU error flag'
   units=''
   call compress3D(tmp,nx,ny,1,name,long_name,units, &
        save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
      do i=1,nx
         tmp(i,j,1)=cbmf_xy(i,j)/float(nsave2D)
         cbmf_xy(i,j) = 0.0
      end do
   end do
   name='cbmf'
   long_name='Cloud base mass flux'
   units='J/kg'
   call compress3D(tmp,nx,ny,1,name,long_name,units, &
        save2Dbin,dompi,rank,nsubdomains)
end if

if(doprecip) then
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j,1)=prec_xy(i,j)*dz/dt*86400./float(nsave2D)
      prec_xy(i,j) = 0.
    end do
   end do
  name='Prec'
  long_name='Surface Precip. Rate'
  units='mm/day'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
end if

if(.not.SFC_FLX_FXD) then
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j,1)=shf_xy(i,j)*rhow(1)*cp/float(nsave2D)
      shf_xy(i,j) = 0.
    end do
   end do
  name='SHF'
  long_name='Sensible Heat Flux'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j,1)=lhf_xy(i,j)*rhow(1)*lcond/float(nsave2D)
      lhf_xy(i,j) = 0.
    end do
   end do
  name='LHF'
  long_name='Latent Heat Flux'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
        
end if

if(dolongwave) then
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=lwns_xy(i,j)/float(nsave2D)
       lwns_xy(i,j) = 0.
     end do
   end do
  name='LWNS'
  long_name='Net LW at the surface'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=lwnsc_xy(i,j)/float(nsave2D)
       lwnsc_xy(i,j) = 0.
     end do
   end do
  name='LWNSC'
  long_name='Net clear-sky LW at the surface'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=lwnt_xy(i,j)/float(nsave2D)
       lwnt_xy(i,j) = 0.
     end do
   end do
  name='LWNT'
  long_name='Net LW at TOA'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=lwntc_xy(i,j)/float(nsave2D)
       lwntc_xy(i,j) = 0.
     end do
   end do
  name='LWNTC'
  long_name='Clear-Sky Net LW at TOA'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

if(doshortwave) then
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=solin_xy(i,j)/float(nsave2D)
       solin_xy(i,j) = 0.
     end do
   end do
  name='SOLIN'
  long_name='Solar TOA insolation'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=swns_xy(i,j)/float(nsave2D)
       swns_xy(i,j) = 0.
     end do
   end do
  name='SWNS'
  long_name='Net SW at the surface'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=swnsc_xy(i,j)/float(nsave2D)
       swnsc_xy(i,j) = 0.
     end do
   end do
  name='SWNSC'
  long_name='Net Clear-sky SW at the surface'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=swnt_xy(i,j)/float(nsave2D)
       swnt_xy(i,j) = 0.
     end do
   end do
  name='SWNT'
  long_name='Net SW at TOA'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=swntc_xy(i,j)/float(nsave2D)
       swntc_xy(i,j) = 0.
     end do
   end do
  name='SWNTC'
  long_name='Net Clear-Sky SW at TOA'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! top-of-model net downward sw flux, useful for energy budgets
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=swntm_xy(i,j)/float(nsave2D)
       swntm_xy(i,j) = 0.
     end do
   end do
  name='SWNTM'
  long_name='top-of-model net SW_down'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

if(docloud) then
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=cw_xy(i,j)/float(nsave2D)
       cw_xy(i,j) = 0.
     end do
   end do
  name='CWP'
  long_name='Cloud Water Path'
  units='mm'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=iw_xy(i,j)/float(nsave2D)
       iw_xy(i,j) = 0.
     end do
   end do
  name='IWP'
  long_name='Ice Path'
  units='mm'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=pw_xy(i,j)/float(nsave2D)
     pw_xy(i,j) = 0.
     end do
   end do
  name='PW'
  long_name='Precipitable Water'
  units='mm'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=usfc_xy(i,j)/float(nsave2D)
       usfc_xy(i,j) = 0.
     end do
   end do
  name='USFC'
  long_name='U at the surface'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)


   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
      tmp(i,j,1)=vsfc_xy(i,j)/float(nsave2D)
      vsfc_xy(i,j) = 0.
     end do
   end do
  name='VSFC'
  long_name='V at the surface'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=u200_xy(i,j)/float(nsave2D)
       u200_xy(i,j) = 0.
     end do
   end do
  name='U200'
  long_name='U at 200 mb'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=v200_xy(i,j)/float(nsave2D)
       v200_xy(i,j) = 0.
     end do
   end do
  name='V200'
  long_name='V at 200 mb'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

if(pres(nzm).lt.500.) then
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=w500_xy(i,j)/float(nsave2D)
       w500_xy(i,j) = 0.
     end do
   end do
  name='W500'
  long_name='OMEGA at 500 mb'
  units='Pa/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
end if

if(dodynamicocean) then
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j,1)=qocean_xy(i,j)/float(nsave2D)
      qocean_xy(i,j) = 0.
    end do
   end do
  name='QOCN'
  long_name='Deep Ocean Cooling'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
end if

if((doxy).or.(dodynamicocean).or.(ocean_type.ne.0)) then
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j,1)=sstxy(i,j)
    end do
   end do
  name='SST'
  long_name='Sea Surface Temperature'
  units='K'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
end if

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=psfc_xy(i,j)/float(nsave2D)/100.
       psfc_xy(i,j) = 0.
     end do
   end do
  name='PSFC'
  long_name='P at the surface'
  units='mbar'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

if (dosaturatedwvp) then
  nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=swvp_xy(i,j)/float(nsave2D)
       swvp_xy(i,j) = 0.
     end do
   end do
  name='SWVP'
  long_name='Saturated Water Vapor Path'
  units='mm'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
       save2Dbin,dompi,rank,nsubdomains)
end if

if(doreflectivity) then
   nfields1=nfields1+1
   do j=1,ny
      do i=1,nx
         tmp(i,j,1)=dbZe1km_xy(i,j)
      end do
   end do
   name='dbZe1km'
   long_name='Equiv Radar Refl at z = 1km'
   units='dbZe'
   call compress3D(tmp,nx,ny,1,name,long_name,units, &
        save2Dbin,dompi,rank,nsubdomains)
        
   nfields1=nfields1+1
   do j=1,ny
      do i=1,nx
         tmp(i,j,1)=dbZe3km_xy(i,j)
      end do
   end do
   name='dbZe3km'
   long_name='Equiv Radar Refl at z = 3km'
   units='dbZe'
   call compress3D(tmp,nx,ny,1,name,long_name,units, &
        save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
      do i=1,nx
         tmp(i,j,1)=dbZe6km_xy(i,j)
      end do
   end do
   name='dbZe6km'
   long_name='Equiv Radar Refl at z = 6km'
   units='dbZe'
   call compress3D(tmp,nx,ny,1,name,long_name,units, &
        save2Dbin,dompi,rank,nsubdomains)
end if

if(doinstantolr) then 

   ! instantaneous top-of-atmosphere net upward lw flux
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=lwntxy(i,j)
     end do
   end do
  name='LWNTI'
  long_name='Instant. top-of-atmos net LW_up'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! instantaneous top-of-atmosphere net upward clearsky lw flux
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=lwntcxy(i,j)
     end do
   end do
  name='LWNTCI'
  long_name='Instant. clearsky top-of-atmos net LW_up'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

if(docolumnbudgets) then 

   ! save storage and advection terms for column-integrated frozen
   ! moist static energy and liquid water ice static energy budgets.
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=hstor_xy(i,j)/1000.
       hstor_xy(i,j) = 0.
     end do
   end do
  name='HSTOR'
  long_name='Column-int. Frozen MSE Storage'
  units='kW/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=hadv_xy(i,j)/1000.
       hadv_xy(i,j) = 0.
     end do
   end do
  name='HADV'
  long_name='Column-int. Frozen MSE Advection'
  units='kW/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=sstor_xy(i,j)/1000.
       sstor_xy(i,j) = 0.
     end do
   end do
  name='SSTOR'
  long_name='Column-int. S-LI Storage'
  units='kW/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=sadv_xy(i,j)/1000.
       sadv_xy(i,j) = 0.
     end do
   end do
  name='SADV'
  long_name='Column-int. S-LI Advection'
  units='kW/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
end if

if(do850mbarwinds) then 

   ! 850 mbar zonal velocity
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=u850_xy(i,j)/float(nsave2D)
       u850_xy(i,j) = 0.
     end do
   end do
  name='U850'
  long_name='850 mbar zonal velocity'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! meridional wind at 850 mbar
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=v850_xy(i,j)/float(nsave2D)
       v850_xy(i,j) = 0.
     end do
   end do
  name='V850'
  long_name='850 mbar meridional velocity'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

if(dooutputheights) then 

   ! 850 mbar geopotential height
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=h850_xy(i,j)/float(nsave2D)
       h850_xy(i,j) = 0.
     end do
   end do
  name='H850'
  long_name='850 mbar geopotential height'
  units='m'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! 200 mbar geopotential height
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=h200_xy(i,j)/float(nsave2D)
       h200_xy(i,j) = 0.
     end do
   end do
  name='H200'
  long_name='200 mbar geopotential height'
  units='m'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

if(docloudechoheights) then 

   ! cloud top height
   nfields1=nfields1+1
   do j=1,ny
      do i=1,nx
         tmp(i,j,1)=cloudtopheight(i,j)/1000.
         cloudtopheight(i,j) = 0.
      end do
   end do
   name='ZC'
   long_name='Cloud top height'
   units='km'
   call compress3D(tmp,nx,ny,1,name,long_name,units, &
        save2Dbin,dompi,rank,nsubdomains)

   ! cloud top temperature
   nfields1=nfields1+1
   do j=1,ny
      do i=1,nx
         tmp(i,j,1)=cloudtoptemp(i,j)
         cloudtoptemp(i,j) = 0.
      end do
   end do
   name='TB'
   long_name='Cloud top temperature'
   units='K'
   call compress3D(tmp,nx,ny,1,name,long_name,units, &
        save2Dbin,dompi,rank,nsubdomains)

   ! echo top height
   nfields1=nfields1+1
   do j=1,ny
      do i=1,nx
         tmp(i,j,1)=echotopheight(i,j)/1000.
         echotopheight(i,j) = 0.
      end do
   end do
   name='ZE'
   long_name='Echo top height'
   units='km'
   call compress3D(tmp,nx,ny,1,name,long_name,units, &
        save2Dbin,dompi,rank,nsubdomains)

end if

if(dolowlevelthetae) then 

   ! 850 mbar equivalent potential temperature
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=the850_xy(i,j)/float(nsave2D)
       the850_xy(i,j) = 0.
     end do
   end do
  name='THE850'
  long_name='850 mbar equivalent potential temperature'
  units='K'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! 850 mbar saturated equivalent potential temperature
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=thes850_xy(i,j)/float(nsave2D)
       thes850_xy(i,j) = 0.
     end do
   end do
  name='THES850'
  long_name='850 mbar saturated equivalent potential temperature'
  units='K'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! 1000 mbar equivalent potential temperature
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=the1000_xy(i,j)/float(nsave2D)
       the1000_xy(i,j) = 0.
     end do
   end do
  name='THE1000'
  long_name='1000 mbar equivalent potential temperature'
  units='K'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! 1000 mbar saturated equivalent potential temperature
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=thes1000_xy(i,j)/float(nsave2D)
       thes1000_xy(i,j) = 0.
     end do
   end do
  name='THES1000'
  long_name='1000 mbar saturated equivalent potential temperature'
  units='K'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

if(dowaveoutput) then 

   ! First mode vertical velocity, time-averaged
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=wmode1_xy(i,j)/float(nsave2D)
       wmode1_xy(i,j) = 0.
     end do
   end do
  name='WMODE1'
  long_name='First mode vertical velocity'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! Second mode vertical velocity, time-averaged
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=wmode2_xy(i,j)/float(nsave2D)
       wmode2_xy(i,j) = 0.
     end do
   end do
  name='WMODE2'
  long_name='Second mode vertical velocity'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! First mode water vapor perturbation, time-averaged
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=qmode1_xy(i,j)/float(nsave2D)
       qmode1_xy(i,j) = 0.
     end do
   end do
  name='QMODE1'
  long_name='First mode water vapor perturbation'
  units='g/kg'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! Second mode water vapor perturbation, time-averaged
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=qmode2_xy(i,j)/float(nsave2D)
       qmode2_xy(i,j) = 0.
     end do
   end do
  name='QMODE2'
  long_name='Second mode water vapor perturbation'
  units='g/kg'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! First mode theta perturbation, time-averaged
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=thmode1_xy(i,j)/float(nsave2D)
       thmode1_xy(i,j) = 0.
     end do
   end do
  name='THMODE1'
  long_name='First mode theta perturbation'
  units='K'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! Second mode theta perturbation, time-averaged
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=thmode2_xy(i,j)/float(nsave2D)
       thmode2_xy(i,j) = 0.
     end do
   end do
  name='THMODE2'
  long_name='Second mode theta perturbation'
  units='K'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! First mode vertical velocity, instantaneous
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=wmode1i_xy(i,j)
       wmode1_xy(i,j) = 0.
     end do
   end do
  name='WMODE1I'
  long_name='instantaneous first mode vertical velocity'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! Second mode vertical velocity, instantaneous
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=wmode2i_xy(i,j)
       wmode2_xy(i,j) = 0.
     end do
   end do
  name='WMODE2I'
  long_name='instantaneous second mode vertical velocity'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

if(dowaveenergetics) then

   ! First mode vertical velocity
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=wdwdt1_xy(i,j)/float(nsave2D)*3600.
       wdwdt1_xy(i,j) = 0.
     end do
   end do
  name='WDWDT1'
  long_name='First mode vertical velocity energetics'
  units='m2/s2/hour'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! Second mode vertical velocity
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=wdwdt2_xy(i,j)/float(nsave2D)*3600.
       wdwdt2_xy(i,j) = 0.
     end do
   end do
  name='WDWDT2'
  long_name='Second mode vertical velocity energetics'
  units='m2/s2/hour'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! First mode liquid-ice static energy
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=sdsdt1_xy(i,j)/float(nsave2D)*3600.
       sdsdt1_xy(i,j) = 0.
     end do
   end do
  name='SDSDT1'
  long_name='First mode liquid-ice static energy energetics'
  units='K2/hour'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! Second mode liquid-ice static energy
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=sdsdt2_xy(i,j)/float(nsave2D)*3600.
       sdsdt2_xy(i,j) = 0.
     end do
   end do
  name='SDSDT2'
  long_name='Second mode liquid-ice static energy energetics'
  units='K2/hour'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! First mode frozen moist static energy
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=hdhdt1_xy(i,j)/float(nsave2D)*3600.
       hdhdt1_xy(i,j) = 0.
     end do
   end do
  name='HDHDT1'
  long_name='First mode frozen moist static energy energetics'
  units='K2/hour'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! Second mode frozen moist static energy
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=hdhdt2_xy(i,j)/float(nsave2D)*3600.
       hdhdt2_xy(i,j) = 0.
     end do
   end do
  name='HDHDT2'
  long_name='Second mode frozen moist static energy energetics'
  units='K2/hour'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

if(doinstantoutput) then 

   ! instantaneous precipitation flux
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=prec_inst(i,j)*86400.*1000./rhor
     end do
   end do
  name='PRECI'
  long_name='Inst. Surface Precipitation Rate'
  units='mm/day'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! instantaneous surface evaporation rate
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=fluxbq(i,j)*rhow(1)*86400.*1000./rhor
     end do
   end do
  name='EVAPI'
  long_name='Inst. Surface Evaporation Rate'
  units='mm/day'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! instantaneous precipitable water
   nfields1=nfields1+1
   tmp = 0.
   do k = 1,nzm
      do j=1,ny
         do i=1,nx
            tmp(i,j,1)=tmp(i,j,1) + dz*adz(k)*rho(k)*q(i,j,k)
         end do
      end do
   end do
  name='PWI'
  long_name='Inst. Precipitable Water'
  units='kg/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! instantaneous column integrated condensate (qc+qi+qs)
   nfields1=nfields1+1
   if(.not.dokruegermicro) then
   tmp = 0.
   do k = 1,nzm
      do j=1,ny
         do i=1,nx
            tmp(i,j,1)=tmp(i,j,1) + dz*adz(k)*rho(k)*(qn(i,j,k) &
                 + (1.-omegap(tabs(i,j,k)))*(1.-omegag(tabs(i,j,k)))*qp(i,j,k))
         end do
      end do
   end do
   elseif(dokruegermicro) then
   do k = 1,nzm
      do j=1,ny
         do i=1,nx
            tmp(i,j,1)=tmp(i,j,1) + dz*adz(k)*rho(k)*(qn(i,j,k) + qs(i,j,k))
         end do
      end do
   end do
   end if
  name='CONDI'
  long_name='Inst. Column-integrated Condensate (cloud+snow)'
  units='kg/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! instantaneous top-of-atmosphere net upward lw flux
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=lwntxy(i,j)
     end do
   end do
  name='LWNTI'
  long_name='Inst. top-of-atmos net LW_up'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! instantaneous top-of-atmosphere net upward sw flux
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=swntxy(i,j)
     end do
   end do
  name='SWNTI'
  long_name='Inst. top-of-atmos net SW_up'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! instantaneous top-of-atmosphere net upward clearsky lw flux
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=lwntcxy(i,j)
     end do
   end do
  name='LWNTCI'
  long_name='Inst. clearsky top-of-atmos net LW_up'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! instantaneous top-of-atmosphere net upward clearsky sw flux
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=swntcxy(i,j)
     end do
   end do
  name='SWNTCI'
  long_name='Inst. clearsky top-of-atmos net SW_up'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! instantaneous lowest-level temperature
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=tabs(i,j,1)
     end do
   end do
  name='TSFC'
  long_name='Inst. lowest-level temperature'
  units='K'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! instantaneous lowest-level specific humidity
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=1000.*(q(i,j,1)-qn(i,j,1))
     end do
   end do
  name='QSFC'
  long_name='Inst. lowest-level specific humidity'
  units='g/kg'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! instantaneous lowest-level wind speed
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=sqrt(u(i,j,1)*u(i,j,1) + v(i,j,1)*v(i,j,1))
     end do
   end do
  name='WSSFC'
  long_name='Inst. lowest-level wind speed'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

! peters.  output soil wetness if needed
if(doxysoil_wet) then
   nfields1=nfields1+1
   tmp(:,:,1) = xysoil_wet(1:nx,1:ny)
   name = 'SOILWET'
   long_name = 'Soil Wetness'
   units = ''
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
end if


call task_barrier()


if(nfields.ne.nfields1) then
  if(masterproc) print*,'write_fields2D: error in nfields!!'
  call task_abort()
end if

if(masterproc) then
     close(46)
     if(save2Dsep.and.dogzip2D) call systemf('gzip -f '//filename)
     print*, 'Appending 2D data. file:'//filename
endif


end
