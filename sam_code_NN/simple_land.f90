
!  peters module simple_land
!
!  A collection of subroutines for implementing a simple land scheme.
!
!  These read in an evolving realistic distribution of evolving soil
!  wetness, a land/sea mask and a horizontally varying distribution of 
!  surface roughness.
!
!  This module also contains some utility subroutines needed for extentions
!  made to simple_ocean.
!
!  Coded by Matthew Peters (Apr-May 2006).

module simple_land

use vars
implicit none

private

! public subroutines
public set_landmask            ! set the land/sea mask
public set_surfrough           ! set the xy varying roughness
public checkLANDinput      ! check the input params; called from setparm

! public vars
public :: dolandseamask, land_filename, land_ncfldname, &
     doxysurf_rough, xysurf_rough, surfrough_filename, surfrough_ncfldname, &
     doxysoil_wet, xysoil_wet, evolvesoilwet, &
     dointeractiveland, tau_landdamping, hml_land, &
     evolveSST, maxsurfrough, nstart_surfrough, landmask_frac, sstxyequil

! public utility subroutines needed by simple_ocean.
public :: getfld, print_nc_error, set_evolveforcing, setdefault_evolveparms,&
          checkevolveparms, set_sstevolveforcing, write_fields_invariant, &
          initialize_simple_land
public evolveparms


!------------------------------------------
! type declarations
!------------------
type evolveparms

  ! parameters set in namelist
  character(400) :: filename       ! full path needed!
  character(400) :: ncfldname       ! field name in nc file
  integer        ::  evolvetype       ! default=1
  integer        ::  nstepevolve ! read in field this many steps
                                  ! need if evolvetype=2 or 3
  integer        ::  nevolve     ! number of flds to cycle through;
                                  ! need if evolvetype=2
  integer        ::  nctimestepin

  ! internal variables used by SAM
  integer        ::  nctimestep 
  integer        ::  ncid
  real           ::  readflds(3,nx,ny)   ! fields read in from file to interp

end type evolveparms
!------------------------------------------


!-------------------------------------------
!  VARIABLES
!-----------

! land/sea mask stuff
logical                   ::  dolandseamask       ! default=false
real           ::  landmask_frac(0:nxp1,1-YES3D:nyp1)! frac. ional land (0-1)
character(160) ::  land_filename       ! full path needed!
character(160) ::  land_ncfldname

! parameters atmos/land/sea interactions
logical   dointeractiveland  ! compute land temp interactively? default=false
real tau_landdamping ! damping time for land to equil. temp; default = -1
        ! if negative, then don't do land damping (default=no damping)
real hml_land        ! mixed layer depth for land (detault=0.1 m)

! surface roughness stuff
logical           ::  doxysurf_rough   ! default=.false.
real              ::  xysurf_rough(0:nxp1,1-YES3D:nyp1) ! [m]surface rough
character(160) ::  surfrough_filename       ! full path needed!
character(160) ::  surfrough_ncfldname
real  maxsurfrough        ! [m] maximum allowable value of surface roughness
integer nstart_surfrough  ! gradually increase surface roughness over this
                          ! many time steps to maximum value

! soil wetness stuff
logical                   ::  doxysoil_wet      ! default=.false.
type(evolveparms) evolvesoilwet
real                   ::  xysoil_wet(0:nxp1,1-YES3D:nyp1) ! surface wetness


! evolving SST variable, needed here because it's used in set_sstevolveforcing
type(evolveparms) evolveSST
real sstxyequil(nx,ny)




CONTAINS
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------


!-------------------------------------------------------------------
! sets evolving surface temp forcing 
!  modelled after sst_evolve; modified for land/sea mask and to use
!  presribed evolving SST
subroutine set_sstevolveforcing()

  use vars
  use params
  use rad, only : swnsxy, lwnsxy
  implicit none

  real, parameter :: rhocp = 2.e6   ! rho*cp for soil, ref: Hartmann p 85
  integer i,j
  real coef

  ! interpolate to get the new equilibrum temp if icycle.eq.1
  if(icycle.eq.1) then
    call set_evolveforcing(evolveSST,sstxyequil)
  end if

  ! initialize sstxy_dble.  sstxy_dble in this implementation the unsmoothed
  ! version of the SST.  it is smoothed below to give the sstxy that SAM
  ! uses for surface/radiative calculations
  if(firststep.and.nrestart.eq.0) sstxy_dble(1:nx,1:ny) = sstxyequil

  ! now compute the sst from surface energy balance if needed
  if(dointeractiveland) then
    do j=1,ny
      do i=1,nx
         ! Use forward Euler to integrate the differential equation
         ! for the ocean mixed layer temperature: dT/dt = S - E.
         ! The source: CPT?GCSS WG4 idealized Walker-circulation
         ! RCE Intercomparison proposed by C. Bretherton.
         sstxy_dble(i,j) = sstxy_dble(i,j) &
               + dtn*(swnsxy(i,j)          & ! SW Radiative Heating
               - lwnsxy(i,j)               & ! LW Radiative Heating
               - cp*rhow(1)*fluxbt(i,j)     & ! Sensible Heat Flux
               - lcond*rhow(1)*fluxbq(i,j))     & ! Latent Heat Flux
               /(rhocp*hml_land)        ! Convert W/m^2 Heating to K/s

         ! Add in nudging to specified temp dT/dt = (T-Te)/tau_landdamping
         if(tau_landdamping.gt.0) then
           sstxy_dble(i,j) = sstxy_dble(i,j) &
              - dtn*(sstxy_dble(i,j) - sstxyequil(i,j))/tau_landdamping
         end if

      end do
    end do

    ! compute final surface temp.  need to blend land temp with fixed ocean
    ! temp for total surface temp
    do j=1,ny
      do i=1,nx
         sstxy(i,j) = landmask_frac(i,j)*sstxy_dble(i,j) + &
                        (1.-landmask_frac(i,j))*sstxyequil(i,j)
         !sstxy(i,j) = sstxy_dble(i,j)
      end do
    end do

  else                        ! using prescribed surface temp
    sstxy(1:nx,1:ny) = sstxyequil
  end if       ! if(dointeractiveland)


end subroutine set_sstevolveforcing


!-------------------------------------------------------------------
!  subroutine to set evolving field.  use the type definition above
!
!  help below was originally written for SST.  everything is identical now,
!  except field is more general than sst.
!  ASSUMES:  sst is in a netcdf file input through the namelist
!            into sst_filename.  The actual field in the netcdf file
!            has name equal to sst_ncfldname and is dimensioned as
!            sst_ncfldname(time,ny_gl,nx_gl)
!
!            The x,y dimensions increase as the dimension increase (ie,
!            x(1) is the eastern extent of the domain and y(1) the southern
!
!            evolveSSTtype determines the nature of the evolving SST
!            if evolveSSTtype == 1 then fixes SST to be the field in the first
!              position in the netcdf file for entire run (time==1)
!            if evolveSSTtype == 2 then cycles through the first nevolveSST
!              in the netcdf file every nstepsevolveSST steps (ie,
!              nevolveSST=12 and nstepsevolveSST=30 days will read in
!              monthly climo)
!            if evolveSSTtype == 3 then will read in SST flds every
!             nstepsevolveSST until the end of the run (eg, force w/ 23 year
!             climo)
!
!
!  Interpolation follows the following scheme:  input in the netcdf file
!  is assumed to be the mean value over the time period of interest,
!  ie, averaging over every time step in january of SST will give value in
!  netcdf for January.  Interpolation is then done piecewise from values at
!  Jan 1 to Jan 15 and from Jan 15 to Feb 1, where the value at Jan 15 is
!  chosen to give the specified mean value in the netcdf file.
!
!  The values at Jan 1, Feb 1 are fixed to be the mean value of the neighboring
!  months, ie, value at Jan 1 = 0.5*(Jan mean + Dec mean), ect.  The value
!  at day 15 is then determined by solving:
!   Let  b' = determined value at jan 15. b = specified jan mean, a=specified
!   Dec mean, c=specified march mean.
!  SST in jan is linear interpolation from (a+b)/2 to b' (from Jan 1->15) and
!  from b' to (b+c)/2 (Jan 15->Feb 1)
!
!  We then have
!  b = 0.5*( (a+b)/2+b')/2 + ((b+c)/2+b')/2) ==> b' = 1.5*b - 0.25*(a+c)
!
!  Below, a occupies readflds(1,:,:), b occupies readflds(2,:,:), ect.
subroutine set_evolveforcing(fldparms,fldxy)

  use vars, only : nx, ny, firststep, nrestart, nstep
  implicit none
  include 'netcdf.inc'

  type(evolveparms) fldparms           ! input.  parameters for evolving field
  real fldxy(nx,ny)                   ! output.  SAM variable to store in

  integer err, tmpint
  real timeinterp
  real sstbp(nx,ny)             ! in notation, sstbp = b'
  real slopeinterp(nx,ny)

  ! if the first time step then open the netcdf file.  need to do this
  !  even after restart
  if(firststep) then
    err = NF_OPEN(trim(fldparms%filename), NF_NOWRITE, fldparms%ncid)
    if(err .ne. NF_NOERR) call print_nc_error(fldparms%filename)
  end if

  ! on restart, need to set fldxy if evolvetype.eq.1 since it won't be set
  ! below
  if(firststep .and. nrestart.ge.1) then
    fldxy = fldparms%readflds(1,:,:)
  end if

  ! if the first time step, read in fld from file and setup readflds
  if(nstep.eq.1) then

    ! read in the fld for position 1
    call getfld(fldparms%ncid,fldparms%ncfldname,fldxy,fldparms%nctimestep)
    fldparms%readflds(1,:,:) = fldxy

    ! setup readflds if necessary
    if(fldparms%evolvetype .ge. 2) then
      fldparms%readflds(2,:,:) = fldxy
      call getfld(fldparms%ncid,fldparms%ncfldname,fldparms%readflds(3,:,:),&
                      fldparms%nctimestep+1)
      fldparms%nctimestep = fldparms%nctimestep+1
    end if
  end if

  !---------------------------------------------------------
  ! if not the first time step, do interpolation as necessary
  !  on restart, this will set sst=0 when calling from setdata, but will
  !  overwrite sst and readflds with reading in restart file
  if(fldparms%evolvetype.ge.2 .and. nstep.ge.2) then

    ! update the readflds if stepping forward to the next time
    if(mod(nstep,fldparms%nstepevolve).eq.1) then
      fldparms%readflds(1,:,:) = fldparms%readflds(2,:,:)
      fldparms%readflds(2,:,:) = fldparms%readflds(3,:,:)

      ! update nctimestep
      fldparms%nctimestep = fldparms%nctimestep+1
      if(fldparms%evolvetype.eq.2 .and. &
                  fldparms%nctimestep.gt.fldparms%nevolve) then
        fldparms%nctimestep = fldparms%nctimestep-fldparms%nevolve
      end if

      call getfld(fldparms%ncid,fldparms%ncfldname,&
                  fldparms%readflds(3,:,:),fldparms%nctimestep)

    end if

    ! do the interpolation.  Use piece wise linear from mean value at
    !  monthly interface (Feb 1, Mar 1) to some determined value at day 15
    !  such that the monthly mean is equal to prescribed value
    tmpint = mod(nstep,fldparms%nstepevolve)
    if(tmpint .eq. 0) then
      tmpint = tmpint+fldparms%nstepevolve !1<=tmpint<=nstepevolve
    end if

    timeinterp=2.*(float(tmpint)-1.)/float(fldparms%nstepevolve-1)-1.

    ! compute value at day 15
    sstbp = 1.5*fldparms%readflds(2,:,:)&
            - 0.25*(fldparms%readflds(1,:,:)+fldparms%readflds(3,:,:))
    ! compute the slope of the interpolation
    if(timeinterp.le.0) then
      slopeinterp = sstbp - &
              0.5*(fldparms%readflds(1,:,:)+fldparms%readflds(2,:,:))
    else
      slopeinterp = 0.5*(fldparms%readflds(2,:,:)+fldparms%readflds(3,:,:)) - &
              sstbp
    end if

    ! compute the interpolated SST
    fldxy = slopeinterp*timeinterp + sstbp

  end if


end subroutine set_evolveforcing

!-------------------------------------------------------------------
!  sets the default evolving parameters
subroutine setdefault_evolveparms(fldparms)

  implicit none

  type(evolveparms) fldparms

  fldparms%filename = ''
  fldparms%ncfldname = ''
  fldparms%evolvetype = -1
  fldparms%nstepevolve=-1
  fldparms%nevolve=-1
  fldparms%nctimestepin=1

end subroutine setdefault_evolveparms

!-------------------------------------------------------------------
!  initializes the simple land stuff to default values
subroutine initialize_simple_land()

  use vars
  use params
  implicit none

  if(nrestart.eq.0) then
    evolvesoilwet%readflds=0.
    evolveSST%readflds=0.
  end if

  ! set the matrix values of xysoil_wet and xysurf_rough to
  ! the values input in namelist so surface.f90 works
  !if doxysoil_wet and doxysurf_rough are set then these are overwritten
  xysoil_wet(:,:) = soil_wetness
  xysurf_rough(:,:) = z0

  ! initialize landmask
  if(OCEAN) then
    landmask_frac=0.
  else
    landmask_frac=1.
  end if

end subroutine initialize_simple_land


!-------------------------------------------------------------------
!  checks the input parameters
subroutine checkevolveparms(fldparms,fldname)

  implicit none

  type(evolveparms) fldparms
  character(*) fldname

  if(fldparms%filename.eq.'') then
    print*,'Error with filename for evolving field ',fldname
    call task_abort()
  end if
  if(fldparms%ncfldname.eq.'') then
    print*,'Error with ncfldname for evolving field ',fldname
    call task_abort()
  end if
  if(fldparms%evolvetype.ge.2.and.fldparms%nstepevolve.eq.-1) then
    print*,'Error with nstepevolve for evolving field ',fldname
    call task_abort()
  end if
  if(fldparms%evolvetype.eq.2.and.fldparms%nevolve.eq.-1) then
    print*,'Error with nevolve for evolving field ',fldname
    call task_abort()
  end if
  if(nrestart.eq.0) then
    fldparms%nctimestep = fldparms%nctimestepin
  end if

end subroutine checkevolveparms


!-----------------------------------------------------------------------
!  subroutine to read in surface roughness file
!
!  ASSUMES:  surface roughness is in a netcdf file input through the namelist
!            into surfrough_filename.  The actual field in the netcdf file
!            has name equal to surfrough_ncfldname and is dimensioned as
!            surfrough_ncfldname(ny_gl,nx_gl)
!
!            The x,y dimensions increase as the dimension increase (ie,
!            x(1) is the eastern extent of the domain and y(1) the southern
!
subroutine set_surfrough()

  use vars
  implicit none
  include 'netcdf.inc'

  integer ncid, err, i,j
  real maxs, maxs2

  ! open the netcdf file
  err = NF_OPEN(trim(surfrough_filename), NF_NOWRITE, ncid)
  if(err .ne. NF_NOERR) call print_nc_error(land_filename)

  ! read in the surface roughness
  call getfld(ncid,surfrough_ncfldname,xysurf_rough(1:nx,1:ny),0)

  err = NF_CLOSE(ncid)

  ! now make sure surface roughness isn't bigger than max specified value
  maxs = maxval(xysurf_rough(1:nx,1:ny))
  call task_max_real(maxs,maxs2,1)

  ! now maxs2 is max input surface roughness in entire domain.  limit
  ! this by the specified input value so maxsurfrough is the max allowable
  ! value
  maxsurfrough = min(maxs2,maxsurfrough)

  ! limit surface roughness to maxsurfrough
  do j=1,ny
    do i=1,nx
      xysurf_rough(i,j) = min(xysurf_rough(i,j),maxsurfrough)
    end do
  end do

  ! swap boundaries for surface roughness
  call boundaries(13)


end subroutine set_surfrough


!-----------------------------------------------------------------------
!  subroutine to read in the land-sea mask from a file
!  and initialize it
!
!  ASSUMES:  land/sea mask is in a netcdf file input through the namelist
!            into land_filename.  The actual field in the netcdf file
!            has name equal to land_ncfldname and is dimensioned as
!            land_ncfldname(ny_gl,nx_gl)
!
!            The x,y dimensions increase as the dimension increase (ie,
!            x(1) is the eastern extent of the domain and y(1) the southern
!
subroutine set_landmask()

  use vars
  implicit none
  include 'netcdf.inc'

  integer ncid, varid, err
  integer i, j, it, jt, start(2), countvara(2)


  ! open the netcdf file containing the land/sea mask
  err = NF_OPEN(trim(land_filename), NF_NOWRITE, ncid)
  if(err .ne. NF_NOERR) call print_nc_error(land_filename)

  ! read in the land/sea mask for this subdomain
  call task_rank_to_index(rank,it,jt)
  start(1) = it+1
  start(2) = jt+1
  countvara(1) = nx
  countvara(2) = ny

  err = NF_INQ_VARID(ncid, trim(land_ncfldname), varid)
   if(err .ne. NF_NOERR) call print_nc_error('Getting varid for'//trim(land_ncfldname))

  err=NF_GET_VARA_REAL(ncid, varid, start, countvara, landmask_frac(1:nx,1:ny))
   if(err .ne. NF_NOERR) call print_nc_error('Reading in landsea mask')

  ! close the netcdf file
  err = NF_CLOSE(ncid)

  ! swap boundaries for landmask_frac
  call boundaries(12)


end subroutine set_landmask


!-------------------------------------------------------------------
!  generic subroutine for reading in a field from netcdf file
subroutine getfld(ncid,fldname,fld,timenum)

  implicit none
  include 'netcdf.inc'

  ! input
  integer ncid                   ! open netcdf file
  character(*) fldname           ! name of field in netcdf file
  real fld(:,:)                    ! returned field
  integer timenum                ! time number to fetch (0=2D field)

  ! local
  integer err, varid, start(3), countvara(3)
  integer it,jt

  !-------------------------------------------------
  call task_rank_to_index(rank,it,jt)
  start(1) = it+1
  start(2) = jt+1
  countvara(1) = nx
  countvara(2) = ny

  if(timenum.ne.0) then
    start(3) = timenum
    countvara(3) = 1
  end if

  !  get variable id
  err = NF_INQ_VARID(ncid, trim(fldname), varid)
  if(err .ne. NF_NOERR) call print_nc_error('Getting varid for'//trim(fldname))

  ! fetch field from netcdf file
  if(timenum.ne.0) then
    err = NF_GET_VARA_REAL(ncid, varid, start, countvara, fld)
  else
    err = NF_GET_VARA_REAL(ncid, varid, start(1:2), countvara(1:2), fld)
  end if

  if(err .ne. NF_NOERR) then
     call print_nc_error('Reading in data for'//trim(fldname))
  end if

end subroutine getfld


!-----------------------------------------------------------------------
!  generic print message and die routine for netcdf stuff
subroutine print_nc_error(fname)

  implicit none
  character(*)  ::  fname

  print*,'Error with netcdf file: ',trim(fname)
  print*,'Quitting.'
  call task_abort()

end subroutine print_nc_error


!-----------------------------------------------------------------------
!  check input parameters for land to make sure don't input incompatable ones
subroutine checkLANDinput()

  use vars
  use params
  implicit none


  ! options for land-sea mask
  if(dolandseamask .and. (OCEAN.or.LAND) ) then
     print*,'Both dolandseamask and either OCEAN or LAND was set!'
     print*,'Quitting.'
     call task_abort()
  end if
  if(dolandseamask .and. LES) then
     print*,'dolandseamask and LES both set!'
     print*,'Quitting.'
     call task_abort()
  end if
  if(doxysurf_rough .and. (SFC_FLX_FXD.or.LES)) then
     print*,'doxysurf_rough requires SFC_FLX_FXD and LES=.false.!'
     print*,'Quitting.'
     call task_abort()
  end if
  if(doxysoil_wet .and. (SFC_FLX_FXD.or.LES)) then
     print*,'doxysoil_wet requires SFC_FLX_FXD and LES=.false.!'
     print*,'Quitting.'
     call task_abort()
  end if

  if(doxysoil_wet) call checkevolveparms(evolvesoilwet,'evolvesoilwet')

end subroutine checkLANDinput


!-------------------------------------------------------------------------
! writes the time invariant fields at beginning of run
!  writes the field in 2Dcom format so can use 2Dcom2nc to convert
subroutine write_fields_invariant()

  use vars
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
  logical old2Dfile

  nfields = 0
  nfields1 = 0
  if(dolandseamask) nfields = nfields+1        ! output land/sea mask
  if(doxysurf_rough) nfields = nfields+1       ! output surface rough

  if(masterproc) then

    write(rankchar,'(i4)') nsubdomains
    write(timechar,'(i10)') nstep
    do i=1,11-lenstr(timechar)-1
      timechar(i:i)='0'
    end do

    ! Make sure that the new run doesn't overwrite the file from the old run 
    if(save2Dbin) then
      filetype = '.2Dbin'
    else
      filetype = '.2Dcom'
    end if

    if(save2Dsep) then
       filename='./DATA3D/'//trim(case)//'_'//trim(caseid)//'_'// &
          rankchar(5-lenstr(rankchar):4)//'_INVARIANT_'//timechar(1:10)//filetype 
          open(46,file=filename,status='unknown',form='unformatted')
    else
       filename='./DATA3D/'//trim(case)//'_'//trim(caseid)//'_INVARIANT_'// &
           rankchar(5-lenstr(rankchar):4)//filetype

       !bloss: check for existence of old 2D file
       inquire(FILE=filename,EXIST=old2Dfile)
       if(old2Dfile) then
          filestatus='old'
          open(46,file=filename,status=filestatus,form='unformatted', &
                                                     position='append')
       else
          filestatus='new'
          open(46,file=filename,status=filestatus,form='unformatted')	
       end if
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
       write(c_time,'(f10.2)') day*sqrt(betafactor)
       write(46) long_name(1:32)
       write(46) c_time,c_dx,c_dy

     end if ! save2Dbin

  end if! masterproc

if(dolandseamask) then
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1) = landmask_frac(i,j)
     end do
   end do
  name='LANDMASK'
  long_name='Land-sea mask.  1=land, 0=ocean'
  units=''
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
end if

if(doxysurf_rough) then
   nfields1=nfields1+1
   tmp(:,:,1) = xysurf_rough(1:nx,1:ny)
   name = 'SURF_RGH'
   long_name = 'Surface roughness'
   units = 'm'
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
endif


end subroutine write_fields_invariant




end module simple_land


