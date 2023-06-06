module simple_ocean

!------------------------------------------------------------
! Purpose:
!
! A collection of routines used to specify fixed 
! or compute interactive SSTs, like slab-ocean model, etc.
!
! Author: Marat Khairoutdinov
! Based on dynamic ocean impelemntation from the UW version of SAM.
!
! Modified by Peter Blossey (January 2005)
! Implemented a greater variety of SST/ocean cooling profiles.
!
! Modified by Matthew Peters (Apr 2006) to read in evolving SST from
!  netcdf file.
!------------------------------------------------------------
use grid, only: nx,ny,nx_gl,dx,nsubdomains_x,rank,masterproc,dtn, &
     sst_target, delta_sst, ssty,&
     Szero, deltaS, mlohflx, hml, &
     nstep, firststep, ny_gl, dtfactor
use vars, only: latitude,rho,fluxbq,fluxbt
use simple_land              ! peters

implicit none

private 

real sfctflx(ny_gl),sfcqflx(ny_gl)
! Public methods:

public set_sst     ! set SST 
public sst_evolve ! evolve SST according to a model set by the ocean_type
public readmlohflx ! read slab ocean heating profile (fn of y) from file

public checkSSTinput      ! peters error check input parms in setparm
public doevolveSSTforcing ! peters public var

! peters.  stuff for evolving SST.  note:  evolveSST in simple_land
logical                   ::  doevolveSSTforcing  ! default=false


CONTAINS


SUBROUTINE set_sst()
 use grid, only: dodynamicocean,SFC_FLX_FXD,betafactor,SFC_TAU_MASK,sfctaumask
 use vars, only: sstxy, sstxy_dble, rho
 use params, only: tabs_s, ocean_type, lcond, cp

 real tmpx(nx), pii, lx
 integer i,j,it,jt

 !bloss: Move definition of x coordinate to top of routine.
 lx = float(nx_gl)*dx
 do i = 1,nx
    tmpx(i) = float(mod(rank,nsubdomains_x)*nx+i-1)*dx
 end do

 !bloss: Define parameters for initial sst if using dynamic ocean.
 if (dodynamicocean) then
    tabs_s = sst_target
    delta_sst = deltaS/35.
 end if

 select case (ocean_type)

   case(0) ! Uniform SST/Uniform Cooling
           !   Fixed SST: SST set to value from sfc file.
           !   Dynamic ocean: Initial SST set to sst_target (from prm file),
           !                  Ocean cooling uniform, set to Szero (prm file).

      sstxy = tabs_s

   case(1) !  Mock walker circulation
           !    Fixed SST: sinusoidal SST profile.
           !    Dynamic ocean: Initial SST profile sinusoidal,
           !                   Ocean cooling piecewise linear w/min in center.

     pii = atan2(0.d0,-1.d0)
     do j=1,ny
       do i=1,nx
         sstxy(i,j) = sst_target-delta_sst*cos(2.*pii*tmpx(i)/lx)
       end do
     end do

!================================================================
! UW ADDITIONS

   case(2) !  ENSO-like SST perturbation
           !    Fixed SST: sinusoidal SST profile.
           !    Dynamic ocean: Initial SST profile sinusoidal,
           !                   Ocean cooling piecewise linear w/min in center.

     do j=1,ny
       do i=1,nx
         sstxy(i,j) = tabs_s &
              - delta_sst*max(0.,1.-3.*tmpx(i)/lx,3.*tmpx(i)/lx-2.)
       end do
     end do

   case(3) !  beta plane simulation -- zonally-uniform SST
           !    Fixed SST: latitudinally-varying SST (read from file)
           !    Dynamic ocean: Initial SST profile as above,
           !                   latitudinally-varying ocean cooling profile
           !                    (also read in from file).
     
     call task_rank_to_index(rank,it,jt) 

     ! READ SST PROFILE FROM ssty FILE
     call readssty
     do j=1,ny
        sstxy(:,j) = ssty(j+jt)
     enddo
     !SET MASK SO THAT THERE IS NO SURFACE DRAG WITHIN 1 ROSSBY RADIUS
     if(SFC_TAU_MASK) then
        do j=1,ny
           sfctaumask(j) = max(min((abs(latitude(1,j))-6.)/6,1.),0.)
           write(*,*) 'lati=',latitude(1,j),'mask=', sfctaumask(j)
        enddo
     endif
     if (SFC_FLX_FXD.and.firststep) then
        call readsfcflx
        do j=1,ny
           fluxbq(:,j)=sfcqflx(j+jt)/(rho(1)*lcond)*sqrt(betafactor)
           fluxbt(:,j)=sfctflx(j+jt)/(rho(1)*cp)*sqrt(betafactor)
        enddo
     endif

     
  case(4) !  beta plane simulation with a cold tongue,

     call task_rank_to_index(rank,it,jt) 
     ! READ SST PROFILE FROM ssty FILE
     call readssty
     do j=1,ny
        sstxy(:,j) = ssty(j+jt)&
             - delta_sst*exp(-(latitude(i,j)/2.)**2.) &
             *max(exp(-(tmpx(i)*8./lx)**2), &
             exp(-(tmpx(i)*8./lx-8.)**2))
     enddo

  case(5) !  beta dynamic ocean with a warm pool,

     call task_rank_to_index(rank,it,jt) 
     ! READ SST PROFILE FROM ssty FILE
     call readssty
     do j=1,ny
        sstxy(:,j) = ssty(j+jt) &
             + delta_sst*exp(-(latitude(i,j)/18.)**2.) &
             *exp(-((tmpx(i)-lx/2.)/0.3/lx)**2.)
     enddo

  case(6) !  beta dynamic ocean with a sin wave temperature perturbation,

     call task_rank_to_index(rank,it,jt) 
     ! READ SST PROFILE FROM ssty FILE
     call readssty
     pii = atan2(0.d0,-1.d0)
     do j=1,ny
        sstxy(:,j) = ssty(j+jt) &
             + delta_sst*exp(-(latitude(i,j)/18.)**2.) &
             *sin(2.*pii*tmpx(i)/lx)
     enddo

! END UW ADDITIONS
!================================================================

   case default

     if(masterproc) then
         print*, 'unknown ocean type in set_sst. Exitting...'
         call task_abort
     end if

 end select

end subroutine set_sst



SUBROUTINE sst_evolve
 use grid, only: dofixedoceancooling, dostatis, dompi, resetsst
 use vars, only: sstxy, fluxbt, fluxbq, rhow,qocean_xy, qoceanxy, sstxy_dble, lhobs
 use params, only: rhor, cp, lcond, tabs_s
 use rad, only: swnsxy, lwnsxy

 real, parameter :: cpr = 4000.   ! Liquid Water Cp = 4000 J/kg/K
 real factor_cp, factor_lc, tmp
 integer i,j
 
      if(firststep) then !bloss
         call set_qocean() ! initialize qoceanxy
         if(resetsst) call set_sst() ! reset profile of sst if requested
         sstxy_dble = sstxy(1:nx,1:ny) ! initialize double precision sst.
      end if

      if(.not.dofixedoceancooling) call set_qocean_control()

      ! Define weight factors for the mixed layer heating due to
      ! the model's sensible and latent heat flux.
      factor_cp = rhow(1)*cp
      factor_lc = rhow(1)*lcond

      ! Use forward Euler to integrate the differential equation
      ! for the ocean mixed layer temperature: dT/dt = S - E.
      ! The source: CPT?GCSS WG4 idealized Walker-circulation 
      ! RCE Intercomparison proposed by C. Bretherton.
      do j=1,ny
         do i=1,nx
	   qocean_xy(i,j) = qocean_xy(i,j) + qoceanxy(i,j)*dtfactor

            sstxy_dble(i,j) = sstxy_dble(i,j) &
                 + dtn*(swnsxy(i,j)          & ! SW Radiative Heating
                 - lwnsxy(i,j)               & ! LW Radiative Heating
                 - factor_cp*fluxbt(i,j)     & ! Sensible Heat Flux
                 - factor_lc*fluxbq(i,j)     & ! Latent Heat Flux
                 + qoceanxy(i,j))            & ! Ocean Heating
                 /(rhor*cpr*hml)        ! Convert W/m^2 Heating to K/s

            sstxy(i,j) = sstxy_dble(i,j)
         end do
      end do

      if (dostatis) then
         lhobs = lhobs + SUM(qoceanxy)/float(nx*ny)
         tabs_s = SUM(sstxy(1:nx,1:ny))/float(nx_gl*ny_gl)
         if (dompi) then
            call task_sum_real(tabs_s,tmp,1)
            tabs_s = tmp
         end if
      end if

end subroutine sst_evolve

!================================================================
! UW ADDITIONS

SUBROUTINE set_qocean()
  use params, only: ocean_type
  use vars, only: qoceanxy, betafactor

  real tmpx(nx), pii, lx
  integer i,j,it,jt

  !bloss: Move definition of x coordinate to top of routine.
  lx = float(nx_gl)*dx
  do i = 1,nx
     tmpx(i) = float(mod(rank,nsubdomains_x)*nx+i-1)*dx
  end do

  select case (ocean_type)

  case(0) ! Uniform SST/Uniform Cooling
     !   Fixed SST: SST set to value from sfc file.
     !   Dynamic ocean: Initial SST set to sst_target (from prm file),
     !                  Ocean cooling uniform, set to Szero (prm file).

     qoceanxy = Szero

  case(1) !  Mock walker circulation
     !    Fixed SST: sinusoidal SST profile.
     !    Dynamic ocean: Initial SST profile sinusoidal,
     !                   Ocean cooling piecewise linear w/min in center.

     do j=1,ny
        do i=1,nx
           qoceanxy(i,j) = Szero + deltaS*(0.5 - abs(2.*tmpx(i)/lx - 1))
        end do
     end do

  case(2) !  ENSO-like SST perturbation
     !    Fixed SST: sinusoidal SST profile.
     !    Dynamic ocean: Initial SST profile sinusoidal,
     !                   Ocean cooling piecewise linear w/min in center.

     do j=1,ny
        do i=1,nx
           qoceanxy(i,j) = Szero &
                + deltaS*(1. - max(0.,3.-9.*tmpx(i)/lx,9.*tmpx(i)/lx-6.))
        end do
     end do

  case(3) !  beta plane simulation -- zonally-uniform SST
     !    Fixed SST: latitudinally-varying SST (read from file)
     !    Dynamic ocean: Initial SST profile as above,
     !                   latitudinally-varying ocean cooling profile
     !                    (also read in from file).

     call task_rank_to_index(rank,it,jt) 

     ! READ OCEAN COOLING PROFILE FROM mlohflx FILE
     if (firststep.and.(nstep.eq.1)) call readmlohflx
     do j=1,ny
        qoceanxy(:,j)=mlohflx(j+jt)
     enddo

  case(4) !  beta plane simulation with a cold tongue,

     call task_rank_to_index(rank,it,jt) 
     ! READ OCEAN COOLING PROFILE FROM mlohflx FILE
     if (firststep.and.(nstep.eq.1)) call readmlohflx

     ! SET OCEAN COOLING PROFILE TO THAT READ IN FROM FILE AND
     ! SUBTRACT GAUSSIAN HUMP ON EQUATOR TO INDUCE COLD POOL
     do j = 1,ny
        do i = 1,nx
           qoceanxy(i,j) = mlohflx(j+jt) &
                - deltaS*exp(-(latitude(i,j)/2.)**2.) &
                *max(exp(-(tmpx(i)*8./lx)**2), &
                exp(-(tmpx(i)*8./lx-8.)**2))
        end do
     end do

  case(5) !  beta dynamic ocean with a warm pool,

     call task_rank_to_index(rank,it,jt) 
     ! READ OCEAN COOLING PROFILE FROM mlohflx FILE
     if (firststep.and.(nstep.eq.1)) call readmlohflx

     do j = 1,ny
        do i = 1,nx
           qoceanxy(i,j) = mlohflx(j+jt) &
                + deltaS*exp(-(latitude(i,j)/18.)**2.) &
                *exp(-((tmpx(i)-lx/2.)/0.3/lx)**2.)
        end do
     end do

  case(6) !  beta dynamic ocean with a warm pool,

     call task_rank_to_index(rank,it,jt) 
     ! READ OCEAN COOLING PROFILE FROM mlohflx FILE
     if (firststep.and.(nstep.eq.1)) call readmlohflx
     pii = atan2(0.d0,-1.d0)
     do j = 1,ny
        do i = 1,nx
           qoceanxy(i,j) = sqrt(betafactor)*(mlohflx(j+jt) &
                + deltaS*exp(-(latitude(i,j)/18.)**2.) &
                *sin(2.*pii*tmpx(i)/lx))
        end do
     end do

  case default

     if(masterproc) then
        print*, 'unknown ocean type in set_qocean. Exitting...'
        call task_abort
     end if

  end select

end SUBROUTINE set_qocean

SUBROUTINE set_qocean_control()
  use grid, only: dompi
 use vars, only: sstxy, fluxbt, fluxbq, rhow,qocean_xy, qoceanxy
 use params, only: rhor, cp, lcond, ocean_type
 use rad, only: swnsxy, lwnsxy

  double precision tmp_sum(ny_gl), tmp_sum2(ny_gl)
  real, parameter :: cpr = 4000.   ! Liquid Water Cp = 4000 J/kg/K
  real factor_cp, factor_lc
  real tmpx(nx), pii, lx
  integer i,j,it,jt

  ! Define weight factors for the mixed layer heating due to
  ! the model's sensible and latent heat flux.  
  factor_cp = rhow(1)*cp
  factor_lc = rhow(1)*lcond

  !bloss: Move definition of x coordinate to top of routine.
  lx = float(nx_gl)*dx
  do i = 1,nx
     tmpx(i) = float(mod(rank,nsubdomains_x)*nx+i-1)*dx
  end do

  tmp_sum = 0.

  select case (ocean_type)

  case(0:2) ! Modify Szero to maintain domain-average sst equal to sst_target

     do j = 1,ny
        do i = 1,nx
           tmp_sum(1) = tmp_sum(1) &
                - (swnsxy(i,j)              & ! SW Radiative Heating
                - lwnsxy(i,j)               & ! LW Radiative Heating
                - factor_cp*fluxbt(i,j)     & ! Sensible Heat Flux
                - factor_lc*fluxbq(i,j)     & ! Latent Heat Flux
                + qoceanxy(i,j))              ! Ocean Heating
        end do
     end do
     tmp_sum(1) = tmp_sum(1)/dble(nx_gl*ny_gl)
     if (dompi) then
        call task_sum_real8(tmp_sum,tmp_sum2,1)
        tmp_sum(1) = tmp_sum2(1)
     end if

     qoceanxy = qoceanxy + tmp_sum(1)

  case(3:5) ! Modify qoceanxy to maintain zonally-averaged sst
     ! profile equal to ssty (read in from file)

     call task_rank_to_index(rank,it,jt) 

     do j = 1,ny
        do i = 1,nx
           tmp_sum(j+jt) = tmp_sum(j+jt) &
                - (swnsxy(i,j)              & ! SW Radiative Heating
                - lwnsxy(i,j)               & ! LW Radiative Heating
                - factor_cp*fluxbt(i,j)     & ! Sensible Heat Flux
                - factor_lc*fluxbq(i,j)     & ! Latent Heat Flux
                + qoceanxy(i,j))              ! Ocean Heating
        end do
     end do
     tmp_sum = tmp_sum/dble(nx_gl)
     if (dompi) then
        call task_sum_real8(tmp_sum,tmp_sum2,ny_gl)
        tmp_sum = tmp_sum2
     end if

     do j = 1,ny
        do i = 1,nx
           qoceanxy(i,j) = qoceanxy(i,j) + tmp_sum(j+jt)
        end do
     end do

  case default

     if(masterproc) then
        print*, 'unknown ocean type in set_qocean_control. Exitting...'
        call task_abort
     end if

  end select

end SUBROUTINE set_qocean_control

subroutine readssty
  use vars
  integer ioerror, j

  open(8,file='./'//trim(case)//'/ssty',status='old', &
       form='formatted', iostat=ioerror)
  if (ioerror.eq.0) then
     do j=1,ny_gl      
        read(8,*) ssty(j)
     end do
     close (8)
  else
     write(*,*) 'FAILED TO FIND ssty FILE'
     write(*,*) 'UNABLE TO INITIALIZE SST PROFILE'
     call task_abort()
  end if

end subroutine readssty

subroutine readmlohflx
  use vars
  integer ioerror, j

  open(8,file='./'//trim(case)//'/mlohflx',status='old', &
       form='formatted', iostat=ioerror)
  if (ioerror.eq.0) then
     do j=1,ny_gl      
        read(8,*) mlohflx(j)
     end do
     close (8)
  else
     write(*,*) 'FAILED TO FIND mlohflx FILE'
     write(*,*) 'UNABLE TO INITIALIZE OCEAN COOLING PROFILE'
     call task_abort()
  end if

end subroutine readmlohflx

subroutine readsfcflx
  use vars
  integer ioerror, j

  open(8,file='./'//trim(case)//'/sfcflx',status='old', &
       form='formatted', iostat=ioerror)
  if (ioerror.eq.0) then
     do j=1,ny_gl      
        read(8,*) sfcqflx(j),sfctflx(j)
     end do
     close (8)
  else
     write(*,*) 'FAILED TO FIND mlohflx FILE'
     write(*,*) 'UNABLE TO INITIALIZE OCEAN COOLING PROFILE'
     call task_abort()
  end if

end subroutine readsfcflx

!  check input parameters for SST to make sure don't input incompatable ones
subroutine checkSSTinput()

  use vars
  use params
  use simple_land
  implicit none

  ! first check SSTevolve stuff
  if(doevolveSSTforcing .and. dosfcforcing) then
     print*,'doevolveSSTforcing and doscfcorcing both set.  Evolving'
     print*,'SST will not affect simulation'
     call task_abort()
  end if
  if(doevolveSSTforcing .and. dodynamicocean) then
     print*,'doevolveSSTforcing and dodynamicocean both set'
     print*,'Cannot have both mixed layer and fixed SST!!'
     call task_abort()
  end if
  if(doevolveSSTforcing .and. SFC_FLX_FXD) then
     print*,'doevolveSSTforcing and SFC_FLX_FXD both set.'
     call task_abort()
  end if

  ! check the SST parms
  if(doevolveSSTforcing) call checkevolveparms(evolveSST,'evolveSST')

  ! set resetsst = false if doevolveSSTforcing is set
  if(doevolveSSTforcing) resetsst = .false.

  ! set ocean_type = 1 if evolving SST
  ! this is important for the LW radiative scheme.  otherwise, rad_full
  ! will use the value sst(1,1) for the entire subdomain!
  if(doevolveSSTforcing) ocean_type = 1


end subroutine checkSSTinput


 ! END UW ADDITIONS
 !================================================================

end module simple_ocean
