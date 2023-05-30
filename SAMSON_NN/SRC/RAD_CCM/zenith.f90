      subroutine zenith(calday  ,lat  ,lon   ,coszrs  ,eccf  )
      implicit none	
!-----------------------------------------------------------------------
!
! Compute cosine of solar zenith angle for albedo and radiation 
!   computations.
!
!
! Input arguments
!
      real calday              ! Calendar day, including fraction
      real lat                ! Current centered latitude (degrees)
      real lon          ! Centered longitude (degrees)
!
! Output arguments
!
      real coszrs   ! Cosine solar zenith angle
      real eccf     ! Earth orbit eccentricity factor
!
!---------------------------Local variables----------------------------
!
      real sin_lat, cos_lat, sun_azim, sun_lat, pi	

      pi = 3.141593	
      eccf = 1.+0.034*cos(2.*pi*calday/365.)
      sin_lat = sin(pi/180.*lat)
      cos_lat = cos(pi/180.*lat)
      sun_lat = (pi*23.5/180.)*cos((calday - 172.)*pi/183.)
      sun_azim = 2.*pi*(calday + lon/360.)
      coszrs = -cos(sun_lat)*cos(sun_azim)*cos_lat+sin(sun_lat)*sin_lat
      return
      end


real function perpetual_factor(day, lat, lon)
use grid, ONLY: dt, nrad
implicit none

!  estimate the factor to multiply the solar constant
!  so that the sun hanging perpetually right above
!  the head (zenith angle=0) would produce the same
!  total input the TOA as the sun subgect to diurnal cycle.
!  coded by Marat Khairoutdinov, 2004
!
! Input arguments
!
real day             ! Calendar day, without fraction
real lat                ! Current centered latitude (degrees)
real lon          ! Centered longitude (degrees)

! Local:
real :: dtime 
real :: coszrs, eccf, time

time = day
dtime = dt*float(nrad)/86400.
perpetual_factor = 0.

do while (time.lt.day+1.)
  call zenith(time, lat, lon, coszrs, eccf)
  perpetual_factor = perpetual_factor &
       + min(dtime,day+1.-time)*max(0.,eccf*coszrs)
  time = time+dtime
end do

end


real function perpetual_equinox(lat, lon)
use grid, ONLY: dt, nrad
implicit none

!  estimate the factor to multiply the solar constant
!  so that the sun hanging perpetually right above
!  the head (zenith angle=0) would produce the same
!  total input the TOA as the sun subgect to diurnal cycle.
!  coded by Marat Khairoutdinov, 2004
!
! Input arguments
!
real lat                ! Current centered latitude (degrees)
real lon          ! Centered longitude (degrees)

! Local:
real :: dtime 
real :: coszrs, eccf, time, pii

pii = 3.141593	
time = 0.
dtime = dt*float(nrad)/86400.
perpetual_equinox = 0.

do while (time.lt.1.)
   eccf = 1.
   coszrs = -cos(2*pii*time)*cos(pii*lat/180.)
   perpetual_equinox = perpetual_equinox &
        + min(dtime,1.-time)*max(0.,eccf*coszrs)
   time = time+dtime
end do
end

