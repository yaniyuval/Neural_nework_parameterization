      subroutine albedo(ocean, coszrs,asdir,aldir,asdif,aldif )
C-----------------------------------------------------------------------
C Computes surface albedos over ocean for Slab Ocean Model (SOM)

C and the surface (added by Marat Khairoutdinov)

C Two spectral surface albedos for direct (dir) and diffuse (dif)
C incident radiation are calculated. The spectral intervals are:
C   s (shortwave)  = 0.2-0.7 micro-meters
C   l (longwave)   = 0.7-5.0 micro-meters
C
C Uses knowledge of surface type to specify albedo, as follows:
C
C Ocean           Uses solar zenith angle to compute albedo for direct
C                 radiation; diffuse radiation values constant; albedo
C                 independent of spectral interval and other physical
C                 factors such as ocean surface wind speed.
C
C For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
C Approximation for Solar Radiation in the NCAR Community Climate Model,
C Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
C
C---------------------------Code history--------------------------------
C
C Original version:         CCM1
C Standardized:             J. Rosinski, June 1992
C Reviewed:                 J. Kiehl, B. Briegleb, August 1992
C Rewritten for ocean only: J. Rosinski, May 1994
C Modified for SOM:         B. Briegleb  April 1995
C Reviewed:                 B. Briegleb  March 1996
C Reviewed:                 J. Kiehl     April 1996
C Reviewed:                 B. Briegleb  May   1996
c Added Land albedo	    M. Kairoutdinov Sep. 1999
C
C-----------------------------------------------------------------------
C
C $Id: albedo.f,v 1.1.1.1 2007/01/11 22:58:18 cw Exp $
C $Author: cw $
C
      use domain, only: nz_gl
	implicit none
C------------------------------Parameters-------------------------------
	include 'prgrid.h'
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real coszrs(plond)    ! Cosine solar zenith angle
C
C Output arguments
C
      real asdir(plond)     ! Srf alb for direct rad   0.2-0.7 micro-ms
      real aldir(plond)     ! Srf alb for direct rad   0.7-5.0 micro-ms
      real asdif(plond)     ! Srf alb for diffuse rad  0.2-0.7 micro-ms
      real aldif(plond)     ! Srf alb for diffuse rad  0.7-5.0 micro-ms
      logical ocean	! logical switch for ocean
C
C---------------------------Local variables-----------------------------
C
      integer i             ! Longitude index
      real adif
      data adif    /  0.06 /
C
C
C Initialize all ocean surface albedos to zero
C
      do i=1,plon
          asdir(i) = 0.
          aldir(i) = 0.
          asdif(i) = 0.
          aldif(i) = 0.
      enddo

      if(ocean) then
C
C
C Ice-free ocean albedos function of solar zenith angle only, and
C independent of spectral interval:
C
      do i=1,plon
          if (coszrs(i).gt.0.) then
            aldir(i) = (.026/(coszrs(i)**1.7 + .065)) +
     $                 (.15*(coszrs(i) - 0.10)*
     $                      (coszrs(i) - 0.50)*
     $                      (coszrs(i) - 1.00)  )
            asdir(i) = aldir(i)
            aldif(i) = adif
            asdif(i) = adif
          endif
      enddo

      else ! land

 ! Albedos for land type I (Briegleb)

      do i=1,plon
          if (coszrs(i).gt.0.) then
            asdir(i) = 1.4 * 0.06 / ( 1. + 0.8 * coszrs(i))
            asdif(i) = 1.2 * 0.06
            aldir(i) = 1.4 * 0.24 / ( 1. + 0.8 * coszrs(i))
            aldif(i) = 1.2 * 0.24
          endif
      enddo

      end if
C   
      return
      end
 
