      subroutine radinp(pmid    ,pint    ,h2ommr  ,cld     ,
     $                  pmidrd  ,pintrd  ,plco2   ,plh2o   ,tclrsf )
C-----------------------------------------------------------------------
C
C Modified for generalized orbit
C 19 November 1996    Bruce P. Briegleb
C
C Set latitude and time dependent arrays for input to solar
C and longwave radiation.
C
C Convert model pressures to cgs, compute path length arrays needed for the 
C longwave radiation, and compute ozone mixing ratio, needed for the solar
C radiation.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C Reviewed:          J. Kiehl, April 1996
C Reviewed:          B. Briegleb, May 1996
C Bug fixed:         B. Briegleb, September 1996
C
C-----------------------------------------------------------------------
c
c $Id: radinp.f,v 1.1.1.1 2007/01/11 22:58:18 cw Exp $
c
      use domain, only: nz_gl
	implicit none
C------------------------------Parameters-------------------------------
	include 'prgrid.h'
C-----------------------------------------------------------------------
	include 'crdcon.h'
C-----------------------------------------------------------------------
	include 'comvmr.h'
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real pmid(plond,plev)    ! Pressure at model mid-levels (pascals)
      real pint(plond,plevp)   ! Pressure at model interfaces (pascals)
      real h2ommr(plond,plev)  ! H2o mass mixing ratio
      real cld(plond,plevp)    ! Fractional cloud cover
C
C Output arguments
C
      real pmidrd(plond,plev)  ! Pressure at mid-levels (dynes/cm*2)
      real pintrd(plond,plevp) ! Pressure at interfaces (dynes/cm*2)
      real plco2(plond,plevp)  ! Vert. pth lngth of co2 (prs-weighted)
      real plh2o(plond,plevp)  ! Vert. pth lngth h2o vap.(prs-weighted)
      real tclrsf(plond,plevp) ! Product of clr-sky fractions from top
C                                  of atmosphere to level.
C---------------------------Local variables-----------------------------
C
      integer i                ! Longitude loop index
      integer k                ! Vertical loop index

      real p0                  ! Standard pressure (dynes/cm**2)
      real amd                 ! Effective molecular weight of dry air (g/mol)
      real amo                 ! Molecular weight of ozone (g/mol)
      real amco2               ! Molecular weight of co2   (g/mol)
      real cpwpl               ! Const in co2 mix ratio to path length conversn
      real vmmr                ! Ozone volume mixing ratio

      save     p0   ,amd   ,amo  ,amco2

      data p0    /  1.01325e6 /
      data amd   /  28.9644   /
      data amo   /  48.0000   /
      data amco2 /  44.0000   /

c bloss: moved definition of co2vmr to radini
c	co2vmr = 3.55e-4

C
C Convert pressure from pascals to dynes/cm2
C
      do k=1,plev
        do i=1,plon
          pmidrd(i,k) = pmid(i,k)*10.0
          pintrd(i,k) = pint(i,k)*10.0
        end do
      end do
      do i=1,plon
        pintrd(i,plevp) = pint(i,plevp)*10.0
      end do
C
C Compute path quantities used in the longwave radiation:
C
      vmmr  = amco2/amd
      cpwpl = vmmr*0.5/(gravit*p0)
      do i=1,plon
        plh2o(i,1)  = rgsslp*h2ommr(i,1)*pintrd(i,1)*pintrd(i,1)
        plco2(i,1)  = co2vmr*cpwpl*pintrd(i,1)*pintrd(i,1)
        tclrsf(i,1) = 1.
      end do
      do k=1,plev
        do i=1,plon

          plh2o(i,k+1)  = plh2o(i,k) + rgsslp*
     $             (pintrd(i,k+1)**2 - pintrd(i,k)**2)*h2ommr(i,k)
          plco2(i,k+1)  = co2vmr*cpwpl*pintrd(i,k+1)**2
          tclrsf(i,k+1) = tclrsf(i,k)*(1.-cld(i,k+1))
        end do
      end do

C
      return
      end





 
