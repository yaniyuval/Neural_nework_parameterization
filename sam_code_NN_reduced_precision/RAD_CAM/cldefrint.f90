      subroutine cldefrint(m,n,landfrac,tlayer,rel,rei,psurface, player,landm,icefrac,snowh)
!-----------------------------------------------------------------------
!
! interface for cldefr to work with isccp simulator calls and CAM3 radiation
!
!-----------------------------------------------------------------------
  use ppgrid
  use pkg_cldoptics, only: cldefr
  implicit none
!------------------------------Parameters-------------------------------

! Input arguments
!
  integer m,n
  real(r4) landfrac(pcols)
  real(r4) icefrac(pcols)       ! Ice fraction
  real(r4) psurface(pcols)      ! Surface pressure
  real(r4) tlayer(pcols,pver)   ! Temperature
  real(r4) player(pcols,pver)   ! Midpoint pressures
  real(r4) landm(pcols)         ! Land fraction
  real(r4) snowh(pcols)         ! snow depth, water equivalent (meters)

!
! Output arguments
!
  real rel(pcols,plev)      ! Liquid effective drop size (microns)
  real rei(pcols,plev)      ! Ice effective drop size (microns)

  call cldefr(m,n,landfrac,tlayer,rel,rei, & ! CAM3 interface
                psurface,player,landm,icefrac, snowh)

      return
      end
 
