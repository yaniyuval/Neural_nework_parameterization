module rad

use grid

implicit none

!--------------------------------------------------------------------
!
! Variables accumulated between two calls of radiation routines


        real tabs_rad(nx, ny, nzm)       ! accumulated temperature
        real qc_rad  (nx, ny, nzm)       ! accumulated cloud water (g/g)
        real qi_rad  (nx, ny, nzm)       ! accumulated cloud ice (g/g)
        real qs_rad  (nx, ny, nzm)       ! accumulated snow (g/g)
        real qv_rad  (nx, ny, nzm)       ! accumulated water vapor (g/g)
	real qrad    (nx, ny, nzm)	 ! radiative heating(K/s)
	real lwnsxy  (nx, ny)
	real swnsxy  (nx, ny)
	real lwntxy  (nx, ny)
	real swntxy  (nx, ny)
	real swntmxy  (nx, ny)
	real lwnscxy  (nx, ny)
	real swnscxy  (nx, ny)
	real lwntcxy  (nx, ny)
	real swntcxy  (nx, ny)
	real lwdsxy  (nx, ny)
	real swdsxy  (nx, ny)
	real solinxy  (nx, ny)

        real lwutxy  (nx, ny)
        real lwutcxy  (nx, ny)

	logical initrad		! flag to initialize profiles of traces
	integer nradsteps	! curent number of steps done before
				!   calling radiation routines
	data initrad/.true./

!       Gas traces (mass mixing ratios):  
        
        real o3(nzm)            ! Ozone
        real n2o(nzm)           ! N2O
        real ch4(nzm)           ! CH4
        real cfc11(nzm)         ! CFC11
        real cfc12(nzm)         ! CFC12

	integer ndiv, nxdiv, nydiv
!	parameter (ndiv=4)
! pog mod
	parameter (ndiv=1)
	parameter (nxdiv=max(1,nx/ndiv))
	parameter (nydiv=max(1,ny/ndiv))

        real absnxt(nz, 4, nxdiv, nydiv) ! Nearest layer absorptivities
        real abstot(nz, nz,nxdiv, nydiv)! Non-adjacent layer absorptivites
        real emstot(nz, nxdiv, nydiv)! Total emissivity
	real p_factor(nx, ny) ! perpetual-sun factor

end module rad
