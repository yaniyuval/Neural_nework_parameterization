module rad

use grid

implicit none

!--------------------------------------------------------------------
!
! Variables accumulated between two calls of radiation routines


	real lwnsxy  (nx, ny)
	real swnsxy  (nx, ny)
	real lwntxy  (nx, ny)
	real swntxy  (nx, ny)
	real lwnscxy  (nx, ny)
	real swnscxy  (nx, ny)
	real lwntcxy  (nx, ny)
	real swntcxy  (nx, ny)
	real lwdsxy  (nx, ny)
	real swdsxy  (nx, ny)
	real solinxy  (nx, ny)

end module rad
