  
subroutine precip_init_krueger

! Initialize precipitation related stuff

use vars
use params

implicit none

real :: pie, gammafff

pie = atan2(0.,-1.)

gamr3 = gammafff(4.80)
gams3 = gammafff(4.25)
gamg3 = gammafff(4.50)

end subroutine precip_init_krueger


