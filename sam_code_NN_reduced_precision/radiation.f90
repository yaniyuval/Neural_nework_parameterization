
subroutine radiation()

!	Radiation interface

use grid
implicit none

	
if(doradsimple) then

!  A simple predefined radiation (longwave only)


    call rad_simple()
	 
else


! Call full radiation package:
 

    call rad_full()	
 
endif

end


