
subroutine boundaries(flag)

use grid, only: dompi

        	
implicit none
integer flag

if(dompi) then
  call task_boundaries(flag)
else
  call periodic(flag)
end if

end subroutine boundaries
