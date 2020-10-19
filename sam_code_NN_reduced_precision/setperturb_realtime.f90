
subroutine setperturb_realtime

!  Random noise added to t and q below 1km every 15 minutes

use vars

implicit none
	
integer i,j,k
real rrr,ranf_

if (mod(nstep*floor(dt*sqrt(betafactor)+0.0001),900).eq.0) then
do k=1,nzm
 do j=1,ny
  do i=1,nx
    if(z(k).le.1000) then
       rrr=1.-2.*ranf_()
       t(i,j,k)=t(i,j,k)+0.1*rrr
       rrr=1.-2.*ranf_()
       q(i,j,k)=q(i,j,k)+1.e-4*rrr
    endif
  end do
 end do
end do
endif



end
