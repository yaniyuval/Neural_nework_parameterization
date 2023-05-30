
subroutine zero
	
use vars

implicit none
	
integer k
	
dudt(:,:,:,na) = 0.
dvdt(:,:,:,na) = 0.
dwdt(:,:,:,na) = 0.
misc(:,:,:) = 0.

total_water_before = 0.
total_water_after = 0.
total_water_evap = 0.
total_water_prec = 0.
total_water_ls = 0.
do k=1,nzm
 total_water_before = total_water_before + &
            (sum(q(1:nx,1:ny,k))+sum(qp(1:nx,1:ny,k)))*adz(k)*dz *rho(k)
end do

end
