
subroutine setperturb

!  Random noise

use vars

implicit none
	
integer i,j,k
real rrr,ranf_

call ranset_(30*rank)

do k=1,nzm
 do j=1,ny
  do i=1,nx
    rrr=1.-2.*ranf_()
!    rrr = 0.0

    if(k.le.5) then
      t(i,j,k)=t(i,j,k)+0.02*rrr*(6-k)
    endif

    !bloss  add noise in regions where cloud is present in initial condition
    if(qc0(k)+qi0(k).gt.0.) then
      t(i,j,k)=t(i,j,k)+0.1*rrr
    endif

    if(k.le.4.and..not.dosmagor) then
      tke(i,j,k)=tke(i,j,k)+0.04*(5-k)
    endif

    if(doscalar) then
      tke(i,j,k) = q(i,j,k)
    end if
  end do
 end do
end do



end
