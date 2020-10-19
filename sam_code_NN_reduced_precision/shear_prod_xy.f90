
subroutine shear_prod_xy(def2)
	
use vars
implicit none
	
real def2(nx,ny,nzm)
	
real rdx, rdy, rdz
integer i,j,k,ib,ic,jb,jc

rdx=1./dx 
rdy=1./dy

do k=1,nzm

 rdz = 1./(dz*adz(k))

 do j=1,ny
   jb=j-YES3D
   jc=j+YES3D
   do i=1,nx
     ib=i-1
     ic=i+1
	 
      def2(i,j,k)= ( &
          ( (u(ic,j,k)-u(i,j,k))*rdx)**2 - &
          ( 2.*( (u(ic,j,k)-u(i,j,k))*rdx * (v(i,jc,k)-v(i,j,k))*rdy ) ) + &
          ( (v(i,jc,k)-v(i,j,k))*rdy)**2 ) + &
          ( (u(i,jc,k)-u(i,j,k))*rdy+(v(ic,j,k)-v(i,j,k))*rdx )**2

! spatial average
!!$        + 0.25 * ( &
!!$          ( (u(ic,jc,k)-u(ic,j ,k))*rdy+(v(ic,jc,k)-v(i ,jc,k))*rdx )**2 +  &
!!$          ( (u(i ,jc,k)-u(i ,j ,k))*rdy+(v(i ,jc,k)-v(ib,jc,k))*rdx )**2 +  &
!!$          ( (u(ic,j ,k)-u(ic,jb,k))*rdy+(v(ic,j ,k)-v(i ,j ,k))*rdx )**2 +  &
!!$          ( (u(i ,j ,k)-u(i ,jb,k))*rdy+(v(i ,j ,k)-v(ib,j ,k))*rdx )**2 )   

    end do
 end do

end do ! k
	
end subroutine shear_prod_xy

