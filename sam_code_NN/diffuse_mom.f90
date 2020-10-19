subroutine diffuse_mom

!  Interface to the diffusion routines

use vars
use mse          ! peters
implicit none
integer i,j,k
real du(nx,ny,nz,3)

real dtu_old(nxp1,ny,nzm), dtv_old(nx,nyp1,nzm)  ! peters
dtu_old = dudt(:,:,:,na)           ! peters
dtv_old = dvdt(:,:,:,na)           ! peters

if(dostatis) then
	
  do k=1,nzm
   do j=1,ny
    do i=1,nx
     du(i,j,k,1)=dudt(i,j,k,na)
     du(i,j,k,2)=dvdt(i,j,k,na)
     du(i,j,k,3)=dwdt(i,j,k,na)
    end do
   end do
  end do

endif

if(RUN3D) then
   call diffuse_mom3D()
else
   call diffuse_mom2D()
endif


! update the diffusive mom fluxes
if(isallocateMSE)  then
  dtu_diff(:,:,:,na) = dudt(:,:,:,na) - dtu_old    ! peters
  dtv_diff(:,:,:,na) = dvdt(:,:,:,na) - dtv_old    ! peters
end if

if(dostatis) then
	
  do k=1,nzm
   do j=1,ny
    do i=1,nx
     du(i,j,k,1)=dudt(i,j,k,na)-du(i,j,k,1)
     du(i,j,k,2)=dvdt(i,j,k,na)-du(i,j,k,2)
     du(i,j,k,3)=dwdt(i,j,k,na)-du(i,j,k,3)
    end do
   end do
  end do

  call stat_tke(du,tkelediff)
  call stat_mom(du,momlediff)
  call setvalue(twlediff,nzm,0.)
  call setvalue(qwlediff,nzm,0.)
  call setvalue(swlediff,nzm,0.)
  call stat_sw1(du,twlediff,qwlediff,swlediff)

endif


end subroutine diffuse_mom

