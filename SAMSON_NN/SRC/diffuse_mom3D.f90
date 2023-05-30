
subroutine diffuse_mom3D
	
!        momentum tendency due to SGS diffusion

use vars
implicit none

real rdx2,rdy2,rdz2,rdz,rdx25,rdy25
real rdx21,rdy21,rdx251,rdy251,rdz25
real dxy,dxz,dyx,dyz,dzx,dzy

integer i,j,k,ic,ib,jb,jc,kc,kcu
real tkx, tky, tkz, rhoi, iadzw, iadz
real fu(0:nx,0:ny,nz),fv(0:nx,0:ny,nz),fw(0:nx,0:ny,nz)

rdx2=1./(dx*dx)
rdy2=1./(dy*dy)

rdx25=0.25*rdx2
rdy25=0.25*rdy2

dxy=dx/dy
dxz=dx/dz
dyx=dy/dx
dyz=dy/dz


do k=1,nzm
 kc=k+1
 kcu=min(kc,nzm)
 dxz=dx/(dz*adzw(kc))
 dyz=dy/(dz*adzw(kc))
  rdx21=rdx2    * grdf_x(k)
  rdy21=rdy2    * grdf_y(k)
  rdx251=rdx25  * grdf_x(k)
  rdy251=rdy25  * grdf_y(k)
  do j=1,ny
   jb=j-1
   if(LES) then !bloss  do geometric averaging of tk in vertical
    do i=0,nx
     ic=i+1
     tkx=rdx21*tk_xy(i,j,k)
     fu(i,j,k)=-2.*tkx*(u(ic,j,k)-u(i,j,k))*ravefactor
     tkx=rdx251*(tk_xy(i,j,k)+tk_xy(i,jb,k)+tk_xy(ic,j,k)+tk_xy(ic,jb,k))
     fv(i,j,k)=-tkx*(v(ic,j,k)-v(i,j,k)+(u(ic,j,k)-u(ic,jb,k))*dxy)*ravefactor
     tkx=2.*rdx251*sqrt((tk_xy(i,j,k)+tk_xy(ic,j,k))*(tk_xy(i,j,kcu)+tk_xy(ic,j,kcu)))
     fw(i,j,k)=-tkx*((w(ic,j,kc)-w(i,j,kc))*ravefactor+(u(ic,j,kcu)-u(ic,j,k))*dxz/ravefactor)
    end do 
   else
    do i=0,nx
     ic=i+1
     tkx=rdx21*tk_xy(i,j,k)
     fu(i,j,k)=-2.*tkx*(u(ic,j,k)-u(i,j,k))*ravefactor
     tkx=rdx251*(tk_xy(i,j,k)+tk_xy(i,jb,k)+tk_xy(ic,j,k)+tk_xy(ic,jb,k))
     fv(i,j,k)=-tkx*(v(ic,j,k)-v(i,j,k)+(u(ic,j,k)-u(ic,jb,k))*dxy)*ravefactor
     tkx=rdx251*(tk_xy(i,j,k)+tk_xy(ic,j,k)+tk_xy(i,j,kcu)+tk_xy(ic,j,kcu)) 	
     fw(i,j,k)=-tkx*((w(ic,j,kc)-w(i,j,kc))*ravefactor+(u(ic,j,kcu)-u(ic,j,k))*dxz/ravefactor)
    end do 
   end if
   do i=1,nx
    ib=i-1
    dudt(i,j,k,na)=dudt(i,j,k,na)-(fu(i,j,k)-fu(ib,j,k))
    dvdt(i,j,k,na)=dvdt(i,j,k,na)-(fv(i,j,k)-fv(ib,j,k))
    dwdt(i,j,kc,na)=dwdt(i,j,kc,na)-(fw(i,j,k)-fw(ib,j,k))
   end do  
  end do 

  do j=0,ny
   jc=j+1
   if(LES) then !bloss  do geometric averaging of tk in vertical
    do i=1,nx
     ib=i-1
     tky=rdy21*tk_xy(i,j,k)
     fv(i,j,k)=-2.*tky*(v(i,jc,k)-v(i,j,k))*ravefactor
     tky=rdy251*(tk_xy(i,j,k)+tk_xy(ib,j,k)+tk_xy(i,jc,k)+tk_xy(ib,jc,k))
     fu(i,j,k)=-tky*(u(i,jc,k)-u(i,j,k)+(v(i,jc,k)-v(ib,jc,k))*dyx)*ravefactor
     tky=2.*rdy251*sqrt((tk_xy(i,j,k)+tk_xy(i,jc,k))*(tk_xy(i,j,kcu)+tk_xy(i,jc,kcu)))
     fw(i,j,k)=-tky*((w(i,jc,kc)-w(i,j,kc))*ravefactor+(v(i,jc,kcu)-v(i,jc,k))*dyz/ravefactor)
    end do 
   else
    do i=1,nx
     ib=i-1
     tky=rdy21*tk_xy(i,j,k)
     fv(i,j,k)=-2.*tky*(v(i,jc,k)-v(i,j,k))*ravefactor
     tky=rdy251*(tk_xy(i,j,k)+tk_xy(ib,j,k)+tk_xy(i,jc,k)+tk_xy(ib,jc,k))
     fu(i,j,k)=-tky*(u(i,jc,k)-u(i,j,k)+(v(i,jc,k)-v(ib,jc,k))*dyx)*ravefactor
     tky=rdy251*(tk_xy(i,j,k)+tk_xy(i,jc,k)+tk_xy(i,j,kcu)+tk_xy(i,jc,kcu)) 	
     fw(i,j,k)=-tky*((w(i,jc,kc)-w(i,j,kc))*ravefactor+(v(i,jc,kcu)-v(i,jc,k))*dyz/ravefactor)
    end do 
   end if
  end do 
  do j=1,ny
    jb=j-1
    do i=1,nx	    
     dudt(i,j,k,na)=dudt(i,j,k,na)-(fu(i,j,k)-fu(i,jb,k))
     dvdt(i,j,k,na)=dvdt(i,j,k,na)-(fv(i,j,k)-fv(i,jb,k))
     dwdt(i,j,kc,na)=dwdt(i,j,kc,na)-(fw(i,j,k)-fw(i,jb,k))
   end do 
  end do 

end do 
 
!-------------------------
rdz=1./dz
dzx=dz/dx
dzy=dz/dy

do k=1,nzm-1
 kc=k+1
 uwsb(kc)=0.
 vwsb(kc)=0.
 iadz = 1./adz(k)
 iadzw= 1./adzw(kc)
 rdz2 = rdz*rdz * grdf_z(k)
 if(LES) then !bloss  do geometric averaging of tk in vertical
   rdz25 = 0.5*rdz2
   do j=1,ny
    jb=j-1
    do i=1,nx
     ib=i-1
     tkz=rdz2*tk_z(i,j,k)
     fw(i,j,kc)=-2.*tkz*(w(i,j,kc)-w(i,j,k))*rho(k)*iadz/ravefactor
     tkz=rdz25*sqrt((tk_z(i,j,k)+tk_z(ib,j,k))*(tk_z(i,j,kc)+tk_z(ib,j,kc)))
     fu(i,j,kc)=-tkz*( (u(i,j,kc)-u(i,j,k))*iadzw/ravefactor + &
                        (w(i,j,kc)-w(ib,j,kc))*dzx*ravefactor)*rhow(kc) 	
     tkz=rdz25*sqrt((tk_z(i,j,k)+tk_z(i,jb,k))*(tk_z(i,j,kc)+tk_z(i,jb,kc)))
     fv(i,j,kc)=-tkz*( (v(i,j,kc)-v(i,j,k))*iadzw/ravefactor + &
                        (w(i,j,kc)-w(i,jb,kc))*dzy*ravefactor)*rhow(kc)
     uwsb(kc)=uwsb(kc)+fu(i,j,kc)
     vwsb(kc)=vwsb(kc)+fv(i,j,kc)
   end do 
  end do
 else
   rdz25 = 0.25*rdz2
   do j=1,ny
    jb=j-1
    do i=1,nx
     ib=i-1
     tkz=rdz2*tk_z(i,j,k)
     fw(i,j,kc)=-2.*tkz*(w(i,j,kc)-w(i,j,k))*rho(k)*iadz/ravefactor
     tkz=rdz25*(tk_z(i,j,k)+tk_z(ib,j,k)+tk_z(i,j,kc)+tk_z(ib,j,kc))
     fu(i,j,kc)=-tkz*( (u(i,j,kc)-u(i,j,k))*iadzw/ravefactor + &
                        (w(i,j,kc)-w(ib,j,kc))*dzx*ravefactor)*rhow(kc) 	
     tkz=rdz25*(tk_z(i,j,k)+tk_z(i,jb,k)+tk_z(i,j,kc)+tk_z(i,jb,kc))
     fv(i,j,kc)=-tkz*( (v(i,j,kc)-v(i,j,k))*iadzw/ravefactor + &
                        (w(i,j,kc)-w(i,jb,kc))*dzy*ravefactor)*rhow(kc)
     uwsb(kc)=uwsb(kc)+fu(i,j,kc)
     vwsb(kc)=vwsb(kc)+fv(i,j,kc)
   end do 
  end do
 end if
end do

uwsb(1) = 0.
vwsb(1) = 0.
	
do j=1,ny
 do i=1,nx
   tkz=rdz2*grdf_z(nzm)*tk_z(i,j,nzm)
   fw(i,j,nz)=-2.*tkz*(w(i,j,nz)-w(i,j,nzm))/adz(nzm)*rho(nzm)/ravefactor
!!$   fu(i,j,1)=(fluxbu(i,j)) * rdz * rhow(1)
!!$   fv(i,j,1)=(fluxbv(i,j)) * rdz * rhow(1)
   fu(i,j,1)=(fluxbu(i,j)-uwpbl_fluxbu(i,j)) * rdz * rhow(1)
   fv(i,j,1)=(fluxbv(i,j)-uwpbl_fluxbv(i,j)) * rdz * rhow(1)
   fu(i,j,nz)=fluxtu(i,j) * rdz * rhow(nz)
   fv(i,j,nz)=fluxtv(i,j) * rdz * rhow(nz)
   uwsb(1) = uwsb(1) + fu(i,j,1)
   vwsb(1) = vwsb(1) + fv(i,j,1)
  end do
 end do
 
 do k=1,nzm
  kc=k+1
  rhoi = 1./(rho(k)*adz(k))
  do j=1,ny	  
   do i=1,nx
    dudt(i,j,k,na)=dudt(i,j,k,na)-(fu(i,j,kc)-fu(i,j,k))*rhoi
    dvdt(i,j,k,na)=dvdt(i,j,k,na)-(fv(i,j,kc)-fv(i,j,k))*rhoi
   end do
  end do
 end do ! k

 do k=2,nzm
  rhoi = 1./(rhow(k)*adzw(k))
  do j=1,ny
   do i=1,nx	 
    dwdt(i,j,k,na)=dwdt(i,j,k,na)-(fw(i,j,k+1)-fw(i,j,k))*rhoi
   end do
  end do
 end do ! k


end subroutine diffuse_mom3D
