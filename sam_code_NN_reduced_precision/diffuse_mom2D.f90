
subroutine diffuse_mom2D
	
!        momentum tendency due to SGS diffusion

use vars
implicit none

real rdx2,rdz2,rdz,rdx25,rdz25,rdx21,rdx251
real dxz,dzx

integer i,j,k,ic,ib,kc,kcu
real tkx, tkz, rhoi, iadzw, iadz
real fu(0:nx,1,nz),fv(0:nx,1,nz),fw(0:nx,1,nz)

rdx2=1./dx/dx
rdx25=0.25*rdx2

dxz=dx/dz

j=1

if(.not.docolumn) then

do k=1,nzm

 kc=k+1
 kcu=min(kc,nzm)
 dxz=dx/(dz*adzw(kc))
 rdx21=rdx2 * grdf_x(k)

 if(LES) then ! do geometric averaging of tk in vertical
   rdx251=2.*rdx25 * grdf_x(k)
   do i=0,nx
    ic=i+1
    tkx=rdx21*tk_xy(i,j,k)
    fu(i,j,k)=-2.*tkx*(u(ic,j,k)-u(i,j,k))*ravefactor
    fv(i,j,k)=-tkx*(v(ic,j,k)-v(i,j,k))*ravefactor
    tkx=rdx251*sqrt((tk_xy(i,j,k)+tk_xy(ic,j,k))*(tk_xy(i,j,kcu)+tk_xy(ic,j,kcu)))
    fw(i,j,k)=-tkx*((w(ic,j,kc)-w(i,j,kc))*ravefactor+(u(ic,j,kcu)-u(ic,j,k))*dxz/ravefactor)
   end do 

 else
   rdx251=rdx25 * grdf_x(k)
   do i=0,nx
    ic=i+1
    tkx=rdx21*tk_xy(i,j,k)
    fu(i,j,k)=-2.*tkx*(u(ic,j,k)-u(i,j,k))*ravefactor
    fv(i,j,k)=-tkx*(v(ic,j,k)-v(i,j,k))*ravefactor
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

end if 

!-------------------------
rdz=1./dz
dzx=dz/dx

do k=1,nzm-1
 kc=k+1
 uwsb(kc)=0.
 vwsb(kc)=0.
 iadz = 1./adz(k)
 iadzw= 1./adzw(kc)
 rdz2=rdz*rdz *grdf_z(k)
 if(LES) then ! do geometric averaging of tk in vertical
   rdz25=0.5*rdz2
   do i=1,nx
    ib=i-1
    tkz=rdz2*tk_z(i,j,k)
    fw(i,j,kc)=-2.*tkz*(w(i,j,kc)-w(i,j,k))*rho(k)*iadz/ravefactor
    tkz=rdz25*sqrt((tk_z(i,j,k)+tk_z(ib,j,k))*(tk_z(i,j,kc)+tk_z(ib,j,kc)))
    fu(i,j,kc)=-tkz*( (u(i,j,kc)-u(i,j,k))*iadzw/ravefactor + &
                       (w(i,j,kc)-w(ib,j,kc))*dzx*ravefactor)*rhow(kc) 	
    fv(i,j,kc)=-tkz*(v(i,j,kc)-v(i,j,k))*iadzw*rhow(kc)/ravefactor
    uwsb(kc)=uwsb(kc)+fu(i,j,kc)
    vwsb(kc)=vwsb(kc)+fv(i,j,kc)
  end do 
 else
   rdz25=0.25*rdz2
   do i=1,nx
    ib=i-1
    tkz=rdz2*tk_z(i,j,k)
    fw(i,j,kc)=-2.*tkz*(w(i,j,kc)-w(i,j,k))*rho(k)*iadz/ravefactor
    tkz=rdz25*(tk_z(i,j,k)+tk_z(ib,j,k)+tk_z(i,j,kc)+tk_z(ib,j,kc))
    fu(i,j,kc)=-tkz*( (u(i,j,kc)-u(i,j,k))*iadzw/ravefactor + &
                       (w(i,j,kc)-w(ib,j,kc))*dzx*ravefactor)*rhow(kc) 	
    fv(i,j,kc)=-tkz*(v(i,j,kc)-v(i,j,k))*iadzw*rhow(kc)/ravefactor
    uwsb(kc)=uwsb(kc)+fu(i,j,kc)
    vwsb(kc)=vwsb(kc)+fv(i,j,kc)
  end do 
 end if
end do

uwsb(1) = 0.
vwsb(1) = 0.
	
do i=1,nx
 tkz=rdz2*grdf_z(nzm)*tk_z(i,j,nzm)
 fw(i,j,nz)=-2.*tkz*(w(i,j,nz)-w(i,j,nzm))/adz(nzm)*rho(nzm)/ravefactor
 fu(i,j,1)=fluxbu(i,j) * rdz * rhow(1)
 fv(i,j,1)=fluxbv(i,j) * rdz * rhow(1)
 fu(i,j,nz)=fluxtu(i,j) * rdz * rhow(nz)
 fv(i,j,nz)=fluxtv(i,j) * rdz * rhow(nz)
 uwsb(1) = uwsb(1) + fu(i,j,1)
 vwsb(1) = vwsb(1) + fv(i,j,1)
end do
	 

do k=1,nzm
  kc=k+1
  rhoi = 1./(rho(k)*adz(k))
  do i=1,nx
    dudt(i,j,k,na)=dudt(i,j,k,na)-(fu(i,j,kc)-fu(i,j,k))*rhoi
    dvdt(i,j,k,na)=dvdt(i,j,k,na)-(fv(i,j,kc)-fv(i,j,k))*rhoi
  end do
end do ! k

do k=2,nzm
  rhoi = 1./(rhow(k)*adzw(k))
  do i=1,nx	 
    dwdt(i,j,k,na)=dwdt(i,j,k,na)-(fw(i,j,k+1)-fw(i,j,k))*rhoi
  end do
end do ! k


end subroutine diffuse_mom2D


