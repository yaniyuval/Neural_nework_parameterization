
subroutine rad_simple
 	
!	Simple Interactive Radiation
!  Coded by Ping Zhu for DycomsII LES intercomparison


use grid
use vars
use params
implicit none
	
real f0, xk, coef, coef1
integer i,k

!pzhu
real qq1,qq2,hfact,flux(nzm),FTHRL(nzm)
integer j,kk,itop,item3
!pzhu-end


if(.not.dolongwave) return

coef=70.
coef1=22.
f0=3.75e-6
xk=85.

do k=1,nzm
radlwdn(k) =0.
radqrlw(k) =0.
enddo

do i=1,nx
do j=1,ny
 do k=nzm,2,-1
  if(qn(i,j,k).gt.0.) then
   itop=k
   item3=1
   go to 101
  else
   item3=0
  endif
 enddo
101 continue
 do k=1,nzm
  qq1=0. 
  do kk=k,nzm-1
  qq1=qq1+xk*rho(kk)*qn(i,j,kk)*(Z(kk+1)-Z(kk))
  enddo
  qq2=0.
  do kk=k,2,-1
  qq2=qq2+xk*rho(kk)*qn(i,j,kk)*(Z(kk)-Z(kk-1))
  enddo
  flux(k)=coef*exp(-qq1)+coef1*exp(-qq2)
  if(item3.eq.1) then
   if(Z(k).gt.Z(itop)) then
    hfact=(Z(k)-Z(itop))**(4./3.)/4.+Z(itop)*(Z(k)-Z(itop))**(1./3.)
   else
    hfact=0.
   endif
   flux(k)=flux(k)+1015.*rho(k)*f0*hfact
  endif
 enddo

 do k=1,nzm-1
 FTHRL(k)=-(flux(k+1)-flux(k))/(Z(k+1)-Z(k))/rho(k)/1015.
 enddo
 FTHRL(nzm)=FTHRL(nzm-1)
 
 do k=1,nzm 
  t(i,j,k) = t(i,j,k) + FTHRL(k) * dt
  radlwdn(k) = radlwdn(k) + flux(k) 
  radqrlw(k) = radqrlw(k) + FTHRL(k)
 enddo
enddo
enddo    

end 



