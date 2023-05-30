
! Monin-Obukhov Similarity
! Coded by Marat Khairoutdinov (C) 2003
! Modified by Zhiming Kuang 2004, following scheme in Garrat's book.

subroutine landflx(p0, th, ts, qh, qs, uh, vh, h, z0, shf, lhf, taux, tauy, xlmo)

implicit none

! Input:

real p0   ! surface pressure, mb
real th   ! pot. temperature at height h
real ts   ! pot. Temperature at z0
real qh   ! vapor at height h
real qs   ! vapor at z0
real uh   ! zonal wind at height h
real vh   ! merid wind at height h
real h    ! height h
real z0   ! friction height
        
! Output:

real shf   ! sensible heat flux (K m/s)
real lhf   ! latent heat flux (m/s)
real taux  ! zonal surface stress (N/m2)
real tauy  ! merid surface stress (N/m2)
real xlmo  ! Monin-Obukhov length

real r, x, pii, zody, vel
real a, b, c, d, ustar, tstar
real xm, xh, xsi, xsi1, xsi2, dxsi, fm, fh
integer iter
real gm1, gh1, fm1, fh1

!Use Garratt (p.52), instead of the Kansas experiment
gm1(x)=(1.-16.*x)**0.25
gh1(x)=sqrt(1.-16.*x)
fm1(x)=2.*alog((1.+x)/2.)+alog((1.+x*x)/2.)-2.*atan(x)+pii
fh1(x)=2.*alog((1.+x)/2.)

pii=acos(-1.)/2.
zody=alog(h/z0)

vel = sqrt(max(0.5,uh**2+vh**2))
r=9.81/ts*(th*(1+0.61*qh)-ts*(1.+0.61*qs))*h/vel**2

!limit Ri to be -10<Ri<0.1999
r=max(-10.,r)
r=min(0.1999,r)

xsi=r*zody

if(xsi .lt.0) then!3 iterations
   xm=gm1(xsi)
   xh=gh1(xsi)
   fm=zody-fm1(xm)
   fh=zody-fh1(xh)
   xsi1=r/fh*fm**2

   xsi=xsi1
   xm=gm1(xsi)
   xh=gh1(xsi)
   fm=zody-fm1(xm)
   fh=zody-fh1(xh)
   xsi1=r/fh*fm**2

   xsi=xsi1
   xm=gm1(xsi)
   xh=gh1(xsi)
   fm=zody-fm1(xm)
   fh=zody-fh1(xh)
   xsi1=r/fh*fm**2
   xsi=xsi1
else 
   !following Garratt's book
   xsi=zody*r/(1.-5.*r)
   fm=zody+5.*xsi
   fh=zody+5.*xsi
endif


vel = sqrt(uh**2+vh**2)
shf=0.4**2/fm/fh*vel*(ts-th)
lhf=0.4**2/fm/fh*vel*(qs-qh)
taux=-0.4**2/fm/fm*vel*uh*(p0*100./287./ts)
tauy=-0.4**2/fm/fm*vel*vh*(p0*100./287./ts)
      
ustar = 0.4/fm*vel
tstar = 0.4/fh*(th-ts)
if(xsi.ge.0.) then
   xsi = max(1.e-5,xsi)
else
   xsi = min(-1.e-5,xsi)
end if
xlmo = h/xsi

return
end
