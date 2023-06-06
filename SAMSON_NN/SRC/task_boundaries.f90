
        subroutine task_boundaries(flag)
        	
!  These routines exchanges overlaping information for various subdomains to
!  be able to calculate various quantities in common /com3d/
!  near the subdomain boundaries.

use vars
use simple_land, only : landmask_frac, xysurf_rough, doxysoil_wet, xysoil_wet
implicit none

integer flag

if(flag.eq.0) then

 call task_exchange(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,1,1,1,1,1)
 call task_exchange(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,1,1,1,1,2)
 call task_exchange(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,1,1,1,1,3)	

endif

if(flag.eq.1) then

 call task_exchange(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,2,3,2,2,1)
 call task_exchange(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,2,2,2,3,2)
 call task_exchange(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,3,3,3,3,3)	

! call task_exchange(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,2,2,2,2,3)	

endif


if(flag.eq.2) then

 call task_exchange(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,4)
 if(.not.dokruegermicro) &
      call task_exchange(q,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,5)
 if(dosgs.and..not.dosmagor.or.doscalar) &
     call task_exchange(tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,6)
 if(docloud.and.doprecip.and.(.not.dokruegermicro)) &
     call task_exchange(qp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,9)

!===========================================================================
!  UW ADDITIONS

if(dokruegermicro) then
   call task_exchange(qv,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,19)
   call task_exchange(qc,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,20)
   call task_exchange(qi,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,21)
   call task_exchange(qr,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,22)
   call task_exchange(qg,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,23)
   call task_exchange(qs,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,24)
end if

!kzm tracer exchange
 if(dotrz) &
      call task_exchange(trz,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,10)
 if(dotrx) &
      call task_exchange(trx,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,11)
 if(dotry) &
      call task_exchange(try,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,12)
 if(dotrzz) &
      call task_exchange(trzz,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,13)
 if(dotro3) &
      call task_exchange(tro3,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,14)

!  END UW ADDITIONS
!===========================================================================
endif

if(flag.eq.3) then

 call task_exchange(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4)
 call task_exchange(q,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,5)
 if(dosgs.and..not.dosmagor.or.doscalar) &
     call task_exchange(tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,6)
 if(docloud.and.doprecip) &
     call task_exchange(qp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,9)

!===========================================================================
!  UW ADDITIONS

if(dokruegermicro) then
   call task_exchange(qv,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,19)
   call task_exchange(qc,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,20)
   call task_exchange(qi,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,21)
   call task_exchange(qr,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,22)
   call task_exchange(qg,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,23)
   call task_exchange(qs,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,24)
end if

 !kzm tracer exchange
 if(dotrz) &
      call task_exchange(trz,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,10)
 if(dotrx) &
      call task_exchange(trx,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,11)
 if(dotry) &
      call task_exchange(try,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,12)
 if(dotrzz) &
      call task_exchange(trzz,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,13)
 if(dotro3) &
      call task_exchange(tro3,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,14)
!  END UW ADDITIONS
!===========================================================================

endif

if(dosgs.and.flag.eq.4) then

 call task_exchange(tk_xy,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,27)
 call task_exchange(tk_z,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,28)        
 call task_exchange(tkh_xy,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,37)
 call task_exchange(tkh_z,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,38)        

endif

! peters
if(flag.eq.11) then
 call task_exchange(sstxy,0,nxp1,1-YES3D,nyp1,1,1,1,1,1,70)
 if(doxysoil_wet) then
   call task_exchange(xysoil_wet,0,nxp1,1-YES3D,nyp1,1,1,1,1,1,71)
 end if
endif

! peters
if(flag.eq.12) then
 call task_exchange(landmask_frac,0,nxp1,1-YES3D,nyp1,1,1,1,1,1,72)
endif

! peters
if(flag.eq.13) then
 call task_exchange(xysurf_rough,0,nxp1,1-YES3D,nyp1,1,1,1,1,1,73)
endif


end
	
	
