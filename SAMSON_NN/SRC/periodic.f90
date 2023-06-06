
subroutine periodic(flag)

use vars
implicit none

integer flag

if(flag.eq.0) then

  call bound_exchange(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,1,1,1,1,1)
  call bound_exchange(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,1,1,1,1,2)
  call bound_exchange(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,1,1,1,1,3)

endif

if(flag.eq.1) then

  call bound_exchange(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,2,3,2,2,1)
  call bound_exchange(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,2,2,2,3,2)
  call bound_exchange(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,2,2,2,2,3)

endif

if(flag.eq.2) then

 call bound_exchange(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,4)
 if(.not.dokruegermicro) &
      call bound_exchange(q,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,5)
 if(dosgs.and..not.dosmagor.or.doscalar) &
     call bound_exchange(tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,6)
 if(docloud.and.doprecip.and.(.not.dokruegermicro)) &
     call bound_exchange(qp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,9)

!====================================================================
! UW ADDITIONS

 if(dokruegermicro) then !bloss
    call bound_exchange(qv,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,19)
    call bound_exchange(qc,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,20)
    call bound_exchange(qi,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,21)
    call bound_exchange(qr,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,22)
    call bound_exchange(qg,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,23)
    call bound_exchange(qs,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,24)
 end if

 !kzm tracer exchange
 if(dotrz) call bound_exchange(trz,dimx1_s,dimx2_s,dimy1_s,dimy2_s&
      &,nzm,3,3,3,3,10)	
 if(dotrx)  call bound_exchange(trx,dimx1_s,dimx2_s,dimy1_s,dimy2_s&
      &,nzm,3,3,3,3,11)	  
 if(dotry)  call bound_exchange(try,dimx1_s,dimx2_s,dimy1_s,dimy2_s&
      &,nzm,3,3,3,3,12)	  
 if(dotrzz)  call bound_exchange(trzz,dimx1_s,dimx2_s,dimy1_s,dimy2_s&
      &,nzm,3,3,3,3,13) 
 if(dotro3) call bound_exchange(tro3,dimx1_s,dimx2_s,dimy1_s,dimy2_s&
      &,nzm,3,3,3,3,14)	  
!====================================================================

endif
        
if(flag.eq.3) then
        
 call bound_exchange(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4)
 if(.not.dokruegermicro) &
      call bound_exchange(q,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,5)
 if(dosgs.and..not.dosmagor.or.doscalar) &
     call bound_exchange(tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,6)
 if(docloud.and.doprecip.and.(.not.dokruegermicro)) &
     call bound_exchange(qp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,9)

!====================================================================
! UW ADDITIONS

 if(dokruegermicro) then !bloss
    call bound_exchange(qv,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,19)
    call bound_exchange(qc,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,20)
    call bound_exchange(qi,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,21)
    call bound_exchange(qr,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,22)
    call bound_exchange(qg,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,23)
    call bound_exchange(qs,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,24)
 end if

 !kzm tracer exchange
 if(dotrz) call bound_exchange(trz,dimx1_s,dimx2_s,dimy1_s,dimy2_s&
      &,nzm,1,1,1,1,10)	  
 if(dotrx) call bound_exchange(trx,dimx1_s,dimx2_s,dimy1_s,dimy2_s&
      &,nzm,1,1,1,1,11)	  
 if(dotry) call bound_exchange(try,dimx1_s,dimx2_s,dimy1_s,dimy2_s&
      &,nzm,1,1,1,1,12)	  
 if(dotrzz) call bound_exchange(trzz,dimx1_s,dimx2_s,dimy1_s,dimy2_s&
      &,nzm,1,1,1,1,13) 
 if(dotro3) call bound_exchange(tro3,dimx1_s,dimx2_s,dimy1_s,dimy2_s&
      &,nzm,1,1,1,1,14)	  
!====================================================================
        
endif
        
if(dosgs.and.flag.eq.4) then
        
 call bound_exchange(tk_xy,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,7) 
 call bound_exchange(tk_z,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,7)
 call bound_exchange(tkh_xy,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,8)
 call bound_exchange(tkh_z,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,8)
        
endif
        
        
end subroutine periodic

