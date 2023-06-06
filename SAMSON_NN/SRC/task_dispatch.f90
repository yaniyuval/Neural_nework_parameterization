
subroutine task_dispatch(buff,tag)
	
! dispathes the messages according to the field sent.

use vars
use simple_land, only : landmask_frac, xysurf_rough, xysoil_wet
implicit none
	
real buff(*)	! buff for sending data
integer tag	! tag of the message
integer field   ! id of field
	
field = tag/100000
	
if(field.eq.1) then        
  call task_assign_bnd(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,buff,tag)
elseif(field.eq.2) then        
  call task_assign_bnd(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,buff,tag)
elseif(field.eq.3) then        
  call task_assign_bnd(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,buff,tag)	
elseif(field.eq.4) then        
  call task_assign_bnd(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)
elseif(field.eq.5) then        
  call task_assign_bnd(q,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)
elseif(field.eq.6) then        
  call task_assign_bnd(tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)
elseif(field.eq.9) then        
  call task_assign_bnd(qp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)	 

!=======================================================================
!  UW ADDITIONS

!bloss: added fields for krueger microphysics
elseif(field.eq.19) then        
  call task_assign_bnd(qv,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)	 
elseif(field.eq.20) then        
  call task_assign_bnd(qc,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)	 
elseif(field.eq.21) then        
  call task_assign_bnd(qi,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)	 
elseif(field.eq.22) then        
  call task_assign_bnd(qr,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)	 
elseif(field.eq.23) then        
  call task_assign_bnd(qg,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)	 
elseif(field.eq.24) then        
  call task_assign_bnd(qs,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)	 

!kzm added some tracers
elseif(field.eq.10) then        
  call task_assign_bnd(trz,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)
elseif(field.eq.11) then        
  call task_assign_bnd(trx,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)
elseif(field.eq.12) then        
  call task_assign_bnd(try,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)
elseif(field.eq.13) then        
  call task_assign_bnd(trzz,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)
elseif(field.eq.14) then        
  call task_assign_bnd(tro3,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)

! cw added for UW PBL scheme
elseif(field.eq.27) then        
  call task_assign_bnd(tk_xy,0,nxp1,1-YES3D,nyp1,nzm,buff,tag)
elseif(field.eq.28) then        
  call task_assign_bnd(tk_z,0,nxp1,1-YES3D,nyp1,nzm,buff,tag)
elseif(field.eq.37) then        
  call task_assign_bnd(tkh_xy,0,nxp1,1-YES3D,nyp1,nzm,buff,tag)        
elseif(field.eq.38) then        
  call task_assign_bnd(tkh_z,0,nxp1,1-YES3D,nyp1,nzm,buff,tag)    
elseif(field.eq.47) then
  call task_assign_bnd(ttend_uwpbl,0,nxp1,1-YES3D,nyp1,nzm,buff,tag)    
elseif(field.eq.48) then
  call task_assign_bnd(utend_uwpbl,0,nxp1,1-YES3D,nyp1,nzm,buff,tag)    
elseif(field.eq.49) then
  call task_assign_bnd(vtend_uwpbl,0,nxp1,1-YES3D,nyp1,nzm,buff,tag)    
elseif(field.eq.50) then
  call task_assign_bnd(qtend_uwpbl,0,nxp1,1-YES3D,nyp1,nzm,buff,tag)    
elseif(field.eq.51) then
  call task_assign_bnd(qptend_uwpbl,0,nxp1,1-YES3D,nyp1,nzm,buff,tag)    
elseif(field.eq.52) then
  call task_assign_bnd(wtend_uwpbl,0,nxp1,1-YES3D,nyp1,nz,buff,tag)    

! peters stuff for surface scheme
elseif(field.eq.70) then
  call task_assign_bnd(sstxy,0,nxp1,1-YES3D,nyp1,1,buff,tag)
elseif(field.eq.71) then
  call task_assign_bnd(xysoil_wet,0,nxp1,1-YES3D,nyp1,1,buff,tag)
elseif(field.eq.72) then
  call task_assign_bnd(landmask_frac,0,nxp1,1-YES3D,nyp1,1,buff,tag)
elseif(field.eq.73) then
  call task_assign_bnd(xysurf_rough,0,nxp1,1-YES3D,nyp1,1,buff,tag)

	
!  UW ADDITIONS
!=======================================================================
end if

end
	     
	     
