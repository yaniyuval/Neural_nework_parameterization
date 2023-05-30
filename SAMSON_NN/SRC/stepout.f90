subroutine stepout(nstatsteps)

use vars
use params
use hbuffer
use isccp, only : isccp_write
implicit none	
	
integer i,j,k,ic,jc,nstatsteps
real div, divmax, divmin
real rdx, rdy, rdz, coef
integer im,jm,km
real wmax, qnmax(1), qnmax1(1)
real buffer(5), buffer1(5)

if(masterproc) then
   write(*,999) nstep, ncycle, maxcfl
999 format('NSTEP = ',i20,'    NCYCLE=',i4,'    MAX CFL=',f7.3)
end if

if(mod(nstep,nprint).eq.0) then
	

 divmin=1.e20
 divmax=-1.e20
	 
 rdx = 1./dx
 rdy = 1./dy

 wmax=0.
 do k=1,nzm
  coef = rho(k)*adz(k)*dz
  rdz = 1./coef
  if(ny.ne.1) then
   do j=1,ny-1*YES3D
    jc = j+1*YES3D
    do i=1,nx-1
     ic = i+1
     div = (u(ic,j,k)-u(i,j,k))*rdx + (v(i,jc,k)-v(i,j,k))*rdy + &
		  (w(i,j,k+1)*rhow(k+1)-w(i,j,k)*rhow(k))*rdz
     divmax = max(divmax,div)
     divmin = min(divmin,div)
     if(w(i,j,k).gt.wmax) then
	wmax=w(i,j,k)
	im=i
	jm=j
	km=k
     endif
    end do
   end do
  else
    j = 1
    do i=1,nx-1
    ic = i+1
     div = (u(ic,j,k)-u(i,j,k))*rdx +(w(i,j,k+1)*rhow(k+1)-w(i,j,k)*rhow(k))*rdz
     divmax = max(divmax,div)
     divmin = min(divmin,div)
     if(w(i,j,k).gt.wmax) then
	wmax=w(i,j,k)
	im=i
	jm=j
	km=k
     endif
    end do
  endif
 end do

 if(dompi) then
   buffer(1) = total_water_before
   buffer(2) = total_water_after
   buffer(3) = total_water_evap
   buffer(4) = total_water_prec
   buffer(5) = total_water_ls
   call task_sum_real(buffer, buffer1,5)
   total_water_before = buffer1(1)
   total_water_after = buffer1(2)
   total_water_evap = buffer1(3)
   total_water_prec = buffer1(4)
   total_water_ls = buffer1(5)
 end if

!write(6,'(16f7.2)')((p(i,1,k),i=1,16),k=nzm,1,-1)
!print*
!write(6,'(16f7.2)')((u(i,1,k),i=1,16),k=nzm,1,-1)
!print*
!write(6,'(16f7.2)')((v(i,1,k),i=1,16),k=nzm,1,-1)
!print*
!write(6,'(16f7.2)')((t(i,1,k),i=1,16),k=nzm,1,-1)

!--------------------------------------------------------
 if(masterproc) then
	
    print*,'DAY = ',day	
    write(6,*) 'NSTEP=',nstep
    write(6,*) 'div:',divmax,divmin
    write(6,*) 'SST=',tabs_s, '  surface pressure=',pres0

 endif

 call fminmax_print('u:',u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm)
 call fminmax_print('v:',v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm-5)
 call fminmax_print('w:',w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz)
 call fminmax_print('p:',p,0,nx,1-YES3D,ny,nzm)
 call fminmax_print('t:',t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 call fminmax_print('tabs:',tabs,1,nx,1,ny,nzm)
 call fminmax_print('q:',q,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 call fminmax_print('tke:',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 call fminmax_print('tk_xy:',tk_xy,0,nxp1,1-YES3D,nyp1,nzm)
 call fminmax_print('tk_z:',tk_z,0,nxp1,1-YES3D,nyp1,nzm)
 call fminmax_print('tkh_xy:',tkh_xy,0,nxp1,1-YES3D,nyp1,nzm)
 call fminmax_print('tkh_z:',tkh_z,0,nxp1,1-YES3D,nyp1,nzm)
 call fminmax_print('qn:',qn,1,nx,1,ny,nzm)
 call fminmax_print('qp:',qp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 if(dokruegermicro) then
    call fminmax_print('qv:',qv,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
    call fminmax_print('qc:',qc,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
    call fminmax_print('qi:',qi,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
    call fminmax_print('qr:',qr,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
    call fminmax_print('qg:',qg,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
    call fminmax_print('qs:',qs,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 endif
 if(dotrx) call fminmax_print('trx:',trx,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 if(dotry) call fminmax_print('try:',try,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 if(dotrz) call fminmax_print('trz:',trz,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 if(dotrzz) call fminmax_print('trzz:',trzz,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 if(dotro3) call fminmax_print('tro3:',tro3,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 if(dolongwave.or.doshortwave) call fminmax_print('qrad(K/day):',misc*86400.,1,nx,1,ny,nzm)
 call fminmax_print('shf:',fluxbt*cp*rho(1),1,nx,1,ny,1)
 call fminmax_print('lhf:',fluxbq*lcond*rho(1),1,nx,1,ny,1)
 call fminmax_print('uw:',fluxbu,1,nx,1,ny,1)
 call fminmax_print('vw:',fluxbv,1,nx,1,ny,1)
 call fminmax_print('sst:',sstxy(1:nx,1:ny),1,nx,1,ny,1)
 if (dodynamicocean) call fminmax_print('qocn:',qoceanxy,1,nx,1,ny,1)

 total_water_before = total_water_before/float(nx_gl*ny_gl)
 total_water_after = total_water_after/float(nx_gl*ny_gl)
 total_water_evap = total_water_evap/float(nx_gl*ny_gl)
 total_water_prec = total_water_prec/float(nx_gl*ny_gl)
 total_water_ls = total_water_ls/float(nx_gl*ny_gl)
 
 if(masterproc) then
   
   print*,'total water budget:'
   print*,'before (mm):',total_water_before
   print*,'after (mm) :',total_water_after
   print*,'evap (mm/day):',total_water_evap/dtn*86400.
   print*,'prec (mm/day):',total_water_prec/dtn*86400.
   print*,'ls (mm/day):',total_water_ls/dtn*86400.
   print*,'disbalance (mm/day)', &
     (total_water_after-(total_water_before+total_water_evap+total_water_ls-total_water_prec))/dtn * 86400.
   print*,'disbalance (mm)', &
     (total_water_after-(total_water_before+total_water_evap+total_water_ls-total_water_prec))

 end if



end if ! (mod(nstep,nprint).eq.0)
	
if(nstep.eq.nsave2Dstart-nsave2D) then
     prec_xy(:,:) = 0.
     shf_xy(:,:) = 0. 
     lhf_xy(:,:) = 0. 
     lwnt_xy(:,:) = 0. 
     swnt_xy(:,:) = 0. 
     lwntc_xy(:,:) = 0. 
     swntc_xy(:,:) = 0. 
     pw_xy(:,:) = 0.
     cw_xy(:,:) = 0.
     iw_xy(:,:) = 0.
     u200_xy(:,:) = 0.
     v200_xy(:,:) = 0.
     usfc_xy(:,:) = 0.
     vsfc_xy(:,:) = 0.

     w500_xy = 0.

     lwns_xy(:,:) = 0.
     swns_xy(:,:) = 0.
     solin_xy(:,:) = 0.
     lwnsc_xy(:,:) = 0.
     swnsc_xy(:,:) = 0.
     qocean_xy(:,:) = 0.

     swntm_xy(:,:) = 0. 
     psfc_xy = 0.
     h200_xy = 0.
     h850_xy = 0.
     u850_xy = 0.
     v850_xy = 0.
     the850_xy = 0.
     thes850_xy = 0.
     the1000_xy = 0.
     thes1000_xy = 0.
     swvp_xy(:,:) = 0.
     wmode1_xy(:,:) = 0.
     wmode2_xy(:,:) = 0.
     qmode1_xy(:,:) = 0.
     qmode2_xy(:,:) = 0.
     thmode1_xy(:,:) = 0.
     thmode2_xy(:,:) = 0.
     wmode1i_xy(:,:) = 0.
     wmode2i_xy(:,:) = 0.
     hstor_xy(:,:) = 0.
     hadv_xy(:,:)  = 0.
     sstor_xy(:,:) = 0.
     sadv_xy(:,:)  = 0.
end if
if(mod(nstep,nsave2D).eq.0.and.nstep.ge.nsave2Dstart &
                                   .and.nstep.le.nsave2Dend) then
  call write_fields2D()
endif
	
if(mod(nstep,nstat).eq.0) then
  if(masterproc) print *,'Writting statistics:nstatsteps=',nstatsteps

  call hbuf_average(nstatsteps)
  call hbuf_write(nstatsteps)
  call hbuf_flush()	  
  nstatsteps = 0
	  
  call isccp_write()

  ! peters.  if nwriterestart is set, only write restart file if time to do so
  if (.not.dorestart_last .and.  &
      (nwriterestart.eq.-1 .or. mod(nstep,nwriterestart).eq.0)) then
    call write_all() ! save restart file
  end if
	  
endif


! determine if the maximum cloud water exceeds the threshold
! value to save 3D fields:

qnmax(1)=0.
do k=1,nzm
  do j=1,ny
    do i=1,nx
       qnmax(1) = max(qnmax(1),qn(i,j,k))
    end do
  enddo
enddo
if(dompi) then
   call task_max_real(qnmax(1),qnmax1(1),1)
   qnmax(1) = qnmax1(1)
end if
if(mod(nstep,nsave3D).eq.0.and.nstep.ge.nsave3Dstart.and.nstep.le.nsave3Dend &
                                  .and.qnmax(1).ge.qnsave3D) then
  call write_fields3D()
endif


end
