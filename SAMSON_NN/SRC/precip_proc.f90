  
subroutine precip_proc

use vars
use params
use microscaling

implicit none

integer i,j,k
real autor, autos, accrr, accris, accrcs, accrig, accrcg
real dq, omn, omp, omg, qsat 
real pows1, pows2, powg1, powg2, powr1, powr2, tmp
real qii, qcc, qrr, qss, qgg
real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real f0(nzm),df0(nzm)
real dtn_old

!if (masterproc) then
!write(*,*) 'dtn', dtn
!write(*,*) 'tabs', tabs(1,1,:)
!write(*,*) 'qt', q(1,1,:)
!write(*,*) 'qp', qp(1,1,:)
!write(*,*) 'qn', qn(1,1,:)
!write(*,*) 'pres', pres
!endif

dtn_old = dtn
if(dobetafactor.and.daremicrophysics.and..not.domicroscaling) dtn=dtn*mpbetafactor**0.5

if(dokruegermicro) then
   call precip_proc_krueger()
   return
end if

powr1 = (3 + b_rain) / 4.
powr2 = (5 + b_rain) / 8.
pows1 = (3 + b_snow) / 4.
pows2 = (5 + b_snow) / 8.
powg1 = (3 + b_grau) / 4.
powg2 = (5 + b_grau) / 8.
      
     
if(dostatis) then
        
  do k=1,nzm
    do j=dimy1_s,dimy2_s
     do i=dimx1_s,dimx2_s

        if(domicroscaling) dtn=dtn_scaled(i,j,k)

      df(i,j,k) = q(i,j,k)
     end do
    end do
  end do
         
endif


do k=1,nzm
 qpsrc(k)=0.
 qpevp(k)=0.
 tlevp(k)=0.
 do j=1,ny
  do i=1,nx	  
     
                if(domicroscaling) dtn=dtn_scaled(i,j,k)

!-------     Autoconversion/accretion 

   if(qn(i,j,k)+qp(i,j,k).gt.0.) then


         omn = max(0.,min(1.,(tabs(i,j,k)-tbgmin)*a_bg))
         omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
         omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))
         qrr = qp(i,j,k) * omp
         qss = qp(i,j,k) * (1.-omp)*(1.-omg)
         qgg = qp(i,j,k) * (1.-omp)*omg

	 if(qn(i,j,k).gt.0.) then
     
           qcc = qn(i,j,k) * omn
           qii = qn(i,j,k) * (1.-omn)

           if(qcc .gt. qcw0) then
            autor = alphaelq
           else
            autor = 0.
           endif 

           if(qii .gt. qci0) then
            autos = betaelq*coefice(k)
           else
            autos = 0.
           endif 

           accrr = accrrc(k) * qrr ** powr1
           tmp = qss ** pows1
           accrcs = accrsc(k) * tmp
           accris = accrsi(k) * tmp
           tmp = qgg ** powg1
           accrcg = accrgc(k) * tmp
           accrig = accrgi(k) * tmp
           qcc = (qcc+dtn*autor*qcw0)/(1.+dtn*(accrr+accrcs+accrcg+autor))
           qii = (qii+dtn*autos*qci0)/(1.+dtn*(accris+accrig+autos))
           dq = dtn *(accrr*qcc + autor*(qcc-qcw0)+ &
             (accris+accrig)*qii + (accrcs+accrcg)*qcc + autos*(qii-qci0))
           dq = min(dq,qn(i,j,k))
           qp(i,j,k) = qp(i,j,k) + dq
           q(i,j,k) = q(i,j,k) - dq
           qn(i,j,k) = qn(i,j,k) - dq
	   qpsrc(k) = qpsrc(k) + dq

!if (masterproc) then
!write(*,*) 'a', -dq
!endif

         elseif(qp(i,j,k).gt.qp_threshold.and.qn(i,j,k).eq.0.) then
            if(tabs(i,j,k).gt.tmin_evap) then !kzm limit evaporation to temperatures above tmin_evap

           qsat = omn*qsatw(tabs(i,j,k),pres(k))+ &
                 (1.-omn)*qsati(tabs(i,j,k),pres(k))
           dq = dtn * (evapr1(k) * sqrt(qrr) + &
                		evapr2(k) * qrr**powr2 + &
				evaps1(k) * sqrt(qss) + &
				evaps2(k) * qss**pows2 + &
				evapg1(k) * sqrt(qgg) + &
				evapg2(k) * qgg**powg2)* &
			 	(q(i,j,k) /qsat-1.) 
           dq = max(-0.5*qp(i,j,k),dq) 
           qp(i,j,k) = qp(i,j,k) + dq
           q(i,j,k) = q(i,j,k) - dq
	   qpevp(k) = qpevp(k) + dq/dtn
	   tlevp(k) = tlevp(k) - dq*(fac_cond + fac_fus*(1.-omp))/dtn
!if (masterproc) then
!write(*,*) 'b', -dq
!endif
           endif
	 else
	
           q(i,j,k) = q(i,j,k) + qp(i,j,k)
	   qpevp(k) = qpevp(k) + qp(i,j,k)/dtn
	   tlevp(k) = tlevp(k) - qp(i,j,k)*(fac_cond + fac_fus*(1.-omp))/dtn
!if (masterproc) then
!write(*,*) 'c', -qp(1,1,k)
!endif
           qp(i,j,k) = 0.

         endif

    endif

    qp(i,j,k)=max(0.,qp(i,j,k))

  enddo
 enddo

enddo


    
if(dostatis) then
                  
  call stat_varscalar(q,df,f0,df0,q2leprec)
  call setvalue(qwleprec,nzm,0.)
  call stat_sw2(q,df,qwleprec)

endif

dtn = dtn_old

end subroutine precip_proc

