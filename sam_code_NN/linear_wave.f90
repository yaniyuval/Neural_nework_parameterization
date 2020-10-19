
subroutine linear_wave

!Zhiming Kuang, Nov. 27, 2006
!parameterized large-scale linear wave dynamics
use vars
use params
implicit none
integer k,i

real dampingprofile(nzm)
real wn2,wn,N2top,aa(nzm),bb(nzm),cc(nzm),rhs(nzm)
real ztop
integer kztop
ztop=14000.

!profile of wave damping coefficients
do k=1,nz
   dampingprofile(k)=1.+(1-max(min(1.,(zi(k)/1000.-1.)/4),0.))*(pbldampingfactor-1)
enddo

if(nstep.gt.nstartlinearwave) then
   if(doevolvewavebg) then
      !first compute the background, which evolves with time
      t_wavebg_hist(:,wavecounter)=t_wavebg_hist(:,wavecounter)+tabs0/ncycle
      q_wavebg_hist(:,wavecounter)=q_wavebg_hist(:,wavecounter)+q0/ncycle
     if(mod(nstep-nstartlinearwave,nsteplinearwavebg).eq.0.and.icycle.eq.ncycle) then
         
         t_wavebg_hist(:,wavecounter)=t_wavebg_hist(:,wavecounter)/nsteplinearwavebg
         q_wavebg_hist(:,wavecounter)=q_wavebg_hist(:,wavecounter)/nsteplinearwavebg
 
         if (nstep.eq.nstartlinearwave+nsteplinearwavebg) then !first time
            do i=1,10 
               t_wavebg_hist(:,i)=t_wavebg_hist(:,wavecounter)
               q_wavebg_hist(:,i)=q_wavebg_hist(:,wavecounter)
             enddo 
         endif
         t_wavebg=0.
         q_wavebg=0.
         do i=1,10 
            t_wavebg=t_wavebg+t_wavebg_hist(:,i)/10
            q_wavebg=q_wavebg+q_wavebg_hist(:,i)/10
         enddo
         
         wavecounter=wavecounter+1
         if (wavecounter.gt.10) wavecounter=1
         t_wavebg_hist(:,wavecounter)=0.
         q_wavebg_hist(:,wavecounter)=0.
       endif
   else
      !fixed background
      if(nstep.le.nstartlinearwave+nsteplinearwavebg) then !first compute the background
         t_wavebg=t_wavebg+tabs0/ncycle
         q_wavebg=q_wavebg+q0/ncycle
          if (nstep.eq.nstartlinearwave+nsteplinearwavebg.and.icycle.eq.ncycle) then
            t_wavebg=t_wavebg/nsteplinearwavebg
            q_wavebg=q_wavebg/nsteplinearwavebg
          endif
      endif
   endif

   if(nstep.gt.nstartlinearwave+nsteplinearwavebg) then 
      t_wave_local=t_wave_local+tabs0/ncycle
      q_wave_local=q_wave_local+q0/ncycle
       if(mod(nstep-nstartlinearwave-nsteplinearwavebg,nsteplinearwave).eq.0.and.icycle.eq.ncycle) then
         t_wave=t_wave_local/nsteplinearwave
         q_wave=q_wave_local/nsteplinearwave
         if(dolinearwavelid) then
            tv_wavebg=t_wavebg*(1+0.61*q_wavebg)
            tv_wave=t_wave*(1+0.61*q_wave)
         else
            tv_wavebg=t_wavebg*(1+0.61*q_wavebg)
            tv_wave=t_wave*(1+0.61*q_wave)            
         endif
        
         wn=sqrt(betafactor)*wavenumber_factor*6.283e-6 !2pi/1000km	
         wn2=wn*wn
         N2top=4.e-4

         rhs=0.
         !Neglect the top layer
         do k=1,nzm-1 !assume constant damping time with height
            rhs(k)=-rho(k)*ggr*wn2*(tv_wave(k)-tv_wavebg(k))/tv_wavebg(k)*adzw(k)*adzw(k-1)*dz*dz/2.
         end do
         !set up the tridiagonal matrix
         do k=1,nzm-1
            aa(k)=adzw(k+1)/(adzw(k)+adzw(k+1))
            bb(k)=-1.
            cc(k)=adzw(k)/(adzw(k)+adzw(k+1))
         end do
         !symmetric lower BC
         aa(1)=0.
         bb(1)=-(2*adzw(2)+adzw(1))/(adzw(1)+adzw(2))
         
         if(dolinearwavelid) then
            !symmetric upper BC
            bb(nzm-1)=-(2*adzw(nzm-1)+adzw(nzm))/(adzw(nzm-1)+adzw(nzm))
            cc(nzm-1)=0.
         else
            !Radiating upper BC
            aa(nzm-1)=adzw(nzm)/adzw(nzm-1)
            bb(nzm-1)=-1.*adzw(nzm)/adzw(nzm-1)
            rhs(nzm-1)=rhs(nzm-1)+sqrt(N2top)*wn*w_wave(nzm-1)*adzw(nzm)
         endif
         
         !Gaussian Elimination with no pivoting
         do k=1,nzm-2
            bb(k+1)=bb(k+1)-aa(k+1)/bb(k)*cc(k)
            rhs(k+1)=rhs(k+1)-aa(k+1)/bb(k)*rhs(k)
         end do
         !backward substitution
         rhs(nzm-1)=rhs(nzm-1)/bb(nzm-1);
         do k= nzm-2,1, -1
            rhs(k)=(rhs(k)-cc(k)*rhs(k+1))/bb(k)
         end do
         
          if(dosteadylinearwave) then
            w_wave=rhs/rho*wavedampingtime*86400. !WTG
         else
            w_wave=w_wave*(1.-float(nsteplinearwave)*dt*sqrt(betafactor)/86400./wavedampingtime*dampingprofile)+rhs/rho*dt*nsteplinearwave !Evolving wave
         end if

      t_wave_local=0.
      q_wave_local=0.
   endif
end if
end if


end subroutine linear_wave
