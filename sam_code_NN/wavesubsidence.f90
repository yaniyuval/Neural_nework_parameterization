
subroutine wavesubsidence()
	
use vars
implicit none

integer i,j,k,k1,k2,k3,k4
real rdz,coef,rdz2

real t_tend(nzm), q_tend(nzm)

!As in subsidence.f90, here the advection is upwind 
if(nstep.gt.nstartlinearwave+nsteplinearwavebg) then

      ttend_wave=0.
      qtend_wave=0.
      
      q_tend=0.
      t_tend=0.

      do k=1,nzm-1
         rdz=w_wave(k)/(dz*adzw(k)) 
         if (rdz.ge.0) then
            k2 = max(1,k-1) ! deals with the end points
            k1 = k2+1
            rdz=w_wave(k)/(dz*adzw(k1))
         else
            k1 = min(nzm-1,k+1)  ! deals with the end points
            k2 = k1-1
            rdz=w_wave(k)/(dz*adzw(k2))
         end if
            
         if (.not.dowavesubsidence_tonly) then
            ! apply subsidence velocity to all fields (q,T)
            if(dosubsidenceonbackground) then 
               t_tend(k) =  - rdz * max((t_wavebg(k1)+gamaz(k1)-t_wavebg(k2)-gamaz(k2)),5.e-4*dz*adzw(k))
               q_tend(k) =  - rdz * (q_wavebg(k1)-q_wavebg(k2))
            else
               t_tend(k) =  - rdz * max((t_wave(k1)+gamaz(k1)-t_wave(k2)-gamaz(k2)),5.e-4*dz*adzw(k))
               q_tend(k) =   - rdz *(q_wave(k1)-q_wave(k2))
            end if
         else
            !kzm limit subsidence to temperature so it is solely a dynamical
            ! heating term
            t_tend(k) =  - rdz * max((t_wavebg(k1)+gamaz(k1)-t_wavebg(k2)-gamaz(k2)),5.e-4*dz*adzw(k))
         endif
      end do
      
      t_tend=t_tend*sqrt(betafactor)
      q_tend=q_tend*sqrt(betafactor)
      
      !   endif
      if(dowavedampingt) then 
         do k=1,nzm
            do j=1,ny
               do i=1,nx
                  t(i,j,k) = t(i,j,k) + dtn * (t_tend(k)-(t_wave(k)-t_wavebg(k))/86400./wavedampingtime)
                  q(i,j,k) = q(i,j,k) + dtn * (q_tend(k)-(q_wave(k)-q_wavebg(k))/86400./wavedampingtime)
                  q(i,j,k) = max(0.,q(i,j,k))
               end do
            end do
         end do
      else
         do k=1,nzm
            do j=1,ny
               do i=1,nx
                  t(i,j,k) = t(i,j,k) + dtn * t_tend(k)
                  q(i,j,k) = q(i,j,k) + dtn * q_tend(k) 
                  q(i,j,k) = max(0.,q(i,j,k))
               end do
            end do
         end do
      end if
      if(icycle.eq.ncycle) then
         ttend_wave = t_tend 
         qtend_wave = q_tend 
      endif

endif

end
