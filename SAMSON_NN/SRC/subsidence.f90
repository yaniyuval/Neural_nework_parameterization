
subroutine subsidence()
	
use vars
implicit none

integer i,j,k,k1,k2
real rdz,coef
real t_vtend, q_vtend
real t_tend(nx,ny,nzm), q_tend(nx,ny,nzm)

! Initialize large-scale vertical advective tendencies.
do k = 1,nzm
   qlsvadv(k) = 0.
   tlsvadv(k) = 0.
end do

if (doparameterizedwave.and.nstep.gt.nstartlinearwave+nsteplinearwavebg.and.doadvectbg) then
   if(dozerowsub) wsub=0.
   do k=2,nzm-1
      if(wsub(k).ge.0) then
         rdz=wsub(k)
         k1 = k
         k2 = k-1 
      else
         rdz=wsub(k)
         k1 = k+1
         k2 = k
      end if
      if (.not.dosubsidence_tonly) then
         ! apply subsidence velocity to all fields (u,v,q,T)
         if (dosubsidenceonbackground) then
            do j=1,ny
               do i=1,nx
                  t_tend(i,j,k) =  - rdz * max((t_wavebg(k1)+gamaz(k1)-t_wavebg(k2)-gamaz(k2))/(dz*adzw(k1)),5.e-4)
                  q_tend(i,j,k) =  - rdz * (q_wavebg(k1)-q_wavebg(k2))/(dz*adzw(k1)) 
               end do
            end do
         else ! JAA advect the actual (mean) state, rather than background profiles
            do j=1,ny
               do i=1,nx
                  t_tend(i,j,k) =  - rdz * max((t_wave(k1)+gamaz(k1)-t_wave(k2)-gamaz(k2))/(dz*adzw(k1)),5.e-4)
                  q_tend(i,j,k) =  - rdz * (q_wave(k1)-q_wave(k2))/(dz*adzw(k1)) 
               end do
            end do
         end if
      else
         !kzm limit subsidence to temperature so it is solely a dynamical
         ! heating term
         do j=1,ny
            do i=1,nx
               t_tend(i,j,k) =  - rdz * max((t_wavebg(k1)+gamaz(k1)-t_wavebg(k2)-gamaz(k2))/(dz*adzw(k1)),5.e-4)
            end do
         end do
      end if
   end do
else
   coef=1.
   
   if(doperiodic_forcing) &
        coef=sqrt(betafactor)*(1-cos(6.28319*max(0.,(day-day0_forcing/sqrt(betafactor)))/period_forcing*sqrt(betafactor)))/2.

   do k=2,nzm-1
      if(wsub(k).ge.0) then
         rdz=coef*wsub(k)/(dz*adzw(k))	
         k1 = k
         k2 = k-1 
      else
         rdz=coef*wsub(k)/(dz*adzw(k+1))       
         k1 = k+1
         k2 = k
      end if
      if (.not.dosubsidence_tonly) then
         ! apply subsidence velocity to all fields (u,v,q,T)
         do j=1,ny
            do i=1,nx
               dudt(i,j,k,na) = dudt(i,j,k,na) - rdz*(u(i,j,k1)-u(i,j,k2)) 
               dvdt(i,j,k,na) = dvdt(i,j,k,na) - rdz*(v(i,j,k1)-v(i,j,k2)) 
               t_tend(i,j,k) =  - rdz * (t(i,j,k1)-t(i,j,k2))
               q_tend(i,j,k) =  - rdz * (q(i,j,k1)-q(i,j,k2))
            end do
         end do
      else
         !kzm limit subsidence to temperature so it is solely a dynamical
     ! heating term
         do j=1,ny
            do i=1,nx
               t_tend(i,j,k) =  - rdz * (t(i,j,k1)-t(i,j,k2))
            end do
         end do
      end if
   end do
endif

do k=2,nzm-1
  t_vtend = 0.
  q_vtend = 0.
  do j=1,ny
    do i=1,nx
      t(i,j,k) = t(i,j,k) + dtn * t_tend(i,j,k)
      q(i,j,k) = q(i,j,k) + dtn * q_tend(i,j,k)
      q(i,j,k) = max(0.,q(i,j,k))
      t_vtend = t_vtend + t_tend(i,j,k)
      q_vtend = q_vtend + q_tend(i,j,k)
    end do
  end do
  t_vtend = t_vtend / float(nx*ny) 
  q_vtend = q_vtend / float(nx*ny) 
  ttend(k) = ttend(k) + t_vtend
  qtend(k) = qtend(k) + q_vtend
  tlsvadv(k) = t_vtend
  qlsvadv(k) = q_vtend
end do 
	
end
