subroutine nudging
	
use vars
use params
implicit none

real coef
integer i,j,k
	
coef = 1./tauls

if(donudging_uv) then

    coef = 1. / tauls
    do k=1,nzm
      do j=1,ny
        do i=1,nx
           dudt(i,j,k,na)=dudt(i,j,k,na)-(u0(k)-ug0(k))*coef
           dvdt(i,j,k,na)=dvdt(i,j,k,na)-(v0(k)-vg0(k))*coef
        end do
      end do
    end do

endif


if((donudging_tq).and.(.not.domoistadiabaticnudging).and.(.not.dotauz)) then

    coef = dtn / tauls * sqrt(betafactor)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
           t(i,j,k)=t(i,j,k)-(t0(k)-tg0(k)-gamaz(k))*coef
           q(i,j,k)=q(i,j,k)-(q0(k)-qg0(k))*coef
           if(dokruegermicro) qv(i,j,k)=qv(i,j,k)-(q0(k)-qg0(k))*coef
        end do
      end do
      tnudge(k) = (t0(k)-tg0(k)-gamaz(k))*coef/dtn
      qnudge(k) = (q0(k)-qg0(k))*coef/dtn
    end do

endif

!Nov.28,05, nudging q only to a factor times the original profile
if((donudging_q)) then

   do k=1,nzm
      coef = dtn * sqrt(betafactor) / tauls
!      coef = dtn * sqrt(betafactor) / (tauls + &
!           1e10 * min(1.,1.-(min( (z(k) - 2000.) / (3000. - 2000.),1.))))
      do j=1,ny
        do i=1,nx
            q(i,j,k)=q(i,j,k)-(q0(k)-qg0(k)*nudging_q_factor)*coef
           if(dokruegermicro) qv(i,j,k)=qv(i,j,k)-(q0(k)-qg0(k)*nudging_q_factor)*coef
        end do
      end do
      qnudge(k) = (q0(k)-qg0(k)*nudging_q_factor)*coef/dtn
    end do

endif

!=========================================================================
! UW ADDITIONS

if((donudging_tq).and.((domoistadiabaticnudging).or.(dotauz))) then

   ! Initialize nudging tendencies:
   do k=1,nzm
      tnudge(k)=0.
      qnudge(k)=0.
   end do

   if (domoistadiabaticnudging) then

      call moist_adiabatic_nudging()

   elseif (dotauz) then

      do k=1,nzm
         tnudge(k)= -(t0(k)-tg0(k)-gamaz(k)) * itauz(k) * sqrt(betafactor)
         qnudge(k)= -(q0(k)-qg0(k)) * itauz(k) * sqrt(betafactor)
      end do

   end if

   do k=1,nzm
      do j=1,ny
         do i=1,nx
            t(i,j,k)=t(i,j,k) + tnudge(k) * dtn 
            q(i,j,k)=q(i,j,k) + qnudge(k) * dtn
           if(dokruegermicro) qv(i,j,k)=qv(i,j,k) + qnudge(k) * dtn
         end do
      end do
   end do

end if

!kzm: nudge scalar profiles

!kzm: always do nudging for tro3 if dotro3 ( it is more as a decay term)
if (dotro3) then 
    coef = dtn * itauo3
    do k=1,nzm
      do j=1,ny
        do i=1,nx
           tro3(i,j,k)=tro3(i,j,k)-(tro3(i,j,k)-o3g0(k))*coef
        end do
      end do
    end do
endif

if (dotrz) then 
    coef = dtn * itautrz
    do k=1,nzm
      do j=1,ny
        do i=1,nx
           trz(i,j,k)=trz(i,j,k)-trz(i,j,k)*coef
        end do
      end do
    end do
endif

if (dotrx) then 
    coef = dtn * itautrx
    do k=1,nzm
      do j=1,ny
        do i=1,nx
           trx(i,j,k)=trx(i,j,k)-trx(i,j,k)*coef
        end do
      end do
    end do
endif

if (dotry) then 
    coef = dtn * itautry
    do k=1,nzm
      do j=1,ny
        do i=1,nx
           try(i,j,k)=try(i,j,k)-try(i,j,k)*coef
        end do
      end do
    end do
endif

if (dotrzz) then 
    coef = dtn * itautrzz
    do k=1,nzm
      do j=1,ny
        do i=1,nx
           trzz(i,j,k)=trzz(i,j,k)-trzz(i,j,k)*coef
        end do
      end do
    end do
endif

! END UW ADDITIONS
!=========================================================================

end subroutine nudging

!=========================================================================
! UW ADDITION

subroutine moist_adiabatic_nudging()

use vars
use params
implicit none

real, parameter :: press_nudge = 150. 
real tdrift, tdrift_mean, tnorm, tau_drift, tau_nudge, eta, dhsdT, dhdT
real tmp1, tmp2
integer i,j,k

! Initialize the nudging to balance the mean drift in the temperature
! profile.  The nudging will attempt to introduce a uniform-with
!-height perturbation in the moist static energy while leaving the
! horizontally-averaged relative humidity fixed (or at least
! unaffected by the nudging).

if (pres(1).ne.0.) then
   ! Compute the drift of the temperature profile relative to the
   ! sounding averaged over the domain (both horizontally and
   ! vertically up to press_nudge)
   tdrift = 0.0
   tnorm = 0.0
   do k = 1,nzm
      if (pres(k).ge.press_nudge) then
         tdrift = tdrift + dz*adz(k)*rho(k) &
              *(tg0(k) - tabs0(k) + fac_cond*(qg0(k) - qv0(k)))
         tnorm  = tnorm  + dz*adz(k)*rho(k)
         dhsdT  = dhsdT  + dz*adz(k)*rho(k)*(cp + lcond*lcond &
              *qsatw(tabs0(k),pres(k))/(rv*tabs0(k)**2))
         dhdT   = dhdT   + dz*adz(k)*rho(k)*(cp + lcond*lcond &
              *qv0(k)/(rv*tabs0(k)**2))
      end if
   end do
   ! Compute mass-weighted, column-averaged temperature drift
   ! below press_nudge (nominally 150 millibar).
   tdrift = tdrift/tnorm
   ! Compute the ratio of the partial derivatives wrt temperature
   ! of saturated moist static energy to moist static energy.
   eta    = dhsdT/dhdT 

   ! The relaxation time scale for the forcing to counteract drift
   ! in the mean temperature profile is two days.  The nudging
   ! timescale applied above press_nudge is tauls, the large-scale
   ! forcing timescale specified in the prm file.
   tau_drift = 48.*3600.
   tau_nudge = tauls

   do k=1,nzm
      if (pres(k).ge.press_nudge) then
         ! Moist static energy nudging is mean drift divided by
         ! tau_drift, i.e. cp*tdrift/tau_drift.  Temperature
         ! nudging is computed to give a uniform nudging to the
         ! moist static energy. 
         tnudge(k)=(cp*tdrift/tau_drift)*eta/(cp+lcond*lcond&
              *qsatw(tabs0(k),pres(k))/(rv*tabs0(k)**2))
         ! Apply qnudge so that it doesn't change relative humidity 
         qnudge(k)=qv0(k)*tnudge(k)*lcond/(rv*tabs0(k)**2) 
      else
         ! Above press_nudge:
         !   - Nudge t profiles back to sounding
         !   - Do not nudge q profile to avoid damping high cloud
         tnudge(k) = (tg0(k) - tabs0(k))/tau_nudge
         qnudge(k) = 0.
      end if
   end do

end if

end subroutine moist_adiabatic_nudging

! END UW ADDITION
!=========================================================================
