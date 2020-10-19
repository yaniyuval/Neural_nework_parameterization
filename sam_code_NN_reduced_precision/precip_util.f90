subroutine precip_iter(vterm,nk,dt_in,deltaz,rhofac,tabs,rq,totalflux)

  ! subroutine which computes precipitation of hydrometeors in two
  ! steps.  In the first step, the precipitation velocity is computed
  ! by calling the function that is passed in by the argument vterm.
  ! This function should take three real arguments:
  !
  !    vterm(rq(k),rhofac(k),tabs(k))
  !
  ! where rq = rho*qm with qm the hydrometeor mass mixture fraction, 
  !       rhofac = sqrt(rho_reference/rho), and
  !       tabs is a dummy argument (here temperature) which can be
  !               used to help compute the terminal velocity.
  !
  ! In the second step, the wave speed of the precipitation is
  ! computed and the time step is reduced if the
  ! maximum cfl of the precipitation velocity exceeds 0.9.  The
  ! precipitation fluxes are computed and rhoq are updated in a
  ! separate routine.  If the time step was reduced, these two steps
  ! are repeated until the full duration of the time step has been
  ! completed.
  !
  ! written by Peter Blossey (UW) September--October 2004.
  !
  ! Note: The form of this subroutine was inspired by the WRF module 
  ! which implements the Lin scheme.

  implicit none

  ! in/outputs
  real, external                     :: vterm
  integer, intent(in)                :: nk
  real, intent(in)                   :: dt_in
  real, dimension(nk), intent(in)    :: deltaz, rhofac, tabs
  real, dimension(nk), intent(inout) :: rq
  real, dimension(nk+1), intent(out) :: totalflux

  !local variables
  integer k, kb, mink, maxk
  real    :: dt_step, prec_cfl, dt_remaining
  real, dimension(nk)   :: wp, s
  real, dimension(nk+1) :: flux
  real, parameter :: cfl_limit = 0.9
  logical :: not_done

  not_done = .true.
  dt_remaining = dt_in
  totalflux = 0.

  do while (not_done)

     do k = 1,nk
        kb = max(1,k-1)
        ! compute precipitation velocity
        wp(k) = vterm(rq(k),rhofac(k),tabs(k))
        ! compute wave speed 
        if (rq(kb).eq.rq(k)) then
           s(k) = abs(wp(k))
        else
           s(k) = abs((rq(k)*wp(k) - rq(kb)*wp(kb))/(rq(k) - rq(kb)))
        end if
     end do
     prec_cfl = MAXVAL(s/deltaz)*dt_remaining
     if (prec_cfl.gt.cfl_limit) then
        dt_step = (cfl_limit/prec_cfl)*dt_remaining
        dt_remaining = dt_remaining - dt_step
     else
        dt_step = dt_remaining
        not_done = .false.
     end if

     call advect_precip1D(nk,rq,-wp,s,flux,deltaz,dt_step)

     totalflux = totalflux + flux*(dt_step/dt_in) !accumulate fluxes for stats
  end do

end subroutine precip_iter

subroutine advect_precip1D(nk,rq,w,s,flux,dz,dt) 

  ! positive-definite advection scheme to compute falling of precipitation 
  ! in a single vertical column.  Relies on the input of 
  !    - a CELL-CENTERED precipitation fall velocity 
  !         (which should not exceed a cfl of 1.0) and
  !    - the hydrometeor density (rho*qm) where m indicates the type of precip.
  ! This subroutine returns the hydrometeor density profile after a single time
  ! step of length dt.

  implicit none

  ! inputs
  integer, intent(in) :: nk       ! number of points in vertical direction
  real, intent(in)    :: w(nk)    ! cell-centered precipitation fall velocity
                                  ! ASSUMED NEGATIVE DEFINITE
  real, intent(in)    :: s(nk)    ! wave speed at bottom interface of cell k
  real, intent(in)    :: dz(nk)   ! grid size (in same length units as in w)
  real, intent(in)    :: dt       ! time step (in same time units as in w)

  ! in/outputs
  real, intent(inout) :: rq(nk) ! cell-centered hydrometeor density 
                                  ! (e.g. rho*qr for rain)

  real, intent(out)   :: flux(nk+1) ! Fluxes across cell faces

  ! local variables
  real, dimension(nk)    :: tmp_rq, mx, mn, cfl
  real, dimension(nk+1)  :: adiff, fcorr
  real    :: rqmx, rqmn, rqu, rqc, rqd, tmp_phi, tmp_theta, waveup, wave
  integer :: k, kb, kc

!!$  if (MAXVAL(w).gt.0.) then
!!$     write(*,*) 'positive precipitation velocity: aborting in advect_precip'
!!$!     call task_abort()
!!$     STOP 'in advect_precip'
!!$  end if

  ! compute upwind flux -- assumes all entries of wp are negative or zero.
  flux(1:nk) = rq(1:nk)*w(1:nk)

  ! apply second-order correction to flux using a cfl based on the
  ! wavespeeds computed above.  Note that there is no flux correction
  ! at the top or bottom of the domain since we choose to use
  ! extrapolation boundary conditions for outflow/inflow: rq(0)=rq(1)
  ! and rq(nk+1)=rq(nk).
  do k = 2,nk-1
     if (rq(k-1).ne.rq(k)) then
        waveup = rq(k+1)-rq(k)
        wave   = rq(k)-rq(k-1)
        tmp_theta = (wave*waveup)/(wave*wave)
        tmp_phi = max(0.,min(0.5*(1.+tmp_theta),2.,2.*tmp_theta))
        flux(k) = flux(k) &
             + 0.5*s(k)*(1.-s(k)*dt/dz(k))*tmp_phi*(rq(k) - rq(k-1))
     end if
  end do
  ! flux is zero at top of domain.
  flux(nk+1) = 0.

  ! update solution using new flux
  rq = rq - (flux(2:nk+1) - flux(1:nk))*dt/dz

end subroutine advect_precip1D
