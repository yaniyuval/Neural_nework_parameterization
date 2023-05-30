
subroutine precip_proc_krueger

  use vars
  use params

  implicit none

  integer i,j,k, ii, kk
  real qmixz(6,nx*ny), tabsz(nx*ny), pssz(35), rhoz(nx*ny), presz(nx*ny)
  real tssz(7,nx*ny), dtscale, dq

  real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
  real f0(nzm),df0(nzm)

  real, parameter :: rhoref = 1.226

  if(dostatis) then

     do k=1,nzm
        do j=dimy1_s,dimy2_s
           do i=dimx1_s,dimx2_s
              df(i,j,k) = qv(i,j,k) + qc(i,j,k) + qi(i,j,k)
           end do
        end do
     end do

  endif

  if(dobetafactor) dtn=dtn*betafactor**0.5

  call setupm(dtn,rhoref)

  ! Outer loop over vertical layers on this processor.
  do k = 1,nzm

     ii = 0
     do j = 1,ny
        do i = 1,nx
           ii = ii + 1
           qmixz(1,ii) = qv(i,j,k)
           qmixz(2,ii) = qc(i,j,k)
           qmixz(3,ii) = qi(i,j,k)
           qmixz(4,ii) = qs(i,j,k)
           qmixz(5,ii) = qg(i,j,k)
           qmixz(6,ii) = qr(i,j,k)
           tabsz(ii)   = tabs(i,j,k)
           presz(ii)   = pres(k)
           rhoz(ii)    = rho(k)
        end do
     end do

     call imicro_driver(nx*ny,tabsz,presz,qmixz,tssz,pssz,rhoz,dostatis)

     ! update microphysical quantities using simple forward Euler scheme.
     ! TODO(???): Use a second-order rk scheme here, alternating
     !       calls to imicro() and sat() within each substep.
     !       Would require re-setting smax within imicro() since the
     !       maximum change in qmix over one substep would depend on the
     !       rhs from the previous substep.
     ii = 0
     pladjz(k) = 0.
     piadjz(k) = 0.
     pvaporz(k)   = 0.
     pclwz(k)     = 0.
     pcliz(k)     = 0.
     psnowz(k)    = 0.
     pgraupelz(k) = 0.
     prainz(k)    = 0.

     do j = 1,ny
        do i = 1,nx
           ii = ii + 1
           ! account for changes to vapor and individual hydrometeors
           ! Note that the microphysics subroutine has limited the
           ! size of tssz, so that qc cannot become smaller than zero.
           ! We set the microphysical quantities to be stricly non
           !-negative here because rounding error could give negative
           ! values.
           qv(i,j,k) = max(qv(i,j,k) + dtn*tssz(1,ii),0.)
           qc(i,j,k) = max(qc(i,j,k) + dtn*tssz(2,ii),0.)
           qi(i,j,k) = max(qi(i,j,k) + dtn*tssz(3,ii),0.)
           qs(i,j,k) = max(qs(i,j,k) + dtn*tssz(4,ii),0.)
           qg(i,j,k) = max(qg(i,j,k) + dtn*tssz(5,ii),0.)
           qr(i,j,k) = max(qr(i,j,k) + dtn*tssz(6,ii),0.)
           
           ! evaporate remaining precipitation if less then qp_threshold
           if ((qr(i,j,k).gt.0).and.(qr(i,j,k).lt.qp_threshold)) then
              dq = qr(i,j,k)
              qv(i,j,k) = qv(i,j,k) + dq
              qr(i,j,k) = 0.
              pssz(9) = pssz(9) + dq/dtn ! stats: rain evaporation
           end if
           if ((qg(i,j,k).gt.0).and.(qg(i,j,k).lt.qp_threshold)) then
              dq = qg(i,j,k)
              qv(i,j,k) = qv(i,j,k) + dq
              qg(i,j,k) = 0.
              pssz(22) = pssz(22) + dq/dtn ! stats: graupel sublimation
           end if
           if ((qs(i,j,k).gt.0).and.(qs(i,j,k).lt.qp_threshold)) then
              dq = qs(i,j,k)
              qv(i,j,k) = qv(i,j,k) + dq
              qs(i,j,k) = 0.
              pssz(21) = pssz(21) + dq/dtn ! stats: snow sublimation
           end if
              

           if (min(qv(i,j,k),qc(i,j,k),qi(i,j,k),qr(i,j,k), &
                qg(i,j,k),qs(i,j,k)).lt.0.) then
              write(*,*) nstep, i,j,k, tabs(i,j,k), q(i,j,k), &
                   qv(i,j,k),qc(i,j,k),qi(i,j,k),qr(i,j,k), &
                   qg(i,j,k),qs(i,j,k)
              write(*,*) 'error after imicro()'
              call task_abort()
           end if

           ! re-compute absolute temperature and water vapor mixture fraction
           ! using new hydrometeor levels
!!$           q(i,j,k) = qv(i,j,k) + qc(i,j,k) + qi(i,j,k) ! now in main.f90
           tabs(i,j,k) = t(i,j,k) - gamaz(k) &
                + fac_cond*(qc(i,j,k)+qr(i,j,k)) &
                + fac_sub*(qi(i,j,k)+qg(i,j,k)+qs(i,j,k))

           pladjz(k) = pladjz(k) - qc(i,j,k)
           piadjz(k) = piadjz(k) - qi(i,j,k)
           pvaporz(k)   = pvaporz(k)   + tssz(1,ii)
           pclwz(k)     = pclwz(k)     + tssz(2,ii)
           pcliz(k)     = pcliz(k)     + tssz(3,ii)
           psnowz(k)    = psnowz(k)    + tssz(4,ii)
           pgraupelz(k) = pgraupelz(k) + tssz(5,ii)
           prainz(k)    = prainz(k)    + tssz(6,ii)

           ! call saturation adjustment routine
           call sat(i,j,tabs(i,j,k),pres(k),qv(i,j,k),qc(i,j,k),qi(i,j,k))

           if (min(qv(i,j,k),qc(i,j,k),qi(i,j,k)).lt.0.) then
              write(*,*) nstep, i,j,k, tabs(i,j,k), q(i,j,k), &
                   qv(i,j,k),qc(i,j,k),qi(i,j,k),qr(i,j,k), &
                   qg(i,j,k),qs(i,j,k)
              write(*,*) 'error after SAT()'
              call task_abort()
           end if

           pladjz(k) = pladjz(k) + qc(i,j,k)
           piadjz(k) = piadjz(k) + qi(i,j,k)

           ! update absolute temperature based on new qc, qi
           tabs(i,j,k) = t(i,j,k) - gamaz(k) &
                + fac_cond*(qc(i,j,k)+qr(i,j,k)) &
                + fac_sub*(qi(i,j,k)+qg(i,j,k)+qs(i,j,k))
!!$           NOW IN main.f90
!!$           ! update total cloud condensate (qn) and total precipitate (qp)
!!$           qn(i,j,k) = qc(i,j,k) + qi(i,j,k)
!!$           qp(i,j,k) = qr(i,j,k) + qs(i,j,k) + qg(i,j,k)
!!$           ! Note that total water doesn't change through sat adjustment.
        end do
     end do


     ! accumulate statistics for microphysical quantities.
     if (dostatis) then
        dtscale = sqrt(betafactor)

        pladjz(k) = pladjz(k)*dtscale/dtn
        piadjz(k) = piadjz(k)*dtscale/dtn

        pvaporz(k)   = pvaporz(k)   *dtscale/dtn
        pclwz(k)     = pclwz(k)     *dtscale/dtn
        pcliz(k)     = pcliz(k)     *dtscale/dtn
        psnowz(k)    = psnowz(k)    *dtscale/dtn
        pgraupelz(k) = pgraupelz(k) *dtscale/dtn
        prainz(k)    = prainz(k)    *dtscale/dtn

        prautz(k) = pssz(1)*dtscale
        pracwz(k) = pssz(2)*dtscale
        pracsz(k) = pssz(3)*dtscale
        psacwz(k) = pssz(4)*dtscale
        !           pwacsz(k) = pssz(5)*dtscale ! set to zero in IMICRO()
        pgacwz(k) = pssz(6)*dtscale
        pgmltz(k) = pssz(7)*dtscale
        psmltz(k) = pssz(8)*dtscale
        prevpz(k) = pssz(9)*dtscale
        piacrz(k) = pssz(10)*dtscale
        psacrz(k) = pssz(11)*dtscale
        pgacrz(k) = pssz(12)*dtscale
        pgfrz(k)  = pssz(13)*dtscale
        pgacsz(k) = pssz(14)*dtscale
        pgautz(k) = pssz(15)*dtscale
        pgaciz(k) = pssz(16)*dtscale
        pgwordz(k) = pssz(17)*dtscale ! wet or dry graupel processes
        praciz(k) = pssz(18)*dtscale
        psautz(k) = pssz(19)*dtscale
        psaciz(k) = pssz(20)*dtscale
        pssubz(k) = pssz(21)*dtscale
        pgsubz(k) = pssz(22)*dtscale
        psfwz(k)  = pssz(23)*dtscale
        psfiz(k)  = pssz(24)*dtscale
        pidwz(k)  = pssz(25)*dtscale
        pgwetz(k) = pssz(26)*dtscale

        vtrz(k)     = pssz(27) ! rain terminal velocity avg. over whole domain
        rainfrac(k) = pssz(28) ! rain area fraction

        vtgz(k)     = pssz(29) ! same for graupel
        graufrac(k) = pssz(30)

        vtsz(k)     = pssz(31) ! same for snow
        snowfrac(k) = pssz(32)

        vtiz(k)     = pssz(33) ! same for ice
        clifrac(k)  = pssz(34)

        clwfrac(k)  = pssz(35) ! cloud water area fraction
        
        ! bloss: precipitation process diagnostics output by WRF Lin scheme
        !        but not by Kreuger's version of the Lord scheme.
        !           pimltz = pimltz + pimlt
        !           pihomz = pihomz + pihom
        !           psdepz = psdepz + psdep
        !           psmltevpz = psmltevpz + psmltevp
        !           pgacipz = pgacipz + pgacip
        !           pgacrpz = pgacrpz + pgacrp
        !           pgacspz = pgacspz + pgacsp
        !           pdryz = pdryz + pdry
        !           pgdepz = pgdepz + pgdep
        !           pgmltevpz = pgmltevpz + pgmltevp

     end if

  end do

  if(dobetafactor) dtn=dtn/betafactor**0.5

  if(dostatis) then

     call stat_varscalar(qv+qc+qi,df,f0,df0,q2leprec)
     call setvalue(qwleprec,nzm,0.)
     call stat_sw2(qv+qc+qi,df,qwleprec)

  endif

end subroutine precip_proc_krueger

