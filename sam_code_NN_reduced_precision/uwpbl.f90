module uwpbl
  !---------------------------------------------------------------------------------
  ! Module to compute eddy diffusion coefficients associated with turbulence in the 
  ! planetary boundary layer and elsewhere.
  !
  ! Public interfaces
  !
  !    init_eddy_diff       initialize time independent coefficients
  !    compute_eddy_diff    compute diffusivities based on current state
  !    
  !    vdiff_eddy
  !      trbintr        interface for vertical diffusion and pbl scheme and write output
  !          trbintd       initializes time dependent variables
  !          vd_k_freeatm  computes mixing coefficients for free atmosphere
  !          vd_k_pbl      computes mixing coefficients for pbl
  !          caleddy       computes all turbulent diffusivities
  !             exacol     classifies turbulent layers in each atmospheric column 
  !                        using moist Ri
  !             zisocl     determines breadth and layer-mean TKE of convective layers.
  !---------------------------Code history------------------------------------------
  ! updated:           J. McCaa, December 2003
  !---------------------------------------------------------------------------------
  use uwcu, only: mcshallow

  use diffusion_solver, only: vdiff_selector
  
!!$  added by cw 8/15/05
  use wv_saturation, only: gestbl, fqsatd

!!$  added by cw 2/16/05
  use vars
  use params
  use grid, only : maxtkz, maxcfl, fixedndiff

  implicit none
  private
  save

  public compute_eddy_diff
  public uwpbldrv

!  integer, parameter :: r8 = selected_real_kind(6,30) ! 4 byte real
  integer, parameter :: r8 = selected_real_kind(6,70) ! 8 byte real
  ! added by cw 
  integer, parameter :: ncnst = 3   ! number of moisture species to be diffused

  real(r8), dimension(:, :), allocatable ::                                   &
       cldn_flip               ! cloud fraction

  real(r8), dimension(:, :), allocatable ::                                   &
       qrl_flip,           &   ! LW heating rate (K/s) * cp * dp (W/kg*Pa)
       s_flip,             &   ! dry static energy
       rpdel_flip,         &   ! 1/dp
       z_flip,             &   ! half-level heights above surface (m)
       qv_sh_flip,         &   ! water vapor specific humidity (kg/kg)
       qc_sh_flip,         &   ! cloud water specific humidity (kg/kg)
       qi_sh_flip,         &   ! cloud ice specific humidity (kg/kg)
       qp_sh_flip,         &   ! precipitating water specific humidity (kg/kg)
       q_sh_flip,          &   ! total water specific humidity (kg/kg)
       tabs_flip,          &   ! temperature with index kts at surface
       u_flip,             &   ! zonal wind with index kts at surface
       v_flip,             &   ! meridional wind with index kts at surface
       w_flip,             &   ! vertical wind with index kts at surface
       p_flip,             &   ! half-level pressure
       dtk_flip,           &   ! T tendency from KE dissipation (not used)
       p8w_flip_temp,      &   ! full-level pressures (Pa)
       p8w_flip,           &   ! full-level pressures (Pa)
       z8w_flip,           &   ! full-level heights above surface (m)
       bprod_flip,         &   ! buoyancy production
       sprod_flip,         &   ! shear production
       sfi_flip,           &   ! interfacial saturation fraction
       tke_flip,           &   ! turbulent kinetic energy
       cgh_flip,           &   ! counter-gradient term for heat 
       cgs_flip,           &   ! counter-gradient star
       kvm_in_flip,        &   ! eddy diffusivity (mom), last timestep (m2/s)
       kvh_in_flip,        &   ! eddy diffusivity (heat), last timestep (m2/s)
       kvm_out_flip,       &   ! eddy diffusivity (mom)  (m2/s)
       kvh_out_flip,       &   ! eddy diffusivity (heat) (m2/s)
       kvq_out_flip            ! eddy diffusivity (constituents)  (m2/s)

  real(r8), dimension(:), allocatable ::                                      &
       taux_local,         &   ! zonal stress (N/m2)
       tauy_local,         &   ! meridional stress (N/m2)
       ksrftms,            &   ! sfc "drag" coeff. for mountains (set to zero)
       tautmsx,            &   ! zonal mountain stress (not used)
       tautmsy,            &   ! meridional mountain stress (not used)
       topflx,             &   ! molecular heat flux at TOA (not used)
       tpert,              &   ! convective temperature excess
       qpert,              &   ! convective humidity excess
       ustar_local,        &   ! one-dimensional ustar 
       sfx_local,          &   ! one-dimensional surface heat flux
       pblh_local              ! one-dimensional PBL height

  real(r8), dimension(:, :, :), allocatable ::                                &
       q_local                 ! variable with all moisture species

  real(r8), dimension(:, :), allocatable ::                                   &
       qfx_local               ! surface fluxes of water species

  real(r8), dimension(ncnst) ::                                               &
       qmincg                  ! minimum water vapor mixing ratio from cg flux

  integer ::                                                                  &
       nx_p                    ! number of patches in x direction
  !
  ! PBL limits
  !
  real(r8), parameter :: ustar_min = 0.01        ! min permitted value of ustar
  real(r8), parameter :: pblmaxp   = 4.e4        ! pbl max depth in pressure units
  real(r8), parameter :: zkmin     = 0.01        ! Minimum kneutral*f(ri)
  !
  ! PBL Parameters
  !
  real(r8), parameter :: onet  = 1./3. ! 1/3 power in wind gradient expression
  real(r8), parameter :: betam = 15.0  ! Constant in wind gradient expression
  real(r8), parameter :: betas =  5.0  ! Constant in surface layer gradient expression
  real(r8), parameter :: betah = 15.0  ! Constant in temperature gradient expression 
  real(r8), parameter :: fakn  =  7.2  ! Constant in turbulent prandtl number
  real(r8), parameter :: fak   =  8.5  ! Constant in surface temperature excess         
  real(r8), parameter :: ricr  =  0.3  ! Critical richardson number
  real(r8), parameter :: sffrac=  0.1  ! Surface layer fraction of boundary layer
  real(r8), parameter :: binm  = betam*sffrac       ! betam * sffrac
  real(r8), parameter :: binh  = betah*sffrac       ! betah * sffrac
  !
  ! Pbl constants set using values from other parts of code
  !
  real(r8) :: zvir       ! rh2o/rair - 1

  integer  :: ntop_turb  ! Top level to which turbulent vertical diffusion is applied.
  integer  :: nbot_turb  ! Bottom level to which turbulent vertical diff is applied.

  real(r8), allocatable :: ml2(:) ! Mixing lengths squared

  logical, parameter :: use_kvf = .false.  ! T => initialize kvh/kvm =  kvf
                                           ! F => initialize kvh/kvm =  0.

  !cb caleddy parameters ---------- (note: vk, fak from above are also used)
  character, parameter :: sftype  = 'z'   ! method for calculating sat frac
  integer, parameter :: ncvmax = 100      ! Max no. of separate convective layers

  real(r8), parameter :: qmin  = 1.e-5    ! Min liq water (kg/kg) counted as cld
  real(r8), parameter :: ntzero = 1.e-12  ! not zero (small positive number)
  real(r8), parameter :: b1  = 5.8        ! TKE dissipation D = e^3/(b1*leng)
  real(r8) :: b123                        ! b1**(2/3)

  real(r8), parameter :: tunl = 0.085     ! Asympt leng = tunl*(turb lay depth)
  real(r8), parameter :: alph1 = 0.5562   ! Galperin stability fn params
  real(r8), parameter :: alph2 = -4.3640
  real(r8), parameter :: alph3 = -34.6764
  real(r8), parameter :: alph4 = -6.1272
  real(r8), parameter :: alph5 = 0.6986
  real(r8), parameter :: ricrit = 0.19    ! Critical Richardson no. for turbulence
  real(r8), parameter :: mu = 70.         ! used in finding e/ebar in conv. layers
  real(r8), parameter :: rinc = -0.5      ! Min W/<W> for incorp into conv layer

  ! params governing entr. efficiency A = a1l*evhc, evhc=1+a2l*a3l*L*ql/jt2slv
  ! where ql is cloud-top liquid water and jt2slv is the jump in slv across
  ! the cloud-top entrainment zone.

  real(r8), parameter :: a2l = 15.       ! Moist entrainment enhancement param (recommended range : 10~20 )
  ! ...tunable within limits
  real(r8), parameter :: a3l = 0.8       ! Approximation to a complicated
  ! thermo. parameter (not for tuning!)
  real(r8), parameter :: jbumin = .001   ! Min buoyancy jump at an entrainment
  ! interface (m/s2), (~ jump in deg K /30)
  real(r8), parameter :: evhcmax = 10.   ! Max entrainment efficiency

  ! parameters specific to TKE-based entrainment closure

  real(r8), parameter :: a1l = 0.10   ! Dry entr. efficiency for TKE closure
  ! a1l = 0.2*tunl*erat^-1.5, where 
  ! erat = <e>/wstar^2 for dry CBL =  0.3.

  ! parameters specific to wstar-based entrainment closure

  real(r8), parameter :: a1i = 0.2            ! Dry entr. efficiency for wstar closure
  real(r8), parameter :: ccrit = 0.5          ! Min allowable sqrt(tke)/wstar
  real(r8), parameter :: wstar3factcrit = 0.5 ! 1/wstar3factcrit is the 
  ! max allowed enhancement of wstar3
  ! due to contribution from entrainment
  ! interfaces

  !  Parameters governing perturbation amplitudes

  real(r8), parameter :: wpertmin =1.e-6 ! min PBL eddy vert vel perturbation
  real(r8), parameter :: wfac = 0.707    ! ratio of wpert to sqrt(tke)
  real(r8), parameter :: tfac = 1.       ! ratio of tpert to (w't')/wpert
  ! (same ratio also used for q)

  !  Parameters limiting TKE

  real(r8), parameter :: rcapmin = 0.1   ! Min allowable e/<e> in a CL
  real(r8), parameter :: rcapmax = 2.0   ! Max allowable e/<e> in a CL
  real(r8), parameter :: tkemax =20.     ! tke capped at tkemax (m2/s2)
  type(vdiff_selector) :: fieldlist      ! Logical switches for moist vs dry mixing ratio diffusion

  !  Parameter controlling iteration in the subroutine, 'compute_eddy_diff'

  ! changed by cw 10/21/05 for consistency with latest UW codes 
   real(r8), parameter :: lambda = 0.7 ! Under-relaxation factor (0<lambda=<1)
!!$   real(r8), parameter :: lambda = 1.0 ! Under-relaxation factor (0<lambda=<1)
   logical :: initialized = .false.
   logical :: kvinit = .false.            ! whether to initialize kv

   real(r8) ::                                                                &
        an,                &   ! microphysics constant
        bn,                &   ! microphysics constant
        ap,                &   ! microphysics constant
        bp,                &   ! microphysics constant
        fac1,              &   ! microphysics factor
        fac2,              &   ! microphysics factor
        avg_fact               ! 1/(uwpbl_ndiv*uwpbl_ndiv)

  !cb end of caleddy params ---------

CONTAINS
  !
  !===============================================================================
  !
  subroutine init_eddy_diff(kind, pver, ntop_eddy, nbot_eddy )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose:  
    ! Initialize time independent variables of turbulence/pbl package.
    ! 
    !-----------------------------------------------------------------------
    use diffusion_solver, only: init_vdiff, vdiff_select
    implicit none
    !------------------------------Arguments--------------------------------
    integer, intent(in)  :: kind      ! kind of reals being passed in
    integer, intent(in)  :: pver     
    integer, intent(in)  :: ntop_eddy ! Top level to which eddy vert diff is applied.
    integer, intent(in)  :: nbot_eddy ! Bottom level to which eddy vert diff is applied.
    character(128) :: errstring                    ! error status for init_vdiff
    !---------------------------Local workspace-----------------------------
    integer :: k                     ! vertical loop index
    !-----------------------------------------------------------------------

    if ( kind .ne. r8 ) then
       write(6,*) 'wrong KIND of reals passed to init_diffusvity -- exiting.'
       call task_abort()
    endif

    if(mod(nx, uwpbl_ndiv).ne.0.or.(RUN3D.and.mod(ny, uwpbl_ndiv).ne.0)) then
       if(masterproc) print*,'nx or ny is not divisible by uwpbl_ndiv'
       if(masterproc) print*,'set in uwpbl.f90'
       if(masterproc) print*,'Stop.'
       call task_abort()
    end if

!    nx_p = nx / uwpbl_ndiv
    nx_p = 1

    allocate(cldn_flip   (nx_p, nzm) )
    allocate(qrl_flip    (nx_p, nzm) )
    allocate(s_flip      (nx_p, nzm) )
    allocate(rpdel_flip  (nx_p, nzm) )
    allocate(z_flip      (nx_p, nzm) )
    allocate(qv_sh_flip  (nx_p, nzm) )
    allocate(qc_sh_flip  (nx_p, nzm) )
    allocate(qi_sh_flip  (nx_p, nzm) )
    allocate(qp_sh_flip  (nx_p, nzm) )
    allocate(q_sh_flip   (nx_p, nzm) )
    allocate(tabs_flip   (nx_p, nzm) )
    allocate(u_flip      (nx_p, nzm) )
    allocate(v_flip      (nx_p, nzm) )
    allocate(w_flip      (nx_p, nz) )
    allocate(p_flip      (nx_p, nzm) )
    allocate(dtk_flip    (nx_p, nzm) )

    allocate(p8w_flip    (nx_p, nz) )
    allocate(z8w_flip    (nx_p, nz) )
    allocate(bprod_flip  (nx_p, nz) )
    allocate(sprod_flip  (nx_p, nz) )
    allocate(sfi_flip    (nx_p, nz) )
    allocate(tke_flip    (nx_p, nz) )
    allocate(cgh_flip    (nx_p, nz) )
    allocate(cgs_flip    (nx_p, nz) )
    allocate(kvm_in_flip (nx_p, nz) )
    allocate(kvh_in_flip (nx_p, nz) )
    allocate(kvm_out_flip(nx_p, nz) )
    allocate(kvh_out_flip(nx_p, nz) )
    allocate(kvq_out_flip(nx_p, nz) )

    allocate(taux_local  (nx_p) )
    allocate(tauy_local  (nx_p) )
    allocate(ksrftms     (nx_p) )
    allocate(tautmsx     (nx_p) )
    allocate(tautmsy     (nx_p) )
    allocate(topflx      (nx_p) )
    allocate(tpert       (nx_p) )
    allocate(qpert       (nx_p) )
    allocate(ustar_local (nx_p) )
    allocate(sfx_local   (nx_p) )
    allocate(pblh_local  (nx_p) )

    allocate(q_local     (nx_p, nzm, ncnst) )

    allocate(qfx_local   (nx_p, ncnst) )

    !
    ! Basic constants
    !
    zvir  = rv/rgas - 1.0
    ntop_turb = ntop_eddy
    nbot_turb = nbot_eddy
    b123 = b1**(2./3.) ! b1**(2/3)
    !
    ! Set the square of the mixing lengths.
    !
    allocate(ml2(pver+1))
    ml2(1:ntop_turb) = 0.
    do k = ntop_turb+1, nbot_turb
       ml2(k) = 30.0**2
    end do
    ml2(nbot_turb+1:pver+1) = 0.
    

    ! Initialize diffusion solver module
    call init_vdiff(r8, 1, dble(rv), dble(ggr), fieldlist, errstring)
!    print *, 'init_vdiff errstring is ',errstring
    ! Select the fields which will be diffused 
    if(vdiff_select(fieldlist,'s').ne.'') print *, 'error: ', vdiff_select(fieldlist,'s')
    if(vdiff_select(fieldlist,'q',1).ne.'') print *, 'error: ', vdiff_select(fieldlist,'q',1)
    if(vdiff_select(fieldlist,'u').ne.'') print *, 'error: ', vdiff_select(fieldlist,'u')
    if(vdiff_select(fieldlist,'v').ne.'') print *, 'error: ', vdiff_select(fieldlist,'v')

    ! added by cw 9/8/05
    call gestbl(            173.16d0,                        375.16d0,        &
                              20.0d0,                          .true.,        &
                             0.622d0,                     dble(lcond),        &
                          dble(lfus),                        dble(rv),        &
                            dble(cp),                      dble(tice) )
  
    kvinit = .false.
    if (nstep .eq. 1) then
       kvinit = .true.
       kvh_out_flip_3d = 0.0
       kvm_out_flip_3d = 0.0
    endif

    initialized = .true.

    ! microphysics initialization
    an = 1./dble(tbgmax-tbgmin)	
    bn = dble(tbgmin) * an
    ap = 1./dble(tprmax-tprmin)	
    bp = dble(tprmin) * ap
    fac1 = dble(fac_cond)+(1.0d0+bp)*dble(fac_fus)
    fac2 = dble(fac_fus)*ap

    avg_fact = 1.0 / dble(uwpbl_ndiv * uwpbl_ndiv)

    ! end cw addition


    return
  end subroutine init_eddy_diff
  !
  !===============================================================================
  !
  ! added by cw 7/26/05


    subroutine uwpbldrv

! ----------------------------------------------------------------------------

    use mse, only  : columnint, doMSE, navgMSE,                               &
                     huwpbl_mse, suwpbl_mse, u_uwpbl, v_uwpbl   ! peters
    implicit none

! ------------------------------ local variables ------------------------------

    integer ::                                                                &
         i,                &   ! x-direction loop index
         j,                &   ! y-direction loop index
         ip,               &   ! longitude loop counter within patch
         jp,               &   ! latitude loop counter within patch
         ig,               &   ! global longitude loop counter 
         jg,               &   ! global latitude loop counter 
         ic,               &
         jc,               &
         k,                &   ! z-direction loop index
         k_flip,           &   ! z-direction index for flipping upside down
         ncol,             &   ! number of columns passed to compute_eddy_diff
         nlev,             &   ! number of full (WRF parlance) layers
         nturb,            &   ! number of predictor-corrector timesteps
         ndiff,            &   ! number of diffusion timesteps
         idiff                 ! diffusion timestep

    logical ::                                                                &
         wstarent              ! whether to use wstar entrainment

    real ::                                                                   &
         omn                   ! partitioning function for non-precip. water

    real(r8) ::                                                               &
         ratio,            &   !
         diff_cfl              ! CFL from diffusion

    real(r8), dimension(nz) ::                                                &
         p8w_flip_temp

    real(r8), dimension(1,nzm) ::                                             &
         new_q,            &
         new_qp,           &
         new_qv,           &
         new_qi,           &
         new_qc,           &
         s_diff,           &   ! diffused-state static energy
         tabs_diff,        &   ! diffused-state temperature energy
         q_diff,           &   ! diffused-state total humidity
         qv_diff,          &   ! diffused-state specific humidity
         qc_diff,          &   ! diffused-state specific humidity
         qi_diff,          &   ! diffused-state specific humidity
         qp_diff,          &   ! diffused-state precip. water
         u_diff,           &   ! diffused-state zonal wind
         v_diff                ! diffused-state meridional wind

    real(r8), dimension(nzm) ::                                               &
         uwcu_du,          &   ! zonal wind tendency from UWCU
         uwcu_dv,          &   ! meridional wind tendency from UWCU
         uwcu_ds,          &   ! s tendency from UWCU
         uwcu_dq,          &   ! q tendency from UWCU
         uwcu_dqv,         &   ! qv tendency from UWCU
         uwcu_dqc,         &   ! qc tendency from UWCU
         uwcu_dqi              ! qi tendency from UWCU
    real(r8) wstar
!kzm stuff for perturb the PBL
    real, dimension(nzm) ::                                               &
         tv_local_mean, &
         tv_local_min, &
         tv_local_max, &
         tv_local_range
    real ::  wrk
    real  ranf_

    real(r8), dimension(1,nz) ::                                               &
         w_diff                ! diffused-state w

    real(r8) :: dtn_frac, ddtn

    real :: dtn_save

    real(r8), dimension(nx, ny, nz) ::  tke_local_flip

!for bilinear interpolation output    
    real wrktend_uwpbl(nx,ny,nzm)

    logical :: printout

    real :: lmaxtkz   ! maximum vertical diffusivity for this patch
    real :: lmaxndiff ! maximum number of PBL timsteps for this patch

    real(r8) :: cush_local! convective scale height for the patch

    real :: tmpcol(nzm), tmp, coefp  ! peters

! ------------------------------ executable code ------------------------------

    dtn_save=dtn
    ddtn=dble(dtn*sqrt(betafactor))
    dtn=dtn*sqrt(betafactor)

    if(.not.initialized) call init_eddy_diff(kind(qmincg), nzm, 1, nzm )

    nturb = 2 ! original value is 2
    wstarent = .true.

    ! no minimum water vapor mixing ratio
    qmincg = 0.d0

    ! cloud fraction (not used)
    cldn_flip(1,:) = 0.0d0

    ! set mountain drag to zero
    ksrftms = 0.d0

    ! domain variables
    ncol = 1
    nlev = nzm

    tk_z_uw = 0.0
    tkh_z_uw = 0.0
    pblh = 0.0
    pblh_local = 0.0
    tke_local_flip = 0.0

    lmaxndiff = 0.0 ! maximum ndiff for this patch

    lmaxtkz = 0.0 ! maximum tk_z for this patch

    do j = 1, ny/uwpbl_ndiv

       do i = 1, nx/uwpbl_ndiv
       
          p8w_flip = 0.0
          z8w_flip = 0.0
          z_flip = 0.0
          p_flip = 0.0
          qv_sh_flip = 0.0
          qc_sh_flip = 0.0
          qi_sh_flip = 0.0
          qp_sh_flip = 0.0
          q_sh_flip = 0.0
          p_flip = 0.0
          rpdel_flip = 0.0
          qrl_flip = 0.0
          z8w_flip = 0.0
          z_flip = 0.0
          tabs_flip = 0.0
          u_flip = 0.0
          v_flip = 0.0
          w_flip = 0.0
          s_flip = 0.0
          qfx_local = 0.0
          taux_local = 0.0
          tauy_local = 0.0
          sfx_local = 0.0
          cush_local = 0.0

          do jp = 1, uwpbl_ndiv
             jg = (j - 1) * uwpbl_ndiv + jp
             jc = jg + YES3D
             do ip = 1, uwpbl_ndiv
                ig = (i - 1) * uwpbl_ndiv + ip
                ic = ig + 1
                
                p8w_flip_temp(:) = 0.0
                
                p8w_flip_temp(nz) = dble(pres(1) +                            &
                     0.01 * (p(ig, jg, 1) + rho(1) * ggr * z(1)))

                p8w_flip(1, nz) = p8w_flip(1, nz) + p8w_flip_temp(nz)
                     
                
                do k = 1, nz-1
                   k_flip = nz - k + 1 
                   p8w_flip_temp(k_flip - 1) =                                &
                        exp( log( dble(pres(k) + 0.01 * p(ig, jg, k)) /       &
                        p8w_flip_temp(k_flip)) * dble(zi(k+1) - zi(k)) /      &
                        dble(z(k) - zi(k)) + log(p8w_flip_temp(k_flip)))
                   
                   p8w_flip(1, k_flip - 1) = p8w_flip(1, k_flip - 1) +        &
                        p8w_flip_temp(k_flip - 1)
                enddo
                
                z8w_flip(1, nz) = 0.0 ! full-level heights (m)
                
                do k = 1, nz
                   k_flip = nz - k + 1
                   kvm_in_flip(1, k) = dble(kvm_out_flip_3d(ig, jg, k))
                   kvh_in_flip(1, k) = dble(kvh_out_flip_3d(ig, jg, k))
                   w_flip(1, k_flip) = w_flip(1, k_flip) +                    &
                        dble(w(ig, jg, k))
                enddo
                
                do k = 1, nzm
                   
                   k_flip = nzm - k + 1 
                   ! k_flip is used to change z-orientation from WRF/SAM (index 1 
                   ! is near the surface) to CAM (index 1 is near the TOA)
                   
                   omn = omegan(tabs(ig, jg, k))
                   
                   qv_sh_flip(1, k_flip) = qv_sh_flip(1, k_flip) +            &
                        dble(q(ig, jg, k) - qn(ig, jg, k))
                   qc_sh_flip(1, k_flip) = qc_sh_flip(1, k_flip) +            &
                        dble(qn(ig, jg, k) * omn)
                   qi_sh_flip(1, k_flip) = qi_sh_flip(1, k_flip) +            &
                        dble(qn(ig, jg, k) * (1.-omn))
                   qp_sh_flip(1, k_flip) = qp_sh_flip(1, k_flip) +            &
                        dble(qp(ig, jg, k))
                   
                   p_flip(1, k_flip) = p_flip(1, k_flip) +                    &
                        100.0d0 * dble(pres(k)) + dble(p(ig, jg, k))
                   
                   ! 1/pdel
                   rpdel_flip(1, k_flip) = rpdel_flip(1, k_flip) +            &
                        1. / ( 100.0d0 * (p8w_flip_temp(k_flip+1) -           &
                        p8w_flip_temp(k_flip) ))
                   
                   ! LW heating rate (K/s) * cp * dp (W/kg*Pa)
                   qrl_flip(1, k_flip) = qrl_flip(1, k_flip) +                &
                        dble(radqrclw(k)) * dble(cp) / (100.0d0 *             &
                        (p8w_flip_temp(k_flip+1) - p8w_flip_temp(k_flip) ))

                   ! half-level heights
                   z8w_flip(1, k_flip) = z8w_flip(1, k_flip) + dble(zi(k+1))
                   
                   ! full-level heights
                   z_flip(1, k_flip) = z_flip(1, k_flip) + dble(z(k))
                   
                   tabs_flip(1, k_flip) = tabs_flip(1, k_flip) +              &
                        dble(tabs(ig, jg, k))
                   u_flip(1, k_flip) = u_flip(1, k_flip) +                    &
                        0.5d0*dble(u(ig, jg, k) + u(ic, jg, k))
                   v_flip(1, k_flip) = v_flip(1, k_flip) +                    &
                        0.5d0*dble(v(ig, jg, k) + v(ig, jc, k))
                  
                enddo
                
                ! units are corrected in the call to compute_eddy_diff
                qfx_local(1, 1) = qfx_local(1, 1) + dble(fluxbq(ig, jg))
                taux_local(1) = taux_local(1) + dble(fluxbu(ig, jg))
                tauy_local(1) = tauy_local(1) + dble(fluxbv(ig, jg))
                sfx_local(1) = sfx_local(1) + dble(fluxbt(ig, jg))

                cush_local = cush_local + dble(cush(ig, jg))
             enddo
          enddo
          
          p8w_flip = p8w_flip * 100.0d0 * avg_fact! from mb to Pa
          z8w_flip = z8w_flip * avg_fact
          qv_sh_flip = qv_sh_flip * avg_fact
          qc_sh_flip = qc_sh_flip * avg_fact
          qi_sh_flip = qi_sh_flip * avg_fact
          qp_sh_flip = qp_sh_flip * avg_fact
          p_flip = p_flip * avg_fact
          rpdel_flip = rpdel_flip * avg_fact
          qrl_flip = qrl_flip * avg_fact
          z_flip = z_flip * avg_fact
          
          tabs_flip = tabs_flip * avg_fact

          ! dry static energy (why isn't this vector?)
          do k = 1, nzm
             k_flip = nzm - k + 1
             s_flip(1, k_flip) = dble(cp)*tabs_flip(1, k_flip) + dble(ggr*z(k))
          enddo

          u_flip = u_flip * avg_fact
          v_flip = v_flip * avg_fact
          w_flip = w_flip * avg_fact
          cush_local = cush_local * avg_fact
          q_sh_flip(1, :) = qv_sh_flip(1, :) + qc_sh_flip(1, :) + qi_sh_flip(1, :)

          qfx_local = qfx_local * avg_fact / sqrt(betafactor)
          taux_local = taux_local * avg_fact / sqrt(betafactor)
          tauy_local = tauy_local * avg_fact / sqrt(betafactor)
          sfx_local = sfx_local * avg_fact / sqrt(betafactor)

          s_diff(1, :) = s_flip(1, :)
          qv_diff(1, :) = qv_sh_flip(1, :)
          qc_diff(1, :) = qc_sh_flip(1, :)
          qi_diff(1, :) = qi_sh_flip(1, :)
          qp_diff(1, :) = qp_sh_flip(1, :)
          u_diff(1, :) = u_flip(1, :)
          v_diff(1, :) = v_flip(1, :)
          w_diff(1, :) = w_flip(1, :)
          q_diff(1, :) = q_sh_flip(1, :)
          tabs_diff(1, :) = tabs_flip(1, :)

          call compute_eddy_diff(                                   1,        &
                                nlev,                               1,        &
                     tabs_diff(:, :),                   qv_diff(:, :),        &
                            2.0*ddtn,                   qc_diff(:, :),        &
                       qi_diff(:, :),                    s_diff(:, :),        &
                    rpdel_flip(:, :),                 cldn_flip(:, :),        &
                      qrl_flip(:, :),                    z_flip(:, :),        &
                      z8w_flip(:, :),                    p_flip(:, :),        &
                      p8w_flip(:, :),                    u_diff(:, :),        &
                        v_diff(:, :),                   taux_local(:),        &
                       tauy_local(:),   sfx_local(:)*dble(cp*rhow(1)),        &
       qfx_local(:, :)*dble(rhow(1)),                        wstarent,        &
                               nturb,                     ustar_local,        &
                          pblh_local,               kvm_in_flip(:, :),        &
                   kvh_in_flip(:, :),              kvm_out_flip(:, :),        &
                  kvh_out_flip(:, :),              kvq_out_flip(:, :),        &
                      cgh_flip(:, :),                  cgs_flip(:, :),        &
                            tpert(:),                        qpert(:),        &
                      tke_flip(:, :),                bprod_flip(:, :),        &
                    sprod_flip(:, :),                  sfi_flip(:, :),        &
                              fqsatd,                          kvinit,wstar )

          do k = 1, nz
             ratio = min(1.0d0, 500.0d0 / max(kvm_out_flip(1, k),             &
                  kvh_out_flip(1, k)))
             kvm_out_flip(1, k) = kvm_out_flip(1, k) * ratio
             kvh_out_flip(1, k) = kvh_out_flip(1, k) * ratio
             
             lmaxtkz = max( lmaxtkz, kvm_out_flip(1, k) )
             
          enddo
          
          diff_cfl = 0.0d0
          do k = 1, nz
             k_flip = nz - k + 1
             diff_cfl = max(diff_cfl,                                         &
                  ddtn / (dble(dz * adzw(k))**2) * kvm_out_flip(1, k_flip))
          end do
          
          ndiff = ceiling(diff_cfl / 0.25d0)
          if(fixedndiff.gt.0) ndiff = fixedndiff
          dtn_frac = ddtn / dble(max(1, ndiff))
          lmaxndiff = max(lmaxndiff, real(ndiff))
          
          do idiff = 1, ndiff
             if(ndiff .gt. 1) then
                call compute_eddy_diff(                             1,        &
                                nlev,                               1,        &
                     tabs_diff(:, :),                   qv_diff(:, :),        &
                        2.0*dtn_frac,                   qc_diff(:, :),        &
                       qi_diff(:, :),                    s_diff(:, :),        &
                    rpdel_flip(:, :),                 cldn_flip(:, :),        &
                      qrl_flip(:, :),                    z_flip(:, :),        &
                      z8w_flip(:, :),                    p_flip(:, :),        &
                      p8w_flip(:, :),                    u_diff(:, :),        &
                        v_diff(:, :),                   taux_local(:),        &
                       tauy_local(:),   sfx_local(:)*dble(cp*rhow(1)),        &
       qfx_local(:, :)*dble(rhow(1)),                        wstarent,        &
                               nturb,                     ustar_local,        &
                          pblh_local,               kvm_in_flip(:, :),        &
                   kvh_in_flip(:, :),              kvm_out_flip(:, :),        &
                  kvh_out_flip(:, :),              kvq_out_flip(:, :),        &
                      cgh_flip(:, :),                  cgs_flip(:, :),        &
                            tpert(:),                        qpert(:),        &
                      tke_flip(:, :),                bprod_flip(:, :),        &
                     sprod_flip(:, :),                  sfi_flip(:, :),       &
                              fqsatd,                          kvinit,wstar )

                do k = 1, nz
                   ratio = min(1.0d0, 500.0 / max(kvm_out_flip(1, k),         &
                        kvh_out_flip(1, k)))
                   kvm_out_flip(1, k) = kvm_out_flip(1, k) * ratio
                   kvh_out_flip(1, k) = kvh_out_flip(1, k) * ratio

                   lmaxtkz = max( lmaxtkz, kvm_out_flip(1, k) )

                enddo

             endif

             ! call uwcu with pre-diffused state
             if(douwcu) then
                call mcshallow(                                               &
                            dtn_frac,                    p_flip(1, :),        &
                      p8w_flip(1, :),                    u_diff(1, :),        &
                        v_diff(1, :),                 tabs_diff(1, :),        &
                        s_diff(1, :),                   qv_diff(1, :),        &
                       qc_diff(1, :),                   qi_diff(1, :),        &
                      tke_flip(1, :),                   pblh_local(1),        &
                          cush_local,                         uwcu_du,        &
                             uwcu_dv,                         uwcu_ds,        &
                             uwcu_dq,                        uwcu_dqv,        &
                            uwcu_dqc,                        uwcu_dqi )
 
                 do jp = 1, uwpbl_ndiv
                    jg = (j - 1) * uwpbl_ndiv + jp
                    do ip = 1, uwpbl_ndiv
                       ig = (i - 1) * uwpbl_ndiv + ip                        
                       do k = 1, nzm
                          k_flip = nz - k + 1
                          qvtend_uwcu(ig, jg, k) = qvtend_uwcu(ig, jg, k) +    &
                               uwcu_dqv(k_flip) / dtn
                          qltend_uwcu(ig, jg, k) = qltend_uwcu(ig, jg, k) +    &
                               uwcu_dqc(k_flip) / dtn
                          qitend_uwcu(ig, jg, k) = qitend_uwcu(ig, jg, k) +    &
                               uwcu_dqi(k_flip) / dtn
                          ttend_uwcu(ig, jg, k) = ttend_uwcu(ig, jg, k) +      &
                               uwcu_ds(k_flip) / cp / dtn
                       end do
                    enddo
                 end do
             endif

             ! mult sfx by cp -- this is diffusing DSE and sfx has units of t
             call diffuse_all(                                                &
                        s_diff(1, :),                    q_diff(1, :),        &
                       qp_diff(1, :),                    u_diff(1, :),        &
                        v_diff(1, :),                    w_diff(1, :),        &
             sfx_local(1) * dble(cp),                 qfx_local(1, 1),        &
                       taux_local(1),                   tauy_local(1),        &
                  kvm_out_flip(1, :),              kvh_out_flip(1, :),        &
                            dtn_frac )
             
             if(ndiff .gt. 1) then

             do k = 1, nzm
                k_flip = nzm - k + 1
                tabs_diff(1, k_flip) = (s_diff(1, k_flip) - dble(ggr*z(k))) / dble(cp)
             enddo

             ! use diffusivities from this iteration for next iteration
             kvm_in_flip = kvm_out_flip
             kvh_in_flip = kvh_out_flip

             ! update _diff values of q, qp, qv, qc, qi, tabs, and s 
             ! to take into account changes in q and qp from diffusion
             call microphysics(    i,                               j,        &
                        s_diff(1, :),                 tabs_diff(1, :),        &
                     tabs_flip(1, :),                    q_diff(1, :),        &
                     q_sh_flip(1, :),                   qp_diff(1, :),        &
                    qp_sh_flip(1, :),                   qv_diff(1, :),        &
                       qc_diff(1, :),                   qi_diff(1, :) )
             endif

             if(douwcu) then
                u_diff(1, :) = u_diff(1, :) + uwcu_du
                v_diff(1, :) = v_diff(1, :) + uwcu_dv
                s_diff(1, :) = s_diff(1, :) + uwcu_ds
                q_diff(1, :) = q_diff(1, :) + uwcu_dq
                qv_diff(1, :) = qv_diff(1, :) + uwcu_dqv
                qc_diff(1, :) = qc_diff(1, :) + uwcu_dqc
                qi_diff(1, :) = qi_diff(1, :) + uwcu_dqi
             endif

             do jp = 1, uwpbl_ndiv
                jg = (j - 1) * uwpbl_ndiv + jp
                do ip = 1, uwpbl_ndiv
                   ig = (i - 1) * uwpbl_ndiv + ip 
                   
                   do k = 1, nzm
                      k_flip = nz - k + 1
                      
                      tk_z_uw(ig, jg, k) = tk_z_uw(ig, jg, k) +               &
                           (real(kvm_out_flip(1, k_flip)) * (zi(k+1) - z(k)) +&
                           real(kvm_out_flip(1, k_flip-1)) * (z(k) - zi(k))) /&
                           ( zi(k+1) - zi(k) )
                      tkh_z_uw(ig, jg, k) = tkh_z_uw(ig, jg, k) +             &
                           (real(kvh_out_flip(1, k_flip)) * (zi(k+1) - z(k)) +&
                           real(kvh_out_flip(1, k_flip-1)) * (z(k) - zi(k))) /&
                           ( zi(k+1) - zi(k) )
                   enddo
                   
                   cush(ig, jg) = real(cush_local)
                   pblh(ig, jg) = pblh(ig, jg) + real(pblh_local(1))
                   tke_local_flip(ig, jg, :) = tke_local_flip(ig, jg, :) + tke_flip(1, :)

                enddo
             enddo
          enddo

!kzm Oct. 13, 06, perturb the PBL tendencies
          tv_local_mean=0.
          tv_local_min=10000.
          tv_local_max=-10000.
          do k = 1, nzm
             do jp = 1, uwpbl_ndiv
                jg = (j - 1) * uwpbl_ndiv + jp
                do ip = 1, uwpbl_ndiv
                   ig = (i - 1) * uwpbl_ndiv + ip 
                   wrk=tabs(ig,jg,k)*(1+0.61*q(ig,jg,k)-1.61*qn(ig,jg,k)-qp(ig,jg,k))
                   tv_local_mean(k)=tv_local_mean(k)+wrk
                   tv_local_min(k)=min(tv_local_min(k),wrk)
                   tv_local_max(k)=max(tv_local_max(k),wrk)                   
                enddo
             enddo
             tv_local_mean(k)=tv_local_mean(k)/float(uwpbl_ndiv)/float(uwpbl_ndiv)
             wrk=max(tv_local_mean(k)-tv_local_min(k),tv_local_max(k)-tv_local_mean(k))/3
             tv_local_range(k)=max(wrk,0.001)
          enddo
             k_flip = nzm - k + 1
!             wrk=((30.*pblh_local(1))**0.333*max(real((s_diff(1, k_flip) - s_flip(1, k_flip)) / dtn/cp*sqrt(betafactor)),0.001)**0.667)
             wrk=min(0.1,wstar*wstar*30/pblh_local(1)*betafactor**0.333)!The range is taken to be 4 scale
          do k=1,nzm
             !note here uses dtn not ddtn
             if(z(k).gt.pblh_local(1)) then
                tv_local_range(k)=10000.
             else
                if (tv_local_range(k).gt.wrk) then
                   tv_local_range(k)=10000.
                else
!                   write(*,*) 'kzm',tv_local_range(k),wrk,k,pblh_local(1)
                   tv_local_range(k)=tv_local_range(k)
                endif
             endif
!             tv_local_range(k)=max(5*wrk*max(exp(wrk-(30./pblh_local(1))**0.333*max(ravefactor*pblh_local(1)*real((s_diff(1, k_flip) - s_flip(1, k_flip)) / ddtn/cp),0.001)**0.667)/0.01),0.05),0.001)
!factor 5 gives 100mm/day, 20 gives 1mm/day
!             if(z(k).gt.pblh_local(1)) tv_local_range(k)=10000.
             tv_local_range(k)=10000.
!             if(k.gt.12) tv_local_range(k)=10000.
          enddo
          
!          write(*,*) 'kzm', tv_local_range(1:20)
          do jp = 1, uwpbl_ndiv
             jg = (j - 1) * uwpbl_ndiv + jp
             do ip = 1, uwpbl_ndiv
                ig = (i - 1) * uwpbl_ndiv + ip 
               
                tk_z_uw(ig, jg, :) = tk_z_uw(ig, jg, :) / float(max(1, ndiff))
                tkh_z_uw(ig, jg, :) = tkh_z_uw(ig, jg, :) / float(max(1, ndiff))
                pblh(ig, jg) = pblh(ig, jg) / float(max(1, ndiff))

                uwpbl_ndiff(ig, jg) = ndiff

                
                ! save surface fluxes so they can be subtracted from SGS
                ! diffusion
               uwpbl_fluxbu(ig, jg) = real(taux_local(1)) * sqrt(betafactor)
               uwpbl_fluxbv(ig, jg) = real(tauy_local(1)) * sqrt(betafactor)
               uwpbl_fluxbt(ig, jg) = real(sfx_local(1)) * sqrt(betafactor)
               uwpbl_fluxbq(ig, jg) = real(qfx_local(1, 1)) * sqrt(betafactor)

                do k = 1, nzm
                   k_flip = nzm - k + 1
                   wrk=tabs(ig,jg,k)*(1+0.61*q(ig,jg,k)-1.61*qn(ig,jg,k)-qp(ig,jg,k))
                   
                   qtend_uwpbl(ig, jg, k) =                                   &
                        real((q_diff(1, k_flip) - q_sh_flip(1, k_flip)) / ddtn)*(1+(wrk-tv_local_mean(k))/tv_local_range(k))
!!$                   qtend_uwpbl(ig, jg, k) =                                   &
!!$                        real((q_diff(1, k_flip) - q_sh_flip(1, k_flip)) / ddtn)*(1+uwcu_pert_cldsrc(ig,jg))

                   qptend_uwpbl(ig, jg, k) =                                  &
                        real((qp_diff(1, k_flip) - qp_sh_flip(1, k_flip)) / ddtn)

                   qvtend_uwpbl(ig, jg, k) =                                  &
                        real((qv_diff(1, k_flip) - qv_sh_flip(1, k_flip)) / ddtn)

                   qctend_uwpbl(ig, jg, k) =                                  &
                        real((qc_diff(1, k_flip) - qc_sh_flip(1, k_flip)) / ddtn)

                   qitend_uwpbl(ig, jg, k) =                                  &
                        real((qi_diff(1, k_flip) - qi_sh_flip(1, k_flip)) / ddtn)
                   
                   utend_uwpbl(ig, jg, k) =                                   &
                        real((u_diff(1, k_flip) - u_flip(1, k_flip)) / ddtn)

                   vtend_uwpbl(ig, jg, k) =                                   &
                        real((v_diff(1, k_flip) - v_flip(1, k_flip)) / ddtn)

                   ttend_uwpbl(ig, jg, k) =                                   &
                        real( (s_diff(1, k_flip) - s_flip(1, k_flip)) / ddtn ) / cp*(1+(wrk-tv_local_mean(k))/tv_local_range(k)) &
                         - fac_cond * qctend_uwpbl(ig, jg, k) -     &
                        fac_fus * qitend_uwpbl(ig, jg, k)
!!$                   ttend_uwpbl(ig, jg, k) =                                   &
!!$                        real( (s_diff(1, k_flip) - s_flip(1, k_flip)) / ddtn ) / cp &
!!$                         - fac_cond * qctend_uwpbl(ig, jg, k) -     &
!!$                        fac_fus * qitend_uwpbl(ig, jg, k)
!!$                   ttend_uwpbl(ig, jg, k) =                                   &
!!$                        real( (s_diff(1, k_flip) - s_flip(1, k_flip)) / ddtn ) / cp*(1+uwcu_pert_cldsrc(ig,jg)) &
!!$                         - fac_cond * qctend_uwpbl(ig, jg, k) -     &
!!$                        fac_fus * qitend_uwpbl(ig, jg, k)
!!$               
                   kvm_out_flip_3d(ig, jg, k) = real(kvm_out_flip(1, k))
                   kvh_out_flip_3d(ig, jg, k) = real(kvh_out_flip(1, k))
                enddo

                do k = 1, nz
                   k_flip = nz - k + 1
                   tke_uw(ig, jg, k) = real(tke_local_flip(ig, jg, k_flip)) / float(ndiff)
                   wtend_uwpbl(ig, jg, k) =                                   &
                     real((w_diff(1, k_flip) - w_flip(1, k_flip)) / ddtn)
                enddo

             enddo
          enddo
       enddo
    enddo

    if(dowally.and.RUN3D.and.rank.lt.nsubdomains_x) then
       do k=1,nzm
          do i=1,nx
             vtend_uwpbl(i,1,k) = 0.
          end do
       end do
    end if

    if(dowallx.and.mod(rank,nsubdomains_x).eq.0) then
       do k=1,nzm
          do j=1,ny
             utend_uwpbl(1,j,k) = 0.
          end do
       end do
    end if

!kzm =========bilinear interpolation begin
 call task_exchange(ttend_uwpbl,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,47)
 call task_exchange(utend_uwpbl,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,48)
 call task_exchange(vtend_uwpbl,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,49)
 call task_exchange(qtend_uwpbl,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,50)
 call task_exchange(qptend_uwpbl,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,51)
 call task_exchange(wtend_uwpbl,0,nxp1,1-YES3D,nyp1,nz,1,1,1,1,52)

call bilinear_uwpbl(ttend_uwpbl,wrktend_uwpbl)
ttend_uwpbl(1:nx,1:ny,1:nzm)=wrktend_uwpbl
call bilinear_uwpbl(utend_uwpbl,wrktend_uwpbl)
utend_uwpbl(1:nx,1:ny,1:nzm)=wrktend_uwpbl
call bilinear_uwpbl(vtend_uwpbl,wrktend_uwpbl)
vtend_uwpbl(1:nx,1:ny,1:nzm)=wrktend_uwpbl
call bilinear_uwpbl(qtend_uwpbl,wrktend_uwpbl)
qtend_uwpbl(1:nx,1:ny,1:nzm)=wrktend_uwpbl
call bilinear_uwpbl(qptend_uwpbl,wrktend_uwpbl)
qptend_uwpbl(1:nx,1:ny,1:nzm)=wrktend_uwpbl
call bilinear_uwpbl(wtend_uwpbl(:,:,1:nzm),wrktend_uwpbl)
wtend_uwpbl(1:nx,1:ny,1:nzm)=wrktend_uwpbl
!kzm =============bilinear interpolation end

       t(1:nx, 1:ny, 1:nzm) = t(1:nx, 1:ny, 1:nzm) + ttend_uwpbl(1:nx,1:ny,1:nzm) * dtn
       
       ! update U and V here with forward Euler 
       u(1:nx, 1:ny, 1:nzm) = u(1:nx, 1:ny, 1:nzm) + utend_uwpbl(1:nx,1:ny,1:nzm) * dtn
       v(1:nx, 1:ny, 1:nzm) = v(1:nx, 1:ny, 1:nzm) + vtend_uwpbl(1:nx,1:ny,1:nzm) * dtn
       
       ! update U and V in pressure with Adams-Bashforth
!!$    dudt(1:nx, 1:ny, 1:nzm, na) = dudt(1:nx, 1:ny, 1:nzm, na) +               &
!!$         utend_uwpbl * sqrt(betafactor)
!!$    dvdt(1:nx, 1:ny, 1:nzm, na) = dvdt(1:nx, 1:ny, 1:nzm, na) +               &
!!$         vtend_uwpbl * sqrt(betafactor)
    
    q(1:nx, 1:ny, 1:nzm) = q(1:nx, 1:ny, 1:nzm) + qtend_uwpbl(1:nx,1:ny,1:nzm) * dtn
    qp(1:nx, 1:ny, 1:nzm) = qp(1:nx, 1:ny, 1:nzm) + qptend_uwpbl(1:nx,1:ny,1:nzm) * dtn
          
    ! update W here with forward Euler                   
    w(1:nx, 1:ny, 1:nz) = w(1:nx, 1:ny, 1:nz) + wtend_uwpbl(1:nx,1:ny,1:nz) * dtn

    ! update W in pressure with Adams-Bashforth
!!$    dwdt(1:nx, 1:ny, 1:nz, na) = dwdt(1:nx, 1:ny, 1:nz, na) +                 &
!!$         wtend_uwpbl * sqrt(betafactor)


    maxndiff = lmaxndiff
    maxtkz = lmaxtkz
    if(dompi) then
       call task_max_real(lmaxndiff, maxndiff, 1)
       call task_max_real(lmaxtkz, maxtkz, 1)
    end if

!!$    if(masterproc) then
!!$       write(*,999) nstep, maxndiff, maxtkz
!!$999    format('DIFFSUBNSTEP = ',i10,'    MAXNDIFF=',f7.0,'    MAXTKZ=',f9.2)
!!$    end if


    dtn=dtn_save

    ! peters accumulate budget stuff
    if(doMSE) then
      coefp = 1./dtfactor/float(navgMSE)
      u_uwpbl = u_uwpbl + utend_uwpbl(1:nx,1:ny,1:nzm)*coefp
      v_uwpbl = v_uwpbl + vtend_uwpbl(1:nx,1:ny,1:nzm)*coefp
      do i=1,nx
        do j=1,ny
         tmpcol = ttend_uwpbl(i,j,:)*cp*coefp
         call columnint(tmpcol,tmp)
         suwpbl_mse(i,j) = suwpbl_mse(i,j) + tmp                ! W/m^2

         tmpcol = lcond*coefp*(qtend_uwpbl(i,j,:) + qptend_uwpbl(i,j,:))
         call columnint(tmpcol,tmp)
         huwpbl_mse(i,j) = huwpbl_mse(i,j) + tmp                    ! W/m^2
        end do
      end do
    end if



  end subroutine uwpbldrv


    !----------------------------------------------------------------------- 

  subroutine microphysics(         i,                               j,        &
                              s_diff,                       tabs_diff,        &
                           tabs_flip,                          q_diff,        &
                              q_flip,                         qp_diff,        &
                             qp_flip,                         qv_diff,        &
                             qc_diff,                         qi_diff )

! calculates updated global values of q, qp, and tabs by adding 
! the block-averaged tendency (_diff - _flip) to each global point.  These
! are then used to calculate new global values of qv, qi, and qc, which are 
! block averaged and put into qv_diff, qi_diff, and qc_diff
    integer, intent(in) :: i, j

    real(r8), dimension(:), intent(in) ::                                     &
         tabs_flip             ! original block-averaged tabs

    real(r8), dimension(:), intent(in) ::                                     &
         qp_flip,          &   ! original block-averaged qp
         q_flip                ! original block-averaged qp


    real(r8), dimension(:), intent(inout) ::                                  &
         tabs_diff,        &   ! diffused block-averaged tabs
         s_diff

    real(r8), dimension(:), intent(inout) ::                                  &
         q_diff,           &   ! diffused block-averaged q
         qp_diff,          &   ! diffused block-averaged qp
         qv_diff,          &
         qc_diff,          &   !
         qi_diff


    integer :: ip, jp, ig, jg, k, niter, k_flip

    real(r8) ::                                                               &
         global_tabs,      &   ! updated tabs (PBL timestep)
         dtabs,            &
         tabs1                 ! temporary microphysics tabs

    real(r8) ::                                                               &
         coef_cl,          &
         fff,              &
         dfff,             &
         qsat,             &
         dqsat,            &
         dlstarp,          &
         dlstarn,          &
         lstarp,           &
         lstarn,           &
         om,               &
         omn

    real(r8) ::                                                               &
         global_q,         &   ! 
         global_qn,        &   !
         global_qp

    real(r8), dimension(size(q_diff)) ::                                      &
         dblock_tabs

    real(r8), dimension(size(q_diff)) ::                                      &
         block_q,          &   ! 
         block_qp,         &   ! 
         block_qv,         &   ! 
         block_qc,         &   ! 
         block_qi


    dblock_tabs = 0.0
    block_qc = 0.0
    block_qi = 0.0
    block_qv = 0.0
    block_q =  0.0
    block_qp = 0.0
    
    do jp = 1, uwpbl_ndiv
       jg = (j - 1) * uwpbl_ndiv + jp
       do ip = 1, uwpbl_ndiv
          ig = (i - 1) * uwpbl_ndiv + ip

          do k = 1, nzm
             
             k_flip = nzm - k + 1
             
             global_tabs = dble(tabs(ig, jg, k)) +                            &
                  (tabs_diff(k_flip) - tabs_flip(k_flip)) 
             global_qp = dble(qp(ig, jg, k)) +                                &
                  (qp_diff(k_flip) - qp_flip(k_flip))
             global_q = dble(q(ig, jg, k)) +                                  &
                  (q_diff(k_flip) - q_flip(k_flip))
             
             ! Initial guess for temperature assuming no cloud water/ice:
             tabs1 = ( global_tabs + fac1 * global_qp) /                      &
                  (1.d0 + fac2 * global_qp)
             
             ! Warm cloud:
             
             if(tabs1 .ge. dble(tbgmax)) then
                
                tabs1 = global_tabs + dble(fac_cond)*global_qp
                qsat = dble(qsatw(real(tabs1), pres(k)))
                
                ! Ice cloud:
                
             elseif(tabs1 .le. dble(tbgmin)) then
                
                tabs1 = global_tabs + dble(fac_sub)*global_qp
                qsat = dble(qsati(real(tabs1), pres(k)))
                
                ! Mixed-phase cloud:
                
             else
                
                om = an*tabs1 - bn
                qsat = om*dble(qsatw(real(tabs1), pres(k))) + (1. - om)*dble(qsati(real(tabs1), pres(k)))
                
             endif
             
             
             !  Test if condensation is possible:
             
             
             if(global_q .gt. qsat) then
                
                niter=0
                dtabs = 100.d0
                do while(abs(dtabs).gt.0.01.and.niter.lt.10)
                   if(tabs1.ge.dble(tbgmax)) then
                      om = 1.
                      lstarn = dble(fac_cond)
                      dlstarn = 0.d0
                      qsat = dble(qsatw(real(tabs1), pres(k)))
                      dqsat = dble(dtqsatw(real(tabs1), pres(k)))
                   else if(tabs1.le.dble(tbgmin)) then
                      om = 0.
                      lstarn = dble(fac_sub)
                      dlstarn = 0.
                      qsat = dble(qsati(real(tabs1), pres(k)))
                      dqsat = dble(dtqsati(real(tabs1), pres(k)))
                   else
                      om = an*tabs1 - bn
                      lstarn = dble(fac_cond) + (1. - om)*dble(fac_fus)
                      dlstarn = an
                      qsat = om*dble(qsatw(real(tabs1), pres(k))) +            &
                           (1. - om)*dble(qsati(real(tabs1), pres(k)))
                      dqsat = om*dble(dtqsatw(real(tabs1), pres(k))) +         &
                           (1. - om)*dble(dtqsati(real(tabs1), pres(k)))
                   endif
                   if(tabs1 .ge. dble(tprmax)) then
                      lstarp = dble(fac_cond)
                      dlstarp = 0.
                   else if(tabs1 .le. dble(tprmin)) then
                      lstarp = dble(fac_sub)
                      dlstarp = 0.
                   else
                      lstarp = dble(fac_cond)+(1.-om)*dble(fac_fus)
                      dlstarp = ap
                   endif
                   fff = global_tabs - tabs1 + lstarn * (global_q - qsat) +     &
                        lstarp * global_qp
                   
                   dfff = dlstarn * (global_q - qsat) +                   &
                        dlstarp*global_qp - lstarn*dqsat - 1.d0
                   dtabs = -fff/dfff
                   niter = niter + 1
                   tabs1 = tabs1 + dtabs
                end do
                
                qsat = qsat + dqsat * dtabs
                global_qn = max(0.d0, global_q - qsat)
                
             else
                
                global_qn = 0.
                
             endif

             global_tabs = tabs1
             dblock_tabs(k_flip) = dblock_tabs(k_flip) + global_tabs

             global_qp = max(0., global_qp) ! just in case

             omn = omegan(real(global_tabs))
             block_q(k_flip) = block_q(k_flip) + global_q
             block_qp(k_flip) = block_qp(k_flip) + global_qp
             block_qv(k_flip) = block_qv(k_flip) + global_q - global_qn
             block_qc(k_flip) = block_qc(k_flip) + global_qn * omn
             block_qi(k_flip) = block_qi(k_flip) + global_qn * (1. - omn)
             
          end do

       end do
    end do

    tabs_diff = dblock_tabs * avg_fact    
    q_diff  = block_q  * avg_fact    
    qp_diff = block_qp * avg_fact    
    qv_diff = block_qv * avg_fact    
    qc_diff = block_qc * avg_fact    
    qi_diff = block_qi * avg_fact    
    
    do k = 1, nzm
       k_flip = nzm - k + 1
       s_diff(k_flip) = tabs_diff(k_flip) * dble(cp) + dble(ggr * z(k))
    enddo

  end subroutine microphysics
  

! ----------------------------------------------------------------------------


  subroutine diffuse_all(     s_diff,                          q_diff,        &
                             qp_diff,                          u_diff,        &
                              v_diff,                          w_diff,        &
                           sfx_local,                       qfx_local,        &
                          taux_local,                      tauy_local,        &
                                 tkz,                            tkhz,        &
                           dtn_local )


    real(r8), dimension(:), intent(inout) ::                                  &
         tkhz,             &   !
         tkz                   !

    real(r8), intent(in) ::                                                   &
         qfx_local,        &   ! surface moisture flux
         sfx_local,        &   ! surface sensible heat flux
         taux_local,       &   ! 
         tauy_local,       &   ! 
         dtn_local


    real(r8), dimension(:), intent(inout) ::                                  &
         s_diff

    real(r8), dimension(:), intent(inout) ::                                  &
         q_diff,           &   ! 
         qp_diff,          &   ! 
         u_diff,           &   ! 
         v_diff,           &   ! 
         w_diff                ! 

    integer :: k, k_flip
    real(r8), dimension(nzm) :: tkz_full

    real, dimension(nzm) :: tkhz_full
    real, dimension(nz) :: tkz_half
    real, dimension(nz) :: tkhz_half
    real, dimension(nzm) :: flip

!!$
!!$    do k = 1, nzm
!!$       k_flip = nz - k + 1
!!$       tkz_full(k) = (tkz(k_flip) * (zi(k+1) - z(k)) + tkz(k_flip-1) * (z(k) - zi(k))) /&
!!$            ( zi(k+1) - zi(k) )
!!$       tkhz_full(k) = (tkhz(k_flip) * (zi(k+1) - z(k)) + tkhz(k_flip-1) * (z(k) - zi(k))) /&
!!$            ( zi(k+1) - zi(k) )
!!$    enddo
!!$
!!$
!!$    do k = 1, nzm
!!$       k_flip = nzm - k + 1
!!$       flip(k) = s_diff(k_flip)
!!$    enddo
!!$
!!$    call diffuse_scalar1D_rev_full(                                           &
!!$                                flip,                       sfx_local,        &
!!$                           tkhz_full,                       dtn_local )
!!$
!!$    do k = 1, nzm
!!$       k_flip = nzm - k + 1
!!$       s_diff(k_flip) = flip(k)
!!$       flip(k) = q_diff(k_flip)
!!$    enddo
!!$
!!$    call diffuse_scalar1D_rev_full(                                           &
!!$                                flip,                       qfx_local,        &
!!$                           tkhz_full,                       dtn_local )
!!$
!!$
!!$    do k = 1, nzm
!!$       k_flip = nzm - k + 1
!!$       q_diff(k_flip) = flip(k)
!!$       flip(k) = qp_diff(k_flip)
!!$    enddo
!!$
!!$    call diffuse_scalar1D_rev_full(                                           &
!!$                                flip,                             0.0,        &
!!$                           tkhz_full,                       dtn_local )
!!$
!!$    do k = 1, nzm
!!$       k_flip = nzm - k + 1
!!$       qp_diff(k_flip) = flip(k)
!!$       flip(k) = u_diff(k_flip)
!!$    enddo
!!$
!!$    call diffuse_scalar1D_rev_full(                                           &
!!$                                flip,                      taux_local,        &
!!$                            tkz_full,                       dtn_local )
!!$
!!$    do k = 1, nzm
!!$       k_flip = nzm - k + 1
!!$       u_diff(k_flip) = flip(k)
!!$       flip(k) = v_diff(k_flip)
!!$    enddo
!!$
!!$    call diffuse_scalar1D_rev_full(                                           &
!!$                                flip,                       tauy_local,        &
!!$                             tkz_full,                       dtn_local )
!!$
!!$    do k = 1, nzm
!!$       k_flip = nzm - k + 1
!!$       v_diff(k_flip) = flip(k)
!!$    enddo

! ------

!!$    do k = 1, nz
!!$       k_flip=nz-k+1
!!$       tkhz_half(k) = tkhz(k_flip)
!!$       tkz_half(k) = tkz(k_flip)
!!$    enddo
!!$
!!$    do k = 1, nzm
!!$       k_flip = nzm - k + 1
!!$       flip(k) = s_diff(k_flip)
!!$    enddo
!!$
!!$    call diffuse_scalar1D_rev_half(                                           &
!!$                                flip,                       sfx_local,        &
!!$                                tkhz_half,                       dtn_local )
!!$
!!$    do k = 1, nzm
!!$       k_flip = nzm - k + 1
!!$       s_diff(k_flip) = flip(k)
!!$       flip(k) = q_diff(k_flip)
!!$    enddo
!!$
!!$    call diffuse_scalar1D_rev_half(                                           &
!!$                                flip,                       qfx_local,        &
!!$                                tkhz_half,                       dtn_local )
!!$
!!$
!!$    do k = 1, nzm
!!$       k_flip = nzm - k + 1
!!$       q_diff(k_flip) = flip(k)
!!$       flip(k) = qp_diff(k_flip)
!!$    enddo
!!$
!!$    call diffuse_scalar1D_rev_half(                                           &
!!$                                flip,                             0.0,        &
!!$                                tkhz_half,                       dtn_local )
!!$
!!$    do k = 1, nzm
!!$       k_flip = nzm - k + 1
!!$       qp_diff(k_flip) = flip(k)
!!$       flip(k) = u_diff(k_flip)
!!$    enddo
!!$
!!$    call diffuse_scalar1D_rev_half(                                           &
!!$                                flip,                      taux_local,        &
!!$                                 tkz_half,                       dtn_local )
!!$
!!$    do k = 1, nzm
!!$       k_flip = nzm - k + 1
!!$       u_diff(k_flip) = flip(k)
!!$       flip(k) = v_diff(k_flip)
!!$    enddo
!!$
!!$    call diffuse_scalar1D_rev_half(                                           &
!!$                                flip,                       tauy_local,        &
!!$                                 tkz_half,                       dtn_local )
!!$
!!$    do k = 1, nzm
!!$       k_flip = nzm - k + 1
!!$       v_diff(k_flip) = flip(k)
!!$    enddo
    
    call ddiffuse_scalar1D(                                                   &
                              s_diff,                       sfx_local,        &
                                tkhz,                       dtn_local )

!!$    call diffuse_scalar1D(                                                    &
!!$                              s_diff,                       sfx_local,        &
!!$                                tkhz,                       dtn_local )

    call ddiffuse_scalar1D(                                                   &
                              q_diff,                       qfx_local,        &
                                tkhz,                       dtn_local )

    call ddiffuse_scalar1D(                                                   &
                             qp_diff,                            0.d0,        &
                                tkhz,                       dtn_local )

    call ddiffuse_scalar1D(                                                   &
                              u_diff,                      taux_local,        &
                                 tkz,                       dtn_local )

    call ddiffuse_scalar1D(                                                   &
                              v_diff,                      tauy_local,        &
                                 tkz,                       dtn_local )


    do k = 1, nzm
       k_flip = nz - k + 1
       tkz_full(k_flip-1) = (tkz(k_flip) * (zi(k+1) - z(k)) + tkz(k_flip-1) * (z(k) - zi(k))) /&
            ( zi(k+1) - zi(k) )
    enddo

!!$    tkz_full(1) = 2.0 * tkz(1)
!!$
!!$    do k = 1, nzm-1
!!$       tkz_full(k+1) =  max(0., 2.0*tkz(k) - tkz_full(k))
!!$    enddo
!!$    
!!$
!!$    if(maxval(tkz_full).gt.0.0) write(*,*) 'two second', tkz_full


    call ddiffuse_scalar1Dw(                                                  &
                              w_diff,                        tkz_full,        &
                           dtn_local )

  end subroutine diffuse_all

!  ----------------------------------------------------------------------------

  subroutine compute_eddy_diff( &
       pcols , pver   , ncol     , t      , qv       , ztodt , &
       ql    , qi     , s        , rpdel  , cldn     , qrl   , &
       z     , zi     , pmid     , pi     , u        , v     , &
       taux  , tauy   , shflx    , qflx   , wstarent , nturb , &
       ustar , pblh   , kvm_in   , kvh_in , kvm_out  , kvh_out , kvq   , & 
       cgh   , cgs    , tpert    , qpert  , tke      , bprod , &
       sprod , sfi    , qsat     , kvinit, wstar )

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    !  Interface routines for calcualtion and diatnostics of turbulence related
    !  coefficients
    !
    ! Author: B. Stevens (rewrite August 2000)
    ! 
    !-----------------------------------------------------------------------
    use diffusion_solver, only: compute_vdiff
    implicit none
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: pcols
    integer, intent(in) :: pver                       ! number of atmospheric layers
    integer, intent(in) :: ncol                       ! number of atmospheric columns

    real(r8), intent(in)  :: t(pcols,pver)            ! temperature (used for density)
    real(r8), intent(in)  :: qv(pcols,pver)! specific humidity [kg/kg]
    real(r8), intent(in)  :: ztodt                    ! 2 delta-t
    real(r8), intent(in)  :: ql(pcols,pver)! specific humidity [kg/kg]
    real(r8), intent(in)  :: qi(pcols,pver)! specific humidity [kg/kg]
    real(r8), intent(in)  :: s(pcols,pver)            ! dry static energy
    real(r8), intent(in)  :: rpdel(pcols,pver)        ! 1./pdel
    real(r8), intent(in)  :: cldn(pcols,pver)         ! cloud fraction
    real(r8), intent(in)  :: qrl(pcols,pver)          ! input qrl
    real(r8), intent(in)  :: z(pcols,pver)            ! midpoint height above surface [m]
    real(r8), intent(in)  :: zi(pcols,pver+1)         ! interface height above surface [m]
    real(r8), intent(in)  :: pmid(pcols,pver)         ! midpoint pressures
    real(r8), intent(in)  :: pi(pcols,pver+1)         ! interface pressures
    real(r8), intent(in)  :: u(pcols,pver)            ! zonal velocity
    real(r8), intent(in)  :: v(pcols,pver)            ! meridional velocity
    real(r8), intent(in)  :: taux(pcols)              ! zonal stress [N/m2]
    real(r8), intent(in)  :: tauy(pcols)              ! meridional stress [N/m2]
    real(r8), intent(in)  :: shflx(pcols)             ! sensible heat flux
    real(r8), intent(in)  :: qflx(pcols)              ! constituent flux
    logical,  intent(in)  :: wstarent
    integer,  intent(in)  :: nturb                    ! number of predictor-corrector 'time steps'
    logical,  intent(in)  :: kvinit                   ! 'true' means time step = 1 : need to initialize kvh, kvm (uses kvf or zero)
    real(r8), intent(in)  :: kvm_in(pcols,pver+1)     ! kvm saved from last timestep
    real(r8), intent(in)  :: kvh_in(pcols,pver+1)     ! kvh saved from last timestep
    real(r8), intent(out) :: kvm_out(pcols,pver+1)    ! eddy diffusivity for momentum [m2/s]
    real(r8), intent(out) :: kvh_out(pcols,pver+1)    ! eddy diffusivity for heat [m2/s]
    real(r8), intent(out) :: kvq(pcols,pver+1)        ! eddy diffusivity for constituents [m2/s] (note not having '_out')
    real(r8), intent(out) :: ustar(pcols)             ! surface friction velocity [m/s]
    real(r8), intent(out) :: pblh(pcols)              ! boundary-layer height [m]
    real(r8), intent(out) :: cgh(pcols,pver+1)        ! counter-gradient term for heat [J/kg/m]
    real(r8), intent(out) :: cgs(pcols,pver+1)        ! counter-gradient star (cg/flux)
    real(r8), intent(out) :: tpert(pcols)             ! convective temperature excess
    real(r8), intent(out) :: qpert(pcols)             ! convective humidity excess
    real(r8), intent(out) :: tke(pcols,pver+1)        ! Turbulent kinetic energy (m2s-2)
    real(r8), intent(out) :: bprod(pcols,pver+1) 
    real(r8), intent(out) :: sprod(pcols,pver+1) 
    real(r8), intent(out) :: sfi(pcols,pver+1)        ! interfacial saturation fraction
    real(r8), intent(out) :: wstar 
    integer, external     :: qsat
    !
    !---------------------------Local workspace-----------------------------
    !
    real(r8) :: kvf(pcols,pver+1)       ! free atmospheric eddy diffsvty [m2/s]
    real(r8) :: kvm(pcols,pver+1)       ! eddy diffusivity for momentum [m2/s]
    real(r8) :: kvh(pcols,pver+1)       ! eddy diffusivity for heat [m2/s]
    real(r8) :: s2(pcols,pver)          ! shear squared, defined at interfaces except surface
    real(r8) :: n2(pcols,pver)          ! brunt vaisaila frequency, defined at interfaces except surface
    real(r8) :: ri(pcols,pver)          ! richardson number: n2/s2, defined at interfaces except surface

    !cb The following variables are used only by gbturb.

    real(r8) :: qt(pcols,pver)          ! total specific humidity [kg/kg]
    real(r8) :: sfuh(pcols,pver)        ! sat frac in upper half-layer
    real(r8) :: sflh(pcols,pver)        ! sat frac in lower half-layer
    real(r8) :: sl(pcols,pver)          ! liquid water static energy
    real(r8) :: slv(pcols,pver)         ! virtual liquid water static energy
    real(r8) :: slslope(pcols,pver)     ! sl slope in thermo layer
    real(r8) :: qtslope(pcols,pver)     ! qt slope in thermo layer
    real(r8) :: rrho(pcols)             ! density at the lowest grid point
    real(r8) :: qvfd(pcols,pver)        ! specific humidity for diffusion
    real(r8) :: tfd(pcols,pver)         ! temperature for diffusion [kg/kg]
    real(r8) :: slfd(pcols,pver)        ! liquid static energy (conserved during diffusion)
    real(r8) :: qtfd(pcols,pver)        ! total water (conserved during diffusion)
    real(r8) :: qlfd(pcols,pver)        ! liquid water for diffusion
    real(r8) :: ufd(pcols,pver)         ! u-wind for diffusion [kg/kg]
    real(r8) :: vfd(pcols,pver)         ! v-wind for diffusion [kg/kg]

    ! Coefficients in buoyancy flux calculation

    real(r8) :: chu(pcols,pver+1)       ! heat var. coef for dry states, defined at layer midpoints & surface
    real(r8) :: chs(pcols,pver+1)       ! heat var. coef for sat states, defined at layer midpoints & surface 
    real(r8) :: cmu(pcols,pver+1)       ! mstr var. coef for dry states, defined at layer midpoints & surface
    real(r8) :: cms(pcols,pver+1)       ! mstr var. coef for sat states, defined at layer midpoints & surface
    real(r8) :: jnk1d(pcols)
    real(r8) :: jnk2d(pcols,pver+1)  
    integer iturb, i, k, status
    real(r8) :: zero(pcols)
    real(r8) :: zero2d(pcols,pver+1)
    real(r8) :: es(1)                   ! saturation vapor pressure
    real(r8) :: qs(1)                   ! saturation spec. humidity
    real(r8) :: gam(1)                  ! (L/cp)*dqs/dT
    real(r8) :: ep2, templ, temps
    character(128) :: errstring         ! error status for compute_vdiff
    real(r8) :: minpblh(pcols)          ! minimum pbl height based on surface stress
    
    zero(:) = 0.
    zero2d(:,:) = 0.

    ufd(:ncol,:)  = u(:ncol,:)
    vfd(:ncol,:)  = v(:ncol,:)
    tfd(:ncol,:)  = t(:ncol,:)
    qvfd(:ncol,:) = qv(:ncol,:) !
    
    qlfd(:ncol,:) = ql(:ncol,:)
    do iturb = 1, nturb
       !
       ! Calculate (qt,sl,n2,s2,ri) from a given set of (t,qv,ql,qi,u,v)
       ! 
       call trbintd( &
            pcols , pver , ncol  , z       , ufd     , vfd     , tfd   , pmid    , &
            taux  , tauy , ustar , rrho    , s2      , n2      , ri    , zi      , &
            pi    , cldn , qtfd  , qvfd    , qlfd    , qi      , sfi   , sfuh    , &
            sflh  , slfd , slv   , slslope , qtslope , chs     , chu   , cms     , &
            cmu   , minpblh, qsat )

       ! Save initial (i.e., before iterative diffusion) profile of (qt,sl) at each time step         
       ! Only necessary for (qt,sl) not (u,v) because (qt,sl) are newly calculated variables. 
       if(iturb.eq.1) then
          qt(:ncol,:) = qtfd(:ncol,:)
          sl(:ncol,:) = slfd(:ncol,:)
       endif
       !
       ! Get free atmosphere exchange coefficients
       !
       call austausch_atm(pcols, pver, ncol    ,ri      ,s2      ,kvf)
       !
       ! Initialize kvh/kvm to send to caleddy, depending on model timestep and iteration number
       !
       if ( iturb .eq. 1 ) then
          if ( kvinit ) then
             ! first iteration of first model timestep, using free tropospheric value
             if ( use_kvf ) then
                kvh(:ncol,:) = kvf(:ncol,:)
                kvm(:ncol,:) = kvf(:ncol,:)
             else
                kvh(:ncol,:) = 0._r8
                kvm(:ncol,:) = 0._r8
             endif
          else
             ! first iteration on any model timestep except 1, use value from previous timestep
             kvh(:ncol,:) = kvh_in(:ncol,:)
             kvm(:ncol,:) = kvm_in(:ncol,:)
          endif
       else
          ! use value from previous iteration
          kvh(:ncol,:) = kvh_out(:ncol,:)
          kvm(:ncol,:) = kvm_out(:ncol,:)
       endif
       !
       ! Calculate eddy diffusivity (kvh_out,kvm_out) and (tke,bprod,sprod) using
       ! a given (kvh,kvm) which are used only for initializing (bprod,sprod)  at
       ! the first part of caleddy. (bprod,sprod) are fully updated at the end of
       ! caleddy after calculating (kvh_out,kvm_out) 
       !  
       call caleddy(pcols, pver, ncol    , &
            slfd    ,qtfd    ,qlfd    ,slv     ,ufd     , &
            vfd     ,pi      ,z       ,zi      , &
            qflx    ,shflx   ,slslope ,qtslope , &
            chu     ,chs     ,cmu     ,cms     ,sfuh    , &
            sflh    ,n2      ,s2      ,ri      ,    rrho, &
            pblh    ,ustar   , &
            kvh     ,kvm     ,kvh_out ,kvm_out , &
            tpert   ,qpert   ,qrl     ,kvf, tke, wstarent,bprod, sprod, minpblh,wstar)


!      Check if iteration process at each time step produces convergent solution.
!      do k = 1, pver + 1
!!$      do k = pver - 10, pver + 1
!!$        write(6,*)  iturb, k,kvh(1,k), kvh_out(1,k)
!!$      end do


       ! Eddy diffusivities which will be used for the initialization of (bprod,
       ! sprod) in 'caleddy' at the next iteration step. 
        if(iturb.gt.1) then
           kvm_out(:ncol,:) = lambda*kvm_out(:ncol,:)+(1-lambda)*kvm(:ncol,:)
           kvh_out(:ncol,:) = lambda*kvh_out(:ncol,:)+(1-lambda)*kvh(:ncol,:)     
        endif

       ! Set nonlocal terms to zero for flux diagnostics, since not used by caleddy.
       cgh(:ncol,:) = 0.
       cgs(:ncol,:) = 0.      

       if (iturb.lt.nturb) then
          ! Each time we diffuse the original state
          slfd(:ncol,:)  = sl(:ncol,:)
          qtfd(:ncol,:)  = qt(:ncol,:)
          ufd(:ncol,:)   = u(:ncol,:)
          vfd(:ncol,:)   = v(:ncol,:)
          !
          ! Diffuse initial profile of each time step using a given (kvh_out,kvm_out)
          ! (slfd,qtfd,ufd,vfd) : in and out variables.
          !
          call compute_vdiff(                                                   &
               pcols      , pver     , 1               , ncol      , pmid      , &
               pi         , rpdel    , t               , ztodt     , taux      , &
               tauy       , shflx    , qflx            , ntop_turb , nbot_turb , &
               kvh_out    , kvm_out  , kvh_out         , cgs       , cgh       , &
               zi         , zero     , zero            , fieldlist , &
               ufd        , vfd      , qtfd            , slfd      , &
               jnk1d      , jnk1d    , jnk2d           , jnk1d     , errstring , &
               .false.    )

          ! Retrieve (tfd,qvfd,qlfd) from (slfd,qtfd) in order to use 'trbintd'
          do k = 1,pver
             do i = 1, ncol
                templ = ( slfd(i,k) - ggr*z(i,k) ) / cp
                status = qsat(templ,pmid(i,k),es(1),qs(1),gam(1),1)
                ep2    = .622 ! shr_const_mwwv/shr_const_mwdair
                temps=templ + (qtfd(i,k)-qs(1)) / ( cp/lcond + &
                     ep2 * lcond * qs(1) / (rgas * templ**2))
                status = qsat(temps,pmid(i,k),es(1),qs(1),gam(1),1)
! debug
                qlfd(i,k) = max( qtfd(i,k) - qi(i,k) - qs(1) ,0.)
                qvfd(i,k) = qtfd(i,k) - qi(i,k) - qlfd(i,k)
                tfd(i,k) = (slfd(i,k) + lcond*qlfd(i,k) + lsub*qi(i,k) - ggr*z(i,k)) / cp
             end do
          end do
       endif
    end do  ! iturb

    kvq(:ncol,:) = kvh_out(:ncol,:)


!   Check if previous' time step's kvh is successfully transfered 
!                           to the current time step   
!   do k = 1, pver + 1
!     write(6,*)  k, kvh_in(:ncol,k), kvh_out(:ncol,k)   
!   end do

   
    return
  end subroutine compute_eddy_diff
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !==============================================================================
  subroutine sfdiag(pcols, pver, ncol    ,qt      ,ql      ,sl      ,pi      , &
       pm      ,zi      ,cld     ,sfi     ,sfuh    , &
       sflh    ,slslope ,qtslope , qsat)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    !  Interface for computing saturation fraction at upper and lower-half layers,
    !  and interfaces for use by turbulence scheme
    ! 
    ! Method: 
    !
    ! Author: B. Stevens and C. Bretherton (August 2000)
    ! 
    !-----------------------------------------------------------------------
    implicit none
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: pcols             ! number of atmospheric columns   
    integer, intent(in) :: pver              ! number of atmospheric layers   
    integer, intent(in) :: ncol              ! number of atmospheric columns   

    real(r8), intent(in)  :: sl(pcols,pver)      ! liquid water statis energy
    real(r8), intent(in)  :: qt(pcols,pver)      ! total water spec hum.
    real(r8), intent(in)  :: ql(pcols,pver)      ! liquid water spec hum.
    real(r8), intent(in)  :: pi(pcols,pver+1)    ! Interface pressures
    real(r8), intent(in)  :: pm(pcols,pver)      ! Midpoint pressures
    real(r8), intent(in)  :: zi(pcols,pver+1)    ! Interface heights
    real(r8), intent(in)  :: cld(pcols,pver)     ! cloud fraction
    real(r8), intent(in)  :: slslope(pcols,pver) ! sl slope in thermo layer
    real(r8), intent(in)  :: qtslope(pcols,pver) ! qt slope in thermo layer
    integer, external     :: qsat
    !
    ! Output arguments
    !
    real(r8), intent(out) :: sfi(pcols,pver+1)   ! interfacial layer sat frac
    real(r8), intent(out) :: sfuh(pcols,pver)    ! sat frac in upper half-layer
    real(r8), intent(out) :: sflh(pcols,pver)    ! sat frac in lower half-layer
    !
    !---------------------------Local workspace-----------------------------
    !
    real(r8) :: sltop, slbot           ! sl at top/bot of grid layer
    real(r8) :: qttop, qtbot           ! qt at top/bot of grid layer
    real(r8) :: tltop(1), tlbot(1)     ! liq wat temp at top/bot of grid layer
    real(r8) :: qxtop, qxbot           ! sat excess at top/bot of grid layer
    real(r8) :: qxm                    ! sat excess at midpoint
    real(r8) :: es(1)                  ! saturation vapor pressure
    real(r8) :: qs(1)                  ! saturation spec. humidity
    real(r8) :: gam(1)                 ! (L/cp)*dqs/dT

    integer  :: i                      ! longitude index
    integer  :: k                      ! vertical index
    integer  :: km1                    ! k-1
    integer  :: status                 ! status returned by function calls

    !
    !-----------------------------------------------------------------------
    !
    sfi(1:ncol,:) = 0.
    sfuh(1:ncol,:) = 0.
    sflh(1:ncol,:) = 0.
    select case (sftype)
    case default
       !
       ! use cloud fraction ('horizontal' cloud partitioning)
       !
       do k = ntop_turb+1, nbot_turb
          km1 = k-1
          do i=1,ncol
             sfuh(i,k) = cld(i,k)
             sflh(i,k) = cld(i,k)
             sfi(i,k) = 0.5*(sfuh(i,k)+sflh(i,km1))
          end do
       end do
    case ('l')
       !
       ! use modified cloud fraction partitioning
       !
       do k = ntop_turb+1, nbot_turb
          km1 = k-1
          do i=1,ncol
             sfuh(i,k) = cld(i,k)
             sflh(i,k) = cld(i,k)
             if(ql(i,k) .lt. qmin) then
                sfuh(i,k) = 0.
                sflh(i,k) = 0.
             end if
             sfi(i,k) = 0.5*(sflh(i,km1) + min(sfuh(i,k),sflh(i,km1)))
          end do
       end do
    case ('u')
       !
       ! use unsaturated buoyancy - since sfi, sfuh, sflh have already been zeroed
       ! nothing more need be done for this case.
       !
    case ('z')

       ! Calculate saturation fraction based on whether the air just above or
       ! just below the interface is saturated, i.e. with vertical cld partitioning
       ! The saturation fraction of the interfacial layer between midpoints k and k+1
       ! is computed by averaging the
       ! saturation fraction of the half-layers above and below the interface, with
       ! a special provision for cloud tops (more cloud in the half-layer below than 
       ! in the half-layer above).  In each half-layer, vertical partitioning of
       ! cloud based on the slopes diagnosed above is used.

       ! Loop down through the layers, computing the saturation fraction in each 
       ! half-layer (sfuh for upper half, sflh for lower half). Once sfuh(i,k) is 
       ! computed, use with sflh(i,k-1) to determine the 
       ! saturation fraction sfi(i,k) for interfacial layer k-0.5.

       do k = ntop_turb+1, nbot_turb
          km1 = k-1
          do i = 1, ncol

             ! Compute saturation excess at top and bottom of layer k

             sltop = sl(i,k) + slslope(i,k)*(pi(i,k)-pm(i,k))      
             qttop = qt(i,k) + qtslope(i,k)*(pi(i,k)-pm(i,k))
             tltop(1) = (sltop - ggr*zi(i,k))/cp  ! liquid water temp.
             status = qsat(tltop(1),pi(i,k),es(1),qs(1),gam(1),1)
             qxtop = qttop - qs(1)   ! saturation excess at top of layer k

             slbot = sl(i,k) + slslope(i,k)*(pi(i,k+1)-pm(i,k))      
             qtbot = qt(i,k) + qtslope(i,k)*(pi(i,k+1)-pm(i,k))
             tlbot(1) = (slbot - ggr*zi(i,k+1))/cp  ! liquid water temp.
             status = qsat(tlbot(1),pi(i,k+1), es(1), qs(1), gam(1),1)
             qxbot = qtbot - qs(1)   ! saturation excess at bot of layer k
             ! saturation excess at midpt of layer k 
             qxm = qxtop + (qxbot-qxtop)*(pm(i,k)-pi(i,k)) / (pi(i,k+1)-pi(i,k))

             ! Next, find the saturation fraction sfuh(i,k) of the upper half of layer k.

             if ((qxtop.lt.0.) .and. (qxm.lt.0.)) then
                sfuh(i,k) = 0   ! Upper half-layer unsaturated
             else if ((qxtop.gt.0.) .and. (qxm.gt.0.)) then
                sfuh(i,k) = 1   ! Upper half-layer fully saturated
             else ! Either qxm < 0 and qxtop > 0 or vice versa
                sfuh(i,k) = max(qxtop,qxm)/abs(qxtop - qxm)
             end if

             ! Combine with sflh(i) (still for layer k-1) to get interfac layer sat frac

             sfi(i,k) = 0.5*(sflh(i,k-1)+min(sflh(i,k-1),sfuh(i,k)))

             ! Update sflh to be for the lower half of layer k.             

             if ((qxbot.lt.0.) .and. (qxm.lt.0.)) then
                sflh(i,k) = 0   ! Upper half-layer unsaturated
             else if ((qxbot.gt.0.) .and. (qxm.gt.0.)) then
                sflh(i,k) = 1   ! Upper half-layer fully saturated
             else ! Either qxm < 0 and qxbot > 0 or vice versa
                sflh(i,k) = max(qxbot,qxm)/abs(qxbot - qxm)
             end if
          end do  ! i
       end do ! k
       do i = 1,ncol
          sfi(i,pver+1) = sflh(i,pver)  ! Sat frac in lowest half-layer. 
       end do
    end select

    return
  end subroutine sfdiag
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !===============================================================================
  subroutine trbintd(pcols, pver, ncol    ,                            &
       z       ,u       ,v       , &
       t       ,pmid    ,taux    , &
       tauy    ,ustar   , rrho , &
       s2      ,n2      ,ri      ,          &
       zi      ,pi      ,cld     ,                   &
       qt      ,qv      ,ql,qi      ,sfi     ,sfuh    , &
       sflh    ,sl      ,slv     ,slslope ,qtslope , &
       chs     ,chu     ,cms     ,cmu     , minpblh, qsat)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    !  Time dependent initialization
    ! 
    ! Method: 
    !  Diagnosis of variables that do not depend on mixing assumptions or
    !  PBL depth
    !
    ! Author: B. Stevens (extracted from pbldiff, August, 2000)
    ! 
    !-----------------------------------------------------------------------
    implicit none
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: pcols                 ! number of atmospheric columns   
    integer, intent(in) :: pver                  ! number of atmospheric layers   
    integer, intent(in) :: ncol                  ! number of atmospheric columns

    real(r8), intent(in)  :: z(pcols,pver)       ! midpoint height above surface [m]
    real(r8), intent(in)  :: u(pcols,pver)       ! midpoint windspeed x-direction [m/s]
    real(r8), intent(in)  :: v(pcols,pver)       ! midpoint windspeed y-direction [m/s]
    real(r8), intent(in)  :: t(pcols,pver)       ! midpoint temperature (used for density)
    real(r8), intent(in)  :: pmid(pcols,pver)    ! midpoint pressures
    real(r8), intent(in)  :: taux(pcols)         ! surface u stress [N/m2]
    real(r8), intent(in)  :: tauy(pcols)         ! surface v stress [N/m2]
    real(r8), intent(in)  :: zi(pcols,pver+1)    ! interface height [m]
    real(r8), intent(in)  :: pi(pcols,pver+1)    ! interface pressure
    real(r8), intent(in)  :: cld(pcols,pver)     ! cloud fraction
    integer, external     :: qsat

    !
    ! Output arguments
    !

    real(r8), intent(out) :: ustar(pcols)      ! surface friction velocity [m/s]
    real(r8), intent(out) :: s2(pcols,pver)    ! interfacial (except surface) shear squared
    real(r8), intent(out) :: n2(pcols,pver)    ! interfacial (except surface) brunt vaisaila frequency
    real(r8), intent(out) :: ri(pcols,pver)    ! interfacial (except surface) richardson number: n2/s2
    !cb gbturb variables
    real(r8), intent(in) :: qv(pcols,pver)     ! wtr vapor specific humidity
    real(r8), intent(in) :: ql(pcols,pver)     ! liq water specific humidity
    real(r8), intent(in) :: qi(pcols,pver)     ! liq water specific humidity
    real(r8), intent(out) :: qt(pcols,pver)    ! tot water specific humidity
    real(r8), intent(out) :: sfi(pcols,pver+1) ! interfacial sat frac
    real(r8), intent(out) :: sfuh(pcols,pver)  ! sat frac in upper half-layer
    real(r8), intent(out) :: sflh(pcols,pver)  ! sat frac in lower half-layer
    real(r8), intent(out) :: sl(pcols,pver)    ! liquid water static energy
    real(r8), intent(out) :: slv(pcols,pver)   ! virtual liq water static energy
    ! coefficients in b-flux calc
    real(r8), intent(out) :: chu(pcols,pver+1)   ! heat var. coef for dry states at midpoints and surface
    real(r8), intent(out) :: chs(pcols,pver+1)   ! heat var. coef for sat states at midpoints and surface
    real(r8), intent(out) :: cmu(pcols,pver+1)   ! mstr var. coef for dry states at midpoints and surface
    real(r8), intent(out) :: cms(pcols,pver+1)   ! mstr var. coef for sat states at midpoints and surface
    real(r8), intent(out) :: slslope(pcols,pver) ! sl slope in thermo layer
    real(r8), intent(out) :: qtslope(pcols,pver) ! qt slope in thermo layer
    real(r8), intent(out) :: rrho(pcols)         ! 1./bottom level density (temporary)
    real(r8), intent(out) :: minpblh(pcols)      ! minimum pbl height based on surface stress
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i                        ! longitude index
    integer  :: k, km1                   ! level index
    integer  :: status                   ! status returned by function calls
    !
    ! gbturb local variables
    !
    real(r8) :: qs(pcols,pver)           ! saturation specific humidity
    real(r8) :: es(pcols,pver)           ! saturation vapor pressure
    real(r8) :: gam(pcols,pver)          ! (l/cp)*(d(qs)/dT)
    real(r8) :: rdz                      ! 1 / (delta z) between midpoints
    real(r8) :: dsldz                    ! delta sl / delta z at interface
    real(r8) :: dqtdz                    ! delta qt / delta z at interface
    real(r8) :: ch                       ! sfi weighted ch at the interface
    real(r8) :: cm                       ! sfi weighted cm at the interface
    real(r8) :: bfact                    ! buoyancy factor in n2 calculations
    real(r8) :: product                  ! intermediate vars used to find slopes
    real(r8) :: dsldp_a, dqtdp_a         ! slopes across interface above 
    real(r8) :: dsldp_b(pcols), dqtdp_b(pcols) ! slopes across interface below

    !
    ! Compute ustar, and kinematic surface fluxes from surface energy fluxes
    !
    do i=1,ncol
       rrho(i)   = rgas*t(i,pver)/pmid(i,pver)
       ustar(i)  = max(sqrt(sqrt(taux(i)**2 + tauy(i)**2)*rrho(i)),ustar_min)
       minpblh(i) = 100.0 * ustar(i) ! same as HB eddy diff code
    end do

    !
    ! Compute shear squared (s2), brunt vaisaila frequency (n2) and related richardson
    ! number (ri).  For the n2 calcualtion use the dry theta_v derived from virtem
    !

    ! calculate moist thermodynamic variables at midpoints

    do k=ntop_turb,nbot_turb
       status = qsat(t(1,k),pmid(1,k),es(1,k),qs(1,k),gam(1,k),ncol)
       do i=1,ncol
          !
          ! moisture variable partitioning
          !
          qt(i,k)  = qv(i,k) + ql(i,k) + qi(i,k) 

          ! Static energies

          sl(i,k)  = cp*t(i,k) + ggr*z(i,k) - lcond*ql(i,k) - lsub*qi(i,k)

          slv(i,k) = sl(i,k)*(1. + zvir*qt(i,k))

          ! Thermodynamic coefficients for buoyancy flux - in this loop these
          ! are calculated at midpoints; later they will be averaged to interfaces,
          ! where they will ultimately be used. At the surface, the coefficients are
          ! taken from the lowest midpoint.

          bfact    = ggr/(t(i,k)*(1.+zvir*qv(i,k) - ql(i,k)-qi(i,k)))
          chu(i,k) = (1. + zvir * qt(i,k)) * bfact/cp
          chs(i,k) = ((1. + (1. + zvir)*gam(i,k)*cp*t(i,k)/lcond) &
               / (1. + gam(i,k)))   * bfact/cp
          cmu(i,k) = zvir * bfact * t(i,k)
          cms(i,k) = lcond*chs(i,k)  -  bfact * t(i,k)
       end do
    end do
    do i=1,ncol
       chu(i,pver+1) = chu(i,pver)
       chs(i,pver+1) = chs(i,pver)
       cmu(i,pver+1) = cmu(i,pver)
       cms(i,pver+1) = cms(i,pver)
    end do

    ! Compute slopes in conserved variables sl, qt within thermo layer k. 
    ! a indicates the 'above' gradient from layer k-1 to layer k and 
    ! b indicates the 'below' gradient from layer k   to layer k+1.
    ! We take the smaller (in absolute value) of these gradients
    ! as the slope within layer k. If they have opposite sign, gradient in 
    ! layer k is taken to be zero.

    do i = 1,ncol

       ! Slopes at endpoints determined by extrapolation

       slslope(i,pver) = (sl(i,pver)-sl(i,pver-1))/(pmid(i,pver)-pmid(i,pver-1))
       qtslope(i,pver) = (qt(i,pver)-qt(i,pver-1))/(pmid(i,pver)-pmid(i,pver-1))
       slslope(i,1) = (sl(i,2)-sl(i,1))/(pmid(i,2)-pmid(i,1))
       qtslope(i,1) = (qt(i,2)-qt(i,1))/(pmid(i,2)-pmid(i,1))
       dsldp_b(i) = slslope(i,1)
       dqtdp_b(i) = qtslope(i,1)
    end do

    do k = 2,pver-1
       do i = 1,ncol
          dsldp_a = dsldp_b(i)
          dqtdp_a = dqtdp_b(i)
          dsldp_b(i) = (sl(i,k+1)-sl(i,k))/(pmid(i,k+1)-pmid(i,k))
          dqtdp_b(i) = (qt(i,k+1)-qt(i,k))/(pmid(i,k+1)-pmid(i,k))
          product = dsldp_a*dsldp_b(i)
          if (product .le. 0.) then 
             slslope(i,k) = 0.
          else if (product.gt.0. .and. dsldp_a.lt.0.) then 
             slslope(i,k) = max(dsldp_a,dsldp_b(i))
          else if (product.gt.0. .and. dsldp_a.gt.0.) then 
             slslope(i,k) = min(dsldp_a,dsldp_b(i))
          end if

          product = dqtdp_a*dqtdp_b(i)
          if (product .le. 0.) then 
             qtslope(i,k) = 0.
          else if (product.gt.0. .and. dqtdp_a.lt.0.) then 
             qtslope(i,k) = max(dqtdp_a,dqtdp_b(i))
          else if (product.gt.0. .and. dqtdp_a.gt.0.) then 
             qtslope(i,k) = min(dqtdp_a,dqtdp_b(i))
          end if
       end do ! i
    end do ! k

    !  Compute saturation fraction in the interfacial layers for use in buoyancy
    !  flux computation.

    call sfdiag(pcols, pver, ncol    ,qt      ,ql      ,sl      ,pi      , &
         pmid      ,zi      ,cld     ,sfi     ,sfuh    , &
         sflh    ,slslope ,qtslope ,qsat)

    do k = nbot_turb,ntop_turb+1,-1
       km1 = k-1
       do i = 1,ncol
          rdz     = 1./(z(i,km1) - z(i,k))
          dsldz = (sl(i,km1) - sl(i,k)) * rdz
          dqtdz = (qt(i,km1) - qt(i,k)) * rdz 
          chu(i,k)  = (chu(i,km1) + chu(i,k))*0.5
          chs(i,k)  = (chs(i,km1) + chs(i,k))*0.5
          cmu(i,k)  = (cmu(i,km1) + cmu(i,k))*0.5
          cms(i,k)  = (cms(i,km1) + cms(i,k))*0.5
          ch = chu(i,k)*(1.-sfi(i,k)) + chs(i,k)*sfi(i,k)
          cm = cmu(i,k)*(1.-sfi(i,k)) + cms(i,k)*sfi(i,k)
          n2(i,k) = ch*dsldz +  cm*dqtdz
          s2(i,k) = ((u(i,km1)-u(i,k))**2 + (v(i,km1)-v(i,k))**2)*rdz**2
          s2(i,k) = max(ntzero,s2(i,k))
          ri(i,k) = n2(i,k)/s2(i,k)
       end do !i
    end do !k

    return
  end subroutine trbintd
  !
  !===============================================================================
  subroutine austausch_atm(pcols,pver, ncol    ,ri      ,s2      ,kvf     )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    !  Computes exchange coefficients for free turbulent flows. 
    ! 
    ! Method: 
    !
    ! The free atmosphere diffusivities are based on standard mixing length
    ! forms for the neutral diffusivity multiplied by functns of Richardson
    ! number. K = l^2 * |dV/dz| * f(Ri). The same functions are used for
    ! momentum, potential temperature, and constitutents. 
    !
    ! The stable Richardson num function (Ri>0) is taken from Holtslag and
    ! Beljaars (1989), ECMWF proceedings. f = 1 / (1 + 10*Ri*(1 + 8*Ri))
    ! The unstable Richardson number function (Ri<0) is taken from  CCM1.
    ! f = sqrt(1 - 18*Ri)
    ! 
    ! Author: B. Stevens (rewrite, August 2000)
    ! 
    !-----------------------------------------------------------------------
    implicit none
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: pcols              ! number of atmospheric columns   
    integer, intent(in) :: pver              ! number of atmospheric columns   
    integer, intent(in) :: ncol                     ! number of atmospheric columns

    real(r8), intent(in)  ::  s2(pcols,pver)        ! shear squared
    real(r8), intent(in)  ::  ri(pcols,pver)        ! richardson no
    !
    ! Output arguments
    !
    real(r8), intent(out) :: kvf(pcols,pver+1)       ! coefficient for heat and tracers
    !
    !---------------------------Local workspace-----------------------------
    !
    real(r8) :: fofri                  ! f(ri)
    real(r8) :: kvn                    ! neutral Kv

    integer  :: i                      ! longitude index
    integer  :: k                      ! vertical index
    !
    !-----------------------------------------------------------------------
    !
    ! Make sure everything is set = 0., changes made below
    !
    kvf(:ncol,:) = 0.0
    !
    ! The surface diffusivity is always zero
    !
    kvf(:ncol,pver+1) = 0.0
    !
    ! Set the vertical diffusion coefficient above the top diffusion level
    ! Note that nbot_turb != pver is not supported
    !
    kvf(:ncol,1:ntop_turb) = 0.0
    !
    ! Compute the free atmosphere vertical diffusion coefficients: kvh = kvq = kvm. 
    !
    do k = ntop_turb+1, nbot_turb
       do i=1,ncol
          if (ri(i,k) < 0.0) then
             fofri = sqrt(max(1. - 18.*ri(i,k),0._r8))
          else 
             fofri = 1.0/(1.0 + 10.0*ri(i,k)*(1.0 + 8.0*ri(i,k)))    
          end if
          kvn = ml2(k)*sqrt(s2(i,k))
          kvf(i,k) = max(zkmin,kvn*fofri)
       end do
    end do

    return
  end subroutine austausch_atm
  !
  !=============================================================================
  !
  ! Grenier-Bretherton turbulence subroutines follow
  !-----------------------------------------------------------------------------
  subroutine caleddy(pcols, pver, ncol    , &
       sl      ,qt      ,ql      ,slv     ,u       , &
       v       ,pi      ,z      ,zi      , &
       qflx    ,shflx    ,slslope ,qtslope , &
       chu     ,chs     ,cmu     ,cms     ,sfuh, &
       sflh    ,n2      ,s2      ,ri      , rrho, &
       pblh    ,ustar, &
       kvh_in  ,kvm_in  ,kvh     ,kvm     , &
       tpert   ,qpert   ,qrl     ,kvf, tke, wstarent,bprod, sprod, minpblh,wstar)
    !-----------------------------------------------------------------------
    ! Driver routine to compute eddy diffusion coefficients for momentum,
    ! moisture, trace constituents and static energy.  Uses first order
    ! closure for stable turbulent layers. For convective layers, an entrainment
    ! closure is used, coupled to a diagnosis of layer-average TKE from the 
    ! instantaneous thermodynamic and velocity profiles. Convective layers are 
    ! diagnosed by
    ! extending layers of moist static instability into adjacent weakly stably
    ! stratified interfaces, stopping if the stability is too strong.  This
    ! allows a realistic depiction of dry convective boundary layers with a 
    ! downgradient approach.
    !
    ! NOTE: This routine currently assumes ntop_turb = 1, nbot_turb = pver
    ! (turbulent diffusivities computed at all interior interfaces) and will 
    ! require modification to handle a different ntop_turb.
    ! 
    ! Authors :  Herve Grenier, 06/2000, Chris Bretherton 09/2000
    ! Debugged : Sungsu Park, 06/2005.
    !----------------------------------------------------------------------------
    ! Inputs
    !
    implicit none
    integer, intent(in) :: pcols             ! number of atmospheric columns   
    integer, intent(in) :: pver              ! number of atmospheric layers   
    integer, intent(in) :: ncol              ! number of atmospheric columns   

    real(r8), intent(in) :: u(pcols,pver)    ! u wind input
    real(r8), intent(in) :: v(pcols,pver)    ! v wind input
    real(r8), intent(in) :: sl(pcols,pver)   ! liq water static energy
    !  cp*T+g*z-L*ql
    real(r8), intent(in) :: slv(pcols,pver)  ! liq water virtual static energy
    !  sl*(1 + .608*qt)
    real(r8), intent(in) :: qt(pcols,pver)   ! total spec. humidity  qv + ql 
    real(r8), intent(in) :: ql(pcols,pver)   ! liq water spec. humidity
    real(r8), intent(in) :: pi(pcols,pver+1) ! interface pressures
    real(r8), intent(in) :: z(pcols,pver)    ! layer midpoint height above sfc
    real(r8), intent(in) :: zi(pcols,pver+1) ! interface height above sfc
    ! Coefficients used in computing buoyancy fluxes
    real(r8), intent(in) :: chu(pcols,pver+1) ! Unsat sl (heat) coef. at midpoints and surface
    real(r8), intent(in) :: chs(pcols,pver+1) ! Sat sl (heat) coef. at midpoints and surface
    real(r8), intent(in) :: cmu(pcols,pver+1) ! Unsat qt (moisture) coef. at midpoints and surface
    real(r8), intent(in) :: cms(pcols,pver+1) ! Sat qt (moisture) coef. at midpoints and surface
    real(r8), intent(in) :: sfuh(pcols,pver)  ! sat frac in upper half-layer
    real(r8), intent(in) :: sflh(pcols,pver)  ! sat frac in lower half-layer
    !
    real(r8), intent(in) :: n2(pcols,pver)       ! interfacial (except sfc) moist squared buoy freq (s-2)
    real(r8), intent(in) :: s2(pcols,pver)       ! interfacial (except sfc) squared deformation (s-2)
    real(r8), intent(in) :: ri(pcols,pver)       ! interfacial (except sfc) gradient Richardson number
    real(r8), intent(in) :: qflx(pcols)          ! kinematic surf constituent flux (kg/m2/s)
    real(r8), intent(in) :: shflx(pcols)         ! kinematic surface heat flux 
    real(r8), intent(in) :: slslope(pcols,pver)  ! sl slope in thermo layer
    real(r8), intent(in) :: qtslope(pcols,pver)  ! qt slope in thermo layer
    real(r8), intent(in) :: qrl(pcols,pver)      ! LW heating rate (K/s) * cp *dp (W/kg*Pa)
    real(r8), intent(in) :: ustar(pcols)         ! surface friction velocity
    real(r8), intent(in) :: rrho(pcols)          ! 1./bottom midpoint density
    real(r8), intent(in) :: kvf(pcols,pver+1)    ! free atmosphere eddy diffusivity
    logical,  intent(in) :: wstarent
    real(r8), intent(in) :: minpblh(pcols)       ! minimum pbl height based on surface stress
    real(r8), intent(in) :: kvh_in(pcols,pver+1) ! kvh saved from last timestep
    real(r8), intent(in) :: kvm_in(pcols,pver+1) ! kvm saved from last timestep

    ! Outputs

    real(r8), intent(out) :: kvh(pcols,pver+1)   ! diffusivity for heat and tracers
    real(r8), intent(out) :: kvm(pcols,pver+1)   ! diffusivity for momentum

    real(r8), intent(out) :: pblh(pcols)         ! planetary boundary layer height
    real(r8), intent(out) :: tpert(pcols)        ! convective temperature excess
    real(r8), intent(out) :: qpert(pcols)        ! convective humidity excess
    real(r8), intent(out) :: tke(pcols,pver+1)   ! Turbulent kinetic energy (m2s-2)
    real(r8), intent(out) :: bprod(pcols,pver+1) ! Buoyancy production
    real(r8), intent(out) :: sprod(pcols,pver+1) ! shear production
    real(r8), intent(out) :: wstar        ! convective velocity for CL.

    ! Internal variables

    logical :: belongcv(pcols,pver+1) ! True for interfaces in a conv. layer (CL)
    logical :: belongst(pcols,pver+1) ! True for stable turbulent interfaces
    logical :: in_CL               ! True if interfaces k,k+1 both in same CL.
    integer :: i                   ! longitude index
    integer :: k                   ! vertical index
    integer :: ks                  ! vertical index
    integer :: ncvfin(pcols)       ! Total number of CL in column
    integer :: ncvf                ! Total number of CL in column prior to
    !  addition of single-layer rad-driven CLs 
    integer :: ncv                 ! index of current CL
    integer :: ncvnew              ! index of added single-layer rad-driven CL

    integer :: ncvsurf             ! if nonzero, index of CL including surface
    integer :: kbase(pcols,ncvmax) ! vertical index for base interface of CL
    integer :: ktop(pcols,ncvmax)  ! vertical index for top interface of CL
    integer :: kb, kt              ! kbase and ktop for current CL
    integer :: ktblw               ! ktop of below CL regime
    integer :: turbtype(pcols,pver+1) ! Interface turb. type (0 = no turb, 1 = 
    ! stable turb, 2 = CL interior, 3 = bottom 
    ! entr. intfc, 4 = top entr. intfc. 


    real(r8) :: bflxs(pcols)       ! Surface buoyancy flux (m2/s3)
    real(r8) :: tkes(pcols)        ! Surface TKE
    real(r8) :: wpert(pcols)       ! convective temperature excess
    real(r8) :: ebrk(pcols,ncvmax)
    real(r8) :: wbrk(pcols,ncvmax)
    real(r8) :: lbrk(pcols,ncvmax)
    real(r8) :: ghcl(pcols,ncvmax)
    real(r8) :: shcl(pcols,ncvmax)
    real(r8) :: smcl(pcols,ncvmax)
    real(r8) :: leng(pcols,pver+1) ! turbulent length scale
    real(r8) :: wcap(pcols,pver+1) ! W (m2/s2)
    real(r8) :: rcap(pcols,pver+1) ! e/<e>
    integer  :: ktopbl(pcols)      ! index of first midpt inside pbl
    real(r8) :: jtzm     ! Interface layer thickness atop conv layer (CL) ncv
    real(r8) :: jtsl     ! Jump in s_l               atop       "
    real(r8) :: jtqt     ! Jump in q_t               atop       "
    real(r8) :: jtbu     ! Jump in buoyancy          atop       "
    real(r8) :: jtu      ! Jump in zonal wind        atop       "
    real(r8) :: jtv      ! Jump in meridional wind   atop       "
    real(r8) :: jt2slv   ! 2-layer Jump in s_lv              atop       "
    real(r8) :: radf     ! buoy flx jump at cloudtop from longwave rad flx div
    real(r8) :: jbzm     ! Interface layer thickness at base of CL ncv
    real(r8) :: jbsl     ! Jump in s_l               at base    "
    real(r8) :: jbqt     ! Jump in qt                at base    "
    real(r8) :: jbbu     ! Jump in buoyancy          at base    "
    real(r8) :: jbu      ! Jump in zonal wind        at base    "
    real(r8) :: jbv      ! Jump in merid. wind       at base    "
    real(r8) :: ch       ! buoy flux coefs for sl, qt in a half-layer 
    real(r8) :: cm       ! 
    real(r8) :: n2ht     ! Moist squared buoy freq for top half-layer (s-2)
    real(r8) :: n2hb     ! Moist squared buoy freq for bottom half-layer (s-2)
    real(r8) :: ckh      ! Galperin stability function for heat
    real(r8) :: ckm      ! Galperin stability function for momentum
    real(r8) :: gh       ! Normalised buoyancy production
    real(r8) :: lbulk    ! Depth of turbulent layer
    real(r8) :: dzht     ! Thickness of top half-layer
    real(r8) :: dzhb     ! Thickness of bot half-layer
    real(r8) :: rootp    ! sqrt(CL-mean TKE) [m/s]     
    real(r8) :: evhc     ! (1+E) with E = evap. cool. efficiency [nd]
    real(r8) :: kentr    ! effective entrainment diffusivity we*dz [m2/s]
    real(r8) :: lwp      ! liquid water path in layer kt
    real(r8) :: opt_depth   ! optical depth of layer kt [nd]
    real(r8) :: radinvfrac  ! frac of lw cooling in layer kt put at inv.[nd]
    real(r8) :: wet      ! CL top entrainment rate
    real(r8) :: web      ! CL bot entrainment rate (for elevated CL)
    real(r8) :: vyt      ! n2ht/n2 at top inversion [nd]
    real(r8) :: vyb      ! n2hb/n2 at bot inversion [nd]
    real(r8) :: vut      ! inverse Ri at top inversion [nd], s2/n2
    real(r8) :: vub      ! inverse Ri at bot inversion [nd], s2/n2
    real(r8) :: fact     ! Factor relating TKE generation to entrainment [nd]
    real(r8) :: trma     ! intermediate variables
    real(r8) :: trmb
    real(r8) :: trmc
    real(r8) :: trmp
    real(r8) :: trmq
    real(r8) :: qq
    real(r8) :: det

    ! wstar-entrainment-specific variables

    real(r8) :: cet          ! proportionality coeff between wet and wstar^3
    real(r8) :: ceb          ! proportionality coeff between web and wstar^3
!    real(r8) :: wstar        ! convective velocity for CL.
    real(r8) :: wstar3       ! cubed convective velocity for CL.
    real(r8) :: wstar3fact   ! 1/(relative change in wstar^3 from entrainment)
    real(r8) :: rmin         ! sqrt(p)
    real(r8) :: fmin         ! f(rmin), where f(r) = r^3 - 3*p*r - 2q
    real(r8) :: rcrit        ! ccrit*wstar
    real(r8) :: fcrit        ! f(rcrit)
    logical noroot           ! true if f(r) has no root r > rcrit

    ! kvh and kvm now stored over timesteps in vertical_diffusion.F90
    ! and passed in as kvh_in, kvm_in
    ! However, on the first timestep they need to be computed and are done
    ! just before calling 'caleddy'

    do k = 1,pver+1
       do i = 1, ncol
          ! Initialize kvh and kvm to zero or kvf
          if ( use_kvf ) then
            kvh(i,k) = kvf(i,k)
            kvm(i,k) = kvf(i,k)
          else
            kvh(i,k) = 0._r8
            kvm(i,k) = 0._r8
          endif
          ! Zero diagnostic quantities for the new diffusion step.
            wcap(i,k) = 0.
            leng(i,k) = 0.
            rcap(i,k) = 0.
            tke(i,k) = 0.
            turbtype(i,k) = 0
       end do
    end do

    do k = 2,pver
       do i = 1, ncol
          ! bprod and sprod are diagnosed from eddy diffusivities from end of previous
          ! diffusion step. These will later be modified at the top and bottom of all
          ! convective layers.
            bprod(i,k) = -kvh_in(i,k)*n2(i,k)
            sprod(i,k) =  kvm_in(i,k)*s2(i,k)
       end do
    end do

    ! Set values at top and bottom of model
    do i = 1, ncol
       bprod(i,1) = 0.0 ! top of model
       sprod(i,1) = 0.0 ! top of model
       ! Calculate 'surface' (actually lowest half-layer) buoyancy flux.
       ch = chu(i,pver+1)*(1-sflh(i,pver)) + chs(i,pver+1)*sflh(i,pver)   
       cm = cmu(i,pver+1)*(1-sflh(i,pver)) + cms(i,pver+1)*sflh(i,pver)   
       bflxs(i) = ch*shflx(i)*rrho(i) + cm*qflx(i)*rrho(i)
       bprod(i,pver+1) = bflxs(i)
       sprod(i,pver+1) = 0.
    end do

    ! Examine each column and flag any layers with ri < 0.
    call exacol(pcols, pver, ncol, ri, bflxs, minpblh, zi, ktop, kbase, ncvfin)

    do i = 1, ncol
       ! Diagnostic surface TKE used in calculating layer-average TKE.

       ! This commented-out line can make tkes small or even negative. The replacement
       ! has at least equal justifiability, trying to represent tke at the surface
       ! rather than at zm.  Near the surface is where averaging W to get <e> is
       ! not really justifiable.
       !       tkes(i)=(b1*(ustar(i)**3+vk*z(i,pver)*bflxs(i)))**2./3.
       tkes(i)= b123*ustar(i)**2
       tkes(i)= min(tkes(i),tkemax)
       ncvsurf = 0  ! If nonzero, will be set to index of CL containing surface

       ! --------Convective layer (CL) computations----------------------------------

       ! If some convective layers have been found, determine their bulk 
       ! properties (<e>, <Sh>, <Sm>, and indices ktop/bkase for upper/lower 
       ! inversion )

       if (ncvfin(i) .gt. 0) then
          call zisocl(pcols, pver, i, z, zi, n2, s2, bflxs, tkes, &
               kbase, ktop, ncvfin, ebrk, wbrk, ghcl, shcl, smcl, lbrk, belongcv)
          !
          ! CLs found by zisocl are in order of height, so if any CL contains 
          ! the surface, it will be CL1.
          !
          if (kbase(i,1) .eq. pver+1) ncvsurf = 1
       else
          belongcv(i,:) = .false.
       end if

       ! Find single-level radiatively-driven cloud-topped convective layers (SRCLs).
       ! SRCLs extend through a single thermo layer k, with entrainment at interfaces
       ! k and k+1 (unless k+1 is the surface, in which case surface shear generation
       ! contributes to the layer-averaged energy). The conditions for an SRCL are:
       !   1. cloud at level k
       !   2. no cloud at level k+1 (else assuming that some fraction of the longwave
       !      flux div in layer k is concentrated at the top interface is invalid.
       !   3. Longwave radiative cooling (shortwave heating is assumed uniformly
       !      distributed through layer k, so not relevant to buoyancy production of
       !      TKE
       !   4. Internal stratification n2ht of half-layer from level k to interface k
       !      is unstable using similar method as in sfdiag, but applied to internal
       !      slopes of sl, qt in layer k.
       !   5. Interfaces k, k+1 not both in the same existing convective layer.
       !   6. k >= ntop_turb + 1 
       !   7. Ri(k+.5) > ricrit (otherwise stable turb mixing will broadly distribute
       !      the cloud top in the vertical, preventing localized radiative 
       !      radiative destabilization at the interface height.)

       ncv = 1
       ncvf = ncvfin(i)
       do k = nbot_turb,ntop_turb+1,-1
          if (ql(i,k).gt.qmin .and. ql(i,k-1).lt.qmin .and. qrl(i,k).lt.0. &
                              .and. ri(i,k).gt.ricrit) then
             ch = (1 -sfuh(i,k))*chu(i,k) + sfuh(i,k)*chs(i,k)
             cm = (1 -sfuh(i,k))*cmu(i,k) + sfuh(i,k)*cms(i,k)
             n2ht = ch*slslope(i,k) + cm*qtslope(i,k)
             if (n2ht.le.0.) then

                ! Test if k and k+1 are part of the same preexisting CL. If
                ! not, find appropriate index for new SRCL. Note that this
                ! calculation makes use of ncv set from prior passes through
                ! the do loop

                in_CL = .false.
                do while (ncv .le. ncvf)
                   if (ktop(i,ncv) .le. k) then

                      !  If kbase > k, k and k+1 are part of same prior CL

                      if (kbase(i,ncv) .gt. k) then 
                         in_CL = .true.
                      end if
                      exit  !exit from do-loop once CL top at/above intfc k.
                   else
                      ncv = ncv + 1  ! Go up one CL
                   end if
                end do ! ncv

                if (.not.in_CL) then

                   !  Add a new SRCL

                   ncvfin(i) = ncvfin(i)+1
                   ncvnew = ncvfin(i)
                   ktop(i,ncvnew) = k
                   kbase(i,ncvnew) = k+1
                   belongcv(i,k) = .true.
                   belongcv(i,k+1) = .true.

                   !                   if (k.lt.pver) then
                   !cb-should really use line below, to not add this energy in for
                   !   a surface stable layer.
                   if (k.lt.pver) then
                      ebrk(i,ncvnew) = 0.
                      lbrk(i,ncvnew) = 0.
                      shcl(i,ncvnew) = 0.
                      smcl(i,ncvnew) = 0.
                   else ! surface radiatively driven fog
                      if (bflxs(i).gt.0.) then !incorporate surface TKE into CL
                         ebrk(i,ncvnew) = tkes(i)
                         lbrk(i,ncvnew) = z(i,k)
                      else   !stable surf,don't incorporate surface TKE into CL
                         ebrk(i,ncvnew) = 0.
                         lbrk(i,ncvnew) = 0.
                      end if
                      shcl(i,ncvnew) = 0.
                      smcl(i,ncvnew) = 0.
                      ncvsurf = ncvnew
                   end if
                end if
             end if
          end if
       end do !k

       ! ---------------finished with adding SRCLs------------------------

       ! For each CL, compute length scale, r^2, e, Kh and Km

       do ncv = 1, ncvfin(i)
          kt = ktop(i,ncv)
          kb = kbase(i,ncv)
          lbulk = zi(i,kt)-zi(i,kb)
          do k = min(kb,pver), kt, -1 
             leng(i,k) = vk*zi(i,k) / (1.+vk*zi(i,k)/(tunl*lbulk))
             wcap(i,k) = leng(i,k)**2 * &
                         (-shcl(i,ncv)*n2(i,k)+smcl(i,ncv)*s2(i,k))
          end do ! k

          ! Calculate jumps at the lower inversion 

          if (kb.lt.pver+1) then 
             jbzm = z(i,kb-1)-z(i,kb)
             jbsl = sl(i,kb-1)-sl(i,kb)
             jbqt = qt(i,kb-1)-qt(i,kb)
             jbbu = n2(i,kb)*jbzm
             jbbu = max(jbbu,jbumin)
             jbu = u(i,kb-1)-u(i,kb)
             jbv = v(i,kb-1)-v(i,kb)
             ch = (1 -sflh(i,kb-1))*chu(i,kb-1) + sflh(i,kb-1)*chs(i,kb-1)
             cm = (1 -sflh(i,kb-1))*cmu(i,kb-1) + sflh(i,kb-1)*cms(i,kb-1)
             n2hb= (ch*jbsl + cm*jbqt)/jbzm
             vyb = n2hb*jbzm/jbbu
             vub = min(1.,(jbu**2+jbv**2)/(jbbu*jbzm) )
          else  !Zero bottom entr. contribution for CL extending down to sfc
             vyb = 0.
             vub = 0.
             web = 0.
          end if

          ! Calculate jumps at the upper inversion

          jtzm = z(i,kt-1)-z(i,kt)
          jtsl = sl(i,kt-1)-sl(i,kt)
          jtqt = qt(i,kt-1)-qt(i,kt)
          jtbu = n2(i,kt)*jtzm !Note - guaranteed positive by defn of entr.int.
          jtbu = max(jtbu,jbumin)  !...but threshold it anyway to be sure
          jtu = u(i,kt-1)-u(i,kt)
          jtv = v(i,kt-1)-v(i,kt)
          ch = (1 -sfuh(i,kt))*chu(i,kt) + sfuh(i,kt)*chs(i,kt)
          cm = (1 -sfuh(i,kt))*cmu(i,kt) + sfuh(i,kt)*cms(i,kt)
          n2ht = (ch*jtsl + cm*jtqt)/jtzm
          vyt = n2ht*jtzm/jtbu ! Ratio of buoy flux to w'(b_l)'
          vut = min(1.,(jtu**2+jtv**2)/(jtbu*jtzm)) ! Inv shear prd/buoy prd
          !write          write(6,*)'CALEDDY: kt, jtsl(K), jtqt, jtzm, n2,n2ht,vyt,vut '
          !write  &                 kt, jtsl,jtqt,jtzm,n2(i,k),n2ht, vyt,vut

          ! Calculate evaporative entrainment enhancement factor evhc. 
          ! We take the full inversion strength to be 
          ! jt2slv = slv(i,kt-2)  - slv(i,kt), (kt-1 is in the ambiguous layer).  
          ! However, for a cloud-topped CL overlain
          ! by another convective layer, it is possible that slv(i,kt-2) < slv(i,kt). 
          ! To avoid negative or excessive evhc, we lower-bound jt2slv and upper-bound
          ! evhc.

          evhc = 1.
          if (ql(i,kt).gt.qmin .and. ql(i,kt-1).lt.qmin) then 
             jt2slv = slv(i,max(kt-2,1)) - slv(i,kt)
             jt2slv = max(jt2slv, jbumin*slv(i,kt-1)/ggr)
             evhc = 1.+a2l*a3l*lcond*ql(i,kt) / jt2slv
             evhc = min(evhc,evhcmax)
          end if

          ! Radiative forcing at the upper inversion if at a cloud top

          ! radf [m2/s3] is the additional buoyancy flux at the inversion base
          ! associated with cloud-top longwave cooling being mainly at the 
          ! cloud layer top rather than being uniformly distributed through
          ! the layer.

          if (ql(i,kt).gt.qmin .and. ql(i,kt-1).lt.qmin) then 
             lwp = ql(i,kt)*(pi(i,kt+1) - pi(i,kt))/ggr
             opt_depth = 156*lwp     ! est. longwave opt depth in layer kt

             ! Approximation to LW cooling frac at inversion
             ! (polynomial approx to exact formula 1 - 2/opt_depth + 2/(exp(opt_depth)-1))

             radinvfrac  = opt_depth*(4.+opt_depth) / &
                  (6.*(4.+opt_depth) + opt_depth**2)
             radf = qrl(i,kt)/(pi(i,kt)-pi(i,kt+1)) ! Cp*radiative COOLING [W/kg] 
             radf = max(radinvfrac*radf*(zi(i,kt)-zi(i,kt+1)),0.) * chs(i,kt)
             ! Can disable cloud LW cooling contribution to turbulence by uncommenting:
             !            radf = 0.
          else
             radf = 0.
          end if

          ! Calculate 1. interior, 2. cloud-top radiative, 
          ! 3. surface buoyancy flux contribution to 
          ! cubed convective velocity wstar3 based 
          ! on incoming buoyancy and shear production.

          dzht = zi(i,kt) - z(i,kt)     ! Thickness of top half-layer
          dzhb = z(i,kb-1) - zi(i,kb)   ! Thickness of bot half-layer
          wstar3 = radf*dzht
          do k = kt+1,kb-1
             wstar3 =  wstar3 + bprod(i,k)*(z(i,k-1) - z(i,k))
          end do
          if(kb.eq.pver+1) wstar3 = wstar3 + bflxs(i)*dzhb
          wstar3 = max(2.5*wstar3,0.)

          ! Now diagnose top and bottom entrainment rates (and the contribution of
          ! top/bottom entrainment to wstar3) using entrainment closures of the form
          !
          !      wet = cet*wstar3, web = ceb*wstar3
          !
          ! where cet and ceb depend on the entrainment interface jumps, ql, etc.
          ! Substituting these into (#), we can solve for wstar3 and then wet,web.
          ! No entrainment is diagnosed unless the interior wstar3 > 0

          if (wstar3 .gt. 0.) then
             cet = a1i * evhc/(jtbu*lbulk)
             if (kb.eq.pver+1) then    ! surface based CL
                wstar3fact = max(1. + 2.5*cet*n2ht*jtzm*dzht, &
                     wstar3factcrit)
             else                     ! elevated CL with basal entrainment
                ceb = a1i/(jbbu*lbulk)
                wstar3fact = max(1. + 2.5*cet*n2ht*jtzm*dzht &
                     + 2.5*ceb*n2hb*jbzm*dzhb,&
                     wstar3factcrit)
             end if   ! kb
             wstar3 = wstar3/wstar3fact       
          else ! wstar3 == 0
             cet = 0.
             ceb = 0.
          end if   ! wstar3

          ! Solve cubic equation (canonical form for analytic solution)
          !   r^3 - 3*trmp*r - 2*trmq = 0,   r = sqrt<e>
          ! to estimate <e> for CL, derived from layer-mean TKE balance:
          !
          !   <e>^(3/2)/(b_1*<l>) \approx <B + S>   (*)
          !   <B+S> = (<B+S>_int * l_int + <B+S>_et * dzt + <B+S>_eb * dzb)/lbulk
          !   <B+S>_int = <e>^(1/2)/(b_1*<l>)*<e>_int
          !   <B+S>_et  = (-vyt+vut)*wet*jtbu + radf 
          !   <B+S>_eb  = (-vyb+vub)*web*jbbu
          !
          ! where:
          !   <> denotes a vertical avg (over the whole CL unless indicated)
          !   l_int (called lbrk below) is aggregate thickness of interior CL flux layers
          !   dzt = zi(i,kt)-z(i,kt)   is thickness of top entrainment layer
          !   dzb = z(i,kb-1)-zi(i,kb) is thickness of bot entrainment layer
          !   <e>_int (called ebrk below) is the CL-mean TKE if only interior
          !                             interfaces contributed.
          !   wet, web                  are top. bottom entrainment rates
          ! For a single-level radiatively-driven convective layer, there are no 
          ! interior interfaces so ebrk = lbrk = 0. If the CL goes to the 
          ! surface, vyb and vub are zero and ebrk and lbrk have already incorporated 
          ! the surface interfacial layer, so the same formulas still apply.  
          !
          ! In the original formulation 
          !   wet*jtbu = a1l*evhc*<e>^3/2/leng(i,kt)
          !   web*jbbu = a1l*<e>^3/2/leng(i,kt)
          !
          ! In the wstar formulation
          !   wet*jtbu = a1i*evhc*wstar3/lbulk
          !   web*jbbu = a1i*wstar3/lbulk, 

          fact = (evhc*(-vyt+vut)*dzht + (-vyb+vub)*dzhb*leng(i,kb)/leng(i,kt))/lbulk

          if (wstarent) then

             !-- wstar entrainment formulation-------------------------------------
             ! Here trmq can have either sign, and will usually be nonzero even for non-
             ! cloud topped CLs.  If trmq > 0, there will be two positive roots r; we take 
             ! the larger one. If necessary, we limit entrainment and wstar to prevent a
             ! solution with r < ccrit*wstar, where we take ccrit = 0.5.
             !
             trma = 1.          
             trmp = ebrk(i,ncv) *(lbrk(i,ncv)/lbulk)/3. + ntzero
             trmq = 0.5*b1*(leng(i,kt)/lbulk)*(radf*dzht + a1i*fact*wstar3)

             ! Check if there is an acceptable root with r > rcrit = ccrit*wstar. 
             ! To do this, first find local minimum fmin of the cubic f(r) at sqrt(p), 
             ! and value fcrit = f(rcrit).

             rmin = sqrt(trmp)
             fmin = rmin*(rmin*rmin - 3.*trmp) - 2.*trmq
             wstar = wstar3**onet
             rcrit = ccrit*wstar
             fcrit = rcrit*(rcrit*rcrit - 3.*trmp) -2.*trmq

             ! No acceptable root exists (noroot = .true.) if either:
             !    1) rmin < rcrit (in which case cubic is monotone increasing for r > rcrit)
             !       and f(rcrit) > 0.
             ! or 2) rmin > rcrit (in which case min of f(r) in r > rcrit is at rmin)
             !       and f(rmin) > 0.  
             ! In this case, we reduce entrainment and wstar3 such that r/wstar = ccrit;
             ! this changes the coefficients of the cubic.

             noroot =      ((rmin.lt.rcrit).and.(fcrit.gt.0.)) &
                  .or. ((rmin.ge.rcrit).and.(fmin.gt.0.))
             if (noroot) then ! solve cubic for r
                trma = 1. - b1*(leng(i,kt)/lbulk)*a1i*fact/ccrit**3
                trma = max(trma,0.5)  ! limit entrainment enhancement of ebrk
                trmp = trmp/trma 
                trmq = 0.5*b1*(leng(i,kt)/lbulk)*radf*dzht/trma
             end if   ! noroot

             ! Solve the cubic

             qq = trmq**2 - trmp**3
             if (qq .ge. 0.) then 
                rootp =   (trmq+sqrt(qq))**(1./3.) + (max(trmq-sqrt(qq),0.))**(1./3.)
             else
                rootp = 2.*sqrt(trmp)*cos( acos(trmq/sqrt(trmp**3)) /3.)
             end if
 
             if (noroot)  wstar3 = (rootp/ccrit)**3   ! Adjust wstar if necessary
             wet = cet*wstar3                         ! Find entrainment rates
             if (kb.lt.pver+1) web = ceb*wstar3

          else !
             ! wstarentr = .false.  Use original entrainment formulation-------
             ! trmp > 0 if there are interior interfaces in CL, trmp = 0 otherwise.
             ! trmq > 0 if there is cloudtop radiative cooling, trmq = 0 otherwise.

             trma = 1. - b1*a1l*fact
             trma = max(trma,0.5)  ! Prevents runaway entrainment instab.
             trmp = ebrk(i,ncv) *(lbrk(i,ncv)/lbulk)/(3*trma)
             trmq = 0.5*b1*(leng(i,kt)/lbulk)*radf*dzht/trma

             qq = trmq**2 - trmp**3
             if (qq .ge. 0.) then 
                rootp = (trmq+sqrt(qq))**(1./3.) + (max(trmq-sqrt(qq),0.))**(1./3.)
             else ! also part of case 3
                rootp = 2.*sqrt(trmp)*cos( acos(trmq/sqrt(trmp**3)) /3.)
             end if   ! qq

             ! Find entrainment rates and limit them by free-entrainment values a1l*sqrt(e)

             wet = a1l*rootp*min(evhc*rootp**2/(leng(i,kt)*jtbu),1.)   
             if (kb.lt.pver+1)   web = a1l*rootp*min(evhc*rootp**2/(leng(i,kb)*jbbu),1.)
          end if ! wstarentr

          ! now go back to common code

          ebrk(i,ncv) = rootp**2
          ebrk(i,ncv) = min(ebrk(i,ncv),tkemax) ! limit CL-avg TKE used for entrainment
          wbrk(i,ncv) = ebrk(i,ncv)/b1

          ! The only way ebrk = 0 is for single-level radiatively-driven convective
          ! layers which are actually radiatively heated at top interface. In this case,
          ! we remove 'convective' label from the interfaces around this layer.
          ! This case should now be impossible, so we flag it
          if (ebrk(i,ncv) .le. 0.) then
             write(6,*)'CALEDDY: Warning, CL with zero TKE, i,kt,kb ',i,kt,kb
             belongcv(i,kt) = .false.
             belongcv(i,kb) = .false. 
          end if

          rcap(i,kt) = 1.   ! We approximate TKE = <e> at entr. interfaces
          rcap(i,kb) = 1.   ! consistent with entrainment closure.

          ! Calculate ratio rcap = e/<e> in convective layer interior. Bound it by
          ! limits rcapmin = 0.1 to rcapmax = 2.0 to take care of some pathological cases.

          do k = kb-1, kt+1, -1
             rcap(i,k) = max( (mu*leng(i,k)/lbulk + wcap(i,k)/wbrk(i,ncv)) / &
                  (mu*leng(i,k)/lbulk + 1.), rcapmin)
             rcap(i,k) = min(rcap(i,k), rcapmax)
          end do

          ! Compute TKE throughout CL, and cap by tkemax.

          do k = kb,kt,-1
             tke(i,k) = ebrk(i,ncv) * rcap(i,k)
             tke(i,k) = min(tke(i,k),tkemax)
          end do

          ! Compute CL interior diffusivities, buoyancy and shear production

          do k = kb-1, kt+1, -1
             kvh(i,k) = leng(i,k)*sqrt(tke(i,k))*shcl(i,ncv)
             kvm(i,k) = leng(i,k)*sqrt(tke(i,k))*smcl(i,ncv)
             bprod(i,k) = -kvh(i,k)*n2(i,k)
             sprod(i,k) =  kvm(i,k)*s2(i,k)
          end do

          ! Compute diffusivity wet*dz and some diagnostics at the upper inversion

             kentr = wet*jtzm
             kvh(i,kt) = kentr
             kvm(i,kt) = kentr
             bprod(i,kt) = -kentr*n2ht+radf
             sprod(i,kt) =  kentr*s2(i,kt)
             turbtype(i,kt) = 4

          ! Compute Kh, Km and some diagnostics at the lower inversion

          if (kb.lt.pver+1) then 
             kentr = web*jbzm
             if(kb.ne.ktblw) then
               kvh(i,kb) = kentr
               kvm(i,kb) = kentr
               bprod(i,kb) = -kvh(i,kb)*n2hb
               sprod(i,kb) =  kvm(i,kb)*s2(i,kb)
             else
             ! Handle case of two CLs entraining into each other by adding entr. diffusivities
             ! This works because we always go from surface upward.
             ! Question : Is it OK to directly sum bprod and sprod as below ?
             !            Should I do thickness weighting average (dzht, dzbt) instead ?  
               kvh(i,kb) = kvh(i,kb) + kentr 
               kvm(i,kb) = kvm(i,kb) + kentr
               bprod(i,kb) = bprod(i,kb)-kentr*n2hb
               sprod(i,kb) = sprod(i,kb)+kentr*s2(i,kb)
             end if
             turbtype(i,kb) = max(turbtype(i,kb),3)
          end if

          !cb Print out the W's at the different heights in the CL. Put min threshold 
          ! on TKE to prevent possible division by zero.

          wcap(i,kt) = (bprod(i,kt)+sprod(i,kt))*leng(i,kt)/sqrt(max(tke(i,kt),1e-6))
          if(kb.lt.pver+1) then
             wcap(i,kb) = (bprod(i,kb)+sprod(i,kb))*leng(i,kb)/sqrt(max(tke(i,kb),1e-6))
          else
             wcap(i,kb) = tkes(i)/b1
          end if

          ! Save the index of upper external interface of current CL-regime in order to
          ! handle the case when this interface is also the lower external interface of 
          ! CL-regime located above. 

          ktblw = kt 

! Writing entrainment rate at Sc-topped inversion layer top                                                                                      
!         ! if (ql(i,kt).gt.qmin .and. ql(i,kt-1).lt.qmin) then
!            write(6,500) ncvfin(i),ncv,zi(i,kt),wet*1000,wstar3**(1./3.),evhc
! 500        format ('ncvfin=',i3,3x,'ncv=',i3,3x,'zi=',f8.2,3x,'wet=',f8.2,3x, &
!                    'wstar=',f7.2,3x,'evhc=',f7.2,3x)
!         ! end if


       end do        ! ncv


       ! If the lowest CL reaches the surface, define the PBL depth as the CL top
       ! Also define PBL char horizontal turbulent perturbations tpert, qpert
       ! (These estimates could use some more refinement).

       if (ncvsurf .gt. 0) then
          ktopbl(i) = ktop(i,ncvsurf)
          pblh(i) = zi(i, ktopbl(i))
          wpert(i)  = max(wfac*sqrt(ebrk(i,ncvsurf)),wpertmin)
          tpert(i)  = max(abs(shflx(i)*rrho(i)/cp)*tfac/wpert(i),0._8)
          qpert(i)  = max(abs(qflx(i)*rrho(i))*tfac/wpert(i),0._8)
          if(bflxs(i).ge.0.) then
             turbtype(i,pver+1) = 2
          else
             turbtype(i,pver+1) = 3
          end if
       end if                 ! End of convective layer calculations

       ! ...Now specify diffusivities in stable turbulent layers......................

       ! Start by finding turbulent lengthscales in all such layers.

       belongst(i,1) = .false.   ! k = 1 assumed nonturbulent
       do k = 2, pver
          belongst(i,k)=(ri(i,k).lt.ricrit).and.(.not.belongcv(i,k))
          if (belongst(i,k).and.(.not.belongst(i,k-1))) then
             kt = k     ! Top of stable turb layer
          else if (.not.belongst(i,k) .and. belongst(i,k-1)) then
             kb = k-1   ! Base of stable turb layer
             lbulk = z(i,kt-1) - z(i,kb)
             do ks = kt,kb
                leng(i,ks)=vk*zi(i,ks) / (1.+vk*zi(i,ks)/(tunl*lbulk))
             end do
          end if
       end do ! k

       ! Now look whether stable turb layer extends to ground. Note that interface
       ! pver+1 is assumed to always be stable-turbulent if it is not convective. Note
       ! that if it is convective, pver will also be convective, so the above loop 
       ! will have finished finding all elevated turbulent layers.

       belongst(i,pver+1)=.not.belongcv(i,pver+1)
       if (belongst(i,pver+1)) then     ! kb = pver+1 (surface stable turb layer)
          turbtype(i,pver+1) = 1
          if (belongst(i,pver)) then    ! surface stable layer includes interior
             ! stable turbulent interface
             ! kt already defined above
             lbulk = z(i,kt-1)          ! (since z(i,kb) = 0)
          else                          ! surface stable BL with no interior turbulence
             kt = pver+1
          end if
          lbulk = z(i,kt-1)
          ktopbl(i) = kt-1
          pblh(i) = zi(i,ktopbl(i))     ! PBL Height <-> lowest stable turb. layer
          do ks = kt,pver
             leng(i,ks)=vk*zi(i,ks) / (1.+vk*zi(i,ks)/(tunl*lbulk))
          end do ! ks
          tpert(i)   = max(shflx(i)*rrho(i)/cp*fak/ustar(i),0._8) ! CCM stable-layer forms
          qpert(i)   = max(qflx(i)*rrho(i)*fak/ustar(i),0._8)
       end if

       do k = 2,pver
          if (belongst(i,k)) then    ! stable turbulent layer
             turbtype(i,k) = 1
             trma = alph3*alph4*ri(i,k) + 2.*b1*(alph2-alph4*alph5*ri(i,k))
             trmb = (alph3+alph4)*ri(i,k) + 2.*b1*(-alph5*ri(i,k)+alph1)
             trmc = ri(i,k)
             det = max(trmb*trmb-4.*trma*trmc,0.)
             gh = (-trmb + sqrt(det))/(2.*trma)
             gh = max(gh,-0.28)

             ! Now deduce stability fns, TKE, and diffusivities knowing Gh

             ckh = alph5 / (1.+alph3*gh)
             ckm = (alph1 + alph2*gh)/(1.+alph3*gh)/(1.+alph4*gh)
             if (ri(i,k) .le. ricrit) then 
                tke(i,k) = b1*leng(i,k)**2*(-ckh*n2(i,k)+ckm*s2(i,k))
                tke(i,k) = min(tke(i,k),tkemax)
                kvh(i,k) = leng(i,k) * sqrt(tke(i,k))*ckh
                kvm(i,k) = leng(i,k) * sqrt(tke(i,k))*ckm
                bprod(i,k) = -kvh(i,k)*n2(i,k)
                sprod(i,k) =  kvm(i,k)*s2(i,k)  
             else
                tke(i,k) = 0.
             end if
          end if
       end do  ! k
    end do   ! i

    return
  end subroutine caleddy
  !==============================================================================
  subroutine exacol(pcols, pver, ncol, ri, bflxs, minpblh, zi, ktop, kbase, ncvfin) 

    ! object : determine whether the column has adjacent regions where 
    !          Ri < 0 (unstable layers or ULs) and determine the indices 
    !          kbase, ktop which delimit these unstable layers : 
    !          ri(kbase) > 0 and ri(ktop) > 0, but ri(k) < 0 for ktop < k < kbase. 
    ! author : H. Grenier    05/2000, 
    !          C. Bretherton 08/2000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    ! input

    integer, intent(in) :: pcols              ! number of atmospheric columns   
    integer, intent(in) :: pver              ! number of atmospheric columns   
    integer, intent(in) :: ncol             ! number of atmospheric columns   

    real(r8), intent(in) :: ri(pcols,pver)   ! Moist gradient Richardson no.
    real(r8), intent(in) :: bflxs(pcols)     ! Surface buoyancy flux
    real(r8), intent(in) :: minpblh(pcols)   ! minimum pbl height based on surface stress
    real(r8), intent(in) :: zi(pcols,pver+1) ! interface heights

    ! output       

    integer, intent(out) :: kbase(pcols,ncvmax) ! vertical index of UL base
    integer, intent(out) :: ktop(pcols,ncvmax)  ! vertical index of UL top
    integer, intent(out) :: ncvfin(pcols)       ! number of ULs

    ! local variables

    integer :: i
    integer :: k

    integer :: ncv

    real(r8) :: rimaxentr
    real(r8) :: riex(pver+1)       ! Column Ri profile extended to surface
    ! by taking riex > rimaxentr for
    ! bflxs < 0, riex < rimaxentr for
    ! bflxs > 0.

    ! Initialize the outputs to zero. 

    do i = 1, ncol
       ncvfin(i) = 0
       do ncv = 1, ncvmax
          ktop(i,ncv) = 0
          kbase(i,ncv) = 0
       end do
    end do
    rimaxentr = 0.            ! Limiting Ri for entraining turb layer
    do i = 1, ncol
       riex(2:pver) = ri(i,2:pver)
       riex(pver+1) = rimaxentr-bflxs(i)    ! Allows consistent treatment
       ! of surface and other interfaces
       ncv = 0
       k = pver+1                           ! Work upward from surf interface
       do while ( k.gt.ntop_turb+1 )
          if (riex(k) .lt. rimaxentr .or. zi(i,k).lt.minpblh(i)) then 

             ! New CL

             ncv = ncv + 1

             ! First define kbase as interface below first unstable one

             kbase(i,ncv) = min(k+1,pver+1)

             ! Decrement k until top unstable level

             do while ((riex(k).lt.rimaxentr .or. zi(i,k).lt.minpblh(i)) &
                  .and. k.gt.ntop_turb+1 )
                k = k-1
             end do

             ! ktop is first interface above unstable layer

             ktop(i,ncv) = k
          else

             ! Search upward for a CL.

             k = k-1
          end if
       end do
       ncvfin(i) = ncv    
    end do  ! i
    return 
  end subroutine exacol
  !==============================================================================
  subroutine zisocl(pcols, pver, long, z, zi, n2, s2, bflxs, tkes, & 
       kbase, ktop, ncvfin, ebrk, wbrk, ghcl, shcl, smcl, &
       lbrk, belongcv)
    ! ------------------------------------------------------------------------
    ! object : find <e>, <W>, <Sh>, <Sm>, ktop(ncv), accounting for the 
    !          presence of stably stratified layers inside the convective layer
    !          (CL, to be defined) but with r2 > ratinv.
    !           Re-arrange the indexing of arrays kbase/ktop if some 
    !           convective layers are found to be coupled such that 
    !           ncv defines the index of each convective layer increasing
    !           with height.

    !           At the routine entrance, kb is the lowermost level of a convective 
    !           layer with N2 < 0.  At the exit of the routine, kb is the level 
    !           caarying the lower inversion (N2 >0) of the convective layer.
    ! author : H. Grenier 05/08/2000
    !------------------------------------------------------------------------------
    implicit none
    ! Input variables

    integer, intent(in) :: long              ! Longitude of the column
    integer, intent(in) :: pcols              ! number of atmospheric columns   
    integer, intent(in) :: pver              ! number of atmospheric columns   
    real(r8), intent(in) :: z(pcols, pver)  ! Layer midpoint height
    real(r8), intent(in) :: zi(pcols, pver+1) ! Interface height
    real(r8), intent(in) :: n2(pcols, pver)  ! Moist squared buoy freq (s-2)
    real(r8), intent(in) :: s2(pcols, pver)  ! Shear deformation (s-2)
    real(r8), intent(in) :: bflxs(pcols)     ! Surface buoyancy flux

    ! Input/output variables

    integer, intent(inout) :: kbase(pcols,ncvmax) ! Vert index of CL base
    integer, intent(inout) :: ktop(pcols,ncvmax)  ! Vert index of CL top
    integer, intent(inout) :: ncvfin(pcols)       ! Total number of CLs

    ! Output variables

    logical, intent(out) :: belongcv(pcols,pver+1) ! True if flux level in a CL
    real(r8), intent(in) :: tkes(pcols)           ! TKE at the surface
    real(r8), intent(out) :: ebrk(pcols, ncvmax)  ! Vert average of TKE over CL
    real(r8), intent(out) :: wbrk(pcols, ncvmax)  !   "          of W^2  "
    real(r8), intent(out) :: ghcl(pcols, ncvmax)  !   "          of Gh   "
    real(r8), intent(out) :: shcl(pcols, ncvmax)  !   "          of Sh   "
    real(r8), intent(out) :: smcl(pcols, ncvmax)  !   "          of Sm   "
    real(r8), intent(out) :: lbrk(pcols, ncvmax)  ! CL depth not within entr
    !  layers

    ! Local variables

    logical :: extend                ! True when CL is extended in zisocl
    logical :: bottom                ! True when CL base at surface(kb = pver+1)

    integer :: i                     ! Local index for the longitude
    integer :: ncv                   ! Index enumerating convective layers in col
    integer :: incv
    integer :: k
    integer :: kb                    ! Local index for kbase
    integer :: kt                    ! Local index for ktop
    integer :: ncvinit               ! Value of ncv at routine entrance 
    integer :: cntu                  ! counts upward no. of merged CLs
    integer :: cntd                  ! counts downward  "          "
    integer :: kbinc                 ! Index for incorporating underlying CL
    integer :: ktinc                 ! Index for incorporating  overlying CL

    real(r8) :: wint
    real(r8) :: dwinc
    real(r8) :: dzinc
    real(r8) :: dwsurf
    real(r8) :: gh
    real(r8) :: sh
    real(r8) :: sm
    real(r8) :: l2n2                 ! Vert. integral of l^2N^2 over CL
    real(r8) :: l2s2                 ! Vert. integral of l^2S^2 over CL
    real(r8) :: dl2n2                ! Vert. int. of l^2N^2 over incorp. layer
    real(r8) :: dl2s2                ! Vert. int. of l^2S^2 over incorp. layer
    real(r8) :: lint                 ! CL depth excluding entrainment layers
    real(r8) :: lbulk                ! Depth of the convective layer
    real(r8) :: lz                   ! Turbulent length scale
    real(r8) :: ricl                 ! Ri Number for the whole convective layer
    real(r8) :: trma
    real(r8) :: trmb
    real(r8) :: trmc
    real(r8) :: det
    real(r8) :: zbot                 ! Height of CL base
    real(r8) :: l2rat                ! Square of ratio of actual to initial CL
    !  depth.

    ! Define parameters

    i = long
    ncv = 1

    ! Loop over convective layers to see if any of them need to be extended

    !    if(i.eq.89) then
    !       write(6,*)'ZISOCL: incoming i, ncvfin(i) ',i,ncvfin(i)
    !       write(6,*)'ZISOCL: kbase(ncv) ',(kbase(i,ncv),ncv=1,ncvfin(i))
    !       write(6,*)'ZISOCL: ktop(ncv) ',(ktop(i,ncv),ncv=1,ncvfin(i))
    !    end if
    do while ( ncv .le. ncvfin(i) )
       ncvinit = ncv
       cntu = 0
       cntd = 0
       kb = kbase(i,ncv) 
       kt = ktop(i,ncv)
       !write       write(6,*)'ZISOCL: Incoming ncv kb kt', ncv,kb,kt
       lbulk = zi(i,kt)-zi(i,kb)

       ! Add contribution (if any) from surface interfacial layer to turbulent 
       ! production and lengthscales.  If there is positive surface buoy flx,
       ! the CL extends to the surface and there is a surface interfacial 
       ! layer contribution to W and the CL 
       ! interior depth. If there is negative buoyancy flux, the surface 
       ! interfacial layer is treated as energetically isolated from the
       ! rest of the CL and does not contribute to the layer-interior W or
       ! depth. This case also requires a redefinition of lbulk.

       bottom = kb .eq. pver+1
       if (bottom .and. (bflxs(i).ge. 0.)) then
          lint = z(i,pver)
          dwsurf = (tkes(i)/b1)*z(i,pver)
       else
          lint = 0.
          dwsurf = 0.
       end if

! Below are Previous Codes
!!       if (bottom .and. (bflxs(i).ge. 0.)) then
!       if(bottom)then
!          lint = z(i,pver)
!          dwsurf = (tkes(i)/b1)*z(i,pver)
!!       else if (bottom .and. (bflxs(i).lt. 0.)) then
!!          lint = 0.
!!          dwsurf = 0.
!!          lbulk = zi(i,kt)-z(i,pver)
!       else
!          lint = 0.
!          dwsurf = 0.
!       end if


       l2n2 = 0.
       l2s2 = 0.

       ! Turbulence contribution from conv layer (CL) interior kt < k < kb,
       ! which at this point contains only unstable stratified interfaces.
       ! Based on the CL interior stratification, initial guesses at the
       ! stability functions are made. If there is no CL interior interface,
       ! neutral stability is assumed for now.

       if (kt .lt. kb-1) then 
          do k = kb-1, kt+1, -1
             lz = vk*zi(i,k)/(1.+vk*zi(i,k)/(tunl*lbulk))
             l2n2 = l2n2 + lz*lz*n2(i,k)*(z(i,k-1)-z(i,k))
             l2s2 = l2s2 + lz*lz*s2(i,k)*(z(i,k-1)-z(i,k))
             lint = lint + (z(i,k-1)-z(i,k))
          enddo

          ! Solve for bulk Sh, Sm, and wint over the CL interior

          ricl = min(l2n2/l2s2,ricrit) !actually we should have ricl < 0 
          trma = alph3*alph4*ricl+2.*b1*(alph2-alph4*alph5*ricl)
          trmb = ricl*(alph3+alph4)+2.*b1*(-alph5*ricl+alph1)
          trmc = ricl
          det = max(trmb*trmb-4.*trma*trmc,0.)
          gh = (-trmb + sqrt(det))/2./trma
          gh = max(gh,-0.28)
          gh = min(gh,0.0233)
          sh = alph5 / (1.+alph3*gh)
          sm = (alph1 + alph2*gh)/(1.+alph3*gh)/(1.+alph4*gh)
          wint = -sh*l2n2 + sm*l2s2 + dwsurf 
       else

          ! There is no CL interior interface. The only way that should happen
          ! at this point is if there is upward surface buoy flux but no
          ! unstable interior interfaces. In that case, the surface interface
          ! turbulent production terms are used as its CL 'interior'.

          if (bottom) then
             !write             write(6,*) 'ZISOCL:single-level surf CL', kb, kt
             wint = dwsurf
             gh = 0.  ! use neutral stability fns for layer extension
             sh = alph5 / (1.+alph3*gh)
             sm = (alph1 + alph2*gh)/(1.+alph3*gh)/(1.+alph4*gh)
             ! follow line added by cw 10/21/05 for consistency with latest UW
             l2n2 = -wint/sh
          else
             write(6,*) 'ZISOCL: impossible case', kb, kt
             call task_abort()
          endif
       endif

       ! Try to extend the top and base of convective layer (indices kt, kb). 

       extend = .false.    ! will become .true. if CL top or base is extended

       ! Extend the CL top if possible.

       ! Compute possible contributions to TKE production and lengthscale,
       ! were the CL top interfacial layer found by exacol incorporated
       ! into the CL interior.

       dzinc = z(i,kt-1)-z(i,kt)
       lz = vk*zi(i,kt)/(1.+vk*zi(i,kt)/(tunl*lbulk))
       dl2n2 = lz*lz*n2(i,kt)*dzinc
       dl2s2 = lz*lz*s2(i,kt)*dzinc
       dwinc = -sh*dl2n2 + sm*dl2s2

       ! Test for incorporation of top layer kt into CL interior. If true,
       ! extend CL by incorporating top layers until test fails.

       !       do while (dwinc .gt. (rinc*dzinc*wint/(lint+(1-rinc)*dzinc)))
       l2n2 = -min(-l2n2,tkemax*lint/(b1*sh))
       do while (-dl2n2 .gt. (-rinc*dzinc*l2n2/(lint+(1-rinc)*dzinc)))

          ! Add contributions from layer kt to interior lengthscale/TKE prod

          lint = lint + dzinc
          wint = wint + dwinc
          l2n2 = l2n2 + dl2n2
          l2n2 = -min(-l2n2,tkemax*lint/(b1*sh))
          l2s2 = l2s2 + dl2s2

          ! Extend top of CL upward a layer.

          kt = kt-1
          extend = .true.
          if (kt .eq. ntop_turb) then
             write(6,*) 'ZISOCL: Abort: Tried to extend convective layer ', &
                  'to top model interface k = 1 at i = ',i
!             write(6,*) 'ZISOCL: kb, ebrk = ',kb,b1*wint/lint
             write(6,*) 'ZISOCL: kb, ebrk = ',kb,b1,wint,lint
             write(*,*) 'n2',n2,'s2',s2,'bflx',bflxs
             save3Dbin = .true.
             call write_fields3D()
             call task_abort()
          end if

          ! Check for existence of an overlying CL which might be merged.
          ! If such exists (ktinc > 1), check for merging by testing for 
          ! incorporation of its top interior interface into current CL.
          ! If no such layer exists, ktop(i,ncv+cntu+1) will equal its
          ! default value of zero, so ktinc will be 1 and the test kt=ktinc
          ! will fail.

          ktinc = ktop(i,ncv+cntu+1)+1
          if (kt .eq. ktinc) then
             ncvfin(i) = ncvfin(i) -1
             cntu = cntu + 1 
          end if

          ! Compute possible lengthscale and TKE production
          ! contributions were layer kt incorporated into CL interior.
          ! Then go back to top of loop to test for incorporation

          dzinc = z(i,kt-1)-z(i,kt)
          lz = vk*zi(i,kt)/(1.+vk*zi(i,kt)/(tunl*lbulk))
          dl2n2 = lz*lz*n2(i,kt)*dzinc
          dl2s2 = lz*lz*s2(i,kt)*dzinc
          dwinc = -sh*dl2n2 + sm*dl2s2
       end do   ! Done with top extension of CL
       !       write(6,*)'ZISOCL i,ncv,kb,kt,cntu,ktex,ncvfin', &
       !                 i,ncv,kb,ktop(i,ncv),cntu,kt,ncvfin(i)

       ! Shift indices appropriately if layers have been merged

       if (cntu .gt. 0) then
          do incv = 1, ncvfin(i)-ncv
             kbase(i,ncv+incv) = kbase(i,ncv+cntu+incv)
             ktop(i,ncv+incv) = ktop(i,ncv+cntu+incv)
          end do
       end if

       ! Extend the CL base if possible.

       if (.not. bottom) then

          ! Compute possible contributions to TKE production and lengthscale,
          ! were the CL base interfacial layer found by exacol incorporated
          ! into the CL interior.

          dzinc = z(i,kb-1)-z(i,kb)
          lz = vk*zi(i,kb)/(1.+vk*zi(i,kb)/(tunl*lbulk))
          dl2n2 = lz*lz*n2(i,kb)*dzinc
          dl2s2 = lz*lz*s2(i,kb)*dzinc
          dwinc = -sh*dl2n2 + sm*dl2s2

          ! Test for incorporation of base layer kb into CL interior. If true,
          ! extend CL by incorporating base layers until test fails.

          !          do while((dwinc.gt.rinc*dzinc*wint/(lint+(1-rinc)*dzinc)) &
          do while((-dl2n2.gt.-rinc*dzinc*l2n2/(lint+(1-rinc)*dzinc)) &
               .and. (.not. bottom) )

             !write             write(6,*) 'ZISOCL: Base incorporation of ', kb

             ! Add contributions from layer kb to interior lengthscale/TKE prod

             lint = lint + dzinc
             wint = wint + dwinc
             l2n2 = l2n2 + dl2n2
             l2n2 = -min(-l2n2,tkemax*lint/(b1*sh))
             l2s2 = l2s2 + dl2s2

             ! Extend base of CL downward a layer

             kb = kb+1
             extend = .true.

             ! Check for existence of an underlying CL which might be merged.
             ! If such exists (kbinc > 1), check for merging by testing for 
             ! incorporation of its top interior interface into current CL.
             ! Note that this top 'interior' interface could be the surface.

             kbinc = 0
             if (ncv .gt. 1) kbinc = ktop(i,ncv-1)+1
             if (kb .eq. kbinc) then

                ! We are incorporating interior of CL ncv-1, so merge
                ! this CL into the current CL.

                ncv = ncv - 1
                ncvfin(i) = ncvfin(i) -1
                cntd = cntd + 1 
             end if

             ! If CL would now reach the surface, check sign of surface
             ! buoyancy flux. If positive, add contributions of surface 
             ! interfacial layer to TKE production and lengthscale. If 
             ! negative, we regard the surface layer as stable and do not
             ! add surface interfacial layer contributions to the CL.
             ! In either case the surface interface is classified as
             ! part of the CL for book-keeping purposes (to ensure no base
             ! entrainment calculation is done). If we are merging with a 
             ! surface-driven CL with no interior unstable interfaces, the
             ! above code will already have handled the merging book-keeping.

             bottom = kb .eq. pver+1
             if (bottom) then 
                if (bflxs(i) .gt. 0.) then 
                   dwsurf = (tkes(i)/b1)*z(i,pver)
                   lint = lint + z(i,pver)
                ! following two lines added by cw 10/21/05 for consistency with UW
                else
                   dwsurf = 0
                end if
             else

                ! Compute possible lengthscale and TKE production
                ! contributions were layer kb incorporated into CL interior.
                ! Then go back to top of loop to test for incorporation

                dzinc = z(i,kb-1)-z(i,kb)
                lz = vk*zi(i,kb) / (1.+vk*zi(i,kb)/(tunl*lbulk))
                dl2n2 = lz*lz*n2(i,kb)*dzinc
                dl2s2 = lz*lz*s2(i,kb)*dzinc
                dwinc = -sh*dl2n2 + sm*dl2s2
             end if
          end do

          !          write(6,*)'ZISOCL kb,kt,cntd,kbex,ncvfin', &
          !                 kbase(i,ncvinit),ktop(i,ncvinit),cntd,kb,ncvfin(i)

          if (bottom .and. ncv .ne. 1) then 
             write(*,*) 'n2',n2
             write(6,*) 'Major mistake ZISOCL: bottom CL not indexed 1'
             save3Dbin = .true.
             call write_fields3D() !bloss add output of 3D fields for diagnostics
             call task_abort()
          end if

       end if   ! Done with bottom extension of CL 

       ! Shift indices if some layers with N2 < 0 have been found

       if (cntd .gt. 0) then
          do incv = 1, ncvfin(i)-ncv
             kbase(i,ncv+incv) = kbase(i,ncvinit+incv)
             ktop(i,ncv+incv) = ktop(i,ncvinit+incv)
          end do
       end if

       ! Sanity check for positive wint.
       if (wint .lt. 0.01) then
!          write(6,*)'ZISOCL: Stop - interior avg TKE < 0.01: wint = ',wint
          wint = 0.01
       end if

       if (extend) then

          ! Recompute base and top indices if necessary

          ktop(i,ncv) = kt
          kbase(i,ncv) = kb

          !! Recompute Ri_cl, Sh, Sm, and <W> after layer extension if necessary
          !! Ideally, we would recompute l2n2 and l2s2 to account for the incorrect
          !! lbulk used in the computation of lz, but we take the simpler approach
          !! of simply multiplying the lz's by the ratio of the actual PBL 
          !! depth to lbulk.

          !zbot = zi(i,kb)
          !if (bottom .and. (bflxs(i).lt.0)) zbot = z(i,pver)
          !l2rat = ((zi(i,kt) - zbot)/lbulk)**2
          !l2n2 = l2n2*l2rat
          !l2s2 = l2s2*l2rat

          ! Recompute Ri_cl, Sh, Sm, and <W> after layer extension if necessary.
          ! l2n2 and l2s2 are recalculated using my methods. Below is s replacement
          ! of the above code.

          lbulk=zi(i,kt)-zi(i,kb)
          l2n2=0.
          l2s2=0.
          do k=kt+1,kb-1
            dzinc=z(i,k-1)-z(i,k)
            lz=vk*zi(i,k)/(1.+vk*zi(i,k)/(tunl*lbulk))
            l2n2=l2n2+lz*lz*n2(i,k)*dzinc
            l2s2=l2s2+lz*lz*s2(i,k)*dzinc
          end do

          ricl = min(l2n2/l2s2,ricrit)
          trma = alph3*alph4*ricl+2.*b1*(alph2-alph4*alph5*ricl)
          trmb = ricl*(alph3+alph4)+2.*b1*(-alph5*ricl+alph1)
          trmc = ricl
          det = max(trmb*trmb-4.*trma*trmc,0.)
          gh = (-trmb + sqrt(det))/2./trma
          gh = max(gh,-0.28)
          gh = min(gh,0.0233)
          sh = alph5 / (1.+alph3*gh)
          sm = (alph1 + alph2*gh)/(1.+alph3*gh)/(1.+alph4*gh)

          ! It is conceivable that even though the original wint was positive, it will
          ! be negative after correction. In this case, correct wint to be a small
          ! positive number
          wint = max(dwsurf + (-sh*l2n2 + sm*l2s2),0.01*wint)
       end if

       lbrk(i,ncv) = lint
       wbrk(i,ncv) = wint/lint
       ebrk(i,ncv) = b1*wbrk(i,ncv)
       ebrk(i,ncv) = min(ebrk(i,ncv),tkemax)
       ghcl(i,ncv) = gh 
       shcl(i,ncv) = sh
       smcl(i,ncv) = sm

       ! Increment counter for next CL

       ncv = ncv + 1

    end do                   ! Loop over convective layers

    !write    write(6,*) 'ZISOCL: Total number of CL', ncvfin(i)

    do k = 1, pver+1
       belongcv(i,k) = .false.
    end do

    do ncv = 1, ncvfin(i)
       do k = ktop(i,ncv), kbase(i,ncv)
          belongcv(i,k) = .true.
       enddo
    enddo

    return
  end subroutine zisocl

! ----------------------------------------------------------------------------------------

subroutine diffuse_scalar1D_rev_full (field, bflux, tkh_z, dtn_local)

! input	
real field(:)	! scalar
real tkh_z(:)	! vertical eddy conductivity
real bflux		! bottom flux
real dtn_local         ! timestep
! local        
real flx(0:nzm)
real dfdt
real rdz2, rdz, rdz5
real tkz,rhoi
integer i,j,k,ib,ic,jb,jc,kc,kb, k_flip, kc_flip, kb_flip

rdz=1./dz
rdz2=1./(dz*dz)

flx(0)=bflux*rdz*rhow(1)
flx(nzm)=0.0

dfdt=0.
do k=1,nzm-1
   kc=k+1	
   rhoi = rhow(kc)/adzw(kc)
   rdz5=0.5*rdz2
   tkz=rdz5*(tkh_z(k)+tkh_z(kc))
   flx(k)=-tkz*(field(kc)-field(k))*rhoi
end do

do k=1,nzm
   kb=k-1
   rhoi = 1./(adz(k)*rho(k))
      field(k)=field(k)-dtn_local*(flx(k)-flx(kb))*rhoi
end do

end subroutine diffuse_scalar1D_rev_full
! ----------------------------------------------------------------------------------------

subroutine diffuse_scalar1D_rev_half (field, bflux, tkh_z, dtn_local)

! input	
real field(:)	! scalar
real tkh_z(:)	! vertical eddy conductivity
real bflux		! bottom flux
real dtn_local         ! timestep
! local        
real flx(0:nzm)
real dfdt
real rdz2, rdz, rdz5
real tkz,rhoi
integer i,j,k,ib,ic,jb,jc,kc,kb, k_flip, kc_flip, kb_flip

rdz=1./dz
rdz2=1./(dz*dz)

flx(0)=bflux*rdz*rhow(1)
flx(nzm)=0.0

dfdt=0.
do k=1,nzm-1
   kc=k+1	
   rhoi = rhow(kc)/adzw(kc)
   rdz5=0.5*rdz2
   tkz=rdz2*tkh_z(k+1)
   flx(k)=-tkz*(field(kc)-field(k))*rhoi
end do

do k=1,nzm
   kb=k-1
   rhoi = 1./(adz(k)*rho(k))
   field(k)=field(k)-dtn_local*(flx(k)-flx(kb))*rhoi
end do

end subroutine diffuse_scalar1D_rev_half

! -------------------------------------------------------------------

subroutine diffuse_scalar1D (field, bflux, tkh_z, dtn_local)

! assumes zero flux at TOA

implicit none
! input	
real field(:)	! scalar
real tkh_z(:)	! vertical eddy conductivity
real bflux		! bottom flux
real dtn_local         ! timestep
! local        
real flx(0:nzm)
real dfdt
real rdz2, rdz
real tkz,rhoi
integer i,j,k,ib,ic,jb,jc,kc,kb, k_flip, kc_flip, kb_flip

rdz=1./dz

dfdt=0.

!-----------------------------------------

!  Vertical diffusion:

rdz2 = 1./(dz*dz)
flx(nzm) = bflux * rdz * rhow(1)
flx(0) = 0.0
do k = 1, nzm-1
   k_flip = nzm - k + 1 
   kc_flip = k_flip - 1
   kc = k + 1	

   rhoi = rhow(kc) / adzw(kc)
   tkz = rdz2 * tkh_z(k_flip)
   flx(kc_flip) = -tkz * (field(kc_flip) - field(k_flip)) * rhoi
end do

do k = 1, nzm

   k_flip = nzm - k
   kb = k - 1
   kb_flip = k_flip + 1

   rhoi = 1. / (adz(k) * rho(k))

   dfdt = -dtn_local * (flx(k_flip) - flx(kb_flip)) * rhoi
   field(k_flip+1) = field(k_flip+1) + dfdt

end do

end subroutine diffuse_scalar1D

! ------------------------------------------

subroutine ddiffuse_scalar1D (field, bflux, tkh_z, dtn_local)

! assumes zero flux at TOA

implicit none
! input	
real(r8) field(:)	! scalar
real(r8) tkh_z(:)	! vertical eddy conductivity
real(r8) bflux		! bottom flux
real(r8) dtn_local         ! timestep
! local        
real(r8) flx(0:nzm)
real(r8) dfdt
real(r8) rdz2, rdz
real(r8) tkz,rhoi
integer i,j,k,ib,ic,jb,jc,kc,kb, k_flip, kc_flip, kb_flip

rdz=1./dble(dz)

dfdt=0.0d0

!-----------------------------------------

!  Vertical diffusion:

rdz2 = 1./dble(dz*dz)
flx(nzm) = bflux * dble(rdz * rhow(1))
flx(0) = 0.0d0
do k = 1, nzm-1
   k_flip = nzm - k + 1 
   kc_flip = k_flip - 1
   kc = k + 1	

   rhoi = dble(rhow(kc) / adzw(kc))
   tkz = rdz2 * tkh_z(k_flip)
   flx(kc_flip) = -tkz * (field(kc_flip) - field(k_flip)) * rhoi
end do

do k = 1, nzm

   k_flip = nzm - k
   kb = k - 1
   kb_flip = k_flip + 1

   rhoi = 1. / dble(adz(k) * rho(k))

   dfdt = -dtn_local * (flx(k_flip) - flx(kb_flip)) * rhoi
   field(k_flip+1) = field(k_flip+1) + dfdt

end do

end subroutine ddiffuse_scalar1D

! ------------------------------------------

subroutine diffuse_scalar1Dw (field, tk_z_full, dtn_local)

! assumes zero flux at TOA

implicit none
! input	
real field(:)	       ! scalar
real tk_z_full(:)      ! vertical eddy conductivity
real dtn_local         ! timestep
! local        
real flx(nz)
real dfdt
real rdz2, rdz
real tkz,rhoi
integer i,j,k,ib,ic,jb,jc,kc,kb, k_flip, kc_flip, kb_flip

rdz=1./dz

dfdt=0.

!-----------------------------------------

!  Vertical diffusion:

rdz2 = 1./(dz*dz)

do k = 1, nzm-1
   k_flip = nzm - k + 1 
   kc_flip = k_flip - 1
   kb_flip = k_flip + 1
   kc = k + 1	

   rhoi = rho(k) / adz(k)
   tkz = rdz2 * tk_z_full(k_flip)
   flx(k_flip) = -2.*tkz * (field(k_flip)-field(kb_flip)) * rhoi
end do

flx(1) = 0.0

do k = 1, nzm-1

   k_flip = nzm - k + 1
   kc_flip = k_flip - 1

   rhoi = 1. / (adzw(k+1) * rhow(k+1))

   dfdt = -dtn_local * (flx(kc_flip) - flx(k_flip)) * rhoi
   field(k_flip) = field(k_flip) + dfdt

end do

end subroutine diffuse_scalar1Dw

! ----------------------------------

subroutine ddiffuse_scalar1Dw (field, tk_z_full, dtn_local)

! assumes zero flux at TOA

implicit none
! input	
real(r8) field(:)	       ! scalar
real(r8) tk_z_full(:)      ! vertical eddy conductivity
real(r8) dtn_local         ! timestep
! local        
real(r8) flx(nz)
real(r8) dfdt
real(r8) rdz2, rdz
real(r8) tkz,rhoi
integer i,j,k,ib,ic,jb,jc,kc,kb, k_flip, kc_flip, kb_flip

rdz=1./dble(dz)

dfdt=0.d0

!-----------------------------------------

!  Vertical diffusion:

rdz2 = 1./dble(dz*dz)

do k = 1, nzm-1
   k_flip = nzm - k + 1 
   kc_flip = k_flip - 1
   kb_flip = k_flip + 1
   kc = k + 1	

   rhoi = dble(rho(k) / adz(k))
   tkz = rdz2 * tk_z_full(k_flip)
   flx(k_flip) = -2.*tkz * (field(k_flip)-field(kb_flip)) * rhoi
end do

flx(1) = 0.0

do k = 1, nzm-1

   k_flip = nzm - k + 1
   kc_flip = k_flip - 1

   rhoi = 1. / dble(adzw(k+1) * rhow(k+1))

   dfdt = -dtn_local * (flx(kc_flip) - flx(k_flip)) * rhoi
   field(k_flip) = field(k_flip) + dfdt

end do

end subroutine ddiffuse_scalar1Dw

subroutine bilinear_uwpbl(indata,outdata)
implicit none
real, intent(in):: indata(0:nxp1,0:nyp1,1:nzm)
real, intent(out):: outdata(nx,ny,nzm)

!---local variables
integer ig,jg,i,ip,j,jp
real wrk(0:nxp1,0:nyp1,1:nzm)

 do ig=1,nx
       i=(ig-1)/uwpbl_ndiv+1
       ip=ig-(i-1)*uwpbl_ndiv
    do jg=0,nyp1
       if(ip.gt.uwpbl_ndiv/2) then
          wrk(ig,jg,:)=indata(ig,jg,1:nzm)*(3*uwpbl_ndiv/2-ip+0.5)+(ip-uwpbl_ndiv/2-0.5)*indata(min(ig+uwpbl_ndiv/2,nxp1),jg,1:nzm)     
       else
     wrk(ig,jg,:)=indata(ig,jg,1:nzm)*(uwpbl_ndiv/2+ip-0.5)+(uwpbl_ndiv/2-ip+0.5)*indata(max(ig-uwpbl_ndiv/2,0),jg,1:nzm)

       endif
    enddo
enddo
wrk=wrk/uwpbl_ndiv
do jg=1,ny
   j=(jg-1)/uwpbl_ndiv+1
   jp=jg-(j-1)*uwpbl_ndiv
   do ig=1,nx
      if(jp.gt.uwpbl_ndiv/2) then
         outdata(ig,jg,:)=wrk(ig,jg,1:nzm)*(3*uwpbl_ndiv/2-jp+0.5)+(jp-uwpbl_ndiv/2-0.5)*wrk(ig,min(jg+uwpbl_ndiv/2,nyp1),1:nzm)
      else
         outdata(ig,jg,:)=wrk(ig,jg,1:nzm)*(uwpbl_ndiv/2+jp-0.5)+(uwpbl_ndiv/2-jp+0.5)*wrk(ig,max(jg-uwpbl_ndiv/2,0),1:nzm)
      endif
   enddo
enddo
outdata=outdata/uwpbl_ndiv
!outdata=indata
end subroutine bilinear_uwpbl
END MODULE uwpbl
