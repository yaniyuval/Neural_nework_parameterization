module diffusion_solver
  !---------------------------------------------------------------------------------
  ! Module to solve vertical diffusion equations using a tri-diagonal solver.
  ! The module will also apply countergradient fluxes, and apply molecular 
  ! diffusion for constituents
  !   
  ! Public interfaces:
  !    init_vdiff       initializes time independent coefficients
  !    compute_vdiff    solves diffusion equations
  !    vd_lu_solve      tridiagonal solver also used by gwd (should be private)
  !    vd_lu_decomp     tridiagonal solver also used by gwd (should be private)
  !    vdiff_selector   type for storing fields selected to be diffused
  !    vdiff_select     selects fields to be diffused
  !    operator(.not.)  extends .not. to operate on type vdiff_selector
  !    any              provides functionality of intrinsic any for type vdiff_selector
  !
  !---------------------------Code history--------------------------------
  ! Initial subroutines:  B. Boville and others, 1991-2004
  ! Modularization:       J. McCaa, September 2004
  !---------------------------------------------------------------------------------

  implicit none
  private          ! Make default type private to the module
  save

!!$  integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real
  ! changed by cw 8/1/05:
!  integer, parameter :: r8 = selected_real_kind(6,30) ! 4 byte real
  integer, parameter :: r8 = selected_real_kind(6,70) ! 8 byte real
  !
  ! Public interfaces
  !
  public init_vdiff      ! Initialization
  public compute_vdiff   ! Full routine
  public vd_lu_solve     ! tridiagonal solver also used by gwd ( should be private! )
  public vd_lu_decomp    ! tridiagonal solver also used by gwd ( should be private! )
  public vdiff_selector  ! type for storing fields selected to be diffused
  public vdiff_select    ! selects fields to be diffused
  public operator(.not.) ! extends .not. to operate on type vdiff_selector
  public any             ! provides functionality of intrinsic any for type vdiff_selector

  type vdiff_selector    ! This stores logical array of fields to be diffused
     private
     logical, pointer, dimension(:) :: fields
  end type vdiff_selector

  interface operator(.not.) ! This extends .not. to operate on type vdiff_selector
     module procedure not
  end interface

  interface any             ! This provides functionality of intrinsic any for type vdiff_selector
     module procedure my_any
  end interface
  !
  ! Private data
  !
  real(r8), private :: cpair      ! Specific heat of dry air
  real(r8), private :: gravit     ! Acceleration due to gravity
  real(r8), private :: rair       ! Gas const for dry air
  real(r8), private :: zvir       ! rh2o/rair - 1
  real(r8), private :: latvap     ! latent heat of vaporization
  real(r8), private :: karman     ! von Karman constant

  ! turbulent mountain stress constants
  real(r8), parameter :: z0fac  = 0.025    ! factor determining z_0 from orographic standard deviation
  real(r8), parameter :: z0max  = 100.     ! max value of z_0 for orography
  real(r8), parameter :: horomin= 10.      ! min value of subgrid orographic height for mountain stress
  real(r8), parameter :: dv2min = 0.01     ! minimum shear squared
  real(r8), private   :: oroconst          ! converts from standard deviation to height

contains

  subroutine init_vdiff(kind, ncnst, rair_in, gravit_in, fieldlist, errstring)

    integer, intent(in) :: kind                    ! Kind used for reals
    integer, intent(in) :: ncnst                   ! number of constituents
    real(r8), intent(in) :: rair_in                ! input gas constant for dry air
    real(r8), intent(in) :: gravit_in              ! input gravititational acceleration
    type(vdiff_selector), intent(out) :: fieldlist ! structure in which to store list of fields to be diffused
    character(128), intent(out) :: errstring       ! output status
    
    errstring = ''
    if ( kind .ne. r8 ) then
       write(6,*) 'KIND of reals passed to init_vdiff -- exiting.'
       errstring = 'init_vdiff'
       return
    endif

    rair = rair_in     
    gravit = gravit_in 
    allocate(fieldlist%fields(3+ncnst))
    fieldlist%fields(:) = .false.
    return

  end subroutine init_vdiff

  subroutine compute_vdiff ( &
       pcols         , pver      , ncnst  , ncol         , pmid      , &
       pint          , rpdel     , t      , ztodt        , taux      , &
       tauy          , shflx     , cflx   , ntop         , nbot      , &
       kvh           , kvm       , kvq    , cgs          , cgh       , &
       zi            , ksrftms   , qmincg , fieldlist    , &
       u             , v         , q      , dse          , &
       tautmsx       , tautmsy   , dtk    , topflx       , errstring , &
       do_molec_diff , compute_molec_diff , vd_lu_qdecomp)
    !-----------------------------------------------------------------------
    ! Driver routine to compute vertical diffusion of momentum, moisture, trace 
    ! constituents and dry static energy. The new temperature is computed from
    ! the diffused dry static energy.

    ! Turbulent diffusivities and boundary layer nonlocal transport terms are 
    ! obtained from the turbulence module.
    !-----------------------------------------------------------------------

    !------------------------------Arguments--------------------------------
    integer, intent(in) :: pcols
    integer, intent(in) :: pver
    integer, intent(in) :: ncnst
    integer, intent(in) :: ncol                    ! number of atmospheric columns
    real(r8), intent(in) :: pmid(pcols,pver)       ! midpoint pressures
    real(r8), intent(in) :: pint(pcols,pver+1)     ! interface pressures
    real(r8), intent(in) :: rpdel(pcols,pver)      ! 1./pdel
    real(r8), intent(in) :: t(pcols,pver)          ! temperature input
    real(r8), intent(in) :: ztodt                  ! 2 delta-t
    real(r8), intent(in) :: taux(pcols)            ! x surface stress (N/m2)
    real(r8), intent(in) :: tauy(pcols)            ! y surface stress (N/m2)
    real(r8), intent(in) :: shflx(pcols)           ! surface sensible heat flux (W/m2)
    real(r8), intent(in) :: cflx(pcols,ncnst)      ! surface constituent flux (kg/m2/s)
    integer,  intent(in) :: ntop                   ! Top level to which vertical diffusion is applied (=1).
    integer,  intent(in) :: nbot                   ! Bottom level to which vertical diffusion is applied (=pver).
    real(r8), intent(in) :: kvh(pcols,pver+1)      ! diffusivity for heat
    real(r8), intent(in) :: kvm(pcols,pver+1)      ! viscosity (diffusivity for momentum)
    real(r8), intent(in) :: kvq(pcols,pver+1)      ! diffusivity for constituents
    real(r8), intent(in) :: cgs(pcols,pver+1)      ! countergradient star (cg/flux)
    real(r8), intent(in) :: cgh(pcols, pver+1)     ! countergradient term for heat
    real(r8), intent(in) :: zi(pcols,pver+1)       ! interface heights
    real(r8), intent(in) :: ksrftms(pcols)         ! surface "drag" coefficient for mountains
    real(r8), intent(in) :: qmincg(ncnst)          ! minimum constituent mixing ratios from cg fluxes
    type(vdiff_selector), intent(in) :: fieldlist  ! array of flags selecting which fields to diffuse

    real(r8), intent(inout) :: u(pcols,pver)       ! u wind
    real(r8), intent(inout) :: v(pcols,pver)       ! v wind
    real(r8), intent(inout) :: q(pcols,pver,ncnst) ! moisture and trace constituents
    real(r8), intent(inout) :: dse(pcols,pver)     ! dry static energy [J/kg]

    real(r8), intent(out) :: dtk(pcols,pver)       ! T tendency from KE dissipation
    real(r8), intent(out) :: tautmsx(pcols)        ! turbulent mountain surface stress (zonal)
    real(r8), intent(out) :: tautmsy(pcols)        ! turbulent mountain surface stress (meridional)
    real(r8), intent(out) :: topflx(pcols)         ! molecular heat flux at top interface
    character(128), intent(out) :: errstring       ! output status

    logical, intent(in) :: do_molec_diff           ! flag indicating multiple constituent diffusivities
    integer, external, optional :: compute_molec_diff ! constituent-independent moleculuar diffusivity routine
    integer, external, optional :: vd_lu_qdecomp      ! constituent-dependent moleculuar diffusivity routine
    !
    !---------------------------Local workspace-----------------------------
    !
    real(r8) :: tmpm(pcols,pver)                   ! potential temperature, ze term in tri-diag sol'n
    real(r8) :: ca(pcols,pver)                     ! -upper diag of tri-diag matrix
    real(r8) :: cc(pcols,pver)                     ! -lower diag of tri-diag matrix
    real(r8) :: dnom(pcols,pver)                   ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))

    real(r8) :: tmp1(pcols)                        ! temporary storage
    real(r8) :: tmpi1(pcols,pver+1)                ! interface KE dissipation
    real(r8) :: tint(pcols,pver+1)                 ! interface temperature
    real(r8) :: rhoi(pcols,pver+1)                 ! rho at interfaces
    real(r8) :: tmpi2(pcols,pver+1)                ! dt*(g*rho)**2/dp at interfaces
    real(r8) :: rrho(pcols)                        ! 1./bottom level density 

    real(r8) :: zero   (pcols)                     ! zero array for surface heat exchange coefficients 
    real(r8) :: tautotx(pcols)                     ! total surface stress (zonal)
    real(r8) :: tautoty(pcols)                     ! total surface stress (meridional)

    real(r8) :: dinp_u(pcols,pver+1)               ! vertical difference at interfaces, input u
    real(r8) :: dinp_v(pcols,pver+1)               ! vertical difference at interfaces, input v
    real(r8) :: dout_u                             ! vertical difference at interfaces, output u
    real(r8) :: dout_v                             ! vertical difference at interfaces, output v
    real(r8) :: dse_top(pcols)                     ! dse on top boundary
    real(r8) :: cc_top(pcols)                      ! lower diagonal at top interface
    real(r8) :: cd_top(pcols)                      ! 
    real(r8) :: rghd(pcols,pver+1)                 ! (1/H_i - 1/H) *(g*rho)^(-1)

    integer  :: i,k,m                              ! longitude,level,constituent indices
    logical  :: lqtst(pcols)                       ! adjust vertical profiles
    real(r8) :: qtm(pcols,pver)                    ! temporary copy of q
    logical  :: need_decomp                        ! whether to compute a new decomposition
    integer  :: status                             ! status indicator
    integer  :: ntop_molec                         ! top level where molecular diffusivity is applied
    integer  :: nbot_molec                         ! bottom level where molecular diffusivity is applied
    real(r8) :: kq_scal(pcols, pver+1)             ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
    real(r8) :: mw_fac(ncnst)                      ! sqrt(1/M_q + 1/M_d) for this constituent
    real(r8) :: cnst_mw(ncnst)                     ! molecular weight (kg/kmole)
    logical  :: cnst_fixed_ubc(ncnst)              ! whether upper bndy condition is fixed
    real(r8) :: ubc_mmr(pcols,ncnst)               ! upper bndy mixing ratios (kg/kg)
    real(r8) :: ubc_t(pcols)                       ! upper bndy temperature (K)

    errstring = ''
    if ((diffuse(fieldlist,'u').or.diffuse(fieldlist,'v')) .and. .not.diffuse(fieldlist,'s')) then
       errstring = 'diffusion_solver.compute_vdiff: must diffuse s if diffusing u or v'
       return
    end if
    zero(:) = 0.

    ! Compute rho at interfaces p/RT,  Ti = (Tm_k + Tm_k-1)/2,  interface k-
    ! Calculate dt * (g*rho)^2/dp at interior interfaces,  interface k-
    tint(:ncol,1) = t(:ncol,1)
    rhoi(:ncol,1) = pint(:ncol,1) / (rair*tint(:ncol,1))
    do k = 2, pver
       do i = 1, ncol
          tint(i,k) = 0.5 * (t(i,k) + t(i,k-1))
          rhoi(i,k) = pint(i,k) / (rair*tint(i,k))
          tmpi2(i,k) = ztodt * (gravit*rhoi(i,k))**2 / (pmid(i,k) - pmid(i,k-1))
       end do
    end do
    tint(:ncol,pver+1) = t(:ncol,pver)
    rhoi(:ncol,pver+1) = pint(:ncol,pver+1) / (rair*tint(:ncol,pver+1))

    rrho(:ncol)   = rair*t(:ncol,pver)/pmid(:ncol,pver)
    tmp1(:ncol)  = ztodt*gravit*rpdel(:ncol,pver)

    !-----------------------------------------------------------------------
    !    Computation of molecular diffusivities
    !-----------------------------------------------------------------------
    if ( do_molec_diff ) then
       if((.not.present(compute_molec_diff)).or.(.not.present(vd_lu_qdecomp))) then
          errstring = 'compute_vdiff: do_molec_diff true but compute_molec_diff or vd_lu_qdecomp missing'
          return
       endif
       ! The next call:
       !     modifies: kvh, kvm, tint, rhoi, and tmpi2
       !     returns:  kq_scal, ubc_mmr, dse_top, cc_top, cd_top, cnst_mw, 
       !     cnst_fixed_ubc , mw_fac , ntop_molec
       status = compute_molec_diff( &
            pcols          , pver   , ncnst      , ncol    , t      , pmid   , pint    , &
            zi             , ztodt  , kvh        , kvm     , tint   , rhoi   , tmpi2   , &
            kq_scal        , ubc_mmr    , dse_top , cc_top , cd_top , cnst_mw , &
            cnst_fixed_ubc , mw_fac , ntop_molec , nbot_molec)
    else
       kq_scal(:,:) = 0.
       cd_top(:) = 0.
       cc_top(:) = 0.
    endif
    !-----------------------------------------------------------------------

    if (diffuse(fieldlist,'u').or.diffuse(fieldlist,'v')) then
       ! Compute the vertical differences of the input u,v for KE dissipation, interface k-
       ! Note, velocity=0 at surface, so then difference at the bottom interface is -u,v(pver)
       do i = 1, ncol
          dinp_u(i,1) = 0.
          dinp_v(i,1) = 0.
          dinp_u(i,pver+1) = -u(i,pver)
          dinp_v(i,pver+1) = -v(i,pver)
       end do
       do k = 2, pver
          do i = 1, ncol
             dinp_u(i,k) = u(i,k) - u(i,k-1)
             dinp_v(i,k) = v(i,k) - v(i,k-1)
          end do
       end do

       ! Add the explicit surface fluxes to the lowest layer (do not include moutain stress)
       u(:ncol,pver) = u(:ncol,pver) + tmp1(:ncol) * taux(:ncol)
       v(:ncol,pver) = v(:ncol,pver) + tmp1(:ncol) * tauy(:ncol)

       ! Diffuse momentum
       call vd_lu_decomp(                                        &
            pcols, pver, ncol     ,                                 &
            ksrftms  ,kvm      ,tmpi2    ,rpdel    ,ztodt    ,zero,   &
            ca       ,cc       ,dnom     ,tmpm     ,ntop     ,nbot     )
       call vd_lu_solve(                                                &
            pcols, pver, ncol     ,                                        &
            u        ,ca       ,tmpm     ,dnom     ,ntop     ,nbot,zero     )
       call vd_lu_solve(                                                &
            pcols, pver, ncol     ,                                        &
            v        ,ca       ,tmpm     ,dnom     ,ntop     ,nbot,zero     )

       ! Calculate the mountain and total stresses
       do i = 1, ncol
          tautmsx(i) = - ksrftms(i) * u(i,pver)
          tautmsy(i) = - ksrftms(i) * v(i,pver)
          tautotx(i) = taux(i) + tautmsx(i)
          tautoty(i) = tauy(i) + tautmsy(i)
       end do

       ! Calculate kinetic energy dissipation
       ! 1. compute dissipation term at interfaces
       k = pver+1
       do i = 1, ncol
          tmpi1(i,1) = 0.
          tmpi1(i,k) = 0.5 * ztodt * gravit * &
               ( (-u(i,k-1) + dinp_u(i,k))*tautotx(i) + (-v(i,k-1)+dinp_v(i,k))*tautoty(i) )
       end do
       do k = 2, pver
          do i = 1, ncol
             dout_u = u(i,k) - u(i,k-1)
             dout_v = v(i,k) - v(i,k-1)
             tmpi1(i,k) = 0.25 * tmpi2(i,k) * kvm(i,k) *   &
                  (dout_u**2 + dout_v**2 + dout_u*dinp_u(i,k) + dout_v*dinp_v(i,k))
          end do
       end do
       ! 2. Compute dissipation term at midpoints, add to dry static energy
       do k = 1, pver
          do i = 1, ncol
             dtk(i,k) = (tmpi1(i,k+1) + tmpi1(i,k)) * rpdel(i,k)
             dse(i,k) = dse(i,k) + dtk(i,k)
          end do
       end do
    end if

    if (diffuse(fieldlist,'s')) then
       ! Add countergradient to input static energy profiles
       do k=1,pver
          dse(:ncol,k) = dse(:ncol,k) + ztodt*rpdel(:ncol,k)*gravit  &
               *(rhoi(:ncol,k+1)*kvh(:ncol,k+1)*cgh(:ncol,k+1)       &
               - rhoi(:ncol,k  )*kvh(:ncol,k  )*cgh(:ncol,k  ))
       end do
       ! Add the explicit surface fluxes to the lowest layer
       dse(:ncol,pver) = dse(:ncol,pver) + tmp1(:ncol) * shflx(:ncol)
       ! Diffuse dry static energy
       call vd_lu_decomp(                                      &
            pcols, pver, ncol     ,                               &
            zero     ,kvh      ,tmpi2    ,rpdel    ,ztodt    ,cc_top, &
            ca       ,cc       ,dnom     ,tmpm     ,ntop     ,nbot    )
       call vd_lu_solve(                                               &
            pcols, pver, ncol     ,                                       &
            dse      ,ca       ,tmpm     ,dnom     ,ntop     ,nbot,cd_top    )
       ! calculate flux at top interface
       if (do_molec_diff) then
          topflx(:ncol) = - kvh(:ncol,ntop_molec) * tmpi2(:ncol,ntop_molec) / (ztodt*gravit) &
               * (dse(:ncol,ntop_molec) - dse_top(:ncol))
       end if
    endif

    ! Loop through constituents
    need_decomp = .true.
    do m = 1, ncnst
       if (diffuse(fieldlist,'q',m)) then
          ! Add the nonlocal transport terms to constituents in the boundary layer.
          ! Check for neg q's in each constituent and put the original vertical
          ! profile back if a neg value is found. A neg value implies that the
          ! quasi-equilibrium conditions assumed for the countergradient term are
          ! strongly violated.
          qtm(:ncol,:pver) = q(:ncol,:pver,m)
          do k = 1, pver
             q(:ncol,k,m) = q(:ncol,k,m) + &
                  ztodt*rpdel(:ncol,k)*gravit*(cflx(:ncol,m)*rrho(:ncol))   &
                  *(rhoi(:ncol,k+1)*kvh(:ncol,k+1)*cgs(:ncol,k+1)           &
                  - rhoi(:ncol,k  )*kvh(:ncol,k  )*cgs(:ncol,k  ))
          end do
          lqtst(:ncol) = all(q(:ncol,1:pver,m) >= qmincg(m), 2)
          do k = 1, pver
             q(:ncol,k,m) = merge (q(:ncol,k,m), qtm(:ncol,k), lqtst(:ncol))
          end do
          ! Add the explicit surface fluxes to the lowest layer
          q(:ncol,pver,m) = q(:ncol,pver,m) + tmp1(:ncol) * cflx(:ncol,m)
          ! Diffuse constituents.
          if ( need_decomp ) then
             call vd_lu_decomp(   &
                  pcols, pver, ncol     ,                               &
                  zero     ,kvq , tmpi2    ,rpdel    ,ztodt, zero, &
                  ca       ,cc       ,dnom     ,tmpm     ,ntop, nbot)
             if(do_molec_diff)then
                ! Update decomposition in molecular diffusion range, include separation velocity term
                status = vd_lu_qdecomp(                                               &
                     pcols, pver, ncol       ,cnst_fixed_ubc(m)      ,cnst_mw(m) ,ubc_mmr(:,m),&
                     kvq        ,kq_scal    ,mw_fac(m)  ,tmpi2      ,rpdel      , &
                     ca         ,cc         ,dnom       ,tmpm       ,rhoi       , &
                     tint       ,ztodt      ,ntop_molec, nbot_molec ,cd_top     )
             else
                need_decomp =  .false.
             endif
          end if
          call vd_lu_solve(                                              &
               pcols, pver, ncol     ,               &
               q(1,1,m) , ca       ,tmpm     ,dnom     ,ntop    ,nbot,cd_top     )
       end if
    end do

    return
  end subroutine compute_vdiff

  !==============================================================================
  subroutine vd_lu_decomp(                                          &
       pcols, pver, ncol       ,                                     &
       ksrf       ,kv         ,tmpi       ,rpdel      ,ztodt      , cc_top, &
       ca         ,cc         ,dnom       ,ze         ,ntop       , &
       nbot       )
    !-----------------------------------------------------------------------
    ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the 
    ! tridiagonal diffusion matrix. 
    ! The diagonal elements (1+ca(k)+cc(k)) are not required by the solver.
    ! Also determine ze factor and denominator for ze and zf (see solver).
    !------------------------------Arguments--------------------------------
    integer, intent(in)  :: pcols                  ! number of allocated atmospheric columns
    integer, intent(in)  :: pver                   ! number of allocated atmospheric levels 
    integer, intent(in)  :: ncol                   ! number of computed atmospheric columns
    integer, intent(in)  :: ntop                   ! top level to operate on
    integer, intent(in)  :: nbot                   ! bottom level to operate on
    real(r8), intent(in) :: ksrf(pcols)            ! surface "drag" coefficient
    real(r8), intent(in) :: kv(pcols,pver+1)       ! vertical diffusion coefficients
    real(r8), intent(in) :: tmpi(pcols,pver+1)     ! dt*(g/R)**2/dp*pi(k+1)/(.5*(tm(k+1)+tm(k))**2
    real(r8), intent(in) :: rpdel(pcols,pver)      ! 1./pdel  (thickness bet interfaces)
    real(r8), intent(in) :: ztodt                  ! 2 delta-t
    real(r8), intent(in) :: cc_top(pcols)          ! lower diagonal on top interface (for fixed ubc only)

    real(r8), intent(out) :: ca(pcols,pver)        ! upper diagonal
    real(r8), intent(out) :: cc(pcols,pver)        ! lower diagonal
    real(r8), intent(out) :: dnom(pcols,pver)      ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
    real(r8), intent(out) :: ze(pcols,pver)        ! term in tri-diag. matrix system

    !---------------------------Local workspace-----------------------------
    integer :: i                                   ! longitude index
    integer :: k                                   ! vertical index
    !-----------------------------------------------------------------------
    ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the 
    ! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are
    ! a combination of ca and cc; they are not required by the solver.
!!$    call t_startf('vd_lu_decomp') ! commented out by cw 7/27/05

    do k = nbot-1, ntop, -1
       do i = 1, ncol
          ca(i,k  ) = kv(i,k+1)*tmpi(i,k+1)*rpdel(i,k  )
          cc(i,k+1) = kv(i,k+1)*tmpi(i,k+1)*rpdel(i,k+1)
       end do
    end do

    ! The bottom element of the upper diagonal (ca) is zero (element not used).
    ! The subdiagonal (cc) is not needed in the solver.

    do i=1,ncol
       ca(i,nbot) = 0.
    end do

    ! Calculate e(k).  This term is 
    ! required in solution of tridiagonal matrix defined by implicit diffusion eqn.

    do i = 1, ncol
       dnom(i,nbot) = 1./(1. + cc(i,nbot) + ksrf(i)*ztodt*gravit*rpdel(i,nbot))
       ze(i,nbot)   = cc(i,nbot)*dnom(i,nbot)
    end do
    do k = nbot-1, ntop+1, -1
       do i=1,ncol
          dnom(i,k) = 1./ (1. + ca(i,k) + cc(i,k) - ca(i,k)*ze(i,k+1))
          ze(i,k) = cc(i,k)*dnom(i,k)
       end do
    end do
    do i=1,ncol
       dnom(i,ntop) = 1./ (1. + ca(i,ntop) + cc_top(i) - ca(i,ntop)*ze(i,ntop+1))
    end do

!!$    call t_stopf('vd_lu_decomp') ! commented out by cw 7/27/05
    return
  end subroutine vd_lu_decomp

  !==============================================================================
  subroutine vd_lu_solve(                                           &
       pcols, pver, ncol       ,                                     &
       q          ,ca         ,ze         ,dnom       ,             &
       ntop       ,nbot       , cd_top)
    !-----------------------------------------------------------------------
    ! Solve the implicit vertical diffusion equation with zero flux boundary conditions.
    ! Actual surface fluxes are explicit (rather than implicit) and are applied separately. 
    ! Procedure for solution of the implicit equation follows Richtmyer and 
    ! Morton (1967,pp 198-200).

    ! The equation solved is
    !
    ! -ca(k)*q(k+1) + cb(k)*q(k) - cc(k)*q(k-1) = d(k), 
    !
    ! where d(k) is the input profile and q(k) is the output profile
    !
    ! The solution has the form
    !
    ! q(k)  = ze(k)*q(k-1) + zf(k)
    !
    ! ze(k) = cc(k) * dnom(k)
    !
    ! zf(k) = [d(k) + ca(k)*zf(k+1)] * dnom(k)
    !
    ! dnom(k) = 1/[cb(k) - ca(k)*ze(k+1)] =  1/[1 + ca(k) + cc(k) - ca(k)*ze(k+1)]

    ! Note that the same routine is used for temperature, momentum and tracers,
    ! and that input variables are replaced.
    !------------------------------Arguments--------------------------------
    integer, intent(in)  :: pcols                  ! number of allocated atmospheric columns
    integer, intent(in)  :: pver                   ! number of allocated atmospheric levels 
    integer, intent(in)  :: ncol                   ! number of computed atmospheric columns
    integer, intent(in)  :: ntop                   ! top level to operate on
    integer, intent(in)  :: nbot                   ! bottom level to operate on
    real(r8), intent(in) :: ca(pcols,pver)         ! -upper diag coeff.of tri-diag matrix
    real(r8), intent(in) :: ze(pcols,pver)         ! term in tri-diag solution
    real(r8), intent(in) :: dnom(pcols,pver)       ! 1./(1. + ca(k) + cc(k) - ca(k)*ze(k+1))
    real(r8), intent(in) :: cd_top(pcols)          ! cc_top * ubc value

    real(r8), intent(inout) :: q(pcols,pver)       ! constituent field

    !---------------------------Local workspace-----------------------------
    real(r8) :: zf(pcols,pver)                     ! term in tri-diag solution
    integer i,k                                    ! longitude, vertical indices
    !-----------------------------------------------------------------------

    ! Calculate zf(k).  Terms zf(k) and ze(k) are required in solution of 
    ! tridiagonal matrix defined by implicit diffusion eqn.
    ! Note that only levels ntop through nbot need be solved for.

    do i = 1, ncol
       zf(i,nbot) = q(i,nbot)*dnom(i,nbot)
    end do
    do k = nbot-1, ntop+1, -1
       do i=1,ncol
          zf(i,k) = (q(i,k) + ca(i,k)*zf(i,k+1))*dnom(i,k)
       end do
    end do

    ! Include boundary condition on top element
    k = ntop
    do i=1, ncol
       zf(i,k) = (q(i,k) + cd_top(i) + ca(i,k)*zf(i,k+1))*dnom(i,k)
    end do

    ! Perform back substitution

    do i=1,ncol
       q(i,ntop) = zf(i,ntop)
    end do
    do k=ntop+1, nbot, +1
       do i=1,ncol
          q(i,k) = zf(i,k) + ze(i,k)*q(i,k-1)
       end do
    end do

    return
  end subroutine vd_lu_solve
  
  character(128) function vdiff_select(fieldlist,name,qindex)
    ! This function sets the field with incoming name as one to be diffused
    type(vdiff_selector), intent(inout) :: fieldlist
    character(*), intent(in) :: name
    integer, intent(in), optional :: qindex
    
    vdiff_select = ''
    select case (name)
    case ('u','U')
       fieldlist%fields(1) = .true.
    case ('v','V')
       fieldlist%fields(2) = .true.
    case ('s','S')
       fieldlist%fields(3) = .true.
    case ('q','Q')
       if ( present(qindex) ) then
          fieldlist%fields(3 + qindex) = .true.
       else
          fieldlist%fields(4) = .true.
       endif
    case default
       write(vdiff_select,*) 'Bad argument to vdiff_index: ',name
    end select
    return
    
  end function vdiff_select

  type(vdiff_selector) function not(a)
    ! This function extends .not. to operate on type vdiff_selector    
    type(vdiff_selector), intent(in)  :: a
    allocate(not%fields(size(a%fields)))
    not%fields(:) = .not. a%fields(:)
  end function not

  logical function my_any(a)
    ! This function extends the intrinsic function 'any' to operate on type vdiff_selector
    type(vdiff_selector), intent(in) :: a
    my_any = any(a%fields)
  end function my_any

  logical function diffuse(fieldlist,name,qindex)
    ! This function reports whether the field with incoming name is to be diffused
    type(vdiff_selector), intent(in) :: fieldlist
    character(*), intent(in) :: name
    integer, intent(in), optional :: qindex
    
    select case (name)
    case ('u','U')
       diffuse = fieldlist%fields(1)
    case ('v','V')
       diffuse = fieldlist%fields(2)
    case ('s','S')
       diffuse = fieldlist%fields(3)
    case ('q','Q')
       if ( present(qindex) ) then
          diffuse = fieldlist%fields(3 + qindex)
       else
          diffuse = fieldlist%fields(4)
       endif
    case default
       diffuse = .false.
    end select
    return
  end function diffuse

end module diffusion_solver
