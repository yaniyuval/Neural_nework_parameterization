module uwcu

  use vars
  use params

  implicit none

  private
  save

  ! added by cw 10/5/05
  public  mcshallow

!  integer, parameter :: r8 = selected_real_kind(6,30) ! 4 byte real
  integer, parameter :: r8 = selected_real_kind(6,70) ! 8 byte real

  real(r8) :: zvir  ! rh2o/rair - 1
  real(r8) :: ep2   ! mol wgt water vapor / mol wgt dry air 
  real(r8) :: p00   ! reference pressure for exner function
  real(r8) :: rovcp ! R/cp

  ! added by cw 10/4/05
  integer, parameter :: ncnst = 1              ! number of scalars

  real(r8), allocatable, dimension(:, :) ::                                   &
       uwcu_umf0,          &   ! upward mass flux
       uwcu_ps0,           &   ! pressure on half sigma levels
       uwcu_zs0,           &   ! height on half sigma levels
       uwcu_tke0,          &   ! turbulent kinetic energy
       uwcu_slflx0,        &   ! updraft liquid static energy flux
       uwcu_qtflx0,        &   ! updraft total water flux
       uwcu_qv0,           &   ! water vapor specific humidity
       uwcu_ql0,           &   ! cloud water mixing ratio
       uwcu_qi0,           &   ! cloud ice mixing ratio
       uwcu_t0,            &   ! temperature
       uwcu_u0,            &   ! zonal wind
       uwcu_v0,            &   ! meridional wind
       uwcu_dp,            &   ! pressure differences between layers
       uwcu_p0,            &   ! pressure on half sigma levels
       uwcu_s0,            &   ! dry static energy
       uwcu_z0,            &   ! height on half sigma levels
       ql00,               &   ! non precip. water
       qi00,               &   ! non precip. water
       sten0,              &   ! tendency of dry static energy
       uten0,              &   ! tendency of zonal wind
       vten0,              &   ! tendency of meridional wind
       qvten0,             &   ! tendency of water vapor specific humidity
       qlten0,             &   ! tendency of cloud water mixing ratio
       qiten0,             &   ! tendency of cloud ice mixing ratio
       qrten0,             &   ! rain tendency
       qsten0,             &   ! snow tendency
       cldfrac0,           &   ! cloud fraction from shallow convection
       qlu0,               &   ! updraft liquid water mixing ratio
       final_qi_sh,        &   ! ice specific humidity post-convection
       final_ql_sh,        &   ! liquid water specific humidity post-convection
       qc_out0,            &   ! recently added by Sungsu (not int. into WRF)
       fer0,               &   ! lateral entrainment rate
       fdr0                    ! lateral detrainment rate
  
  real(r8), dimension(:), allocatable ::                                      &
       cnt,                &   ! recently added by Sungsu
       cnb,                &   ! recently added by Sungsu
       rliq,               &   ! integral of qc
       pblh0,              &   ! PBL height
       precip0,            &   ! precipitation
       cush0,              &   ! convective scale height
       cbmf0,              &   ! cloud base mass flux
       cin0,               &   ! convective inhibition
       uw_err_1d               ! error flag

  integer ::                                                                  &
       nx_p                    ! number of patches in x direction

  logical ::                                                                  &
       initialized = .false.

contains
  
  real(r8) function exnf(pressure)
    real(r8), intent(in) :: pressure
    exnf = (pressure/p00) ** rovcp
    return
  end function exnf

  subroutine init_mcshallow
    !------------------------------------------------------------------------- 
    ! Purpose:  
    ! Initialize time independent variables of the shallow convection package.
    !-------------------------------------------------------------------------
    ! added by cw 10/4/05
    use wv_saturation, only: gestbl

    implicit none

    if(mod(nx,uwcu_ndiv).ne.0.or.(RUN3D.and.mod(ny,uwcu_ndiv).ne.0)) then
       if(masterproc) print*,'nx or ny is not divisible by uwcu_ndiv'
       if(masterproc) print*,'set in uwcu.f90'
       if(masterproc) print*,'Stop.'
       call task_abort()
    end if

    nx_p = 1
    allocate(uwcu_umf0   (nx_p, 0:nzm))
    allocate(uwcu_ps0    (nx_p, 0:nzm))
    allocate(uwcu_zs0    (nx_p, 0:nzm))
    allocate(uwcu_tke0   (nx_p, 0:nzm))
    allocate(uwcu_slflx0 (nx_p, 0:nzm))
    allocate(uwcu_qtflx0 (nx_p, 0:nzm))
    allocate(uwcu_qv0    (nx_p, nzm))
    allocate(uwcu_ql0    (nx_p, nzm))
    allocate(uwcu_qi0    (nx_p, nzm))
    allocate(uwcu_t0     (nx_p, nzm))
    allocate(uwcu_u0     (nx_p, nzm))
    allocate(uwcu_v0     (nx_p, nzm))
    allocate(uwcu_dp     (nx_p, nzm))
    allocate(uwcu_p0     (nx_p, nzm))
    allocate(uwcu_s0     (nx_p, nzm))
    allocate(uwcu_z0     (nx_p, nzm))
    allocate(ql00        (nx_p, nzm))
    allocate(qi00        (nx_p, nzm))
    allocate(sten0       (nx_p, nzm))
    allocate(uten0       (nx_p, nzm))
    allocate(vten0       (nx_p, nzm))
    allocate(qvten0      (nx_p, nzm))
    allocate(qlten0      (nx_p, nzm))
    allocate(qiten0      (nx_p, nzm))
    allocate(qrten0      (nx_p, nzm))
    allocate(qsten0      (nx_p, nzm))
    allocate(cldfrac0    (nx_p, nzm))
    allocate(qlu0        (nx_p, nzm))
    allocate(final_qi_sh (nx_p, nzm))
    allocate(final_ql_sh (nx_p, nzm))
    allocate(qc_out0     (nx_p, nzm))
    allocate(fer0        (nx_p, nzm))
    allocate(fdr0        (nx_p, nzm))

    allocate(cnt         (nx_p))
    allocate(cnb         (nx_p))
    allocate(rliq        (nx_p))
    allocate(pblh0       (nx_p))
    allocate(precip0     (nx_p))
    allocate(cush0       (nx_p))
    allocate(cbmf0       (nx_p))
    allocate(cin0        (nx_p))
    allocate(uw_err_1d   (nx_p))

    zvir = dble(rv/rgas) - 1.0d0
    ep2 = 0.622d0
    p00 = 1.E5
    rovcp = dble(rgas/cp)

    ! added by cw 10/4/05
    call gestbl(            173.16d0,                        375.16d0,        &
                               20.d0,                          .true.,        &
                             0.622d0,                     dble(lcond),        &
                          dble(lfus),                        dble(rv),        &
                            dble(cp),                      dble(tice) )

    initialized = .true.
    ! end cw addition
  end subroutine init_mcshallow

! -----------------------------------------------------------------------------

  subroutine mcshallow(    dtn_local,                          p_flip,        &
                            p8w_flip,                          u_flip,        &
                              v_flip,                       tabs_flip,        &
                              s_flip,                         qv_flip,        &
                             qc_flip,                         qi_flip,        &
                            tke_flip,                      pblh_local,        &
                          cush_local,                         uwcu_du,        &
                             uwcu_dv,                         uwcu_ds,        &
                             uwcu_dq,                        uwcu_dqv,        &
                            uwcu_dql,                        uwcu_dqi )

  use wv_saturation, only: fqsatd

! ------------------------------ input arguments ------------------------------

  real(r8), intent(in) ::                                                     &
       dtn_local,          &   ! timestep
       pblh_local              ! pbl height

  real(r8), dimension(:), intent(in) ::                                       &
       p_flip,             &   ! mid-level pressure
       p8w_flip,           &   ! interface pressure
       u_flip,             &   !
       v_flip,             &   !
       tabs_flip,          &   !
       s_flip,             &   !
       qv_flip,            &   !
       qc_flip,            &   !
       qi_flip,            &   !
       tke_flip

! ------------------------------ output arguments -----------------------------

  real(r8), intent(inout) ::                                                  &
       cush_local    ! convective scale height

  real(r8), dimension(:), intent(out) ::                                      &
       uwcu_du,            &   !
       uwcu_dv,            &   !
       uwcu_ds,            &   !
       uwcu_dq,            &   !        
       uwcu_dqv,           &   !        
       uwcu_dql,           &   !        
       uwcu_dqi

! ------------------------------ local variables ------------------------------

    integer ::                                                                &
         i,                &   ! longitude loop counter
         j,                &   ! latitude loop counter
         ip,               &   ! longitude loop counter within patch
         jp,               &   ! latitude loop counter within patch
         ig,               &   ! global longitude loop counter 
         jg,               &   ! global latitude loop counter 
         k,                &   ! level loop counter
         k_flip,           &
         ic,               &
         jc


! ------------------------------ executable code ------------------------------
    
    if(.not.initialized) call init_mcshallow

    pblh0(1) = pblh_local
    cush0(1) = cush_local
    
    uwcu_zs0(1, 0:nzm) = dble(zi(1:nz))
    uwcu_z0(1, :) = dble(z(1:nzm))
    
    do k = 1, nz
       k_flip = nz - k + 1
       uwcu_ps0(1, k-1) = p8w_flip(k_flip)
       uwcu_tke0(1, k-1) = tke_flip(k_flip)
    enddo

    do k = 1, nzm
       
       k_flip = nzm - k + 1

       uwcu_p0(1, k) = p_flip(k_flip)
       uwcu_dp(1, k) = (uwcu_ps0(1, k-1) - uwcu_ps0(1, k))
       
       uwcu_t0(1, k) = tabs_flip(k_flip)
       uwcu_u0(1, k) = u_flip(k_flip)
       uwcu_v0(1, k) = v_flip(k_flip)
       
       
       ! shallow scheme wants water vapor as specific humidity
       ! cloud water, cloud ice, rain, and snow as mixing ratio
       ! SAM is SH for all, so cloud water, cloud ice, rain, and 
       ! snow are converted
       
       uwcu_qv0(1, k) = qv_flip(k_flip)
       
       ! convert from specific humidity to mixing ratio
       uwcu_ql0(1, k) = qc_flip(k_flip) / (1.d0 - qc_flip(k_flip))              
       uwcu_qi0(1, k) = qi_flip(k_flip) / (1.d0 - qi_flip(k_flip))              

       uwcu_s0(1, k) = s_flip(k_flip)
    enddo

    call compute_mcshallow(                                                   &
                                   1,                             nzm,        &
                                   1,                           ncnst,        &
                           dtn_local,                        uwcu_ps0,        &
                            uwcu_zs0,                         uwcu_p0,        &
                             uwcu_z0,                         uwcu_dp,        &
                             uwcu_u0,                         uwcu_v0,        &
                            uwcu_qv0,                        uwcu_ql0,        &
                            uwcu_qi0,                         uwcu_t0,        &
                             uwcu_s0,                       uwcu_tke0,        &
                               pblh0,                           cush0,        &
                           uwcu_umf0,                     uwcu_slflx0,        &
                         uwcu_qtflx0,                          qvten0,        &
                              qlten0,                          qiten0,        &
                               sten0,                           uten0,        &
                               vten0,                          qrten0,        &
                              qsten0,                         precip0,        &
                            cldfrac0,                            qlu0,        &
                                fer0,                            fdr0,        &
                                cin0,                           cbmf0,        &
                             qc_out0,                            rliq,        &
                                 cnt,                             cnb,        &
                           uw_err_1d,                          fqsatd )
    
       final_ql_sh = (uwcu_ql0 + qlten0*dtn_local) / (1.d0 + uwcu_ql0 + qlten0*dtn_local)
       final_qi_sh = (uwcu_qi0 + qiten0*dtn_local) / (1.d0 + uwcu_qi0 + qiten0*dtn_local)

       do k = 1, nzm
          k_flip = nzm - k + 1
                   
!!$          qvtend_uwcu(k) = qvten0(1, k)
!!$          qltend_uwcu(k) = d_ql / dtn
!!$          qitend_uwcu(k) =  d_qi / dtn
          

          uwcu_dql(k_flip) = final_ql_sh(1, k) - uwcu_ql0(1, k)
          uwcu_dqi(k_flip) = final_qi_sh(1, k) - uwcu_qi0(1, k)
          uwcu_dqv(k_flip) = qvten0(1, k) * dtn_local
          uwcu_du(k_flip) = uten0(1, k) * dtn_local
          uwcu_dv(k_flip) = vten0(1, k) * dtn_local
          uwcu_ds(k_flip) = sten0(1, k) * dtn_local

          uwcu_dq(k) = uwcu_dqv(k) + uwcu_dql(k) + uwcu_dqi(k)
          
!!$          fer(k) = fer0(1, k)
!!$          fdr(k) = fdr0(1, k)
!!$          umf(k) = uwcu_umf0(1, k-1)
!!$          slflx_uwcu( k) = uwcu_slflx0(i, k-1)
!!$          qtflx_uwcu(ig, jg, k) = uwcu_qtflx0(i, k-1)
!!$          
!!$          cldfrc(ig, jg, k) = cldfrac0(i, k)
       enddo
       
!!$       cin(ig, jg) = cin0(i)
!!$       cbmf(ig, jg) = cbmf0(i)
!!$       uwcu_err(ig, jg) = uw_err_1d(i)
       
     end subroutine mcshallow

! ---------------------------------------------------------------------------------------------------

  subroutine compute_mcshallow( mix      , mkx      , iend      , ncnst     , dt     ,           &
                                ps0_in   , zs0_in   , p0_in     , z0_in     , dp0_in ,           &
                                u0_in    , v0_in    , qv0_in    , ql0_in    , qi0_in ,           &
                                t0_in    , s0_in    ,                                            &
                                tke_in   , pblh_in  , cush_inout,                                & 
                                umf_out  , slflx_out, qtflx_out ,                                &
                                qvten_out, qlten_out, qiten_out ,                                & 
                                sten_out , uten_out , vten_out  ,                                &
                                qrten_out, qsten_out, precip_out, cldfrc_out, qlu_out,           &
                                fer_out  , fdr_out  , cinh_out  , cbmf_out, qc_out, rliq_out, cnt_out, cnb_out, err, qsat )

    !======================================================!
    !                                                      ! 
    !     SHALLOW CONVECTION SCHEME                        !
    !     Described in McCaa, Bretherton, and Grenier:     !
    !     (submitted to MWR, December 2001)                !
    !                                                      !
    !     Modified by Sungsu Park. Oct.2005.               !
    !                                                      !
    !======================================================!
    
    !
    ! Input-Output variables
    !

    implicit none

    integer , intent(in)    :: mix
    integer , intent(in)    :: mkx
    integer , intent(in)    :: iend
    integer , intent(in)    :: ncnst
    real(r8), intent(in)    :: dt                             !  time step in seconds
    real(r8), intent(in)    :: ps0_in(mix,0:mkx)              !  environmental pressure at full sigma levels
    real(r8), intent(in)    :: zs0_in(mix,0:mkx)              !  environmental height at full sigma levels
    real(r8), intent(in)    :: p0_in(mix,mkx)                 !  environmental pressure at half sigma levels
    real(r8), intent(in)    :: z0_in(mix,mkx)                 !  environmental height at half sigma levels
    real(r8), intent(in)    :: dp0_in(mix,mkx)                !  environmental layer pressure thickness
    real(r8), intent(in)    :: u0_in(mix,mkx)                 !  environmental zonal wind
    real(r8), intent(in)    :: v0_in(mix,mkx)                 !  environmental meridional wind
    real(r8), intent(in)    :: qv0_in(mix,mkx)                !  environmental specific humidity
    real(r8), intent(in)    :: ql0_in(mix,mkx)                !  environmental liquid water mixing ratio
    real(r8), intent(in)    :: qi0_in(mix,mkx)                !  environmental ice mixing ratio
    real(r8), intent(in)    :: t0_in(mix,mkx)                 !  environmental temperature
    real(r8), intent(in)    :: s0_in(mix,mkx)                 !  environmental dry static energy
    real(r8), intent(in)    :: tke_in(mix,0:mkx)              !  turbulent kinetic energy
    real(r8), intent(in)    :: pblh_in(mix)                   !  height of PBL
    real(r8), intent(inout) :: cush_inout(mix)                !  convective scale height
    real(r8), intent(out)   :: umf_out(mix,0:mkx)             !  updraft mass flux at top of layer
    real(r8), intent(out)   :: qvten_out(mix,mkx)             !  tendency of specific humidity
    real(r8), intent(out)   :: qlten_out(mix,mkx)             !  tendency of liquid water mixing ratio
    real(r8), intent(out)   :: qiten_out(mix,mkx)             !  tendency of ice mixing ratio
    real(r8), intent(out)   :: sten_out(mix,mkx)              !  tendency of static energy
    real(r8), intent(out)   :: uten_out(mix,mkx)              !  tendency of zonal wind
    real(r8), intent(out)   :: vten_out(mix,mkx)              !  tendency of meridional wind
    real(r8), intent(out)   :: qrten_out(mix,mkx)             !  tendency of rain water mixing ratio
    real(r8), intent(out)   :: qsten_out(mix,mkx)             !  tendency of snow mixing ratio
    real(r8), intent(out)   :: precip_out(mix)                !  vertical integral of (qrten+qsten) or precipitation rate at surface ?
    real(r8), intent(out)   :: slflx_out(mix,0:mkx)           !  updraft liquid static energy flux
    real(r8), intent(out)   :: qtflx_out(mix,0:mkx)           !  updraft total water flux
    real(r8), intent(out)   :: cldfrc_out(mix,mkx)            !  shallow cumulus cloud fraction
    real(r8), intent(out)   :: qlu_out(mix,mkx)               !  updraft liquid water mixing ratio
    real(r8), intent(out)   :: fer_out(mix,mkx)               !  lateral entrainment rate
    real(r8), intent(out)   :: fdr_out(mix,mkx)               !  lateral detrainment rate
    real(r8), intent(out)   :: cinh_out(mix)                  !  convective inhibition (CIN) ([J/kg])
    real(r8), intent(out)   :: cbmf_out(mix)                  !  
    real(r8), intent(out)   :: qc_out(mix,mkx)                !  [kg/kg/s]
    real(r8), intent(out)   :: rliq_out(mix)                  !  [m/s]
    real(r8), intent(out)   :: cnt_out(mix)                   !  
    real(r8), intent(out)   :: cnb_out(mix)                   !  
    real(r8), intent(out)   :: err(mix)
    integer , external      :: qsat 

    !
    ! One-dimensional variables at each grid point
    !

    ! 1. Input variables

    real(r8)    ps0(0:mkx)         !  environmental pressure at full sigma levels
    real(r8)    zs0(0:mkx)         !  environmental height at full sigma levels
    real(r8)    p0(mkx)            !  environmental pressure at half sigma levels
    real(r8)    z0(mkx)            !  environmental height at half sigma levels
    real(r8)    dp0(mkx)           !  environmental layer pressure thickness
    real(r8)    u0(mkx)            !  environmental zonal wind
    real(r8)    v0(mkx)            !  environmental meridional wind
    real(r8)    tke(0:mkx)         !  turbulent kinetic energy
    real(r8)    qv0(mkx)           !  environmental specific humidity
    real(r8)    ql0(mkx)           !  environmental liquid water mixing ratio
    real(r8)    qi0(mkx)           !  environmental ice mixing ratio
    real(r8)    t0(mkx)            !  environmental temperature
    real(r8)    s0(mkx)            !  environmental dry static energy
    real(r8)    pblh               !  height of PBL
    real(r8)    cush               !  convective scale height

    ! 2. Environmental variables directly from the input variables

    real(r8)    qt0(mkx)           !  environmental total water mixing ratio
    real(r8)    thl0(mkx)          !  environmental liquid potential temperature
    real(r8)    thvl0(mkx)         !  environmental liquid virtual potential temperature
    real(r8)    ssqt0(mkx)         !  slope of environmental total water mixing ratio
    real(r8)    ssthl0(mkx)        !  slope of environmental liquid potential temperature
    real(r8)    ssu0(mkx)          !  environmental zonal wind speed vertical gradient
    real(r8)    ssv0(mkx)          !  environmental meridional wind speed vertical gradient
    real(r8)    thv0bot(mkx)       !  environmental virtual potential temperature, bottom of layer
    real(r8)    thv0top(mkx)       !  environmental virtual potential temperature, top of layer
    real(r8)    thvl0bot(mkx)      !  environmental liquid virtual potential temperature, bottom of layer
    real(r8)    thvl0top(mkx)      !  environmental liquid virtual potential temperature, top of layer
    real(r8)    exn0(mkx)          !  exner function at midpoints
    real(r8)    exns0(0:mkx)       !  exner function at interfaces

   ! 3. Variables associated with cumulus convection

    real(r8)    umf(0:mkx)              !  updraft mass flux at top of layer
    real(r8)    emf(0:mkx)              !  updraft mass flux at top of layer
    real(r8)    qvten(mkx)              !  tendency of specific humidity
    real(r8)    qlten(mkx)              !  tendency of liquid water mixing ratio
    real(r8)    qiten(mkx)              !  tendency of ice mixing ratio
    real(r8)    sten(mkx)               !  tendency of static energy
    real(r8)    uten(mkx)               !  tendency of zonal wind
    real(r8)    vten(mkx)               !  tendency of meridional wind
    real(r8)    qrten(mkx)              !  tendency of rain water mixing ratio
    real(r8)    qsten(mkx)              !  tendency of snow mixing ratio
    real(r8)    precip                  !  vertical integral of (qrten+qsten)
    real(r8)    slflx(0:mkx)            !  updraft liquid static energy flux
    real(r8)    qtflx(0:mkx)            !  updraft total water flux
    real(r8)    cldfrc(mkx)             !  shallow cumulus cloud fraction
    real(r8)    qlu(mkx)                !  updraft liquid water mixing ratio
    real(r8)    umflx(0:mkx)            !  flux of zonal momentum due to convection
    real(r8)    vmflx(0:mkx)            !  flux of meridional momentum due to convection
    real(r8)    dwten(mkx)              !  detrained water tendency from cumulus updraft
    real(r8)    diten(mkx)              !  detrained ice   tendency from cumulus updraft 
    real(r8)    fer(mkx)                !  fractional lateral entrainment rate
    real(r8)    fdr(mkx)                !  fractional lateral detrainment rate
    real(r8)    uf(mkx)                 !  zonal wind at end of time step
    real(r8)    vf(mkx)                 !  meridional wind at end of time step
    real(r8)    qc(mkx)                 !  [kg/kg/s] reserved 'qlten+qiten' due to detrained 'cloud water + cloud ice' (without rain-snow contribution)
    real(r8)    qc_l
    real(r8)    qc_i
    real(r8)    rliq                    !  [m/s] vertical integral of qc 
    real(r8)    cnt, cnb

    !----- Variables describing cumulus updraft

    real(r8)    wu(0:mkx)                !  updraft vertical velocity at top of layer
    real(r8)    thlu(0:mkx)              !  updraft liquid potential temperature at top of layer
    real(r8)    uu(0:mkx)                !  updraft zonal wind speed
    real(r8)    vu(0:mkx)                !  updraft meridional wind speed
    real(r8)    qtu(0:mkx)               !  updraft total water at top of layer
    real(r8)    thvu(0:mkx)              !  updraft virtual potential temperature at top of layer
    real(r8)    rei(mkx)                 !  updraft mixing rate of with environment
    
    !----- Other internal variables

    integer     id_check
    logical     id_exit   
    real(r8) :: thlsrc,qtsrc,usrc,vsrc,uplus,vplus
    real(r8) :: plcl,plfc,prel,wrel
    real(r8)    frc_rasn, frc_rasn_remain
    real(r8)    ee2,ud2,wtw
    real(r8)    cldhgt,dpsum
    real(r8)    xc                  !  fraction of environmental air in a neutrally buoyant mixture
    real(r8)    cin                 !  convective inhibition (m2/s2)
    integer     k                   !  release level (half-sigma level just below lcl)
    integer     i
    integer     leff
    integer     kmin                !  layer where 'thvl0' is minimum within the PBL
    integer     klcl                !  layer containing LCL of source air
    integer     kinv                !  inversion layer with PBL top interface as a lower interface
    integer     krel                !  release layer where buoyancy sorting mixing occurs for the first time
    integer     klfc                !  LFC layer of cumulus source air
    integer     kbup                !  top layer in which cloud buoyancy is positive both at the top and bottom interfaces
    integer     kpen                !  highest layer with positive updraft vertical velocity - top layer cumulus can reach
    integer     iteration
    integer     kp1,km1,m
    real(r8) :: scaleh,nu,tkeavg,thvusrc,qt0bot,qse
    real(r8) :: thj,qvj, qlj,qij,tscaleh
    real(r8) :: cinlcl,rbuoy,rdrag,pe,dpe,exne,thvebot,thle,qte,ue,ve
    real(r8) :: cbmf,wexp,wcrit,sigmaw,ufrc,thl0bot
    real(r8) :: thvj,exql,exqi,slten,qtten,PGFc,rle,rkm,rpen,thl0top
    real(r8) :: qt0top,thl0lcl,qt0lcl,thv0lcl,thv0rel,rho0inv,thvubot,thvutop
    real(r8) :: rkfre,thv0j,rho0j,tj,rdqsdt
    real(r8) :: rbeta,ths0,aquad,bquad,cquad,xc1,xc2,excessu,excess0,xsat
    real(r8) :: bogbot,bogtop,delbog,expfac,rhos0j, ppen,ufrcbelow,rainflx
    real(r8) :: snowflx,rmaxfrac,drage,qtdef,qpdef(ncnst),qconbelow,rlwp
    real(r8) :: epsvarw                                                     ! vertical velocity variance at inversion base by meso-scale component.  
    real(r8) :: es(1)                                                       !  saturation vapor pressure
    real(r8) :: qs(1)                                                       !  saturation spec. humidity
    real(r8) :: gam(1)                                                      !  (L/cp)*dqs/dT
    integer  :: status                                                      !  return status of qsat call
    real(r8) :: thvlmin                                                     !  minimum theta_vl in the PBL
!    real,external :: erfcfff
    real(r8)    erfarg,qsat_arg,thvlsrc             
    real(r8)    dpi                                                         ! Full(for internal) or half(for external) interface thickness for tkeavg
    real(r8)    x1,x2,f1,f2,xmid,f_thl,fmid,dx,j                            ! Used for 'thlsrc' calculation from ' p, qtsrc, thvmin' 
    integer     kthvmin,iter_scaleh 
    integer     jj
    real(r8) :: xsrc, xmean, xtop, xbot, xflx(0:mkx)
    real(r8) :: qt_dt, thl_dt

    !----- Variables for implicit CIN computation

    real(r8), dimension(mkx)         :: qv0_s  , ql0_s   , qi0_s   , s0_s    , u0_s    ,           & 
                                        v0_s   , t0_s    , qt0_s   , thl0_s  , thvl0_s ,qvten_s ,  &
                                        qlten_s, qiten_s , qrten_s , qsten_s , sten_s  ,           &
                                        uten_s , vten_s  , cldfrc_s, qlu_s   ,                     &
                                        fer_s  ,  fdr_s  , qc_s 
    real(r8), dimension(0:mkx)       :: umf_s  , slflx_s , qtflx_s
    real(r8)                         :: cush_s , precip_s, cin_s  , rliq_s, cbmf_s, cnt_s, cnb_s
    real(r8)                         :: cin_initial,cin_final,del_cin,ke,alpha,thlj
    real(r8)                         :: cinlcl_i,cinlcl_f,del_cinlcl
    integer                          :: iter

    !----- Variables for temporary storages

    real(r8), dimension(mkx)   :: qv0_o, ql0_o, qi0_o, t0_o, s0_o, u0_o, v0_o
    real(r8), dimension(mkx)   :: qt0_o    , thl0_o   , thvl0_o   ,                         &
                                  qvten_o  , qlten_o  , qiten_o   , qrten_o   , qsten_o ,   &
                                  sten_o   , uten_o   , vten_o    , qlu_o     , cldfrc_o,   &
                                  thv0bot_o, thv0top_o, thvl0bot_o, thvl0top_o,             &
                                  ssthl0_o , ssqt0_o  ,  ssu0_o   , ssv0_o    , qc_o
    real(r8), dimension(0:mkx) :: umf_o    , slflx_o  , qtflx_o
    real(r8), dimension(mkx)   :: dwten_o  , diten_o
    real(r8), dimension(mix)   :: cush_o   , precip_o , rliq_o, cbmf_o, cnt_o, cnb_o
    real(r8), dimension(0:mkx) :: umflx_o  , vmflx_o
    real(r8)                   :: tkeavg_o , thvlmin_o, qtsrc_o  , thvlsrc_o, thlsrc_o ,    &
                                  usrc_o   , vsrc_o   , plcl_o   , plfc_o   ,               &
                                  thv0lcl_o, cinlcl_o 
    integer                    :: kmin_o   , kinv_o   , klcl_o   , klfc_o  

    ! ------------------ !
    ! Define Parameters  !
    ! ------------------ !

    !----- For iterative cin calculation

    integer , parameter              :: iter_cin = 2

    !----- For lateral entrainment

    parameter (rle = 0.1)                !  for critical stopping distance for entrainment
    parameter (rkm = 16.0)               !  for fractional mixing rate
    parameter (rpen = 10.0)              !  for entrainment efficiency
    parameter (rkfre = 1.0)              !  vertical velocity variance as fraction of  tke. 
                                         !  '1' rather than '0.5' might be more reasonable. 
    parameter (rmaxfrac = 0.05)          !  maximum allowable updraft fraction
    parameter (rbuoy = 1.0)              !  for nonhydrostatic pressure effects on updraft
    parameter (rdrag = 1.0)

    !----- Vertical velocity variance at inversion base by meso-scale component.  

    parameter (epsvarw = 5.e-4)          

    !----- For momentum transfer

    parameter (PGFc = 0.7)

    !----- Bulk microphysics controlling parameters

    !      frc_rasn: Fraction of rain (snow) in the water detrained from cumulus updraft
    !                The remaining fraction ( 1 - frc_rasn ) is cloud water ( cloud ice)
    !                0 : all detrained water from Cu is in the form of cloud water
    !                1 : all detrained water from Cu is in the form of rain  
    !      frc_rasn_remain : Fraction of rain/snow that remains in the layer where it is
    !                        formed. The remaining rain/snow fraction is detrained into 
    !                        the release layer
    !                1 : all rain/snow remains in the layer where it was formed.
    !                0 : all rain/snow detrained into the release layer    

    parameter ( frc_rasn = 0.0 )
    parameter ( frc_rasn_remain =1.0 )

    !------------------------!
    !                        !
    ! Start Main Calculation !
    !                        !
    !------------------------!

    !----- Initialize Output variables defined for all grid points

    umf_out(:iend,0:mkx)         = 0.0
    slflx_out(:iend,0:mkx)       = 0.0
    qtflx_out(:iend,0:mkx)       = 0.0
    qvten_out(:iend,:mkx)        = 0.0
    qlten_out(:iend,:mkx)        = 0.0
    qiten_out(:iend,:mkx)        = 0.0
    sten_out(:iend,:mkx)         = 0.0
    uten_out(:iend,:mkx)         = 0.0
    vten_out(:iend,:mkx)         = 0.0
    qrten_out(:iend,:mkx)        = 0.0
    qsten_out(:iend,:mkx)        = 0.0
    precip_out(:iend)            = 0.0
    cldfrc_out(:iend,:mkx)       = 0.0
    qlu_out(:iend,:mkx)          = 0.0
    fer_out(:iend,:mkx)          = 0.0
    fdr_out(:iend,:mkx)          = 0.0
    cinh_out(:iend)              = 0.0
    cbmf_out(:iend)              = 0.0
    qc_out(:iend,:mkx)           = 0.0
    rliq_out(:iend)              = 0.0
    cnt_out(:iend)               = float(mkx)
    cnb_out(:iend)               = 0.0

    !---------------------------------------------------------!
    !                                                         !
    ! Start the big i loop where i is a horozontal grid index !
    !                                                         !
    !---------------------------------------------------------!

    do i = 1, iend                                        ! start of big i loop

      id_exit = .false.
      err(i) = 0.0
      !
      ! Define new input variables at each grid point
      !

      ps0(0:mkx)       = ps0_in(i,0:mkx)
      zs0(0:mkx)       = zs0_in(i,0:mkx)
      p0(:mkx)         = p0_in(i,:mkx)
      z0(:mkx)         = z0_in(i,:mkx)
      dp0(:mkx)        = dp0_in(i,:mkx)
      u0(:mkx)         = u0_in(i,:mkx)
      v0(:mkx)         = v0_in(i,:mkx)
      qv0(:mkx)        = qv0_in(i,:mkx)
      ql0(:mkx)        = ql0_in(i,:mkx)
      qi0(:mkx)        = qi0_in(i,:mkx)
      t0(:mkx)         = t0_in(i,:mkx)
      s0(:mkx)         = s0_in(i,:mkx)
      tke(0:mkx)       = tke_in(i,0:mkx)
      pblh             = pblh_in(i)
      cush             = cush_inout(i)

      !
      ! Compute other basic thermodynamic variables directly from the input variables at each grid point
      ! 

      !----- 1. Compute internal environmental variables
      
      exn0(:mkx)   = (p0(:mkx)/p00)**rovcp
      exns0(0:mkx) = (ps0(0:mkx)/p00)**rovcp
      qt0(:mkx)    = (qv0(:mkx) + ql0(:mkx) + qi0(:mkx))
      thl0(:mkx)   = (t0(:mkx) - lcond*ql0(:mkx)/cp - lsub*qi0(:mkx)/cp)/exn0(:mkx)
      thvl0(:mkx)  = (1. + zvir*qt0(:mkx))*thl0(:mkx)
      
      !----- 2. Compute slopes of environmental variables

      ssthl0       = slope(mkx,thl0,p0) ! Dimension of ssthl0(:mkx) is implicit
      ssqt0        = slope(mkx,qt0 ,p0)
      ssu0         = slope(mkx,u0  ,p0)
      ssv0         = slope(mkx,v0  ,p0)

      !----- 3. Compute "thv0" and "thvl0" at top/bottom interfaces

      do k = 1, mkx

        thl0bot = thl0(k) + ssthl0(k)*(ps0(k-1) - p0(k))
        qt0bot  = qt0(k)  + ssqt0(k) *(ps0(k-1) - p0(k))
        call conden(ps0(k-1),thl0bot,qt0bot,thj,qvj,qlj,qij,qse,id_check,qsat)
        if(id_check.eq.1) then
          id_exit = .true.
          err(i) = 1.1
          go to 333
        end if
        thv0bot(k)  = thj*(1. + zvir*qvj - qlj - qij)
        thvl0bot(k) = thl0bot*(1. + zvir*qt0bot)
          
        thl0top       = thl0(k) + ssthl0(k)*(ps0(k) - p0(k))
        qt0top        =  qt0(k) + ssqt0(k) *(ps0(k) - p0(k))
        call conden(ps0(k),thl0top,qt0top,thj,qvj,qlj,qij,qse,id_check,qsat)
        if(id_check.eq.1) then
          id_exit = .true.
          err(i) = 1.2
          go to 333
        end if 
        thv0top(k)  = thj*(1. + zvir*qvj - qlj - qij)
        thvl0top(k) = thl0top*(1. + zvir*qt0top)

      end do

      !
      ! Save input and related thermodynamic variables ( not associated with cumulus convection ) 
      ! for use at "iter_cin=2" when "del_cin >= 0"
      !

      qv0_o(:mkx)          = qv0(:mkx)
      ql0_o(:mkx)          = ql0(:mkx)
      qi0_o(:mkx)          = qi0(:mkx)
      t0_o(:mkx)           = t0(:mkx)
      s0_o(:mkx)           = s0(:mkx)
      u0_o(:mkx)           = u0(:mkx)
      v0_o(:mkx)           = v0(:mkx)
      qt0_o(:mkx)          = qt0(:mkx)
      thl0_o(:mkx)         = thl0(:mkx)
      thvl0_o(:mkx)        = thvl0(:mkx)
      ssthl0_o(:mkx)       = ssthl0(:mkx)
      ssqt0_o(:mkx)        = ssqt0(:mkx)
      thv0bot_o(:mkx)      = thv0bot(:mkx)
      thv0top_o(:mkx)      = thv0top(:mkx)
      thvl0bot_o(:mkx)     = thvl0bot(:mkx)
      thvl0top_o(:mkx)     = thvl0top(:mkx)
      ssu0_o(:mkx)         = ssu0(:mkx) 
      ssv0_o(:mkx)         = ssv0(:mkx) 

      !
      ! Initialize output variables after cumulus convection at each grid point
      !

      umf(0:mkx)          = 0.0
      emf(0:mkx)          = 0.0
      slflx(0:mkx)        = 0.0
      qtflx(0:mkx)        = 0.0
      umflx(0:mkx)        = 0.0
      vmflx(0:mkx)        = 0.0
      qvten(:mkx)         = 0.0
      qlten(:mkx)         = 0.0
      qiten(:mkx)         = 0.0
      sten(:mkx)          = 0.0
      uten(:mkx)          = 0.0
      vten(:mkx)          = 0.0
      qrten(:mkx)         = 0.0
      qsten(:mkx)         = 0.0
      dwten(:mkx)         = 0.0
      diten(:mkx)         = 0.0
      precip              = 0.0
      cldfrc(:mkx)        = 0.0
      qlu(:mkx)           = 0.0
      fer(:mkx)           = 0.0
      fdr(:mkx)           = 0.0
      cin                 = 0.0
      cbmf                = 0.0
      qc(:mkx)            = 0.0
      rliq                = 0.0
      cnt                 = float(mkx)
      cnb                 = 0.0

    !-----------------------------------------------! 
    ! Below 'iter' loop is for implicit CIN closure !
    !-----------------------------------------------!

    do iter = 1, iter_cin

       !---Cumulus scale height
       !   The issue is that if I should locate below two lines 
       !   within or out the cin_iter loop. 

       tscaleh = cush                        
       cush    = -1.

       !----- Find PBL top height interface, 'kinv-1' where 'kinv' is the layer index with PBLH in it.
       !      When PBLH is exactly at interface, 'kinv' is the layer having PBLH as a lower interface.
       !      I set lower limit of 'kinv' 2 for consistency with the other parts of the code later.

       do k = mkx - 1, 1, -1 
         if((pblh + 0.1 - zs0(k))*(pblh + 0.1 - zs0(k+1)) .lt. 0.) then
           kinv = k + 1 
           go to 15
         endif 
       end do
       kinv = 1
15     continue    
       kinv = max(2,kinv)

       !----- Find PBL averaged tke ('tkeavg') and minimum 'thvl' ('thvlmin') in the PBL

       dpsum    = 0.
       tkeavg   = 0.
       thvlmin  = 1000.
       kmin     = 1
       do k = 0, kinv - 1    
         if(k .eq. 0) then
          dpi   = ps0(0) - p0(1)
         elseif(k .eq. (kinv-1)) then 
          dpi   = p0(kinv-1) - ps0(kinv-1)
         else
          dpi   = p0(k) - p0(k+1)
         endif 
         dpsum  = dpsum  + dpi  
         tkeavg = tkeavg + dpi*tke(k) 
!        if( k.ne.0 .and. thvl0(k).lt.thvlmin ) then
!         thvlmin = thvl0(k)
!         kmin = k
!        end if         
! Below commented-out '1~5' lines are old source code set, which should be used as a set.
       if( k.ne.0 ) thvlmin = min(thvlmin,min(thvl0bot(k),thvl0top(k))) !1
       end do
       tkeavg  = tkeavg/dpsum
       thvlmin = min(thvlmin,thvl0bot(kinv)) !2

       !----- Find characteristics of cumulus source air: qtsrc,thlsrc,usrc,vsrc
       !      Note that 'thlsrc' was con-cocked using 'thvlsrc' and 'qtsrc', and
       !      momentum source is defined using the values at inversion interface 

       qtsrc   = qt0(1)  !3                   
       thvlsrc = thvlmin !4
       thlsrc  = thvlsrc / ( 1. + zvir * qtsrc ) !5 

!      qtsrc   = qt0(kmin)
!      thvlsrc = thvl0(kmin) 
!      thlsrc  = thl0(kmin)
       usrc    = u0(kinv-1) + ssu0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )             
       vsrc    = v0(kinv-1) + ssv0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )             

       !----- Find LCL of the source air and a layer index containing LCL (klcl)
       !      When LCL is exactly at the interface, 'klcl' is the layer index  
       !      having 'plcl' as the lower interface.
       !      If LCL is located within the lowest model layer or the top model
       !      layer ( mkx ), no convective adjustment is performed. But I can
       !      choose different exit conditions, e.g., in combination with 'kinv'. 

       plcl = qsinvert(qtsrc,thlsrc,qsat)
       do k = 0, mkx
         if( ps0(k) .lt. plcl ) then
           klcl = k
           go to 25
         endif           
       end do
       klcl = mkx
25     continue
       ! Exit condition by using 'klcl' and potentially 'kinv' 
       if(klcl.le.1.or.klcl.eq.mkx) then          
          id_exit = .true.
          err(i) = 2 ! lcl in lowest or top model layer 
          go to 333
       endif

       !----- Calculate 'thv0lcl' which is repeatedly used solely in cin calculation.

       thl0lcl = thl0(klcl) + ssthl0(klcl) * ( plcl - p0(klcl) )
       qt0lcl  = qt0(klcl)  + ssqt0(klcl)  * ( plcl - p0(klcl) )
       call conden(plcl,thl0lcl,qt0lcl,thj,qvj,qlj,qij,qse,id_check,qsat)
       if(id_check.eq.1) then
          id_exit = .true.
          err(i) = 1.3 ! conden fails to converge
          go to 333
       end if
       thv0lcl = thj * (1. + zvir * qvj - qlj - qij )

       !
       !----- Determine the convective inhibition (CIN)
       !

       !      CIN is computed from the PBL top interface to the LFC using 
       !      standard reconstruction technique of environmental profiles.
       !      1. If LCL is lower  than PBL height, CIN is obtained by summing non-entraining cloud buoyancy 
       !         from the PBL top interface ( NOT from LCL ) to LFC. And 'cinlcl' is set to be zero.
       !      2. If LCL is higher than PBL height, CIN is similarily calculated from PBL top interface to 
       !         to LFC, but with non-zero ( usually, positive due to the definition of thvlsrc ) cinlcl. 
       !      3. If LCL height is higher than PBL interface ( 'pLCL <= ps0(kinv-1)' ) :
       !            'cinlcl' is calculated considering both positive and negative CIN (using 'single_cin' ),
       !            and after calculating 'cinlcl', the remaining CIN are calculated upto LFC counting only
       !            positive value. If either CIN or cinlcl is negative, they are set to be zero.    
       !      4. If LCL height is lower  than PBL interface ( 'pLCL >  ps0(kinv-1)' ) :
       !            'cinlcl' is set to be zero, and CIN is calculated from ps0(kinv-1) upward upto LFC found.
       !            If either CIN or cinlcl is negative, they are set to be zero.

        cin    = 0.
        cinlcl = 0.
        plfc   = 0.
        klfc   = mkx

        !----- Case 1. LCL height is higher than PBL interface ( 'pLCL <= ps0(kinv-1)' )

        if( klcl .ge. kinv ) then

           do k = kinv, mkx - 1

             if( k .lt. klcl ) then

                thvubot = thvlsrc
                thvutop = thvlsrc  
                cin     = cin + single_cin(ps0(k-1),thv0bot(k),ps0(k),thv0top(k),thvubot,thvutop)
             elseif( k .eq. klcl ) then

                !----- Bottom to LCL

                thvubot = thvlsrc
                thvutop = thvlsrc
                cin     = cin + single_cin(ps0(k-1),thv0bot(k),plcl,thv0lcl,thvubot,thvutop)
                cinlcl  = max(cin,0.)
                cin     = cinlcl

                !----- LCL to Top

                thvubot = thvlsrc
                call conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check,qsat)
                if(id_check.eq.1) then
                  id_exit = .true.
                  err(i) = 1.4 ! conden fails to converge
                  go to 333
                end if

                thvutop = thj*(1. + zvir*qvj - qlj - qij )
                call getbuoy(plcl,thv0lcl,ps0(k),thv0top(k),thvubot,thvutop,plfc,cin)
                if(plfc .gt. 0.) then 
                  klfc = k 
                  go to 35
                end if

             else

                thvubot = thvutop
                call conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check,qsat)
                if(id_check.eq.1) then
                   id_exit = .true.
                   err(i) = 1.5 ! conden fails to converge
                   go to 333
                end if
                thvutop = thj*(1. + zvir*qvj - qlj - qij )

                call getbuoy(ps0(k-1),thv0bot(k),ps0(k),thv0top(k),thvubot,thvutop,plfc,cin)
                if(plfc .gt. 0.) then 
                  klfc = k
                  go to 35
                end if 
 
            endif

          end do        

       !----- Case 2. LCL height is lower than PBL interface ( 'pLCL > ps0(kinv-1)' )

       else

          cinlcl = 0. 

          do k = kinv, mkx - 1

             call conden(ps0(k-1),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check,qsat)
             if(id_check.eq.1) then
                id_exit = .true.
                err(i) = 1.6 ! conden fails to converge
                go to 333
             end if
             thvubot = thj*(1. + zvir*qvj - qlj - qij )

             call conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check,qsat)
             if(id_check.eq.1) then
                id_exit = .true.
                err(i) = 1.7 ! conden fails to converge
                go to 333
             end if
             thvutop = thj*(1. + zvir*qvj - qlj - qij )

             call getbuoy(ps0(k-1),thv0bot(k),ps0(k),thv0top(k),thvubot,thvutop,plfc,cin)
             if(plfc .gt. 0.) then 
                klfc = k
                go to 35
             end if 

         end do
  
       endif                                       ! End of CIN case selection

 35    continue
       cin = max(0.,cin)
       ! Exit conditions using 'klfc', 'klcl', 'kinv'
       if( klfc .ge. mkx ) then
         klfc = mkx
         id_exit = .true.
          err(i) =  3 ! no LFC
         go to 333
       endif

       if(iter .eq. 1) then 

          cin_initial = cin
          cinlcl_i    = cinlcl
          ke          = rbuoy / ( rkfre * tkeavg + epsvarw ) 

          kmin_o      = kmin     
          kinv_o      = kinv     
          klcl_o      = klcl     
          klfc_o      = klfc    
          plcl_o      = plcl    
          plfc_o      = plfc     
          tkeavg_o    = tkeavg   
          thvlmin_o   = thvlmin
          qtsrc_o     = qtsrc    
          thvlsrc_o   = thvlsrc  
          thlsrc_o    = thlsrc
          usrc_o      = usrc     
          vsrc_o      = vsrc     
          thv0lcl_o   = thv0lcl  

       endif   

       if(iter .ne. 1) then

          cin_final  = cin
          cinlcl_f   = cinlcl
          del_cin    = cin_final - cin_initial
          del_cinlcl = cinlcl_f - cinlcl_i

          if(del_cin .gt. 0.) then 

             !----- Calculate final 'cin' and 'cinlcl'

             alpha   = compute_alpha( del_cin, ke )
             cin     = cin_initial + alpha * del_cin
             cinlcl  = cinlcl_i + alpha * del_cinlcl

             !----- Restore original values from the previous iter_cin step (1) to compute 
             !      correct tendencies for (n+1) time step using final 'cin' and 'cinlcl'

             kmin      = kmin_o    
             kinv      = kinv_o     
             klcl      = klcl_o     
             klfc      = klfc_o    
             plcl      = plcl_o    
             plfc      = plfc_o     
             tkeavg    = tkeavg_o   
             thvlmin   = thvlmin_o
             qtsrc     = qtsrc_o    
             thvlsrc   = thvlsrc_o  
             thlsrc    = thlsrc_o
             usrc      = usrc_o     
             vsrc      = vsrc_o     
             thv0lcl   = thv0lcl_o  

             qv0(:mkx)            = qv0_o(:mkx)
             ql0(:mkx)            = ql0_o(:mkx)
             qi0(:mkx)            = qi0_o(:mkx)
             t0(:mkx)             = t0_o(:mkx)
             s0(:mkx)             = s0_o(:mkx)
             u0(:mkx)             = u0_o(:mkx)
             v0(:mkx)             = v0_o(:mkx)
             qt0(:mkx)            = qt0_o(:mkx)
             thl0(:mkx)           = thl0_o(:mkx)
             thvl0(:mkx)          = thvl0_o(:mkx)
             ssthl0(:mkx)         = ssthl0_o(:mkx)
             ssqt0(:mkx)          = ssqt0_o(:mkx)
             thv0bot(:mkx)        = thv0bot_o(:mkx)
             thv0top(:mkx)        = thv0top_o(:mkx)
             thvl0bot(:mkx)       = thvl0bot_o(:mkx)
             thvl0top(:mkx)       = thvl0top_o(:mkx)
             ssu0(:mkx)           = ssu0_o(:mkx) 
             ssv0(:mkx)           = ssv0_o(:mkx) 

             ! Initialize all fluxes, tendencies, and other variables resulting from cumulus convection

             umf(0:mkx)          = 0.0
             emf(0:mkx)          = 0.0
             slflx(0:mkx)        = 0.0
             qtflx(0:mkx)        = 0.0
             umflx(0:mkx)        = 0.0
             vmflx(0:mkx)        = 0.0
             qvten(:mkx)         = 0.0
             qlten(:mkx)         = 0.0
             qiten(:mkx)         = 0.0
             sten(:mkx)          = 0.0
             uten(:mkx)          = 0.0
             vten(:mkx)          = 0.0
             qrten(:mkx)         = 0.0
             qsten(:mkx)         = 0.0
             dwten(:mkx)         = 0.0
             diten(:mkx)         = 0.0
             precip              = 0.0
             cldfrc(:mkx)        = 0.0
             qlu(:mkx)           = 0.0
             fer(:mkx)           = 0.0
             fdr(:mkx)           = 0.0
             qc(:mkx)            = 0.0
             rliq                = 0.0
             cbmf                = 0.0
             cnt                 = float(mkx)
             cnb                 = 0.0

          else
            
             ! Restore original output values of "iter_cin=1" and exit

             umf_out(i,0:mkx)         = umf_s(0:mkx)
             qvten_out(i,:mkx)        = qvten_s(:mkx)
             qlten_out(i,:mkx)        = qlten_s(:mkx)  
             qiten_out(i,:mkx)        = qiten_s(:mkx)
             sten_out(i,:mkx)         = sten_s(:mkx)
             uten_out(i,:mkx)         = uten_s(:mkx)  
             vten_out(i,:mkx)         = vten_s(:mkx)
             qrten_out(i,:mkx)        = qrten_s(:mkx)
             qsten_out(i,:mkx)        = qsten_s(:mkx)  
             precip_out(i)            = precip_s
             cush_inout(i)            = cush_s
             cldfrc_out(i,:mkx)       = cldfrc_s(:mkx)  
             slflx_out(i,0:mkx)       = slflx_s(0:mkx)  
             qtflx_out(i,0:mkx)       = qtflx_s(0:mkx)  
             qlu_out(i,:mkx)          = qlu_s(:mkx)  
             fer_out(i,:mkx)          = fer_s(:mkx)  
             fdr_out(i,:mkx)          = fdr_s(:mkx)  
             cinh_out(i)              = cin_s
             cbmf_out(i)              = cbmf_s
             qc_out(i,:mkx)           = qc_s(:mkx)  
             rliq_out(i)              = rliq_s
             cnt_out(i)               = cnt_s
             cnb_out(i)               = cnb_s

             id_exit = .false.
             go to 333

          endif

       endif    

       !----- Define release layer, 'krel'.
       !      When LCL is below inversion base, 'krel' is defined as 'kinv' while
       !      when LCL is above inversion base, 'krel' is defined as 'klcl'.
       !      This ensures only PBL scheme works in the PBL.
       !      'prel' is the lowest level from which buoyancy sorting occurs.  
       !      Taking-out the first 5 lines and last line (endif) in the below block
       !      allows convection starts from 'klcl' always instead of 'kinv' fully consistently.

       if( klcl.lt.kinv ) then
          krel = kinv
          prel = ps0(krel-1)
          thv0rel = thv0bot(krel) 
       else
          krel = klcl
          prel = plcl 
          thv0rel = thv0lcl
       endif  

       !----- Calculate cloud base mass flux (cbmf), fractional area (ufrc), and 
       !      vertical velocity (wexp) of cumulus updraft at inversion interface 
       !----- Also, calculate vertical velocity at 'krel', more specifically, at
       !      LCL (wrel). When LCL is below PBLH, cinlcl = 0 and so wrel = wexp.

       wcrit   = sqrt( 2. * cin * rbuoy )
       sigmaw  = sqrt( rkfre * tkeavg + epsvarw )
       rho0inv = ps0(kinv-1) / ( rgas * thv0bot(kinv) * exns0(kinv-1) )
       cbmf    = rho0inv * sigmaw / 2.5066 * exp(-( ( wcrit / sigmaw )**2 ) / 2)
       ! Set limits to 'cbmf'.
       cbmf    = max( 1.e-10, min( cbmf, dpsum/ggr/dt - 1.e-10 ) ) ! 'dpsum' is the thickness of PBL
                                                                 !  Off-set is necessary to prevent model crash in 'fluxbelowinv'  
       erfarg  = wcrit / 1.4142 / sigmaw
       ! Set upper limit to 'ufrc'            
!       ufrc    = min( rmaxfrac, 0.5*erfcfff(erfarg) )
       ufrc    = min( rmaxfrac, 0.5*erf(erfarg) )
       if(ufrc .gt. 0.001) then
          wexp = cbmf / rho0inv / ufrc 
       else
          id_exit = .true.
          err(i) = 4.0 ! ufrc (fractional area) less than 0.001
          go to 333
       endif

       wtw = wexp * wexp - 2. * cinlcl * rbuoy
       if(wtw .le. 0.) then
         id_exit = .true.
          err(i) = 5.0 ! too much CIN
         go to 333
       endif
       wrel = sqrt(wtw)

       !----- Define updraft properties at the level where buoyancy sorting starts to be
       !      happening, i.e., by definition, at 'prel' level  within the release layer.
       !      Because no lateral entrainment occurs upto 'prel', conservative scalars of 
       !      cumulus updraft at release level is same as those of source air.   However
       !      horizontal momentums of source air are modified by horizontal PGF forcings 
       !      from inversion interface to 'prel'.For this case, we should add additional
       !      horizontal momentum from inversion interface to 'prel' as in below to usrc.
       !      Should I initialize 'thlu','qtu','thvu' at the beginning of the program ?

       umf(krel-1)  = cbmf
       wu(krel-1)   = wrel
       thlu(krel-1) = thlsrc
       qtu(krel-1)  = qtsrc

       call conden(prel,thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check,qsat)
       if(id_check.eq.1) then
          id_exit = .true.
          err(i) = 1.8 ! too much CIN
          go to 333
       end if
       thvu(krel-1) = thj*(1. + zvir*qvj - qlj - qij )       

       uplus = 0.
       vplus = 0.
       if( krel.eq.kinv ) then
           uplus = PGFc * ssu0(kinv) * ( prel - ps0(kinv-1) )
           vplus = PGFc * ssv0(kinv) * ( prel - ps0(kinv-1) )
       else
         do k = kinv, max(krel-1,kinv)
           uplus = uplus + PGFc * ssu0(k) * ( ps0(k) - ps0(k-1) )
           vplus = vplus + PGFc * ssv0(k) * ( ps0(k) - ps0(k-1) )
         end do
           uplus = uplus + PGFc * ssu0(krel) * ( prel - ps0(krel-1) )
           vplus = vplus + PGFc * ssv0(krel) * ( prel - ps0(krel-1) )
       end if

       uu(krel-1) = usrc + uplus
       vu(krel-1) = vsrc + vplus      

       !----- Define the environmental properties at the level ( 'pe', normally, layer midpoint ) 
       !      where mixing occurs. In the first 'krel' layer, however, note that 'pe' is defined
       !      slightly differently because LCL is regarded as lower interface for mixing purpose.

       pe      = 0.5 * ( prel + ps0(krel) )
       dpe     = prel - ps0(krel)
       exne    = exnf(pe)
       thvebot = thv0rel
       thle    = thl0(krel) + ssthl0(krel) * ( pe - p0(krel) )
       qte     = qt0(krel)  + ssqt0(krel)  * ( pe - p0(krel) )
       ue      = u0(krel)   + ssu0(krel)   * ( pe - p0(krel) )
       ve      = v0(krel)   + ssv0(krel)   * ( pe - p0(krel) )

       !-------------------------! 
       ! Buoyancy-Sorting Mixing !
       !-------------------------!

       !----- In order to complete buoyancy-sorting mixing at layer midpoint, and thus to
       !      calculate 'updraft mass flux, updraft w velocity, conservative scalars'  at 
       !      the layer upper interface, we need the following 3 information.  
       !
       !      1. Pressure where mixing occurs ('pe'), and temperature at 'pe' which is
       !         necessary to calculate various thermodynamic coefficients at pe. This
       !         temperature is obtained by undiluted cumulus properties lifted to pe. 
       !      2. Undiluted updraft properties at pe - conservative scalar and vertical
       !         velocity -which are assumed to be the same as the properties at lower 
       !         interface only for calculation of fer(k) and fdr(k).  Actual vertical
       !         variations of cumulus conservative scalars and vertical velocity   at
       !         the top interface from fer(k) and fer(k) are calculated subsequently.
       !      3. Environmental properties at pe.

       !----- Define cumulus scale height for fer(k), fdr(k) calculation
 
       scaleh = tscaleh
       if(tscaleh .lt. 0.0) scaleh = 1000. 

       kbup   = krel
       kpen   = krel

       do iter_scaleh = 1, 3

       kbup    = krel
       kpen    = krel
       wtw     = wexp * wexp - 2. * cinlcl * rbuoy

       pe      = 0.5 * ( prel + ps0(krel) )
       dpe     = prel - ps0(krel)
       exne    = exnf(pe)
       thvebot = thv0rel
       thle    = thl0(krel) + ssthl0(krel) * ( pe - p0(krel) )
       qte     = qt0(krel)  + ssqt0(krel)  * ( pe - p0(krel) )
       ue      = u0(krel)   + ssu0(krel)   * ( pe - p0(krel) )
       ve      = v0(krel)   + ssv0(krel)   * ( pe - p0(krel) )

       do k = krel, mkx - 1

          km1 = k - 1

          !----- Calculate saturation excess.
          !      Note that in order to calculate saturation excess, we should use 
          !      liquid water temperature instead of temperature  as the argument
          !      of "qsat". However, normal argument of "qsat" is temperature.

          call conden(pe,thle,qte,thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             id_exit = .true.
             err(i) = 1.9 ! conden doesn't converge
             go to 333
          end if
          thv0j    = thj * ( 1. + zvir*qvj - qlj - qij )
          rho0j    = pe / ( rgas * thv0j * exne )
          qsat_arg = thle*exne     
          status   = qsat(qsat_arg,pe,es(1),qs(1),gam(1),1)
          excess0  = qte - qs(1)

          call conden(pe,thlu(km1),qtu(km1),thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             id_exit = .true.
             err(i) = 1.10 ! conden doesn't converge
             go to 333
          end if
          thvj     = thj * ( 1. + zvir * qvj - qlj - qij )
          tj       = thj * exne ! Necessary for calculating thermodynamic mixing coefficients below
          qsat_arg = thlu(km1)*exne
          status   = qsat(qsat_arg,pe,es(1),qs(1),gam(1),1)
          excessu  = qtu(km1) - qs(1)

          ! ----- Calculate fractional entrainment rate
          !       Only saturated updrafts with "positive buoyancy" or "negative buoyancy  + 
          !       enough vertical velocity" are kept into the updraft in the below program. 
          !       If the core updraft is unsaturated, we can set xc = 0 and let the cumulus
          !       convection still works or we may exit.

          if( excessu.lt.0 ) then
                       
            xc = 0.
            write(6,*) 'Dry Core Updraft Warning in UW-shallow scheme'
!           id_exit = .true.
!           go to 333 
            err(i) = 6.0 ! dry core updraft

          elseif( excessu.ge.0.and.excess0.ge.0 ) then
                      
            aquad =  wu(km1)**2
            bquad = -2.*rbuoy*ggr*rle*scaleh*(thvj - thv0j)/thv0j - 2.*wu(km1)**2
            cquad =  2.*rbuoy*ggr*rle*scaleh*(thvj - thv0j)/thv0j +    wu(km1)**2

            if( (-2.*aquad-bquad).ge.0. ) then
              xc = 1.
            else
              xc = min(1.,max(0.,-bquad/aquad-1.)) 
            endif
        
          elseif ( excessu.ge.0.and.excess0.lt.0 ) then

            xsat   = excessu / ( excessu - excess0 )

            rdqsdt = ep2 * lcond * qse / rgas / tj**2
            rbeta  = ( 1. + ( 1. + zvir ) * tj * rdqsdt ) / ( 1. + rdqsdt * lcond / cp )
            ths0   = min(thv0j,thvj + rbeta*(thle - thlu(km1)) + (rbeta*lcond/cp/exne - thj)*(qte - qtu(km1)))

            aquad =  wu(km1)**2
            bquad = -2.*rbuoy*ggr*rle*scaleh*(thvj - ths0) /thv0j - 2.*wu(km1)**2
            cquad =  2.*rbuoy*ggr*rle*scaleh*(thvj - thv0j)/thv0j +    wu(km1)**2

            if((bquad**2-4.*aquad*cquad).gt.0.) then
              xc = min(1.,max(0.,min(xsat,(-bquad-sqrt(bquad**2-4.*aquad*cquad))/(2.*aquad))))
            else
              xc = xsat
            endif

          endif

          !----- Calculate lateral entrainment-detrainment rate

          ee2    = xc**2
          ud2    = 1. - 2.*xc + xc**2
          rei(k) = rkm/scaleh/ggr/rho0j
          ! In the below, I should separately use 'dp0(k)' and 'dpe'
          if( xc.gt.0.5 ) rei(k) = min(rei(k),0.9*log(dp0(k)/ggr/dt/umf(km1) + 1.)/dpe/(2.*xc-1.)) 
          fer(k) = rei(k) * ee2
          fdr(k) = rei(k) * ud2
          ! ------------------------------------------------------------------------------ !
          ! Iteration Start due to 'maxufrc' constraint [ ****************************** ] ! 
          ! ------------------------------------------------------------------------------ !

          !----- Calculate cumulus updraft mass flux and entrainment mass flux
          !      Note that entrainment mass flux is set to zero  at interfaces 
          !      below 'kbup' interface ( emf=0 at '=< kbup-1' )
          !      In order to calculate 'umf(k)', we use 'umf(km1)' of lower interface. 
          umf(k) = umf(km1) * exp( dpe * ( fer(k) - fdr(k) ) )
          emf(k) = 0.0    

          !----- Now thermodynamics of the dilute plume

          if( fer(k)*dpe .lt. 1.e-4 ) then

            thlu(k) = thlu(km1) + ( thle + ssthl0(k) * dpe / 2. - thlu(km1) ) * fer(k) * dpe
            qtu(k)  =  qtu(km1) + ( qte  + ssqt0(k) * dpe / 2.  -  qtu(km1) ) * fer(k) * dpe
            uu(k)   =   uu(km1) + ( ue   + ssu0(k) * dpe / 2.   -   uu(km1) ) * fer(k) * dpe - PGFc * ssu0(k) * dpe
            vu(k)   =   vu(km1) + ( ve   + ssv0(k) * dpe / 2.   -   vu(km1) ) * fer(k) * dpe - PGFc * ssv0(k) * dpe

          else

            thlu(k) = ( thle + ssthl0(k) / fer(k) - ssthl0(k) * dpe / 2. ) -          &
                      ( thle + ssthl0(k) * dpe / 2. - thlu(km1) + ssthl0(k) / fer(k) ) * exp(-fer(k) * dpe)
            qtu(k)  = ( qte  + ssqt0(k) / fer(k) - ssqt0(k) * dpe / 2. ) -            &  
                      ( qte  + ssqt0(k) * dpe / 2. - qtu(km1) + ssqt0(k) / fer(k) ) * exp(-fer(k) * dpe)
            uu(k) =   ( ue + ( 1 - PGFc ) * ssu0(k) / fer(k) - ssu0(k) * dpe / 2. ) - &
                      ( ue + ssu0(k) * dpe / 2. - uu(km1) + ( 1 - PGFc ) * ssu0(k) / fer(k) ) * exp(-fer(k) * dpe)
            vu(k) =   ( ve + ( 1 - PGFc ) * ssv0(k) / fer(k) - ssv0(k) * dpe / 2. ) - &
                      ( ve + ssv0(k) * dpe / 2. - vu(km1) + ( 1 - PGFc ) * ssv0(k) / fer(k) ) * exp(-fer(k) * dpe)

          end if

          !----- Detrainment of cloud water and ice frum the upper interface
          !      I should significantly modify below microphysics later.

          call conden(ps0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             id_exit = .true.
             err(i) = 1.11 ! conden doesn't converge
             go to 333
          end if
          if( (qlj + qij) .gt. 5.e-4 ) then
             exql    = ( ( qlj + qij ) - 5.e-4 ) * qlj / ( qlj + qij )
             exqi    = ( ( qlj + qij ) - 5.e-4 ) * qij / ( qlj + qij )
             qtu(k)  = qtu(k) - exql - exqi
             thlu(k) = thlu(k) + ( lcond / cp / exns0(k) ) * exql + ( lsub / cp / exns0(k) ) * exqi 
             ! Detrainment of cloud water and ice into the environment from updraft. 
             ! Correct unit of 'dwten, diten' is [ kg/kg/s ]  and restoration to the 
             ! correct unit by multiplying 'umf(k)*g/dp0(k)' will be performed later 
             ! after updating 'umf' using upper 'maxufrc' constraint near the end of
             ! this updraft loop. 
             dwten(k) = exql   
             diten(k) = exqi
          endif
          ! Update thvu(k) after detraining rain and snow.
          call conden(ps0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             id_exit = .true.
             err(i) = 1.12 ! conden doesn't converge
             go to 333
          end if  
          thvu(k) = thj * ( 1. + zvir * qvj - qlj - qij )

          !----- Calculate vertical velocity at the upper interface 
          !      In order to calculate 'wtw' at the upper interface, we use 'wtw' 
          !      at the lower interface. 
 
          bogbot = rbuoy * ( thvu(km1) / thvebot  - 1. )
          bogtop = rbuoy * ( thvu(k) / thv0top(k) - 1. )

          ! ----- In final, 'kbup' is the upper most layer in which cloud buoyancy 
          !       is positive both at the lower and upper interface.
          if( bogbot .gt. 0. .and. bogtop .gt. 0. ) then 
             kbup = k
          end if

          delbog = bogtop - bogbot
          drage  = fer(k) * ( 1. + rdrag )
          expfac = exp( -2.*drage*dpe )

          if( drage*dpe .gt. 1.e-3 ) then
             wtw = wtw*expfac + (delbog + (1.-expfac)*(bogbot + delbog/(-2.*drage*dpe)))/(rho0j*drage)
          else
             wtw = wtw + dpe * ( bogbot + bogtop ) / rho0j
          endif

          ! ----- 'kpen' is the layer in which updraft vertical velocity becomes zero.
          if(wtw .le. 0.) then
             kpen = k
             go to 45
          end if

          wu(k) = sqrt(wtw)
          if(wu(k) .gt. 100.) then
             go to 333
          endif

          ! ---------------------------------------------------------------------------- !
          ! Iteration End due to 'maxufrc' constraint [ ****************************** ] ! 
          ! ---------------------------------------------------------------------------- !

          !----- Calculate updraft fractional area and set upper limit of 'ufrc' to 'rmaxfrac'. 
          !      It is not clear if this limiting procedure is reasonable or not.   Basically,
          !      this code keeps the consistency among 'ufrc','umf','wu(or wtw)'. If 'ufrc' is
          !      limited, either 'umf' or 'wu' also should be changed. Although both 'umf' and
          !      'wu(wtw)' at the current upper interface are used for updating 'umf' and 'wu'
          !      at the next upper interface, 'umf' is a passive variable not influencing  the
          !      buoyancy mixing process at the next level in contrast to 'wtw'. This may be a
          !      reason why 'umf' is adjusted instead of 'wtw'. In turn, we updated 'fdr' here
          !      instead of 'fer', which guarantees that   all updated thermodynamic variables
          !      at the upper interface before this 'maxufrc' constraint are internally consis
          !      tent without going through the interation loop indicated above.  If we update
          !      'fer' however, then we should go through the above iteration loop.
            
          rhos0j = ps0(k) / ( rgas * 0.5 * ( thv0bot(k+1) + thv0top(k) ) * exns0(k) )
          ufrc   = umf(k) / ( rhos0j * wu(k) )
          if(ufrc .gt. rmaxfrac) then
             ufrc   = rmaxfrac
             umf(k) = rmaxfrac * rhos0j * wu(k)
             fdr(k) = fer(k) - log( rmaxfrac * rhos0j * wu(k) / umf(km1) ) / dpe
          endif
          
          !----- Update environmental properties for next level

          pe      = p0(k+1)
          dpe     = dp0(k+1)
          exne    = exn0(k+1)
          thvebot = thv0bot(k+1)
          thle    = thl0(k+1)
          qte     = qt0(k+1)
          ue      = u0(k+1)
          ve      = v0(k+1) 

       end do       

       !----- End of Updraft Loop

45     cush   = z0(kpen)
       scaleh = cush 

       end do      

       !------ End of iter_scaleh

            
       !-----  Filtering of unerasonable cumulus adjustment here.  Should be careful.
       !       Various ways possible using 'klcl','kinv','krel','klfc','kbup','kpen'.

!      cldhgt = plcl - p0(kpen)         ! Original approach
       cldhgt = p0(kpen)                ! First alternative
!      cldhgt = z0(kpen) - zs0(krel)    ! Second alternative 

       if(kpen .lt. krel + 2 .or. cldhgt .le. 50000 .or. kbup .le. krel + 1) then
!      if(kpen .lt. krel + 2 .or. cldhgt .ge. 40000 .or. kbup .le. krel + 1) then
!      if(kpen .lt. krel + 2 .or. cldhgt .ge. 4500  .or. kbup .le. krel + 1) then
!      if(kpen .lt. krel + 2 .or. cldhgt .ge. 10000 .or. kbup .le. krel + 1) then
!      if(kpen .lt. krel + 1 .or. kbup .le. krel + 1) then
!      if(kpen .eq. krel) then
!      if(cldhgt .ge. 40000) then
          id_exit = .true.
          err(i) = 7.0 ! conden doesn't converge
          go to 333
       end if

 
       !----- Calculate downward penetrative entrainment mass flux, emf(k) < 0.
       !      emf(k) is defined at upper interface from layer 'kbup' to 'kpen-1' 
       !      Similar to MinScCu parameterization, emf(k) is from a similarity theory.

       umf(kpen:mkx) = 0.
       emf(kpen:mkx) = 0.
       dwten(kpen:mkx) = 0.
       diten(kpen:mkx) = 0.
       fer(kpen:mkx) = 0.
       fdr(kpen:mkx) = 0. 
 
       do k = kpen - 1, kbup, -1   ! 'k' is interface index at which 'emf' is calculated. 

          rhos0j = ps0(k) / ( rgas * 0.5 * ( thv0bot(k+1) + thv0top(k) ) * exns0(k) )

          if( k .eq. kpen - 1 ) then

             !-----  Calculate ppen( < 0 ), updarft penetrative distance from the lower 
             !       interface of each layer assuming no lateral entrain.-detrain.   No 
             !       lateral entrain.-detrain. is assumed only for calculation of 'ppen' 
             !       I should check if I should multiply 'rbuoy' here or not 

             bogbot = rbuoy * ( thvu(k)/thv0bot(kpen)    - 1. )
             bogtop = rbuoy * ( thvu(kpen)/thv0top(kpen) - 1. )

             aquad =  ( bogtop - bogbot ) / ( ps0(kpen) - ps0(k) )
             bquad =  2. * bogbot
             cquad = -wu(k)**2 * ps0(k) / ( rgas * thv0bot(kpen) * exns0(k) )
             call roots(aquad,bquad,cquad,xc1,xc2,status)
             if( status .eq. 0 ) then
                if( xc1 .le. 0. .and. xc2 .le. 0. ) then
                   ppen = max( xc1, xc2 )
                else
                   ppen = min( xc1, xc2 )
                endif
                ppen    = min( 0.,max( -dp0(k+1), ppen ) )  
             else
                ppen    = 0.
             endif

             !----- Calculate returning mass flux, emf.
             !      This penetrative entrainment dynamics are not definite at this stage.
             !      Definition of conservative scalars of entrained airs can be refined.
             !      Note that 'emf ~ - umf/ufrc = - w * rho'. Thus, below limit set the
             !      upper limit of |emf| to be ~ 10cm/s, which is very loose constraint.

             emf(k)  = max( umf(k) * ppen * rei(kpen) * rpen, -0.1*rhos0j )
             thlu(k) = thl0(kpen) + ssthl0(kpen) * ( ps0(k) - p0(kpen) )
             qtu(k)  = qt0(kpen)  + ssqt0(kpen)  * ( ps0(k) - p0(kpen) )
             uu(k)   = u0(kpen)   + ssu0(kpen)   * ( ps0(k) - p0(kpen) )     
             vu(k)   = v0(kpen)   + ssv0(kpen)   * ( ps0(k) - p0(kpen) )   

          else

             emf(k)  = max( emf(k+1) - umf(k) * dp0(k+1) * rei(k+1) * rpen, -0.1*rhos0j )
             thlu(k) = ( thlu(k+1) * emf(k+1) + thl0(k+1) * ( emf(k) - emf(k+1) ) ) / emf(k)
             qtu(k)  = ( qtu(k+1)  * emf(k+1) + qt0(k+1)  * ( emf(k) - emf(k+1) ) ) / emf(k)
             uu(k)   = ( uu(k+1)   * emf(k+1) + u0(k+1)   * ( emf(k) - emf(k+1) ) ) / emf(k)
             vu(k)   = ( vu(k+1)   * emf(k+1) + v0(k+1)   * ( emf(k) - emf(k+1) ) ) / emf(k)

          endif

          !----- After calculating emf(k), we set no updraft mass flux at interfaces '>= kbup'
          !      for the purpose of easy calculation of turbulent fluxes. Thus, at interfaces 
          !      '>= kbup', emf(k) is 'NET' downward mass flux. As shown below, at all interf
          !      aces, when emf(k) is non-zero ('>= kbup'), umf(k) is zero,    and vice versa 
          !      (for '< kbup'). Thus emf(k) and umf(k) are mutually independent   each other 
          !      at all interfaces. This does not means that they are physically  independent
          !      each other -  this is just for an easy calculation of turbulent fluxes below
          !      using a single formula.
          !      Important question is : if I don't correctly save 'umf' as a output variable
          !                              it might have a significant influence on calculation
          !                              of the other modules which explictly use 'umf' !     
          !                              Also, the calculation of dwten, diten should be done
          !                              cautiously.       
          !      Turning-off below option might strongly diffuse stratocumulus top interface,
          !      resulting in the reduction of cloud fraction. I should check this later.

          umf(k) = 0.0 

       end do

       !----- Upward turbulent fluxes of conservative scalars by cumulus updraft 
       !      In the calculations below, definition of conservative scalars  for 
       !      entrainment flux (emf(k)) is reasonable, consistent  with a simple 
       !      bulk modelling entrainment flux at inversion base.  


       if( krel.ge.kinv ) then      


       ! 
       ! 1. For interfaces  " 0 <= k <= kinv-1 "
       !  

       xsrc  = qtsrc
       xmean = qt0(kinv)
       xtop  = qt0(kinv+1) + ssqt0(kinv+1) * ( ps0(kinv)   - p0(kinv+1) )
       xbot  = qt0(kinv-1) + ssqt0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )        
       call fluxbelowinv( cbmf, ps0(0:mkx), mkx, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
       qtflx(0:kinv-1) = xflx(0:kinv-1)

       xsrc  = thlsrc
       xmean = thl0(kinv)
       xtop  = thl0(kinv+1) + ssthl0(kinv+1) * ( ps0(kinv)   - p0(kinv+1) )
       xbot  = thl0(kinv-1) + ssthl0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )        
       call fluxbelowinv( cbmf, ps0(0:mkx), mkx, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
       slflx(0:kinv-1) = cp * exns0(0:kinv-1) * xflx(0:kinv-1)

       xsrc  = usrc
       xmean = u0(kinv)
       xtop  = u0(kinv+1) + ssu0(kinv+1) * ( ps0(kinv)   - p0(kinv+1) )
       xbot  = u0(kinv-1) + ssu0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )
       call fluxbelowinv( cbmf, ps0(0:mkx), mkx, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
       umflx(0:kinv-1) = xflx(0:kinv-1)

       xsrc  = vsrc
       xmean = v0(kinv)
       xtop  = v0(kinv+1) + ssv0(kinv+1) * ( ps0(kinv)   - p0(kinv+1) )
       xbot  = v0(kinv-1) + ssv0(kinv+1) * ( ps0(kinv)   - p0(kinv+1) )
       call fluxbelowinv( cbmf, ps0(0:mkx), mkx, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
       vmflx(0:kinv-1) = xflx(0:kinv-1)

       !
       ! 2. For interfaces, " kinv <= k <= krel-1 "
       !

       if( krel.ne.kinv ) then

       uplus = 0.
       vplus = 0.
       do k =  kinv, krel-1 
          qtflx(k) = cbmf * ( qtsrc  - (  qt0(k+1) +  ssqt0(k+1) * ( ps0(k) - p0(k+1) ) ) )          
          slflx(k) = cbmf * ( thlsrc - ( thl0(k+1) + ssthl0(k+1) * ( ps0(k) - p0(k+1) ) ) ) * cp * exns0(k)
          uplus = uplus + PGFc * ssu0(k) * ( ps0(k) - ps0(k-1) )
          vplus = vplus + PGFc * ssv0(k) * ( ps0(k) - ps0(k-1) )
          umflx(k) = cbmf * ( usrc + uplus -  (  u0(k+1)  +   ssu0(k+1) * ( ps0(k) - p0(k+1) ) ) ) 
          vmflx(k) = cbmf * ( vsrc + vplus -  (  v0(k+1)  +   ssv0(k+1) * ( ps0(k) - p0(k+1) ) ) )
       end do
      
       endif 


       else


       ! 
       ! 1-2. For interfaces  " 0 <= k <= krel-1 "
       !  


       xsrc  = qtsrc
       xmean = qt0(krel)
       xtop  = qt0(krel+1) + ssqt0(krel+1) * ( ps0(krel)   - p0(krel+1) )
       xbot  = qt0(krel-1) + ssqt0(krel-1) * ( ps0(krel-1) - p0(krel-1) )        
       call fluxbelowinv( cbmf, ps0(0:mkx), mkx, krel, dt, xsrc, xmean, xtop, xbot, xflx )
       qtflx(0:krel-1) = xflx(0:krel-1)

       xsrc  = thlsrc
       xmean = thl0(krel)
       xtop  = thl0(krel+1) + ssthl0(krel+1) * ( ps0(krel)   - p0(krel+1) )
       xbot  = thl0(krel-1) + ssthl0(krel-1) * ( ps0(krel-1) - p0(krel-1) )        
       call fluxbelowinv( cbmf, ps0(0:mkx), mkx, krel, dt, xsrc, xmean, xtop, xbot, xflx )
       slflx(0:krel-1) = cp * exns0(0:krel-1) * xflx(0:krel-1)
       
       xsrc  = usrc
       xmean = u0(krel)
       xtop  = u0(krel+1) + ssu0(krel+1) * ( ps0(krel)   - p0(krel+1) )
       xbot  = u0(krel-1) + ssu0(krel-1) * ( ps0(krel-1) - p0(krel-1) )
       call fluxbelowinv( cbmf, ps0(0:mkx), mkx, krel, dt, xsrc, xmean, xtop, xbot, xflx )
       umflx(0:krel-1) = xflx(0:krel-1)

       xsrc  = vsrc
       xmean = v0(krel)
       xtop  = v0(krel+1) + ssv0(krel+1) * ( ps0(krel)   - p0(krel+1) )
       xbot  = v0(krel-1) + ssv0(krel+1) * ( ps0(krel)   - p0(krel+1) )
       call fluxbelowinv( cbmf, ps0(0:mkx), mkx, krel, dt, xsrc, xmean, xtop, xbot, xflx )
       vmflx(0:krel-1) = xflx(0:krel-1)


       end if


       !
       ! 3. For interfaces, " krel <= k <= kpen-1 " 
       !

       do k = krel, kpen - 1      
          kp1 = k + 1

          slflx(k) = cp * exns0(k) * umf(k)*(thlu(k) - (thl0(kp1) + ssthl0(kp1)*(ps0(k) - p0(kp1)))) +   &
                     cp * exns0(k) * emf(k)*(thlu(k) - (thl0(k)   + ssthl0(k)  *(ps0(k) - p0(k))))

          qtflx(k) = umf(k)*(qtu(k) - (qt0(kp1) + ssqt0(kp1)*(ps0(k) - p0(kp1)))) +                      &
                     emf(k)*(qtu(k) - (qt0(k)   + ssqt0(k)  *(ps0(k) - p0(k)))) 
          umflx(k) = umf(k)*(uu(k) - (u0(kp1) + ssu0(kp1)*(ps0(k) - p0(kp1)))) +                         &
                     emf(k)*(uu(k) - (u0(k)   + ssu0(k)  *(ps0(k) - p0(k)))) 
          vmflx(k) = umf(k)*(vu(k) - (v0(kp1) + ssv0(kp1)*(ps0(k) - p0(kp1)))) +                         &
                     emf(k)*(vu(k) - (v0(k)   + ssv0(k)  *(ps0(k) - p0(k)))) 
       end do

       !
       !----- Calculate model tendencies. 
       !

       do k = 1, kpen

          km1 = k - 1 

          !----- Horizontal momentum

          uten(k) = ( umflx(km1) - umflx(k) ) * ggr / dp0(k)
          vten(k) = ( vmflx(km1) - vmflx(k) ) * ggr / dp0(k) 
          uf(k)   = u0(k) + uten(k) * dt
          vf(k)   = v0(k) + vten(k) * dt

       end do        

       precip  = 0.
       rliq    = 0.
       rainflx = 0. 
       snowflx = 0.

       do k = 1, kpen

          km1 = k - 1

          !----- Static energy and Water Substance Tendencies
          !
          !      Assumption : 1. Cumulus updraft detrains water at the top interface into  the
          !                      layer just underneath that interface  
          !                   2. Detrained water can be either 'cloud water' or 'rain',  which
          !                      are separate identities, e.g., separated by radius difference
          !                   3. Detrained 'cloud water' must stay in the layer that they were 
          !                      formed. However, detrained rain water can move down to the 
          !                      release level without evaporation

          dwten(k) = dwten(k) * umf(k) * ggr / dp0(k) ! [ kg/kg/s ]
          diten(k) = diten(k) * umf(k) * ggr / dp0(k)

          qrten(k) = frc_rasn_remain * frc_rasn * dwten(k) ! [ kg/kg/s ]
          qsten(k) = frc_rasn_remain * frc_rasn * diten(k) 

          rainflx  = rainflx + ( 1. - frc_rasn_remain ) * frc_rasn * dwten(k) * dp0(k) / ggr ! [kg/m^2/s]
          snowflx  = snowflx + ( 1. - frc_rasn_remain ) * frc_rasn * diten(k) * dp0(k) / ggr

          ! 'qc' are reserved water tendencies, independent of qtten and qrten (or qsten), and
          !  they are always positive.

          qc_l  = ( 1. - frc_rasn ) * dwten(k) ! [ kg/kg/s ]
          qc_i  = ( 1. - frc_rasn ) * diten(k)
          qc(k) = qc_l + qc_i   

          ! Below 'slten' formula assume that 'dwten or diten' is separate independent tendency
          ! not influencing 'qtten' (or qlten), consistent with the treatment of 'qtten' below.
 
          slten   = ( slflx(km1) - slflx(k) ) * ggr / dp0(k) + ( lcond * dwten(k) + lsub * diten(k) )

          if( k.eq.1 ) then
            slten = slten - ggr / 4 / dp0(k) * (                              &
                            umflx(k)*(uf(k+1) - uf(k) + u0(k+1) - u0(k)) +    & 
                            vmflx(k)*(vf(k+1) - vf(k) + v0(k+1) - v0(k)))
          elseif( k.ge.2 .and. k.le.kpen-1 ) then
            slten = slten - ggr / 4 / dp0(k) * (                              &
                            umflx(k)*(uf(k+1) - uf(k) + u0(k+1) - u0(k)) +    &
                            umflx(k-1)*(uf(k) - uf(k-1) + u0(k) - u0(k-1)) +  &
                            vmflx(k)*(vf(k+1) - vf(k) + v0(k+1) - v0(k)) +    &
                            vmflx(k-1)*(vf(k) - vf(k-1) + v0(k) - v0(k-1)))
          elseif( k.eq.kpen ) then
            slten = slten - ggr / 4 / dp0(k) * (                              &
                            umflx(k-1)*(uf(k) - uf(k-1) + u0(k) - u0(k-1)) +  &
                            vmflx(k-1)*(vf(k) - vf(k-1) + v0(k) - v0(k-1)))
          else
            write(6,*) 'Warning : Tendencies should be calculated for 1<= k <= kpen in UW_Shallow'
          endif

          qtten = ( qtflx(km1) - qtflx(k) ) * ggr / dp0(k) - ( dwten(k) + diten(k) )

          ! Option 1.

!         qt_dt    = qt0(k)  +   qtten * dt
!         thl_dt   = thl0(k) +   slten / cp /exn0(k) * dt
!         if( qt_dt.ge.1.e-5 ) then 
!           call conden(p0(k),thl_dt,qt_dt,thj,qvj,qlj,qij,qse,id_check,qsat)
!           if(id_check.eq.1) then
!              id_exit = .true.
!              go to 333
!           end if  
!           qlten(k) = ( qlj - ql0(k) ) / dt + qc_l - qc_l
!           qiten(k) = ( qij - qi0(k) ) / dt + qc_i - qc_i
!           qvten(k) = qtten - qlten(k) - qiten(k)
!         else
!         ! This can occur  at 'kbup' or, less possibly, at 'kinv-1' layers. 
!         ! The former is mainly associated with incomplete penetrative entrainment dynamics.
!           write(6,*) 'Negative qt(t+dt)', i, k, kmin, klcl, kinv, krel, klfc, kbup, kpen
!         ! goto 333
!           qlten(k) = qc_l - qc_l
!           qiten(k) = qc_i - qc_i
!           qvten(k) = qtten
!          endif

          ! Option 2.

!           qlten(k) = qc_l - qc_l
!           qiten(k) = qc_i - qc_i
!           qvten(k) = qtten

          ! Option 3.

          qvten(k) = qtten*qv0(k)/qt0(k)
          qlten(k) = qtten*ql0(k)/qt0(k)
          qiten(k) = qtten*qi0(k)/qt0(k)

          !----- Dry Static Energy        
          !      Note that  dry static energy tendency is not influenced by the conversion
          !      of reserved water tendency (qc(k)) into total water tendency, in contrast  
          !      to qtten (or qlten).

          sten(k) = slten + lcond * qlten(k) + lsub * qiten(k)

          ! Note that 'precip' [m/s] is vertically-integrated total 'rain+snow' formed from 
          ! cumulus updraft, which will fall into the ground in the 'zm_conv_evap' subroutine
          ! Important : 1000 is rhoh2o (water density) [kg/m^3] used for unit conversion from
          !             [kg/m^2/s] to [m/s] for use in stratiform.F90.

          precip  =  precip + frc_rasn * ( dwten(k) + diten(k) ) * dp0(k) / ggr / 1000. ! [m/s]
          rliq    =  rliq   +                            qc(k)   * dp0(k) / ggr / 1000. ! [m/s]

       end do

       !----- Detrain the rest of the rain/snow at the release level

       qrten(krel) = qrten(krel) + rainflx * ggr / dp0(krel) ! [kg/kg/s]
       qsten(krel) = qsten(krel) + snowflx * ggr / dp0(krel)


       !----- Convective cloud fraction ( cldfrc ) & In-cloud 'liquid+ice'
       !      specific humidity ( qlu ) at each layer for radiation scheme

       rhos0j    = ps0(krel-1) / ( rgas * 0.5 * ( thv0bot(krel) + thv0top(krel-1) ) * exns0(krel-1) )
       ufrcbelow = cbmf / ( rhos0j * wu(krel-1) )
       qconbelow = 0.
       rlwp      = 0.

       ! Calculate at each upper/lower interface and do average to get layer value
       do k = krel, kpen - 1
          rhos0j = ps0(k) / ( rgas * 0.5 * ( thv0bot(k+1) + thv0top(k) ) * exns0(k) )
          ! When 'krel>kinv' and CIN>0, below means that 'ufrc' increase with height 
          ! since 'wu' decrease with height, which is contrary to commom sense.   Is
          ! this reasonable ? 
          ufrc   = umf(k) / ( rhos0j * wu(k) )
          call conden(ps0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             id_exit = .true.
             err(i) = 1.13 ! conden doesn't converge
             go to 333
          end if
          rlwp      = rlwp + qlj * ( p0(k) - p0(k+1) ) / ggr * ufrc * 1000.
          cldfrc(k) = 0.5 * ( ufrcbelow + ufrc ) 
          qlu(k)    = 0.5 * ( qconbelow + qlj + qij )
          ufrcbelow = ufrc
          qconbelow = qlj + qij
       end do
 
       cnt = float(kpen)
       cnb = float(krel - 1)
 
       !
       !----- Adjust the original input profiles for iterative CIN calculation
       !

       if(iter .ne. iter_cin) then 

          ! Save the output from "iter_cin=1"  

          qv0_s(:mkx)           = qv0(:mkx) + qvten(:mkx) * dt
          ql0_s(:mkx)           = ql0(:mkx) + qlten(:mkx) * dt
          qi0_s(:mkx)           = qi0(:mkx) + qiten(:mkx) * dt
          s0_s(:mkx)            = s0(:mkx) + sten(:mkx) * dt 
          u0_s(:mkx)            = u0(:mkx) + uten(:mkx) * dt
          v0_s(:mkx)            = v0(:mkx) + vten(:mkx) * dt 
          qt0_s(:mkx)           = qv0_s(:mkx) + ql0_s(:mkx) + qi0_s(:mkx)
          t0_s(:mkx)            = t0(:mkx) + sten(:mkx) * dt / cp

          umf_s(0:mkx)          = umf(0:mkx)
          qvten_s(:mkx)         = qvten(:mkx)
          qlten_s(:mkx)         = qlten(:mkx)  
          qiten_s(:mkx)         = qiten(:mkx)
          sten_s(:mkx)          = sten(:mkx)
          uten_s(:mkx)          = uten(:mkx)  
          vten_s(:mkx)          = vten(:mkx)
          qrten_s(:mkx)         = qrten(:mkx)
          qsten_s(:mkx)         = qsten(:mkx)  
          precip_s              = precip
          cush_s                = cush
          cldfrc_s(:mkx)        = cldfrc(:mkx)  
          slflx_s(0:mkx)        = slflx(0:mkx)  
          qtflx_s(0:mkx)        = qtflx(0:mkx)  
          qlu_s(:mkx)           = qlu(:mkx)  
          fer_s(:mkx)           = fer(:mkx)  
          fdr_s(:mkx)           = fdr(:mkx)  
          cin_s                 = cin
          cbmf_s                = cbmf
          rliq_s                = rliq
          qc_s(:mkx)            = qc(:mkx)
          cnt_s                 = cnt
          cnb_s                 = cnb
 

          ! Recalculate environmental variables for new cin calculation at "iter_cin=2" 
          ! Perform only for variables necessary for the new cin calculation.

          qv0(:mkx) = qv0_s(:mkx)
          ql0(:mkx) = ql0_s(:mkx)
          qi0(:mkx) = qi0_s(:mkx)
          s0(:mkx)  = s0_s(:mkx)
          t0(:mkx)  = t0_s(:mkx)
      
          qt0(:mkx)   = (qv0(:mkx) + ql0(:mkx) + qi0(:mkx))
          thl0(:mkx)  = (t0(:mkx) - lcond*ql0(:mkx)/cp - lsub*qi0(:mkx)/cp)/exn0(:mkx)
          thvl0(:mkx) = (1. + zvir*qt0(:mkx))*thl0(:mkx)

          ssthl0 = slope(mkx,thl0,p0) ! Dimension of ssthl0(:mkx) is implicit
          ssqt0  = slope(mkx,qt0 ,p0)
          ssu0   = slope(mkx,u0  ,p0)
          ssv0   = slope(mkx,v0  ,p0)

          do k = 1, mkx

          thl0bot = thl0(k) + ssthl0(k)*(ps0(k-1) - p0(k))
          qt0bot  = qt0(k)  + ssqt0(k) *(ps0(k-1) - p0(k))
          call conden(ps0(k-1),thl0bot,qt0bot,thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             id_exit = .true.
             err(i) = 1.14 ! conden doesn't converge
             go to 333
          end if
          thv0bot(k)  = thj*(1. + zvir*qvj - qlj - qij)
          thvl0bot(k) = thl0bot*(1. + zvir*qt0bot)
          
          thl0top     = thl0(k) + ssthl0(k)*(ps0(k) - p0(k))
          qt0top      =  qt0(k) + ssqt0(k) *(ps0(k) - p0(k))
          call conden(ps0(k),thl0top,qt0top,thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             err(i) = 1.15 ! conden doesn't converge
             id_exit = .true.
             go to 333
          end if
          thv0top(k)  = thj*(1. + zvir*qvj - qlj - qij)
          thvl0top(k) = thl0top*(1. + zvir*qt0top)

          end do

       endif

     end do                ! end of iter loop (cin_iter)      

     !
     ! Update Output Variables
     ! 

     umf_out(i,0:mkx)         = umf(0:mkx)
     slflx_out(i,0:mkx)       = slflx(0:mkx)
     qtflx_out(i,0:mkx)       = qtflx(0:mkx)
     qvten_out(i,:mkx)        = qvten(:mkx)
     qlten_out(i,:mkx)        = qlten(:mkx)
     qiten_out(i,:mkx)        = qiten(:mkx)
     sten_out(i,:mkx)         = sten(:mkx)
     uten_out(i,:mkx)         = uten(:mkx)
     vten_out(i,:mkx)         = vten(:mkx)
     qrten_out(i,:mkx)        = qrten(:mkx)
     qsten_out(i,:mkx)        = qsten(:mkx)
     precip_out(i)            = precip
     cldfrc_out(i,:mkx)       = cldfrc(:mkx)
     qlu_out(i,:mkx)          = qlu(:mkx)
     fer_out(i,:mkx)          = fer(:mkx)
     fdr_out(i,:mkx)          = fdr(:mkx)
     cush_inout(i)            = cush
     cinh_out(i)              = cin
     cbmf_out(i)              = cbmf
     rliq_out(i)              = rliq
     qc_out(i,:mkx)           = qc(:mkx)
     cnt_out(i)               = cnt
     cnb_out(i)               = cnb
 333 if(id_exit) then ! Exit without cumulus convection

     umf_out(i,0:mkx)         = 0.
     slflx_out(i,0:mkx)       = 0.
     qtflx_out(i,0:mkx)       = 0.
     qvten_out(i,:mkx)        = 0.
     qlten_out(i,:mkx)        = 0.
     qiten_out(i,:mkx)        = 0.
     sten_out(i,:mkx)         = 0.
     uten_out(i,:mkx)         = 0.
     vten_out(i,:mkx)         = 0.
     qrten_out(i,:mkx)        = 0.
     qsten_out(i,:mkx)        = 0.
     precip_out(i)            = 0.
     cldfrc_out(i,:mkx)       = 0.
     qlu_out(i,:mkx)          = 0.
     fer_out(i,:mkx)          = 0.
     fdr_out(i,:mkx)          = 0.
     cush_inout(i)            = -1.
     cinh_out(i)              = 0.
     cbmf_out(i)              = 0.
     rliq_out(i)              = 0.
     qc_out(i,:mkx)           = 0.
     cnt_out(i)               = float(mkx)
     cnb_out(i)               = 0.

     end if
 
     end do                  ! end of big i loop

    return

  end subroutine compute_mcshallow

! ---------------------------------------------------------------------------------------------------


  subroutine getbuoy(pbot,thv0bot,ptop,thv0top,thvubot,thvutop,plfc,cin)
    !
    !----- Subroutine to calculate integrated CIN and 'cinlcl, plfc' if any. 
    !      Assume 'thv' is linear in each layer for both cumulus and enviro.
    !      Note this routine only  calculates positive CIN and negative CIN
    !      is not included into CIN calculation.

    real(r8) pbot,thv0bot,ptop,thv0top,thvubot,thvutop,plfc,cin,frc

    if( thvubot .gt. thv0bot .and. thvutop .gt. thv0top ) then
       plfc = pbot
       return
    elseif( thvubot .le. thv0bot .and. thvutop .le. thv0top ) then 
       cin  = cin - ((thvubot/thv0bot - 1.) + (thvutop/thv0top - 1.)) * (pbot - ptop)/          &
                    (pbot/(rgas*thv0bot*exnf(pbot)) + ptop/(rgas*thv0top*exnf(ptop)))
    elseif( thvubot .gt. thv0bot .and. thvutop .le. thv0top ) then 
       frc  = (thvutop/thv0top - 1.)/((thvutop/thv0top - 1.) - (thvubot/thv0bot - 1.))
       cin  = cin - (thvutop/thv0top - 1.)*((ptop + frc*(pbot - ptop)) - ptop)/                 &
                    (pbot/(rgas*thv0bot*exnf(pbot)) + ptop/(rgas*thv0top*exnf(ptop)))
    else            
       frc  = (thvubot/thv0bot - 1.)/((thvubot/thv0bot - 1.) - (thvutop/thv0top - 1.))
       plfc = pbot - frc*(pbot - ptop)
       cin  = cin - (thvubot/thv0bot - 1.)*(pbot - plfc)/                                       & 
                    (pbot/(rgas*thv0bot*exnf(pbot)) + ptop/(rgas*thv0top * exnf(ptop)))
    endif

    return
  end subroutine getbuoy


  function single_cin(pbot,thv0bot,ptop,thv0top,thvubot,thvutop)
    !
    !----- Function to calculate a single layer CIN summing all positive and negative CIN

    real(r8) :: single_cin
    real(r8)    pbot,thv0bot,ptop,thv0top,thvubot,thvutop 

    single_cin = ((1. - thvubot/thv0bot) + (1. - thvutop/thv0top))*(pbot - ptop)/         &
                 (pbot/(rgas*thv0bot*exnf(pbot)) + ptop/(rgas*thv0top*exnf(ptop)))

    
    return
  end function single_cin   


  subroutine conden(p,thl,qt,th,qv,ql,qi,rvls,id_check,qsat)
    !
    !----- Calculate thermodynamic properties from a given set of ( p, thl, qt )
    !      I should check later if this subroutine is correct or not   

    implicit none
    real(r8), intent(in)  :: p
    real(r8), intent(in)  :: thl
    real(r8), intent(in)  :: qt
    real(r8), intent(out) :: th
    real(r8), intent(out) :: qv
    real(r8), intent(out) :: ql
    real(r8), intent(out) :: qi
    real(r8), intent(out) :: rvls
    integer , intent(out) :: id_check
    integer , external    :: qsat
    real(r8)              :: tc,temps,t
    real(r8)              :: leff, nu, qc
    integer               :: iteration
    real(r8)              :: es(1)                     ! saturation vapor pressure
    real(r8)              :: qs(1)                     ! saturation spec. humidity
    real(r8)              :: gam(1)                    ! (L/cp)*dqs/dT
    integer               :: status                    ! return status of qsat call

    tc     = thl*exnf(p)
    nu     = max(min((268. - tc)/20.,1.0),0.0)         ! Fraction of ice in the condensate. Func. of T.
    leff   = (1. - nu)*lcond + nu*lsub                    ! This is an estimate that hopefully speeds convergence

    ! Below "temps" and "rvls" are just initial guesses for iteration loop below.
    ! Note that the output "temps" from the below iteration loop is "temperature"
    ! NOT "liquid temperature". 

    temps  = tc
    status = qsat(temps,p,es(1),qs(1),gam(1), 1)
    rvls   = qs(1)

    if(qs(1).ge.qt) then  

      id_check = 0
      qv = qt
      qc = 0.
      ql = 0.
      qi = 0.
      th = tc/exnf(p)

    else 

      do iteration = 1, 10
         temps  = temps + (( tc - temps )*cp/leff + qt - rvls)/               & 
                          (cp/leff + ep2*leff*rvls/rgas/temps/temps)
         status = qsat(temps,p,es(1),qs(1),gam(1),1)
         rvls   = qs(1)
      end do
         
      qc = max(qt - qs(1),0.)
      qv = qt - qc
      ql = qc*(1. - nu)
      qi = nu*qc
      th = temps/exnf(p)

      if( abs((temps-(leff/cp)*qc)-tc).ge.1.) then
        id_check = 1
      else
        id_check = 0
      end if

    end if

    return
  end subroutine conden


  subroutine roots(a,b,c,r1,r2,status)
    !
    ! Subroutine to solve the second order polynomial equation
    ! I should check this subroutine later

    real(r8), intent(in)  :: a
    real(r8), intent(in)  :: b
    real(r8), intent(in)  :: c
    real(r8), intent(out) :: r1
    real(r8), intent(out) :: r2
    integer , intent(out) :: status
    real(r8)              :: q

    status = 0
    if(a .eq. 0) then                        ! form b*x + c = 0
       if(b .eq. 0) then                     ! failure: c = 0
          status = 1
       else                                  ! b*x + c = 0
          r1 = -c/b
       endif
       r2 = r1
    else
       if(b .eq. 0.) then                    ! form a*x**2 + c = 0
          if(a*c .gt. 0.) then               ! failure: x**2 = -c/a < 0
             status = 2  
          else                               ! x**2 = -c/a 
             r1 = sqrt(-c/a)
          endif
          r2 = -r1
       else                                  ! form a*x**2 + b*x + c = 0
          if((b**2 - 4.*a*c) .lt. 0.) then   ! failure, no real roots
             status = 3
          else
             q  = -0.5*(b + sign(1.0,b)*sqrt(b**2 - 4.*a*c))
             r1 =  q/a
             r2 =  c/q
          endif
       endif
    endif

    return
  end subroutine roots

  
  function slope(mkx,field,p0)
    !
    !----- Function performing profile reconstruction of conservative scalars in each layer.
    !      Identical to profile reconstruction used in UW-PBL scheme but from bottom to top 
    !      layer here. At the lowest layer near to surface, slope is defined using  the two
    !      lowest layer midpoint values. I checked this subroutine and it is correct. 

    real(r8)             :: slope(mkx)
    integer,  intent(in) :: mkx
    real(r8), intent(in) :: field(mkx)
    real(r8), intent(in) :: p0(mkx)
    
    real(r8)             :: below
    real(r8)             :: above
    integer              :: k

    below = (field(2) - field(1))/(p0(2) - p0(1))
    do k = 2, mkx
       above = (field(k) - field(k-1))/(p0(k) - p0(k-1))
       if (above .gt. 0.) then
          slope(k-1) = max(0.,min(above,below))
       else 
          slope(k-1) = min(0.,max(above,below))
       end if
       below = above
    end do
    slope(mkx) = slope(mkx-1)

    return
  end function slope


  function qsinvert(qt,thl,qsat)
    !
    !----- Function calculating saturation pressure ps (or pLCL) from qt and
    !      thl (liquid potential temperature, NOT liquid virtual potential 
    !      temperature) by inverting Bolton formula. From MinCu.
    !      I should check if the use of 'leff' instead of 'xlv' in this
    !      subroutine is correct or not.

    real(r8)          :: qsinvert    
    real(r8)             qt, thl
    real(r8)             ps, Pis, Ts, err, dlnqsdT, dTdPis
    real(r8)             dPisdps, dlnqsdps, derrdps, dps 
    integer              i
    integer, external :: qsat
    real(r8)          :: es(1)                     ! saturation vapor pressure
    real(r8)          :: qs(1)                     ! saturation spec. humidity
    real(r8)          :: gam(1)                    ! (L/cp)*dqs/dT
    integer           :: status                    ! return status of qsat call
    real(r8)          :: leff, nu

    ps = p00
    do i = 1, 10
      Pis      =  (ps/p00)**rovcp
      Ts       =  thl*Pis
      status   =  qsat(Ts,ps,es(1),qs(1),gam(1),1)
      err      =  qt - qs(1)
      nu       =  max(min((268. - Ts)/20.,1.0),0.0)        
      leff     =  (1. - nu)*lcond + nu*lsub                   
      dlnqsdT  =  gam(1)*(cp/leff)/qs(1)
      dTdPis   =  thl
      dPisdps  =  rovcp*Pis/ps 
      dlnqsdps = -1./(ps - (1. - ep2)*es(1))
      derrdps  = -qs(1)*(dlnqsdT * dTdPis * dPisdps + dlnqsdps)
      dps      = -err/derrdps
      ps       =  ps + dps
    end do  
    qsinvert = ps 

    return
  end function qsinvert


  real(r8) function compute_alpha(del_cin,ke)

    real(r8) :: del_cin, ke
    real(r8) :: x0, x1
    integer  :: iteration

    x0 = 0.
    do iteration = 1, 10
       x1 = x0 - (exp(-x0*ke*del_cin) - x0)/(-ke*del_cin*exp(-x0*ke*del_cin) - 1.)
       x0 = x1
    end do
    compute_alpha = x0

    return
  end function compute_alpha


  subroutine fluxbelowinv(cbmf,ps0,mkx,kinv,dt,xsrc,xmean,xtopin,xbotin,xflx)   
    !
    ! Calculate turbulent fluxes at and below 'kinv-1' interfaces.
    ! Check in the main program : Input 'cbmf' should not be zero.
    ! After calculating turbulent heat fluxes at 'kinv-1' interface by doing appropriate
    ! mass weighting average ( rpinv and rcbmf, or rr = rpinv / rcbmf ) for two possible 
    ! cases, turbulent fluxes at the lower interfaces ( 0 <= k <= 'kinv-2' ) are simply 
    ! linearly interpolated by assuming turbulent flux is zero at surface interface. 
    ! Note that this is a simplified version than the original version, but conceptually
    ! more clear and general, since this method can be applied even when 'cbmf*g*dt'  is
    ! larger than 'dp(kinv-1)'.  In fact, this routine is valid without any modification
    ! as long as 'cbmf*g*dt' is smaller than 'dpi', which is the depth of PBL. Note that
    ! in the main program, the upper limit of 'cbmf*g*dt' is set to be 'dpi'.
    ! This subroutine can be further elaborated for smooth transition across 'kinv-1'.  
        
    integer,  intent(in)                     :: mkx, kinv 
    real(r8), intent(in)                     :: cbmf, dt, xsrc, xmean, xtopin, xbotin
    real(r8), intent(in),  dimension(0:mkx)  :: ps0
    real(r8), intent(out), dimension(0:mkx)  :: xflx  

    integer k
    real(r8) rcbmf, rpeff, dp, rr, pinv_eff, xtop, xbot

    xflx(0:mkx) = 0.
    dp = ps0(kinv-1)-ps0(kinv) 
    
    if( abs(xbotin-xtopin).le.1.e-33 ) then
      xbot = xbotin - 1.e-20
      xtop = xtopin + 1.e-20
    else
      xbot = xbotin
      xtop = xtopin      
    endif

    rcbmf = ( cbmf * ggr * dt ) / dp ! Can be larger than 1 : 'OK'      
    rpeff = ( xmean - xtop ) / ( xbot - xtop ) 
    rpeff = min( max(0.,rpeff), 1. )  ! As of this, 0<= rpeff <= 1   
    if(rpeff.eq.1) xbot = xmean
    if(rpeff.eq.0) xtop = xmean    
    rr = rpeff / rcbmf
    pinv_eff = ps0(kinv-1) + ( rcbmf - rpeff ) * dp ! Effective "pinv" after detraining mass

    ! Below two cases exactly converges when rr = 1.
    if( rr.gt.1 ) then
      do k = 0, kinv - 1
         xflx(k) = cbmf * ( xsrc - xbot ) * ( ps0(0) - ps0(k) ) / ( ps0(0) - pinv_eff )
      end do
    else
         xflx(kinv-1) = cbmf * ( rr * ( xsrc - xbot ) + ( 1 - rr ) * ( xsrc - xtop ) )
      do k = 0, kinv - 2
         xflx(k) =  xflx(kinv-1) * ( ps0(0) - ps0(k) ) / ( ps0(0) - ps0(kinv-1) )
      end do
    endif

    return
  end subroutine fluxbelowinv

            
end module uwcu
