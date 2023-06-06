
! Copyright Steve Klein and Mark Webb 2002 - all rights reserved.
!
! This code is available without charge with the following conditions:
!
!  1. The code is available for scientific purposes and is not for 
!     commercial use.
!  2. Any improvements you make to the code should be made available 
!     to the to the authors for incorporation into a future release.
!  3. The code should not be used in any way that brings the authors 
!     or their employers into disrepute.

! -------------------------------------------------------------------
module isccpSimulator
  implicit none
  private
  public :: icarus, icarus_CRM, computeWaterVaporEmissivity

  !
  ! Overloaded procedures
  !
  interface fluxToTb
    module procedure fluxToTb_1D, fluxToTb_2D, fluxToTb_3d
  end interface ! fluxToTb
  
  interface TbToFlux
    module procedure TbToFlux_1D, TbToFlux_2D, TbToFlux_3d
  end interface ! TbToFlux
  
  interface computeRadiance
    module procedure computeRadiance_1D, computeRadiance_2D, computeRadiance_3D
  end interface ! computeRadiance
  
contains      
! -------------------------------------------------------------------
  subroutine icarus(dtau, pfull, top_height,            & ! Required 
                    dem, dem_wv, at, skt, emsfc_lw, iTrop, & ! Optional
                    boxtau, boxptop)
    !
    ! Required input arguments
    !
    real,    dimension(:, :, :), & ! Dimensions nPoints, nCol, nLev
      intent( in)  :: dtau         !  mean 0.67 micron optical depth of clouds in each model level
                                   !  NOTE:  this the cloud optical depth of only the cloudy part of the grid box, 
                                   !  it is not weighted with the 0 cloud optical depth of the clear part of the grid box
    real,    dimension(:, :),    &
      intent( in) :: pfull         !  pressure of full model levels (Pascals)
                                   !  pfull(npoints,1)    is    top level of model
                                   !  pfull(npoints,nlev) is bottom level of model
    integer,                    &
      intent( in) ::  top_height   !  1 = adjust top height using both a computed
                                   !  infrared brightness temperature and the visible
                                   !  optical depth to adjust cloud top pressure. Note
                                   !  that this calculation is most appropriate to compare
                                   !  to ISCCP data during sunlit hours.
                                   !  2 = do not adjust top height, that is cloud top
                                   !  pressure is the actual cloud top pressure
                                   !  in the model
                                   !  3 = adjust top height using only the computed
                                   !  infrared brightness temperature. Note that this
                                   !  calculation is most appropriate to compare to ISCCP
                                   !  IR only algortihm (i.e. you can compare to nighttime
                                   !  ISCCP data with this option)
    !
    ! Optional input arguments - for computing radiative cloud-top pressure 
    !   All variables except itrop are needed if top_height == 1 .or. top_height == 3
    !
    real,    dimension(:, :, :), optional, & ! Dimensions nPoints, nCol, nLev
     intent( in) ::  dem         !  10.5 micron longwave emissivity of clouds in each model level.  
                                 !  Same note applies as in dtau.
   real,   dimension(:, :),      optional, & ! Dimensions nPoints, nLev
     intent( in) :: dem_wv, &    ! Water vapor emissivity.  
                    at           ! Air temperature in each model level (K)
   real,   dimension(:),         optional, &
     intent( in) :: skt, &       !  skin Temperature (K)
                    emsfc_lw     !  10.5 micron emissivity of surface (fraction)
   ! Users may supply their own tropopause height indicator
   integer, dimension(:),        optional, & ! Dimension nPoints
     intent( in) :: itrop        ! Index of the tropopause location in each column
   
    !     ------
    !     Output
    !     ------
    real, dimension(:, :), &   !  dimension nPoints, nCol
      intent(out) :: boxtau,         &   !  optical thickness in each column
                     boxptop             !  cloud top pressure (mb) in each sub-column
                  
    !     --------------------------------------------
    !     Local variables and parameters
    integer :: nPoints, nCol, nLev
    !     ------
    integer                           :: j, ilev, ibox, icycle
    real,    dimension(size(dtau, 1), &
                       size(dtau, 2)) :: tau, pTop, tb
    real,    dimension(size(dtau, 1), &
                       size(dtau, 2), &
                       size(dtau, 3)) :: dem_local
    
    ! Variables for adjusting cloud top pressure based on TOA IR radiance
    real,    dimension(size(dtau, 1)) :: fluxtop_clrsky
    real,    dimension(size(dtau, 1), &
                       size(dtau, 2)) :: emcld, fluxtop, tauir, taumin, &
                                         fluxtopinit, transmax, btcmin
    logical, dimension(size(dtau, 1), &
                       size(dtau, 2)) :: intermedTrans
  
    ! Historgram quanitities 
    integer, dimension(size(dtau, 1), &
                       size(dtau, 2)) :: levMatch
    
    ! Tropopause 
    integer, dimension(size(dtau, 1)) :: itrop_local    ! index to tropopause level
    real,    dimension(size(dtau, 1)) :: pTrop, atTrop  ! tropopause pressure, temperature
                                   
    !     ------
    !     Local constants
    !     ------
    real,    parameter ::  VisTauToIrTauIce = 1./2.13, VisTauToIrTauLiquid = 1./2.56, &
                           assumedFreezingTb = 260.
    ! -----------------------------------------------------------------------------------
    nPoints = size(dtau, 1); nCol    = size(dtau, 2); nLev    = size(dtau, 3)
    
    ! Error checking
    call check_Lbound(dtau, "cloud optical depth",0.)
    if(present(dem)) &
      call check_bounds(dem, "cloud emissivity",  0., 1.)
    
    ! -----------------------------------------------------------------------------------
    ! The easy part - what's the total cloud optical depth in each column? 
    ! 
    tau(:, :) = sum(dtau, dim = 3)
    
    ! -----------------------------------------------------------------------------------
    ! If we're adjusting cloud-top heights using the VIS and/or IR radiances 
    !    we need to know 1) the TOA radiance, and 2) the position of the tropopause 
    !    (the latter is the pressure level to which really cold clouds are assigned)
    !
    if(top_height == 1 .or. top_height == 3) then
      dem_local(:, :, :) = spread(dem_wv(:, :), dim = 2, nCopies = nCol)
      !
      ! Clear sky calculation - only need one col per GCM grid point
      !
      call computeRadiance(dem_local(:, 1, :), at, skt, emsfc_lw, fluxtop_clrsky(:))
      !
      ! Add contribution to emissivity from cloud
      !
      where(dem(:, :, :) > tiny(dem)) 
        dem_local(:, :, :) = 1. - ( (1. - dem_local(:, :, :)) * (1. -  dem(:, :, :)) )
      end where
      !
      ! And now the all-sky radiance
      !
      call computeRadiance(dem_local, at, skt, emsfc_lw, fluxtop)
    end if
    !     ---------------------------------------------------!
    ! Account for ISCCP procedures to determine cloud top temperature
  
    ! Account for partially transmitting cloud and recompute flux 
    !    ISCCP would see assuming a single layer cloud
    ! Note choice here of VisTauToIrTauIce = 1/2.13, as it is primarily ice
    !    clouds which have partial emissivity and need the adjustment 
    !
    ! If the cloud brightness temperature is greater than 260K,  the liquid cloud conversion
    !   factor is used.
    !
    ! This is discussed on pages 85-87 of the ISCCP D level documentation (Rossow et al. 1996)

    if (top_height == 1 .or. top_height == 3) then
      ! Tropoause height needed if cloud tops are to be adjusted
      !
      if(present(itrop)) then
        itrop_local(:) = itrop(:)
      else 
        call diagnoseTropPressure(pfull, at, itrop_local)
      end if
      do j = 1, nPoints
        ptrop(j)  = pfull(j, itrop_local(j)) 
        attrop(j) =    at(j, itrop_local(j))
      end do

      if (top_height == 1) then
        !compute minimum brightness temperature and optical depth
        btcmin(:, :) = spread(TbToFlux(attrop(:) - 5.), dim = 2, nCopies = nCol)
        
        !note that the initial setting of tauir(:) is needed so that
        !tauir(:) has a realistic value 
        tauir(:, :) = tau(:, :) * VisTauToIrTauIce
        transmax(:, :) = (fluxtop(:,:) - btcmin(:, :)) / &
                         (spread(fluxtop_clrsky(:), dim = 2, nCopies = nCol) - btcmin(:, :))
        taumin(:, :) = -log(max(min(transmax(:, :), 1. - spacing(1.)), tiny(transmax)))

        intermedTrans(:, :) = transmax(:, :) > tiny(transmax) .and. &
                              transmax(:, :) < 1. - spacing(1.)
        where (intermedTrans) 
          fluxtopinit(:, :) = fluxtop(:,:)
          tauir(:, :) = tau(:,:) * VisTauToIrTauIce
        end where
        
        do icycle=1,2
          where (tau(:,:) > tiny(tau) .and. intermedTrans) 
             emcld(:,:) = 1. - exp(-tauir(:, :))
             fluxtop(:,:) = fluxtopinit(:, :) - ( (1. - emcld(:,:)) * &
                                                  spread(fluxtop_clrsky(:), dim = 2, nCopies = nCol) )
             fluxtop(:,:) = max(tiny(fluxtop), fluxtop(:,:)/emcld(:,:))
             tb(:, :) = fluxToTb(fluxtop(:, :))
             where (tb(:, :) .gt. assumedFreezingTb) tauir(:, :) = tau(:,:) * VisTauToIrTauLiquid
          end where
        enddo
      end if  
      
      where(tau(:, :) > tiny(tau)) 
        tb(:, :) = fluxToTb(fluxtop(:, :))
      elsewhere
        ! Clear sky brightness temperature
        tb(:,:) = spread(fluxToTb(fluxtop_clrsky(:)), dim = 2, nCopies = nCol)
      end where 
      
      if(top_height == 1) then
        ! Adjust the brightness temperature and optical depth of very thin clouds
        where (tau(:, :) > tiny(tau) .and. tauir(:, :) < taumin(:, :)) 
           tb(:, :) = spread(attrop(:) - 5., dim= 2, nCopies = nCol)
          tau(:, :) = taumin(:, :) / VisTauToIrTauIce
        end where
      end if
                
    end if
    ! -----------------------------------------------------------------------------------
    ! Determine cloud-top pressure. Three choices:
    !     Radiatively determined cloud top pressure (top_height = 1 for VIS/R;  3 for IR)
    !     Physical cloud top pressure (top_height = 2)
    
    if (top_height .eq. 1 .or. top_height .eq. 3) then  
      !segregate according to optical thickness (??? - RP)
      levmatch(:, :) = 0
      do ibox=1,ncol
        ! find lowest (?) level whose temperature
        ! most closely matches brightness temperature
        !
        do ilev = 1, nlev-1
          where((at(:,ilev) >= tb(:,ibox) .and. at(:,ilev+1) < tb(:,ibox)) .or. &
                (at(:,ilev) <= tb(:,ibox) .and. at(:,ilev+1) > tb(:,ibox)) )
            where(abs(at(:,ilev) - tb(:,ibox)) < abs(at(:,ilev+1) - tb(:,ibox)))
              levmatch(:, ibox) = ilev
            elsewhere 
              levmatch(:, ibox) = ilev + 1
            end where
          end where 
        end do

        ! If we've found a matching level use it, otherwise use the boundary value
        !  
        do j = 1, nPoints
          if(levmatch(j, ibox) >= 1) then
            ptop(j, ibox) = pfull(j, levmatch(j, ibox))
          else if (tb(j, ibox) < minval(at(j, :))) then
            levmatch(j, ibox) = itrop_local(j)
            ptop(j, ibox) = ptrop(j)
          else if (tb(j, ibox) > maxval(at(j, :))) then
            levmatch(j, ibox) = nLev
            ptop(j, ibox) = pFull(j, nLev)
          end if
        end do
      end do
      
    else !  When top_height .eq. 2
      ! The highest cloud top (clouds being where tau > 0). 
      ptop(:, :) = 0.
      do ibox = 1, nCol
        do ilev = 1, nlev
          where(ptop(:, ibox) <= 0. .and. dtau(:, ibox, ilev) > 0.)
            ptop(:, ibox) = pfull(:, ilev)
            levmatch(:,ibox) = ilev
          end where
        end do
      end do
    end if                            
          
    ! No matter how cloud top pressure is determined, 
    !   pTop and levmatch are 0 when the column is clear.   
    ! 
    where (tau(:,:) <= tiny(tau)) 
      ptop(:,:) = 0.
      levmatch(:,:) = 0      
    end where

    !     ---------------------------------------------------!
    ! Output
    !    
    boxptop(:, :) = ptop(:, :) / 100. 
    boxtau(:, :)  = tau(:, :)
    
  end subroutine icarus
  ! -------------------------------------------------------------------
  subroutine icarus_CRM(dtau, pfull, top_height,            & ! Required 
                    dem, dem_wv, at, skt, emsfc_lw, iTrop, & ! Optional
                    boxtau, boxptop)
    !
    ! Required input arguments
    !
    real,    dimension(:, :, :), & ! Dimensions nPoints, nCol, nLev
      intent( in)  :: dtau         !  mean 0.67 micron optical depth of clouds in each model level
                                   !  NOTE:  this the cloud optical depth of only the cloudy part of the grid box, 
                                   !  it is not weighted with the 0 cloud optical depth of the clear part of the grid box
    real,    dimension(:, :),    &
      intent( in) :: pfull         !  pressure of full model levels (Pascals)
                                   !  pfull(npoints,1)    is    top level of model
                                   !  pfull(npoints,nlev) is bottom level of model
    integer,                    &
      intent( in) ::  top_height   !  1 = adjust top height using both a computed
                                   !  infrared brightness temperature and the visible
                                   !  optical depth to adjust cloud top pressure. Note
                                   !  that this calculation is most appropriate to compare
                                   !  to ISCCP data during sunlit hours.
                                   !  2 = do not adjust top height, that is cloud top
                                   !  pressure is the actual cloud top pressure
                                   !  in the model
                                   !  3 = adjust top height using only the computed
                                   !  infrared brightness temperature. Note that this
                                   !  calculation is most appropriate to compare to ISCCP
                                   !  IR only algortihm (i.e. you can compare to nighttime
                                   !  ISCCP data with this option)
    !
    ! Optional input arguments - for computing radiative cloud-top pressure 
    !   All variables except itrop are needed if top_height == 1 .or. top_height == 3
    !
    real,    dimension(:, :, :), optional, & ! Dimensions nPoints, nCol, nLev
     intent( in) ::  dem,   &    !  10.5 micron longwave emissivity of clouds in each model level.  
                                 !  Same note applies as in dtau.
                    dem_wv, &    ! Water vapor emissivity.  
                    at           ! Air temperature in each model level (K)
   real,   dimension(:, :),      optional, & ! Dimensions nPoints, nCol
     intent( in) :: skt, &       !  skin Temperature (K)
                    emsfc_lw     !  10.5 micron emissivity of surface (fraction)
   ! Users may supply their own tropopause height indicator
   integer, dimension(:),        optional, & ! Dimension nPoints
     intent( in) :: itrop        ! Index of the tropopause location in each column
   
    !     ------
    !     Output
    !     ------
    real, dimension(:, :), &   !  dimension nPoints, nCol
      intent(out) :: boxtau,         &   !  optical thickness in each column
                     boxptop             !  cloud top pressure (mb) in each sub-column
                  
    !     --------------------------------------------
    !     Local variables and parameters
    integer :: nPoints, nCol, nLev
    !     ------
    integer                           ::  j, ilev, ibox, icycle
    real,    dimension(size(dtau, 1), &
                       size(dtau, 2)) :: tau, pTop, tb
    real,    dimension(size(dtau, 1), &
                       size(dtau, 2), &
                       size(dtau, 3)) :: dem_local
    real,    dimension(size(dtau, 1), &
                       size(dtau, 3)) :: at_mean
    
    ! Variables for adjusting cloud top pressure based on TOA IR radiance
    real,    dimension(size(dtau, 1), &
                       size(dtau, 2)) :: emcld, fluxtop_clrsky, fluxtop, tauir, taumin, &
                                         fluxtopinit, transmax, btcmin
    logical, dimension(size(dtau, 1), &
                       size(dtau, 2)) :: intermedTrans
  
    ! Historgram quanitities 
    integer, dimension(size(dtau, 1), &
                       size(dtau, 2)) :: levMatch
    
    ! Tropopause 
    integer, dimension(size(dtau, 1)) :: itrop_local    ! index to tropopause level
    real,    dimension(size(dtau, 1)) :: pTrop, atTrop  ! tropopause pressure, temperature
                                   
    !     ------
    !     Local constants
    !     ------
    real,    parameter ::  VisTauToIrTauIce = 1./2.13, VisTauToIrTauLiquid = 1./2.56, &
                           assumedFreezingTb = 260.
    ! -----------------------------------------------------------------------------------
    nPoints = size(dtau, 1); nCol    = size(dtau, 2); nLev    = size(dtau, 3)
    
    ! Error checking
    call check_Lbound(dtau, "cloud optical depth",0.)
    if(present(dem)) &
      call check_bounds(dem, "cloud emissivity",  0., 1.)
    
    ! -----------------------------------------------------------------------------------
    ! The easy part - what's the total cloud optical depth in each column? 
    ! 
    tau(:, :) = sum(dtau, dim = 3)
    
    ! -----------------------------------------------------------------------------------
    ! If we're adjusting cloud-top heights using the VIS and/or IR radiances 
    !    we need to know 1) the TOA radiance, and 2) the position of the tropopause 
    !    (the latter is the pressure level to which really cold clouds are assigned)
    !
    if(top_height == 1 .or. top_height == 3) then
      dem_local(:, :, :) = dem_wv(:, :, :)
      !
      ! Clear sky calculation - only need one col per GCM grid point
      !
      call computeRadiance(dem_local, at, skt, emsfc_lw, fluxtop_clrsky)
      !
      ! Add contribution to emissivity from cloud
      !
      where(dem(:, :, :) > tiny(dem)) 
        dem_local(:, :, :) = 1. - ( (1. - dem_local(:, :, :)) * (1. -  dem(:, :, :)) )
      end where
      !
      ! And now the all-sky radiance
      !
      call computeRadiance(dem_local, at, skt, emsfc_lw, fluxtop)
    end if
    !     ---------------------------------------------------!
    ! Account for ISCCP procedures to determine cloud top temperature
  
    ! Account for partially transmitting cloud and recompute flux 
    !    ISCCP would see assuming a single layer cloud
    ! Note choice here of VisTauToIrTauIce = 1/2.13, as it is primarily ice
    !    clouds which have partial emissivity and need the adjustment 
    !
    ! If the cloud brightness temperature is greater than 260K,  the liquid cloud conversion
    !   factor is used.
    !
    ! This is discussed on pages 85-87 of the ISCCP D level documentation (Rossow et al. 1996)

    if (top_height == 1 .or. top_height == 3) then
      at_mean(:, :) = sum(at, dim = 2)/nCol
      ! Tropoause height needed if cloud tops are to be adjusted
      !
      if(present(itrop)) then
        itrop_local(:) = itrop(:)
      else 
        ! Tropopause height from mean temperature profile
        call diagnoseTropPressure(pfull, at_mean, itrop_local)
      end if
      do j = 1, nPoints
        ptrop(j)  =   pfull(j, itrop_local(j)) 
        attrop(j) = at_mean(j, itrop_local(j))
      end do

      if (top_height == 1) then
        !compute minimum brightness temperature and optical depth
        btcmin(:, :) = spread(TbToFlux(attrop(:) - 5.), dim = 2, nCopies = nCol)
        
        !note that the initial setting of tauir(:) is needed so that
        !tauir(:) has a realistic value 
        tauir(:, :) = tau(:, :) * VisTauToIrTauIce
        transmax(:, :) = (fluxtop(:,:)         - btcmin(:, :)) / &
                         (fluxtop_clrsky(:, :) - btcmin(:, :))
        taumin(:, :) = -log(max(min(transmax(:, :), 1. - spacing(1.)), tiny(transmax)))

        intermedTrans(:, :) = transmax(:, :) > tiny(transmax) .and. &
                              transmax(:, :) < 1. - spacing(1.)
        where (intermedTrans) 
          fluxtopinit(:, :) = fluxtop(:,:)
          tauir(:, :) = tau(:,:) * VisTauToIrTauIce
        end where
        
        tb = 0.
        do icycle=1,2
          where (tau(:,:) > tiny(tau) .and. intermedTrans) 
             emcld(:,:) = 1. - exp(-tauir(:, :))
             fluxtop(:,:) = fluxtopinit(:, :) - ( (1. - emcld(:,:)) * fluxtop_clrsky(:, :) )
             fluxtop(:,:) = max(tiny(fluxtop), fluxtop(:,:)/emcld(:,:))
             tb(:, :) = fluxToTb(fluxtop(:, :))
             where (tb(:, :) .gt. assumedFreezingTb) tauir(:, :) = tau(:,:) * VisTauToIrTauLiquid
          end where
        enddo
      end if  
      
      where(tau(:, :) > tiny(tau)) 
        tb(:, :) = fluxToTb(fluxtop(:, :))
      elsewhere
        ! Clear sky brightness temperature
        tb(:,:) = fluxToTb(fluxtop_clrsky(:, :))
      end where 
      
      if(top_height == 1) then
        ! Adjust the brightness temperature and optical depth of very thin clouds
        where (tau(:, :) > tiny(tau) .and. tauir(:, :) < taumin(:, :)) 
           tb(:, :) = spread(attrop(:) - 5., dim= 2, nCopies = nCol)
          tau(:, :) = taumin(:, :) / VisTauToIrTauIce
        end where
      end if
                
    end if
    ! -----------------------------------------------------------------------------------
    ! Determine cloud-top pressure. Three choices:
    !     Radiatively determined cloud top pressure (top_height = 1 for VIS/R;  3 for IR)
    !     Physical cloud top pressure (top_height = 2)
    
    if (top_height .eq. 1 .or. top_height .eq. 3) then  
      !segregate according to optical thickness (??? - RP)
      levmatch(:, :) = 0
      ! find lowest (?) level whose temperature
      ! most closely matches brightness temperature
      !
      do ilev = 1, nlev-1
        where((at(:, :, ilev) >= tb(:, :) .and. at(:, :, ilev+1) < tb(:, :)) .or. &
              (at(:, :, ilev) <= tb(:, :) .and. at(:, :, ilev+1) > tb(:, :)) )
          where(abs(at(:, :, ilev) - tb(:,:)) < abs(at(:, :, ilev+1) - tb(:,:)))
            levmatch(:, :) = ilev
          elsewhere 
            levmatch(:, :) = ilev + 1
          end where
        end where 
      end do

        ! If we've found a matching level use it, otherwise use the boundary value
        !  
      do ibox = 1, ncol
        do j = 1, nPoints
          if(levmatch(j, ibox) >= 1) then
            ptop(j, ibox) = pfull(j, levmatch(j, ibox))
          else if (tb(j, ibox) < minval(at(j, ibox, :))) then
            levmatch(j, ibox) = itrop_local(j)
            ptop(j, ibox) = ptrop(j)
          else if (tb(j, ibox) > maxval(at(j, ibox, :))) then
            levmatch(j, ibox) = nLev
            ptop(j, ibox) = pFull(j, nLev)
          end if
        end do
      end do
    else !  When top_height .eq. 2
      ! The highest cloud top (clouds being where tau > 0). 
      ptop(:, :) = 0.
      do ibox = 1, nCol
        do ilev = 1, nlev
          where(ptop(:, ibox) <= 0. .and. dtau(:, ibox, ilev) > 0.)
            ptop(:, ibox) = pfull(:, ilev)
            levmatch(:,ibox) = ilev
          end where
        end do
      end do
    end if                            
          
    ! No matter how cloud top pressure is determined, 
    !   pTop and levmatch are 0 when the column is clear.   
    ! 
    where (tau(:,:) <= tiny(tau)) 
      ptop(:,:) = 0.
      levmatch(:,:) = 0      
    end where

    !     ---------------------------------------------------!
    ! Output
    !    
    boxptop(:, :) = ptop(:, :) / 100. 
    boxtau(:, :)  = tau(:, :)
    
  end subroutine icarus_CRM
  ! -------------------------------------------------------------------
  subroutine diagnoseTropPressure(pfull, at, itrop)
    real,    dimension(:, :), intent( in) :: pFull, at
    integer, dimension(:),    intent(out) :: itrop
    
    integer                         :: nPoints, nLev  
    real, dimension(size(pFull, 1)) :: attropmin
    integer                         :: ilev

    nPoints = size(pFull, 1); nLev = size(pFull, 2)
    attropmin(:) = 400.
    itrop(:)     = 1
  
    do  ilev=1,nlev
      where(pfull(:, ilev) < 40000. .and. pfull(:, ilev) > 5000. .and. &
            at(:, ilev) < attropmin(:)) 
        attropmin(:) = at(:, ilev)
        itrop(:)=ilev
      end where
    end do
  end subroutine diagnoseTropPressure 
  ! -------------------------------------------------------------------      
  subroutine computeRadiance_1D(dem, at, skt, emsfc_lw, TOAradiance)
    real,    dimension(:, :), &  ! Dimensions nPoints, nLev
      intent( in) :: dem,     &  !   10.5 micron emissivity of water vapor
                     at          !   air temperature
    real,   dimension(:),     &  ! Dimension nPoints
      intent( in) :: skt,     &  !   skin Temperature (K)
                     emsfc_lw    !   10.5 micron emissivity of surface (fraction) 
    real, dimension(:),       &  ! Dimension nPoint, nCol
      intent(out) :: TOAradiance !   10.5 micron nadir radiance at TOA
    
  !     ------
  !     Local variables and parameters
    integer :: nPoints, nLev
    !     ------
    integer ::  ilev
    real,    dimension(size(dem, 1)) :: trans_layers_above
    real,    dimension(size(dem, 1)) :: bb
                                                          
    !----------------------------------------------------------------------
    ! Computes radiance at TOA from an emitting/absorbing atmosphere
    !     TOAradiance is the 10.5 micron radiance at the top of the
    !              atmosphere
    !     trans_layers_above is the total transmissivity in the layers
    !             above the current layer
    !----------------------------------------------------------------------
  
    !initialize variables
    nPoints = size(dem, 1); nLev    = size(dem, 2)
    TOAradiance(:) = 0.; trans_layers_above(:) = 1.
  
    do ilev=1,nlev
      ! Black body emission at temperature of the layer
      bb(:) = TbToFlux(at(:,ilev))
      !bb(j)= 5.67e-8*at(j,ilev)**4
  
      ! increase TOA flux by flux emitted from layer
      ! times total transmittance in layers above
      TOAradiance(:) = TOAradiance(:) + dem(:, ilev) * trans_layers_above(:) * bb(:)
        
      ! update trans_layers_above with transmissivity
      ! from this layer for next time around loop
      trans_layers_above(:) = trans_layers_above(:) * (1. - dem(:, ilev))
    enddo ! ilev
  
    !surface emission
    bb(:) = TbToFlux(skt(:)) 
    !bb(:)=5.67e-8*skt(:)**4
  
    !add in surface emission
    TOAradiance(:) = TOAradiance(:) + trans_layers_above(:) * emsfc_lw(:) * bb(:)
                       
  end subroutine computeRadiance_1D
 ! -------------------------------------------------------------------      
  subroutine computeRadiance_2D(dem, at, skt, emsfc_lw, TOAradiance)
    real,    dimension(:, :, :), &  ! Dimensions nPoints, nCol, nLev
      intent( in) :: dem            !   10.5 micron emissivity of water vapor
   real,    dimension(:, :),     &  ! Dimensions nPoints, nLev
      intent( in) :: at             !   air temperature
    real,   dimension(:),        &  ! Dimension nPoints
      intent( in) :: skt,        &  !   skin Temperature (K)
                     emsfc_lw       !   10.5 micron emissivity of surface (fraction) 
    real, dimension(:, :),       &  ! Dimension nPoint, nCol
      intent(out) :: TOAradiance    !   10.5 micron nadir radiance at TOA
    
  !     ------
  !     Local variables and parameters
    integer :: nPoints, nCol, nLev
    !     ------
    integer ::  ilev
    real,    dimension(size(dem, 1), &
                       size(dem, 2)) :: trans_layers_above
    real,    dimension(size(dem, 1)) :: bb
 
    !----------------------------------------------------------------------
    ! Computes radiance at TOA from an emitting/absorbing atmosphere
    !     TOAradiance is the 10.5 micron radiance at the top of the
    !              atmosphere
    !     trans_layers_above is the total transmissivity in the layers
    !             above the current layer
    !----------------------------------------------------------------------
  
    !initialize variables
    nPoints = size(dem, 1); nCol    = size(dem, 2); nLev    = size(dem, 3)
    TOAradiance(:, :) = 0.; trans_layers_above(:, :) = 1.
  
    do ilev=1,nlev
      ! Black body emission at temperature of the layer
      bb(:) = TbToFlux(at(:,ilev))
      !bb(j)= 5.67e-8*at(j,ilev)**4
  
      ! increase TOA flux by flux emitted from layer
      ! times total transmittance in layers above
      TOAradiance(:,:) = TOAradiance(:,:) + &
                         dem(:,:, ilev) * trans_layers_above(:,:) * spread(bb(:), dim = 2, nCopies = nCol) 
        
      ! update trans_layers_above with transmissivity
      ! from this layer for next time around loop
      trans_layers_above(:, :) = trans_layers_above(:, :) * (1. - dem(:, :, ilev))
    
    enddo ! ilev
  
  
    !surface emission
    bb(:) = TbToFlux(skt(:)) 
    !bb(:)=5.67e-8*skt(:)**4
  
    !add in surface emission
    TOAradiance(:,:) = TOAradiance(:,:) +  &
                       trans_layers_above(:,:) * spread(emsfc_lw(:) * bb(:), dim = 2, nCopies = nCol)
  
  end subroutine computeRadiance_2D
 ! -------------------------------------------------------------------      
  subroutine computeRadiance_3D(dem, at, skt, emsfc_lw, TOAradiance)
    real,    dimension(:, :, :), &  ! Dimensions nPoints, nCol, nLev
      intent( in) :: dem,        &  !   10.5 micron emissivity of water vapor
                     at             !   air temperature
    real,   dimension(:, :),     &  ! Dimension nPoints, nCol
      intent( in) :: skt,        &  !   skin Temperature (K)
                     emsfc_lw       !   10.5 micron emissivity of surface (fraction) 
    real, dimension(:, :),       &  ! Dimension nPoint, nCol
      intent(out) :: TOAradiance    !   10.5 micron nadir radiance at TOA
    
  !     ------
  !     Local variables and parameters
    integer :: nPoints, nCol, nLev
    !     ------
    integer ::  ilev
    real,    dimension(size(dem, 1),                        &
                       size(dem, 2)) :: trans_layers_above, &
                                        bb
 
    !----------------------------------------------------------------------
    ! Computes radiance at TOA from an emitting/absorbing atmosphere
    !     TOAradiance is the 10.5 micron radiance at the top of the
    !              atmosphere
    !     trans_layers_above is the total transmissivity in the layers
    !             above the current layer
    !----------------------------------------------------------------------
  
    !initialize variables
    nPoints = size(dem, 1); nCol    = size(dem, 2); nLev    = size(dem, 3)
    TOAradiance(:, :) = 0.; trans_layers_above(:, :) = 1.
  
    do ilev=1,nlev
      ! Black body emission at temperature of the layer
      bb(:, :) = TbToFlux(at(:,:, ilev))
  
      ! increase TOA flux by flux emitted from layer
      ! times total transmittance in layers above
      TOAradiance(:,:) = TOAradiance(:,:) + &
                         dem(:,:, ilev) * trans_layers_above(:,:) * bb(:, :)
        
      ! update trans_layers_above with transmissivity
      ! from this layer for next time around loop
      trans_layers_above(:, :) = trans_layers_above(:, :) * (1. - dem(:, :, ilev))
    enddo ! ilev
  
    !surface emission
    bb(:, :) = TbToFlux(skt(:, :)) 
  
    !add in surface emission
    TOAradiance(:,:) = TOAradiance(:,:) +  &
                       trans_layers_above(:,:) * emsfc_lw(:, :) * bb(:, :)
  
  end subroutine computeRadiance_3D
 ! -------------------------------------------------------------------      
   subroutine computeWaterVaporEmissivity(pfull, phalf, qv, at, dem_wv)
     real, dimension(:, :), & ! nPoints, nLev
       intent( in) :: pFull, pHalf, qv, at
     real, dimension(:, :), & ! nPoints, nLev
       intent(out) :: dem_wv
     
     ! Local variables
     integer :: nPoints,  nLev
     integer :: iLev
     real, dimension(size(dem_wv, 1)) :: press, dpress, atmden, rvh20, wk, rhoave, rh20s, &
                                         rfrgn, tmpexp, tauwv
    
     real, parameter :: wtmair = 28.9644, wtmh20 = 18.01534, Navo = 6.023E+23, grav = 9.806650E+02, &
                        pstd = 1.013250E+06, t0 = 296.
    ! -------------------------------------------
    nPoints = size(pFull, 1); nLev    = size(pFull, 2)

    !compute water vapor continuum emissivity
    !this treatment follows Schwarkzopf and Ramasamy
    !JGR 1999,vol 104, pages 9467-9499.
    !the emissivity is calculated at a wavenumber of 955 cm-1, 
    !or 10.47 microns 
    do ilev=1,nlev
      ! press and dpress are dyne/cm2 = Pascals *10
      press(:) = pfull(:,ilev)*10.
      dpress(:) = (phalf(:,ilev+1)-phalf(:,ilev))*10
      !atmden = g/cm2 = kg/m2 / 10 
      atmden(:) = dpress(:)/grav
      rvh20(:) = qv(:,ilev)*wtmair/wtmh20
      wk(:) = rvh20(:)*Navo*atmden(:)/wtmair
      rhoave(:) = (press(:)/pstd)*(t0/at(:,ilev))
      rh20s(:) = rvh20(:)*rhoave(:)
      rfrgn(:) = rhoave(:)-rh20s(:)
      tmpexp(:) = exp(-0.02*(at(:,ilev)-t0))
      tauwv(:) = wk(:)*1.e-20*((0.0224697*rh20s(:)*tmpexp(:)) + (3.41817e-7*rfrgn(:)) )*0.98
      dem_wv(:,ilev) = 1. - exp( -1. * tauwv(:))
    end do
  end subroutine computeWaterVaporEmissivity 
 !     ---------------------------------------------------!
 ! Functions to compute brightness temperature given a 11 micron radiance
 !
 ! -------------------------------------------------------------------      
  function fluxToTb_1D(flux) result(tb)
    real, dimension(:),       intent( in) :: flux
    real, dimension(size(flux))           :: tb
    
    tb(:) = 1307.27 / log(1. + (1./flux(:)))
  end function fluxToTb_1D
 !     ---------------------------------------------------!
   function fluxToTb_2D(flux) result(tb)
     real, dimension(:, :),    intent( in) :: flux
     real, dimension(size(flux, 1), &
                     size(flux, 2))        :: tb
     
     tb(:, :) = 1307.27 / log(1. + (1./flux(:, :)))
   end function fluxToTb_2D
 !     ---------------------------------------------------!
   function fluxToTb_3D(flux) result(tb)
     real, dimension(:, :, :), intent( in) :: flux
     real, dimension(size(flux, 1), &
                     size(flux, 2), &
                     size(flux, 3))        :: tb
     
     tb(:, :, :) = 1307.27 / log(1. + (1./flux(:, :, :)))
   end function fluxToTb_3D
  !     ---------------------------------------------------!
  ! Functions to compute 11 micron radiance given a brightness temperature
  !
  ! -------------------------------------------------------------------      
   function TbToFlux_1D(tb) result(flux)
     real, dimension(:),       intent( in) :: tb
     real, dimension(size(tb))             :: flux
     
     flux(:) = 1. / ( exp(1307.27/tb(:)) - 1. )
   end function TbToFlux_1D
 !     ---------------------------------------------------!
   function TbToFlux_2D(tb) result(flux)
     real, dimension(:, :),    intent( in) :: tb
     real, dimension(size(tb, 1), &  
                     size(tb, 2))          :: flux
     
     flux(:, :) = 1. / ( exp(1307.27/tb(:, :)) - 1. )
   end function TbToFlux_2D
 !     ---------------------------------------------------!
   function TbToFlux_3D(tb) result(flux)
     real, dimension(:, :, :), intent( in) :: tb
     real, dimension(size(tb, 1), &
                     size(tb, 2), &
                     size(tb, 3))          :: flux
     
     flux(:, :, :) = 1. / ( exp(1307.27/tb(:, :, :)) - 1. )
   end function TbToFlux_3D
 !     ---------------------------------------------------!
  ! -------------------------------------------------------------------      
  subroutine check_bounds(array, name, minAllowed, maxAllowed)
    implicit none
    ! Input variables
    real, dimension(:, :, :), intent( in) :: array
    character(len = *),       intent( in) :: name
    real,                     intent( in) :: minAllowed, maxAllowed
    
    ! ---------------------
    if(any(array(:, :, :) < minAllowed .or. array(:, :, :) > maxAllowed)) then
      print *, "Values in arrary ", trim(name), " out of bounds" 
      stop 
    end if
  end subroutine check_bounds
  ! -------------------------------------------------------------------      
  subroutine check_Lbound(array, name, minAllowed)
    implicit none
    ! Input variables
    real, dimension(:, :, :), intent( in) :: array
    character(len = *),       intent( in) :: name
    real,                     intent( in) :: minAllowed
    
    ! ---------------------
    if(any(array(:, :, :) < minAllowed )) then
      print *, "Values in arrary ", trim(name), " less than minimum allowed" 
      stop 
    end if
  end subroutine check_lBound
  ! -------------------------------------------------------------------      
end module isccpSimulator
