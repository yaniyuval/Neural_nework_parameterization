module isccpTables
  use grid, only : case, doisccp
  implicit none
  private
  public :: tautab, invtau, isccp_tables_init
  
  REAL, dimension(0:255) :: tautab
  integer, dimension(-20:45000) :: invtau

CONTAINS

  subroutine isccp_tables_init
    if(.not.doisccp) return
    open (1,file='./'//trim(case)//'/isccp',status="old",form="formatted")
    read(1,*) tautab(0:255)
    read(1,*) invtau(-20:45000)
    close(1)
  end subroutine isccp_tables_init

end module isccpTables


module isccpHistograms
  implicit none
  integer, parameter :: numIsccpPressureIntervals     = 7, &
                        numIsccpOpticalDepthIntervals = 7
  real,    parameter :: minIsccpOpticalDepth = 0.3
  real, dimension(numIsccpPressureIntervals + 1),      parameter ::       &
           isccpPressureBinEdges = (/ tiny(minIsccpOpticalDepth),         &
                                      180., 310., 440., 560., 680., 800., &
                                      huge(minIsccpOpticalDepth) /) 
  real, dimension(numIsccpOpticalDepthIntervals + 1),  parameter ::                    &
          isccpOpticalDepthBinEdges = (/ tiny(minIsccpOpticalDepth),                   &
                                        minIsccpOpticalDepth, 1.3, 3.6, 9.4, 23., 60., & 
                                        huge(minIsccpOpticalDepth) /) 
contains
! ------------------------------------------------------ 
  function computeIsccpJointHistograms(tau, ptop, sunlit) result(isccpJointHistogram)
    ! Dimensions nGridCell, nSubColumn 
    real, dimension(:, :),    intent( in) :: tau, ptop
    ! Dimensions nGridCell 
    integer, dimension(:),    intent( in) :: sunlit
    ! Dimensions nGridCells, nTauLevels (7), nPressureLevels(7)
    real, dimension(size(tau, 1), numIsccpOpticalDepthIntervals, numIsccpPressureIntervals) &
                                          :: isccpJointHistogram
    
    ! Local variables
    integer                     :: i, j
    logical, dimension(size(tau, 1), size(tau, 2)) &
                                :: box_cloudy
    ! --------------------
    box_cloudy(:, :) = tau(:, :) > tiny(tau) .and. ptop(:, :) > 0 
    
    !
    ! Construct the histogram
    !
    do j = 1, numIsccpPressureIntervals 
      do i = 1, numIsccpOpticalDepthIntervals 
        isccpJointHistogram(:, i, j) = count(box_cloudy(:, :) .and.                               &
                                             tau(:, :)  >= isccpOpticalDepthBinEdges(i)     .and. &
                                             tau(:, :)  <  isccpOpticalDepthBinEdges(i + 1) .and. &
                                             pTop(:, :) >= isccpPressureBinEdges(j)         .and. &
                                             pTop(:, :) <  isccpPressureBinEdges(j + 1), dim = 2)
      end do
    end do
    isccpJointHistogram(:, :, :)  = isccpJointHistogram(:, :, :)/size(tau, 2)
  end function computeIsccpJointHistograms
! ------------------------------------------------------ 
  subroutine computeIsccpMeanProperties(tau, ptop, sunlit, top_height, &
                                        totalcldarea, meanptop, meantaucld, &
                         lowcldarea, midcldarea, hghcldarea)
    use isccpTables
    real, dimension(:, :),    intent( in) :: tau, ptop
    integer, dimension(:),    intent( in) :: sunlit
    integer,                  intent( in) :: top_height
    ! Dimensions nGridCell, nSubColumn 
    real, dimension(:), intent(out) :: totalcldarea, meanptop, meantaucld
    real, dimension(:), intent(out) :: lowcldarea, midcldarea, hghcldarea
    ! Dimensions nGridCells
    
    !     Compute grid box mean cloud top pressure and
    !     optical thickness.  The mean cloud top pressure and
    !     optical thickness are averages over the cloudy 
    !     area only. The mean cloud top pressure is a linear
    !     average of the cloud top pressures.  The mean cloud
    !     optical thickness is computed by converting optical
    !     thickness to an albedo, averaging in albedo units,
    !     then converting the average albedo back to a mean
    !     optical thickness.  
    
    integer                  :: nPoints, nCol
    integer                  :: j
    real, dimension(size(tau, 1)) &
                             :: meanalbedocld
    logical, dimension(size(tau, 1), size(tau, 2)) &
                             :: box_cloudy, box_cloudy_low, box_cloudy_mid, box_cloudy_hgh
                             
    ! Optical depth forward and inverse lookup tables for 
    !   averaging optical depth to get the right albedo
                 
    ! ----------------------------------------------------------------------------------    
    nPoints = size(tau, 1); nCol = size(tau, 2)
! Marat: Changed to be able to compare to observations which see tau > some threshold
    box_cloudy(:, :) = tau(:, :) > minIsccpOpticalDepth .and. ptop(:, :) > 0
    box_cloudy_low(:, :) = tau(:, :) > minIsccpOpticalDepth .and. ptop(:, :) >= 700.
    box_cloudy_hgh(:, :) = tau(:, :) > minIsccpOpticalDepth .and. ptop(:, :) <= 400. 
    box_cloudy_mid(:, :) = tau(:, :) > minIsccpOpticalDepth .and. ptop(:, :) < 700. .and. ptop(:, :) > 400.
!    box_cloudy(:, :) = tau(:, :) > tiny(tau) .and. ptop(:, :) > 0
    totalcldarea(:) = count(box_cloudy, dim = 2)/real(nCol)
    lowcldarea(:) = count(box_cloudy_low, dim = 2)/real(nCol)
    midcldarea(:) = count(box_cloudy_mid, dim = 2)/real(nCol)
    hghcldarea(:) = count(box_cloudy_hgh, dim = 2)/real(nCol)
    !
    ! Mean cloudy-sky albedo
    !
    meanalbedocld(:) = 0.
    ! Accumulate albedo over cloudy columns... 
    do j = 1, nCol
      where(box_cloudy(:, j) .and. sunlit(:) == 1) &
        meanalbedocld(:) = meanalbedocld(:) + real(invtau(min(nint(100.*tau(:, j)),45000)))
    end do 
    ! ... and then find the average (sunlit columns only). 
    where(sunlit(:) == 1)
      meanalbedocld(:) = meanalbedocld(:)/ real(max(1,count(box_cloudy(:, :), dim = 2)))
    elsewhere
      meanalbedocld(:) = 0. 
    end where
    
    !
    ! Mean tau (ISCCP definition: the tau that provides the mean albedo)
    !
    where(sunlit == 1)
! it will be computed later after MPI task exchange. - Marat K.
!      meantaucld(:)    = tautab(min(255,max(1,nint(meanalbedocld(:)))))
      meantaucld(:)    = meanalbedocld(:)
    elsewhere
      meantaucld(:)    = 0.
    end where
    
    !
    ! Mean cloud-top height
    !
    where(totalcldarea(:) > tiny(totalcldarea))               &
      meanptop(:) = real(sum(ptop(:, :), mask = box_cloudy(:, :), dim = 2)) / &
                    real(max(1,count(box_cloudy(:, :), dim = 2)))
    ! If we're using a combined VIS/IR determination of cloud-top height
    !   we can't see the tops without sunlit
    if(top_height == 1) &
      where(sunlit(:) /= 1) meanptop(:) = 0.
       
  end subroutine computeIsccpMeanProperties
! ------------------------------------------------------ 
end module isccpHistograms
  
