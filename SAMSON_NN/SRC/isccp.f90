! Interface and output handler for the isccp-simulator
! Marat Khairoutdinov, 2004


module isccp

   use isccpSimulator
   use isccpHistograms
   use params, only : pres0, ggr, coszrs
   use vars

   implicit none

 PRIVATE

   integer, parameter :: ntau = numIsccpOpticalDepthIntervals
   integer, parameter :: npres = numIsccpPressureIntervals

   integer, parameter :: crm_nx = nx, crm_ny = ny, crm_nz = nzm
   integer, parameter :: pver = crm_nz, pverp = pver+1
   real, parameter :: fillvalue = -1.     ! fill value for reduced grid and others

   real :: landm(1), tlayer(crm_nz), qlayer(crm_nz)
   real :: clwp(crm_nz), fracice(crm_nz), emis(crm_nz)
   real :: tmp(8)

   real fq_isccp(ntau*npres)  !accumulated fq_isccp

   real totalcldarea !  the fraction of model grid box columns
      !  with cloud somewhere in them.  This should
               !  equal the sum over all entries of fq_isccp
      ! The following three means are averages over the cloudy areas only.  If no
      ! clouds are in grid box all three quantities should equal zero.
   real meantaucld !  mean optical thickness (dimensionless)
                                        !  linear averaging in albedo performed.
   real meanptop !  mean cloud top pressure (mb) - linear averaging
                                        !  in cloud top pressure.
   integer nisccp ! number of collected samples. Greater than 0 only if at least one
                     ! day-time sample was obtained
   integer nsun ! number of day-time samples.

   real lowcldarea, midcldarea, hghcldarea ! low, mid and high cloud fractions
public :: isccp_get, isccp_write, isccp_zero

CONTAINS

subroutine isccp_zero()
    fq_isccp = 0.
    totalcldarea = 0.
    lowcldarea = 0.
    midcldarea = 0.
    hghcldarea = 0.
    meantaucld = 0.
    meanptop = 0.
    nisccp = 0
    nsun = 0
end subroutine isccp_zero
 


subroutine isccp_get ()

   implicit none

!   local:

   real fq_isccp_s1(ntau*npres)  
   real totalcldarea1 
   real meantaucld1 
   real meanptop1 
   real lowcldarea1 
   real midcldarea1 
   real hghcldarea1 

   if(.not.doisccp) return

   call crm_isccp (fq_isccp_s1, totalcldarea1, meantaucld1, meanptop1, lowcldarea1, midcldarea1, hghcldarea1)

   if(fq_isccp_s1(1).ne.fillvalue) then
    fq_isccp = fq_isccp + fq_isccp_s1
    totalcldarea = totalcldarea + totalcldarea1
    lowcldarea = lowcldarea + lowcldarea1
    midcldarea = midcldarea + midcldarea1
    hghcldarea = hghcldarea + hghcldarea1
    nsun = nsun + 1
   end if
   if(meantaucld1.ne.fillvalue) then
    meantaucld = meantaucld + meantaucld1
    meanptop = meanptop + meanptop1
    nisccp = nisccp + 1
   end if

end subroutine isccp_get


subroutine crm_isccp (fq_isccp_s1, totalcldarea, meantaucld, meanptop, lowcldarea, midcldarea, hghcldarea)

   real    emsfc_lw   !emsfc_lw  longwave emissivity of surface at 10.5 microns
   PARAMETER(emsfc_lw=0.99)
!bloss    real, intent(in) :: pres(crm_nz)  !
!bloss    real, intent(in) :: q   (crm_nx, crm_ny, crm_nz)  ! total water
!bloss    real, intent(in) :: qn  (crm_nx, crm_ny, crm_nz)  ! cloud water/ice
!bloss    real, intent(in) :: t   (crm_nx, crm_ny, crm_nz)  ! temperature
!bloss    real, intent(in) :: tabs0(crm_nz) ! domain-mean temp
!bloss    real, intent(in) :: qv0(crm_nz) ! domain-mean water vapor
  

   real, intent(out) :: fq_isccp_s1(ntau*npres)  !accumulated fq_isccp
   real, intent(out) :: totalcldarea 
   real, intent(out) :: meantaucld 
   real, intent(out) :: meanptop 
   real, intent(out) :: lowcldarea 
   real, intent(out) :: midcldarea 
   real, intent(out) :: hghcldarea 
  
! local
   integer ncolmax
   parameter (ncolmax=crm_nx*crm_ny)
   real pmid(1,pver)  ! Pressure at middle level
   real pint(1,pverp)  ! Pressure at half-levels
   real tau(1,ncolmax, pver)
   real q_down(1,ncolmax, pver)  !
   real t_down(1,ncolmax, pver)  !
   real q0_down(1, pver)  !
   real t0_down(1, pver)  !
   real cld_down  (1,ncolmax, pver)  !cloud cover
   real  rel(pver),rei(pver),irei(pver)
   real  massl, cliqwp,cicewp, qnn, qcc, qii, qss
   real emis_down(1,ncolmax, pver)  ! Cloud longwave emissivity
   real emiswv_down(1, pver)  ! water vapor longwave emissivity
   real cltot  ! Total cloud amount in grid cell
   real cldcol ! Total cloud amount in crm column

   integer top_height
   integer i,j,k,m,ii
   real ts(1,ncolmax)    ! surface temperature
   REAL boxtau(1,ncolmax)
   REAL boxptop(1,ncolmax)
   REAL fq_isccp0(1,ntau,npres)
   real totalcldarea0(1)
   real meantaucld0(1)
   real meanptop0(1)
   integer sunlit(1)
   real lowcldarea0(1)
   real midcldarea0(1)
   real hghcldarea0(1)

   real abarl         ! A coefficient for extinction optical depth
   real bbarl         ! b coefficient for extinction optical depth
   real abari         ! A coefficient for extinction optical depth
   real bbari         ! b coefficient for extinction optical depth
   real cldmin    ! the value taken from radcswmx.F90, note: cldmin much less than cldmin from what on cldnrh
   real cldeps    ! the value taken from radcswmx.F90
   integer itaunpres,it,ip
   real, parameter :: qthresh = 1.e-9 ! minimum cld mixing ratio included
   real, parameter :: eps = 1.e-20

   cldmin = 1.0e-20
!   cldmin = 1.0e-80_r8
   cldeps = 0.0

   abarl = 2.817e-02
   bbarl = 1.305
   abari = 3.448e-03
   bbari = 2.431

   !compute the optical depth
    do m=crm_nz+1,pver
       k = pver-m+1
       cld_down(1,:ncolmax,k)=0.
       emis_down(1,:ncolmax,k)=0.
       tau(1,:ncolmax,k)=0.
       t0_down(1,k) = tabs0(crm_nz)
       q0_down(1,k) = qv0(crm_nz)
       ii=0
       do i=1,crm_nx
       do j=1,crm_ny
          ii=ii+1
          t_down(1,ii,k)=tabs(i,j,crm_nz)
          q_down(1,ii,k)=q(i,j,crm_nz)-qn(i,j,crm_nz)
       end do
       end do
    enddo

    do k=1,crm_nz
       m=crm_nz-k+1
       pmid(1,m) = pres(k)*100.
    end do
    do m=2,crm_nz
       pint(1,m)=0.5*(pmid(1,m)+pmid(1,m-1))
    end do
    pint(1,1)=pmid(1,1)-0.5*(pmid(1,2)-pmid(1,1))
    pint(1,crm_nz+1)=pres0*100.

! set t0_down, q0_down
    do m=1,crm_nz
       k = pver-m+1
       t0_down(1,k) = tabs0(m)
       q0_down(1,k) = qv0(m)
    end do

! Get effective radii:
    landm = 0.
    if (.not.OCEAN) landm = 1.
    tmp = 0.
    tmp(2) = 100.*pres0
    ! compute effective radii using horizontally-avg temperature t0_down
    call cldefrint(1,1,tmp(1),t0_down,rel,rei, tmp(2),pmid,landm,tmp(3),tmp(4))  ! CAM3 interface  

    do m=1,crm_nz
       k = pver-m+1
       ii=0
       irei(k) = 1./rei(k)
       do i=1,crm_nx
       do j=1,crm_ny
       ii=ii+1
          t_down(1,ii,k)=tabs(i,j,m)
          q_down(1,ii,k)=q(i,j,m)-qn(i,j,m)
          qnn = qn(i,j,m)
          if(.not.dokruegermicro) then
          qcc = qn(i,j,m)*omegan(t_down(1,ii,k))
          qii = qn(i,j,m)*(1.-omegan(t_down(1,ii,k)))
          qss = qp(i,j,m)*(1.-omegap(t_down(1,ii,k)))*(1.-omegag(t_down(1,ii,k)))
          elseif(dokruegermicro) then
             qcc = qc(i,j,m)
             qii = qi(i,j,m)
             qss = qs(i,j,m)
          end if

          ! initialize cloud fraction, optical depth and emissivity to zero
          cld_down(1,ii,k) = 0.
          tau(1,ii,k) = 0.0
          emis_down(1,ii,k) = 0.

          ! If cloud is present, set non-zero values for cld, tau, emis
          if (dokruegereffectiveradii) then
             if (qcc+qii+qss.gt.qthresh) then
                ! re-define effective radii
                rel(k) = 10.
                rei(k) = 25.*(qii+qss+eps)/(qii+qss/3.+eps)
                irei(k) = 1./rei(k)
                ! compute cloud thickness/emissivity
                cld_down(1,ii,k) = 0.99
                massl = 1000.*(pint(1,k+1)-pint(1,k))/ggr
                cliqwp = qcc*massl
                cicewp = (qii+qss)*massl
                tau(1,ii,k) = (abarl + bbarl/rel(k)) * cliqwp + &
                     (abari + bbari/rei(k)) * cicewp
                emis_down(1,ii,k) =1.-exp(-(0.15*cliqwp &
                     +1.66*(0.005+irei(k))*cicewp))
             end if
          elseif(qnn.gt.qthresh) then                
             cld_down(1,ii,k) = 0.99
             massl = 1000.*(pint(1,k+1)-pint(1,k))/ggr
             cliqwp = qcc*massl
             cicewp = qii*massl
             tau(1,ii,k) = (abarl + bbarl/rel(k)) * cliqwp + &
                           (abari + bbari/rei(k)) * cicewp
             emis_down(1,ii,k) =1.-exp(-(0.15*cliqwp &
                                        +1.66*(0.005+irei(k))*cicewp))
          end if
       end do
       end do
    end do
    cltot=0.
    do ii=1,ncolmax
       cldcol=0.
       do m=1,crm_nz
          cldcol=max(cldcol,cld_down(1,ii,m))
       end do
       cltot=cltot+cldcol
    end do
    cltot=cltot/(ncolmax)

    top_height = 1

! get surface temperature from sst
   ii = 0
   do i=1,crm_nx
      do j=1,crm_ny
         ii=ii+1
         ts(1,ii) = sstxy(i,j)
      end do
   end do
!
!JR Add coszrs logic to only call the simulator during daytime
!
      if(coszrs > 0.) then        !if cloudy
         if(cltot >= cldmin) then        !if cloudy
           sunlit(1) = 1

         ! Compute water vapor emissivity  
            call computeWaterVaporEmissivity(pmid, pint, q0_down, t0_down, emiswv_down)

            call icarus_CRM(tau, pmid, top_height, emis_down, &
                  spread(emiswv_down, dim = 2, nCopies = ncolmax), &
                  t_down, &
                  ts, &
                  spread((/ (emsfc_lw, j = 1, 1) /), dim = 2, nCopies = ncolmax), &
                  boxtau = boxtau, boxptop = boxptop)

          ! And compute histograms
            fq_isccp0(:, :, :) = computeIsccpJointHistograms(boxtau, boxptop, sunlit)
            call computeIsccpMeanProperties(boxTau, boxPtop, sunlit, top_height, &
                                    totalcldarea0, meanptop0, meantaucld0, &
                                    lowcldarea0, midcldarea0, hghcldarea0)

!           save standard ISCCP type of 7x7 clouds
            do ip=1,npres
               do it=1,ntau
                  itaunpres = (ip-1)*ntau+it
                  fq_isccp_s1(itaunpres) = fq_isccp0(1,it,ip)
               end do
            end do
            totalcldarea = totalcldarea0(1)
            if(totalcldarea.gt.0.) then
              meanptop = meanptop0(1)
              meantaucld = meantaucld0(1)
            else
              meanptop      = fillvalue
              meantaucld    = fillvalue
            endif
            lowcldarea = lowcldarea0(1)
            midcldarea = midcldarea0(1)
            hghcldarea = hghcldarea0(1)
         else                             !cloud free in the (daytime) grid box
            fq_isccp_s1(:ntau*npres) = 0.
            totalcldarea  = 0.
            meanptop      = fillvalue
            meantaucld    = fillvalue
            lowcldarea  = 0.
            midcldarea  = 0.
            hghcldarea  = 0.
         endif
      else                                ! nighttime
         fq_isccp_s1(1:ntau*npres) = fillvalue
         totalcldarea  = fillvalue
         meanptop      = fillvalue
         meantaucld    = fillvalue
         lowcldarea  = fillvalue
         midcldarea  = fillvalue
         hghcldarea  = fillvalue
      end if

   return
   end subroutine crm_isccp


subroutine isccp_write()

   use grid
   use isccpTables, only : tautab
   use isccpHistograms, only:  isccpPressureBinEdges, isccpPressureBinEdges
   implicit none
   integer, parameter :: buf_len = ntau*npres + 6
   real coef, buf(buf_len), buf1(buf_len)
   integer, parameter :: ntape = 88
   integer i,j,nisccp_max,nsun_max
   
   if(.not.doisccp) return

   nisccp_max = nisccp
   nsun_max = nsun
   if(dompi) then
     call task_max_integer (nisccp,nisccp_max,1)
     call task_max_integer (nsun,nsun_max,1)
   end if

   if(nisccp_max.gt.0.and.nsun.eq.nsun_max) then
     if(dompi) then
       if(nisccp.gt.0) then
         buf(1:ntau*npres) = fq_isccp(1:ntau*npres)
         buf(ntau*npres+1) = totalcldarea
         buf(ntau*npres+2) = meantaucld
         buf(ntau*npres+3) = meanptop
         buf(ntau*npres+4) = lowcldarea
         buf(ntau*npres+5) = midcldarea
         buf(ntau*npres+6) = hghcldarea
         coef = 1./float(nsubdomains*nisccp)
         buf = buf * coef
       else
         buf(1:buf_len) = 0.
       endif
       call task_sum_real(buf,buf1,buf_len)
       fq_isccp(1:ntau*npres) = buf1(1:ntau*npres)
       totalcldarea = buf1(ntau*npres+1)
       meantaucld = buf1(ntau*npres+2)
       meanptop = buf1(ntau*npres+3)
       lowcldarea = buf1(ntau*npres+4)
       midcldarea = buf1(ntau*npres+5)
       hghcldarea = buf1(ntau*npres+6)
     else
       coef = 1./float(nisccp)
       fq_isccp(1:ntau*npres) = fq_isccp(1:ntau*npres) * coef
       totalcldarea = totalcldarea * coef
       meantaucld = meantaucld * coef
       meanptop =  meanptop * coef
       lowcldarea = lowcldarea * coef
       midcldarea = midcldarea * coef
       hghcldarea = hghcldarea * coef
     end if
     meantaucld = tautab(min(255,max(1,nint(meantaucld))))
   else if(nsun.gt.0.and.nsun.eq.nsun_max) then
     fq_isccp = 0.
     totalcldarea = 0.
     lowcldarea = 0.
     midcldarea = 0.
     hghcldarea = 0.
     meantaucld = fillvalue
     meanptop = fillvalue
   else
     fq_isccp = fillvalue
     totalcldarea = fillvalue
     lowcldarea = fillvalue
     midcldarea = fillvalue
     hghcldarea = fillvalue
     meantaucld = fillvalue
     meanptop = fillvalue
   end if

   s_acldisccp = totalcldarea
   s_acldlisccp = lowcldarea 
   s_acldmisccp = midcldarea
   s_acldhisccp = hghcldarea


   if(masterproc) then

     open (ntape,file='./'//trim(case)//'/'// &
                  trim(case)//'_'// &
                  trim(caseid)//'.isccp', &
                  status='unknown',form='unformatted')
     if(nstep.ne.nstat) then
       do while(.true.)
          read (ntape,end=222)
       end do
     endif
222  continue
     backspace(88)

     print *,'Writting isccp file ',caseid
     write(ntape) caseid
     write(ntape)  day
     write(ntape)  ntau,npres
     write(ntape) isccpOpticalDepthBinEdges(1:ntau)
     write(ntape) isccpPressureBinEdges(1:npres)
     write(ntape) fq_isccp
     write(ntape) totalcldarea,lowcldarea,midcldarea,hghcldarea
     write(ntape) meantaucld,meanptop
     close(ntape)

     print*,'isccp Table:'
     write(*,'(a,7(3x,f7.1))') '     ',isccpOpticalDepthBinEdges(1:ntau)
     do j=1,npres
       write(*,'(f6.0,7(3x,f7.3))') isccpPressureBinEdges(j),(fq_isccp(i+(j-1)*ntau),i=1,ntau)
     end do
     print*,'lowcldarea (tau > 0.3)=',lowcldarea
     print*,'midcldarea (tau > 0.3)=',midcldarea
     print*,'hghcldarea (tau > 0.3)=',hghcldarea
     print*,'totalcldarea (tau > 0.3)=',totalcldarea
     print*,'meantaucld (tau > 0.3) ',meantaucld
     print*,'meanptop (tau > 0.3)',meanptop
   
   end if

   call isccp_zero()
 
end subroutine isccp_write


end module isccp


