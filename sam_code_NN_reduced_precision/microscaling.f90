module microscaling
  !---------------------------------------------------------------------------------
  ! Module to calculate the appropriate time scale fo microphysics in DARE 
  ! simulations
  !---------------------------Code history------------------------------------------
  ! created:           JA, 8/8/06
  ! Modified:          CW, 5/7/07
  !---------------------------------------------------------------------------------

use vars
use params

implicit none

real, dimension(nx,ny,nz) :: &
	real_stratiformnessc, dtn_scaled	 ! indices of stratiformness/convection

	


CONTAINS
	!
	!=======================================================================
	!


  subroutine testforstratiform
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose
    ! compute the stratiformness of the convection, for the purposes of 
    ! 
    !-----------------------------------------------------------------------

use vars, only: w

implicit none


 
!
!11/7/06 Joe Andersen
!
!Tests for stratiform versus convective motion in netcdf files. To do this,
!it compares the variance (of w, say) in a small area with that in a larger
!encompassing area. Convective motion should show a significantly weaker variance
! after being smoothed slightly on the
!smaller scale compared to the unsmoothed variance, while stratiform motion should
!show similar variance with and without smoothing.

! variables

integer ::									&
	m,		& ! number of points in smaller boxes we smooth over
	n, 		& ! number of points in larger boxes 
	c,a,b,		& ! size of the domain (z,x,y)
	a2,b2,		& ! size of the wrapped domain, to allow for periodicity
	count, count2, count3, count4, i,j,k  
 
real ::										&
	 dtdz, cccc, scaling_number	
	
real, dimension(nz) ::								&
	rhoz

real, dimension(nz-2,-2:nx+3,-2:ny+3) ::							&
	unsmoothedvariance, smoothedvariance,smoothedmean,unsmoothedmean		

real, dimension(nz-2) ::							&
	foo_z2

real, dimension(nz,-2:nx+3,-2:ny+3) ::							&
	ww			 ! vertical velocities
	
real, dimension(nz-2,-2:nx+3,-2:ny+3) ::						&
	slab, 			& ! place to keep the domain, expanded to allow for periodicity
	stratiformnessc,	& ! holding place for index of stratiformness/convection with periodic BCs
	smoothed_slab,   	& ! smoothed data set
	foo5, foo_z
real, dimension(nz-2,4,4) ::								&
	unsmoothedmean_b, unsmootheddiff, unsmoothedsquare, unsmootheddata		

real, dimension(nz-2,2,2) ::							&
	smoothing_data,		 &! holding variable for the generation of the smoothed data set
	smoothedmean_b, smootheddiff, smoothedsquare, smootheddata		



!Set parameters


!ones=1.0
!ones2 = 1.0

m=2
n=4 !Usually n is equal to m*4-3

dtdz = dtn/dz

foo_z = 0.0
smoothedmean = 0.0
unsmoothedvariance = 0.0
smoothedvariance = 0.0
unsmoothedvariance = 0.0
stratiformnessc = 0.0

do k=1,nzm 
  rhoz(k) = rhow(k)*dtdz  
  do j=-2,ny+3
   do i=-2,nx+3
      ww(k,i,j) = w(i,j,k)/(rhoz(k))
   end do 
  end do 
end do 



a=nx
b=ny
c=nz

a2 = a+n
b2 = b+n


slab(:,:,:) = ww(3:c,:,:)

do count2 = -3,a-1
   do count3 = -3,b-1
      
      unsmootheddata = slab(:,1+count2:n+count2,1+count3:n+count3)
      
      
      ! average the unsmootheddata
      
      do count = 1,n
         do count4 = 1,n
            unsmoothedmean(:,count2+1,count3+1) = unsmoothedmean(:,count2+1,count3+1)&
                 + unsmootheddata(:,count,count4)
         end do
      end do
      
      unsmoothedmean(:,count2+1,count3+1) = unsmoothedmean(:,count2+1,count3+1)/((n)*(n))
      
      
      do count=1,nz-2
         do i=1,n
            do j = 1,n
               unsmoothedmean_b(count,i,j) = unsmoothedmean(count,count2+1,count3+1)
            end do
         end do
      end do
      
      
      do k = 1,nz-2
         do i=1,n
            do j = 1,n
               unsmootheddiff(k,i,j) = unsmootheddata(k,i,j)-unsmoothedmean_b(k,i,j)
            end do
         end do
      end do
      
      do k = 1,nz-2
         do i=1,n
            do j = 1,n
               unsmoothedsquare(k,i,j) = unsmootheddiff(k,i,j)*unsmootheddiff(k,i,j)
            end do
         end do
      end do
      
      ! average the unsmoothedsquare
      
      do count = 1,n
         do count4 = 1,n
            unsmoothedvariance(:,count2+1,count3+1) = unsmoothedvariance(:,count2+1,count3+1)&
                 + unsmoothedsquare(:,count,count4)
         end do
      end do
      unsmoothedvariance(:,count2+1,count3+1) = unsmoothedvariance(:,count2+1,count3+1)/((n)*(n))
      
      do count=1,c-2
         do i = 1,n
            do j=1,n
               foo5(count,i,j) = (1.0-erf((1.0*(unsmoothedvariance(count,count2+1,count3+1)*ravefactor**(4./3.)-1.0))))/2.0
            end do
         end do
      end do
      
      do count=1,c-2
         do i = 1,n
            do j=1,n
               stratiformnessc(count,count2+i,count3+j) = stratiformnessc(count,count2+i,count3+j)+foo5(count,i,j)
            end do
         end do
      end do
      
   end do
end do

stratiformnessc = stratiformnessc/float(n*n)
 
do k=1,nzm
   do j=1,ny
      do i=1,nx
         real_stratiformnessc(i,j,k) = 1.0
         
         if (k.eq.2) then
            real_stratiformnessc(i,j,2) = 0.5*stratiformnessc(3,i,j)
         endif
         if (k.gt.2) then
            real_stratiformnessc(i,j,k) = stratiformnessc(k-2,i,j)
         endif
      end do
   end do
end do



do k=1,nzm
   do j=1,ny
      do i=1,nx
         cccc = real_stratiformnessc(i,j,k)
         if(cccc.ge.cmax) then
            scaling_number = 1.0
         elseif(cccc.lt.cmin) then
            scaling_number = 0.3333
         else
            scaling_number = 1.0 + (1.0 - 0.3333) * (cccc-cmax)/(cmax - cmin)
         end if
         dtn_scaled(i,j,k) = dtn*ravefactor**(scaling_number-1.0)
      end do
   end do
end do
end subroutine testforstratiform
!
  !===============================================================================
END MODULE microscaling
