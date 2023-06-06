subroutine diffuse_scalar_vert_flag_surf_flag (f,fluxb,fluxt, &
                          fdiff,flux,f2lediff,f2lediss,fwlediff,doit)

use grid
use vars, only: tkh_xy, tkh_z, rho, rhow, tkh_z_rf, tkh_x_rf, tkh_y_rf
use mse       ! peters
use params    ! peters
implicit none

! input:	
real f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real fluxb(nx,ny)		! bottom flux
real fluxt(nx,ny)		! top flux
real flux(nz)
real fdiff(nz)
real f2lediff(nzm)
real f2lediss(nzm)
real fwlediff(nzm)
logical doit
! Local
real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real f0(nzm),df0(nzm),factor_xy
real r2dx,r2dy,r2dx0,r2dy0,r2dz
integer i,j,k,kb,kc,jb,jc
real hdif(nx,ny), tmpdif(nzm), coef   ! peters

if(dostatis .or. (doMSE .and. (mseflag.ne.0)) ) then  ! peters added mseflag
	
  do k=1,nzm
    do j=dimy1_s,dimy2_s
     do i=dimx1_s,dimx2_s
      df(i,j,k) = f(i,j,k)
     end do
    end do
  end do

endif


if(RUN3D) then
  call diffuse_scalar3D_vert_flag_surf_flag (f,fluxb,fluxt,tkh_xy,tkh_z,rho,rhow,flux, tkh_z_rf, tkh_x_rf, tkh_y_rf)
else  
  call diffuse_scalar2D (f,fluxb,fluxt,tkh_xy,tkh_z,rho,rhow,flux)
endif

! peters split total advection into horizontal and vertical parts
if(doMSE .and. (mseflag.ne.0)) then
  do i=1,nx
    do j=1,ny
      tmpdif = (f(i,j,:)-df(i,j,:))/dtn
      call columnint(tmpdif,hdif(i,j))
    end do
  end do
                                                                                
  ! now store the dif terms in the appropriate place.  need to substract
  ! off the surface heat fluxes and scale by dtfactor
  coef = dtfactor/float(navgMSE)
  select case(mseflag)
    case(1)             ! update difs
      hdif = cp*(hdif - rhow(1)*fluxb)*coef
      hdifs_mse = hdifs_mse + hdif                ! W/m^2
    case(2)             ! update difh
      hdif = lcond*(hdif - rhow(1)*fluxb)*coef
      hdifh_mse = hdifh_mse + hdif                ! W/m^2
  end select
end if   ! end if(doMSE)

if(dostatis) then
	
  do k=1,nzm
    fdiff(k)=0.
    do j=1,ny
     do i=1,nx
      fdiff(k)=fdiff(k)+f(i,j,k)-df(i,j,k)
     end do
    end do
  end do

endif

if(dostatis.and.doit) then
	
  call stat_varscalar(f,df,f0,df0,f2lediff)
  call stat_sw2(f,df,fwlediff)

  factor_xy=1./float(nx*ny)
  r2dx0=1./(2.*dx)
  r2dy0=1./(2.*dy)
  do k=1,nzm
    f2lediss(k)=0.
    kc=min(nzm,k+1)
    kb=max(1,k-1)
    r2dz=2./((kc-kb)*(adzw(k+1)+adzw(k))*dz)
    r2dx=r2dx0*sqrt((kc-kb)*dx*r2dz) ! grid anisotropy correction
    r2dy=r2dy0*sqrt((kc-kb)*dx*r2dz)
    f2lediss(k)=0.
    do j=1,ny
     jc=j+YES3D
     jb=j-YES3D
     do i=1,nx
      f2lediss(k)=f2lediss(k)-tkh_xy(i,j,k)*( &
                       ((f(i+1,j,k)-f(i-1,j,k))*r2dx)**2+ &
                       ((f(i,jc,k)-f(i,jb,k))*r2dy)**2+ &
                       ((f(i,j,kc)-f0(kc)-f(i,j,kb)+f0(kb))*r2dz)**2 )
     end do
    end do
    f2lediss(k)=f2lediss(k)*2.*factor_xy
  end do

endif

end subroutine diffuse_scalar_vert_flag_surf_flag
