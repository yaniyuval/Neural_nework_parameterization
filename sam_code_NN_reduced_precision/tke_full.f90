
subroutine tke_full

!	this subroutine solves the TKE equation

use vars
use params
implicit none

real def2(nx,ny,nzm)	
real grd,betdz,Ck,Ce,Ces,Ce1,Ce2,smix,Pr,Pr1,Cee,Cs,Cs1
real buoy_sgs,ratio,a_prod_sh,a_prod_bu,a_diss
real lstarn, lstarp, bbb, omn, omp
real qsat,dqsat,dimfactor
integer i,j,k,kc,kb


logical smag_xy

Cs = 0.1944
Cs1 = 0.14
Pr = 3.0
Ck=0.1
Ce=Ck**3/Cs**4
Ces=Ce/0.7*3.0	

smag_xy = douwpbl.and.(uwpbl_ndiv.eq.1)

if(RUN3D) then
   if(smag_xy) then
      call shear_prod_xy(def2)
   else
      call shear_prod3D(def2)
      !dimfactor=max(1.,log10(sqrt(dx*dy)/dz/dxdzfactor/ravefactor))
      !Yani hardcoded the dimfactor to be the same as in the high resolution runs.  
      !dimfactor=max(1.,log10(sqrt(dx/16*dy/16)/dz/dxdzfactor/ravefactor))
      dimfactor=max(1.,log10(sqrt(12000.0*12000.0)/dz/dxdzfactor/ravefactor)) 
   endif
else
  call shear_prod2D(def2)
!  dimfactor=sqrt(2.)
  dimfactor=max(1.,log10(sqrt(dx/16*dy/16)/dz/dxdzfactor/ravefactor))
!Yani hardcoded the dimfactor to be the same as in the high resolution runs.
   !dimfactor=max(1.,log10(sqrt(12000*12000)/dz/dxdzfactor/ravefactor))
endif

do k=1,nzm      
  kb=k-1
  kc=k+1

  grd=dz*adz(k)*dimfactor

  ! if xy smagorinsky, mixing length is xy grid spacing -- like WRF
  if(smag_xy) grd = sqrt(dx*dy)

  betdz=bet(k)/dz/(adzw(kc)+adzw(k))
  Ce1=Ce/0.7*0.19
  Ce2=Ce/0.7*0.51
  if(k.eq.1) then
    kb=1
    kc=2
    betdz=bet(k)/dz/adzw(kc)
    Ce1=Ces/0.7*0.19
    Ce2=Ces/0.7*0.51
  end if
  if(k.eq.nzm) then
    kb=nzm-1
    kc=nzm
    betdz=bet(k)/dz/adzw(k)
    Ce1=Ces/0.7*0.19
    Ce2=Ces/0.7*0.51
  end if
  tkelediss(k) = 0.
  tkesbdiss(k) = 0.
  tkesbshear(k)= 0.
  tkesbbuoy(k) = 0.
  
  do j=1,ny
  do i=1,nx

  ! only calculate buoy_sgs for 3D Smagorinsky
  if(.not.smag_xy) then 

!  SGS buoyancy flux

   if(.not.dokruegermicro) then
   omn = max(0.,min(1.,(tabs(i,j,k)-tbgmin)*a_bg)) 
   omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr)) 
   elseif(dokruegermicro) then
      if (qc(i,j,k)+qi(i,j,k).gt.0.) then
         omn = qc(i,j,k)/(qc(i,j,k)+qi(i,j,k))
      else
         omn = max(0.,min(1.,(tabs(i,j,k)-ttfrz)*a_bgkru))
      end if
      if (qr(i,j,k)+qg(i,j,k)+qs(i,j,k).gt.0.) then
         omp = qr(i,j,k)/(qr(i,j,k)+qg(i,j,k)+qs(i,j,k))
      else
         omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr)) 
      end if
   end if
   lstarn = fac_cond+(1.-omn)*fac_fus 
   lstarp = fac_cond+(1.-omp)*fac_fus 

!   if(qn(i,j,kb).gt.0. .and. qn(i,j,k).gt.0. .and. qn(i,j,kc).gt.0.) then
   if(qn(i,j,k) .gt. 0.) then
     
      dqsat = omn*dtqsatw(tabs(i,j,k),pres(k))+ &
                             (1.-omn)*dtqsati(tabs(i,j,k),pres(k))
      qsat = omn*qsatw(tabs(i,j,k),pres(k))+(1.-omn)*qsati(tabs(i,j,k),pres(k))
      bbb = 1. + 0.61*qsat-qn(i,j,k) -qp(i,j,k)+1.61*tabs(i,j,k)*dqsat
      bbb = bbb / (1.+lstarn*dqsat)
      buoy_sgs=betdz*(bbb*(t(i,j,kc)-t(i,j,kb)) &
	+(bbb*lstarn - (1.+lstarn*dqsat)*tabs(i,j,k))*(q(i,j,kc)-q(i,j,kb)) &
	+(bbb*lstarp - (1.+lstarp*dqsat)*tabs(i,j,k))*(qp(i,j,kc)-qp(i,j,kb)) )
   else

      bbb = 1.+0.61*q(i,j,k)-qp(i,j,k)
      buoy_sgs=betdz*( bbb*(t(i,j,kc)-t(i,j,kb)) &
        +0.61*tabs(i,j,k)*(q(i,j,kc)-q(i,j,kb)) &
	+(bbb*lstarp-tabs(i,j,k))*(qp(i,j,kc)-qp(i,j,kb)) )    	     
   end if

   else ! for horizontal Smagorinsky, set buoy_sgs = 0.0
      buoy_sgs = 0.0
   end if

   tke(i,j,k)=max(0.,tke(i,j,k))

   if(buoy_sgs.le.0.) then
     smix=grd
   else
     smix=min(grd,max(0.1*grd,0.76*sqrt(tke(i,j,k)/buoy_sgs+1.e-10)))
   end if
   

   ratio=smix/grd
   Pr1=1.+2.*ratio
   Cee=Ce1+Ce2*ratio

   
   if(dosmagor) then
         
     tk_xy(i,j,k)=sqrt(Ck**3/Cee*max(0.,def2(i,j,k)-Pr1*buoy_sgs))*smix**2
     ! the following is in WRF for xy Smagorinsky -- seems appropriate here too
     if(smag_xy) tk_xy(i,j,k) = min(tk_xy(i,j,k), 10.0*smix)
     if(.not.doscalar) tke(i,j,k) = (tk_xy(i,j,k)/(Ck*smix))**2
     a_prod_sh=(tk_xy(i,j,k)+0.001)*def2(i,j,k)
     a_prod_bu=-(tk_xy(i,j,k)+0.001)*Pr1*buoy_sgs
     a_diss=a_prod_sh+a_prod_bu

    else

     a_prod_sh=(tk_xy(i,j,k)+0.001)*def2(i,j,k)
     a_prod_bu=-(tk_xy(i,j,k)+0.001)*Pr1*buoy_sgs
     a_diss=Cee/smix*tke(i,j,k)**1.5	   
     tke(i,j,k)=max(0.,tke(i,j,k)+dtn*sqrt(betafactor)*(max(0.,a_prod_sh+a_prod_bu)-a_diss))
     !in the DARE framework, speed up the TKE adjustment too.
     tk_xy(i,j,k)=Ck*smix*sqrt(tke(i,j,k))

   end if

   tkh_xy(i,j,k) = Pr1*tk_xy(i,j,k)
   if(.not.smag_xy) then             ! set vertical diffusivity here
      tk_z(i,j,k) = tk_xy(i,j,k)     ! otherwise do it in uwpbldrv
      tkh_z(i,j,k) = Pr1*tk_z(i,j,k) 
   end if

   tkelediss(k) = tkelediss(k) - a_prod_sh
   tkesbdiss(k) = tkesbdiss(k) + a_diss
   tkesbshear(k)= tkesbshear(k)+ a_prod_sh
   tkesbbuoy(k) = tkesbbuoy(k) + a_prod_bu

  end do ! i
  end do ! j

  tkelediss(k) = tkelediss(k)/float(nx*ny)

end do ! k


end


