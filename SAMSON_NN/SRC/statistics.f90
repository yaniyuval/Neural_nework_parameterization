
subroutine statistics()

use vars
use params
use hbuffer
use isccp, only: isccp_get
use rad !bloss
implicit none	
	
	real mse(nzm)
	real dse(nzm)
	real sse(nzm)
	real tpz(nzm)
	real tlz(nzm)
	real tvz(nzm)
	real qcz(nzm)
	real qiz(nzm)
	real tez(nzm)
	real qvz(nzm)
	real qrz(nzm)
	real qsz(nzm)
	real qgz(nzm)
	real relhz(nzm)

	real u2z(nzm)
	real v2z(nzm)
	real w2z(nzm)
	real w22(nzm)
	real w3z(nzm)
	real skw(nzm)
	real t2z(nzm)
	real tqz(nzm)
	real q2z(nzm)
	real s2z(nzm)
	real qc2z(nzm)
	real qi2z(nzm)
	real qs2z(nzm)
	real tkez(nzm)
	real fadv(nz)
	real shear(nz)
	real shearx(nzm)
	real sheary(nzm)
	real presx(nzm)
	real presy(nzm)
	real twgrad(nzm)
	real qwgrad(nzm)
	real swgrad(nzm)
	
	
	real tvwle(nzm)
	real qccwle(nzm)
	real qiiwle(nzm)
	real aup(nzm)
	real wcl(nzm)
	real ucl(nzm)
	real vcl(nzm)
	real tcl(nzm)
	real tacl(nzm)
	real tvcl(nzm)
	real qcl(nzm)
	real qccl(nzm)
	real qicl(nzm)
	real qpcl(nzm)
	real twcl(nzm)
	real swcl(nzm)
	real qwcl(nzm)
	real tvwcl(nzm)
	real qcwcl(nzm)
	real qiwcl(nzm)
	real scl(nzm)
	real wacl(nzm)	
	real cld(nzm)
 	real cldd(nzm)
	real hydro(nzm)
	real qsatwz(nzm)

	real tvirt(nx,ny,nzm)
	
	integer i,j,k,n
	real qcc,qii,qrr,qss,qgg,lstarn,lstarp,coef,coef1
	real factor_xy, factor_n, tmp(4)
        real buffer(nzm,6),buffer1(nzm,6)
	real prof1(nzm),prof2(nzm),prof3(nzm),prof4(nzm)	
	real cwp(nx,ny),cwpl(nx,ny),cwpm(nx,ny),cwph(nx,ny)
        real cwpmax,prec_conv, prec_strat
	real omn, omn1, omp, omg, qsat
	logical condition, condition_cl
	real zero(nzm)

	integer topind(nx,ny)	

!========================================================================
! UW ADDITIONS

        real taz(nzm), threshmin, threshmax, qvv ! added for precip stats.
        real tvcla(nzm)!kzm added Apr. 7,2004 for thetav anomalies
        real tau_isccp(nx,ny)
	real tr0(nzm), tr2(nzm) !kzm tracer

        integer kb,kc

        ! bloss rce intercomparison statistics
        real sfz(nzm), cfz(nzm), ufz(nzm), umfz(nzm), uhz(nzm)
        real hflxz(nzm), rcei(nzm), dfz(nzm), dmfz(nzm), dhz(nzm)
        real tmpmse, tmpw

        !bloss reflectivity computation and histogram data
        real tmp_wgt, tmpr, tmpg, tmps, pii, fac_ice, dbZ(nx,ny)
        real tmp_mu, tmp_exponent
        real cfadm10(nzm), cfadm05(nzm), cfad00(nzm), cfad05(nzm), &
             cfad10(nzm), cfad15(nzm), cfad20(nzm), cfad25(nzm), &
             cfad30(nzm), cfad35(nzm), cfad40(nzm), cfad45(nzm), &
             cfad50(nzm), cfad55(nzm), cfad60(nzm), cfad65(nzm)
        real(4) gammafff
        external gammafff

        !bloss case 5 statistics
        integer :: nts
        real :: fac, thm, thp, qvm, qvp, tkhm
        real, dimension(nx,ny) :: cttemp, cloudxy
        real, dimension(nzm) :: clfracz, hmfracz, cufracz, cumflxz, &
             wufracz, wumflxz, bcufracz, bcumflxz, &
             c5t, c5sst, c5dse, c5qv, c5mse, c5u, c5v, c5shf, c5lhf,&
             & c5tau, c5tav, c5swds, c5lwds, c5swus, c5lwus, c5swdt,&
             & c5swut, c5lwut, c5cld, c5pcp, c5wcf, c5wpcp, c5ifc,&
             & c5ifpcp, c5hcf, c5hpcp, c5lcf, c5lpcp, c5mcf, c5mpcp,&
             & c5pw, c5lwp, c5iwp, c5rp, c5sp, c5gp
        real, dimension(nz) :: thflux, thflxs, qvflux, qvflxs, dthadv, &
             dthdif, dqvadv, dqvdif

        real, parameter :: eps = 1.e-24

! END UW ADDITIONS
!========================================================================



	factor_xy = 1./float(nx*ny)
	factor_n = 1./float(nsubdomains)
	
!-----------------------------------------------
!	Mean thermodynamics profiles:
!-----------------------------------------------	
		
	do k=1,nzm
	 dse(k)=0.
	 mse(k)=0.
	 sse(k)=0.
	 tpz(k) = 0.
	 tlz(k) = 0.
	 tvz(k) = 0.
	 tez(k) = 0.
	 qvz(k) = 0.
	 qcz(k) = 0.
	 qiz(k) = 0.
	 qrz(k) = 0.
	 qsz(k) = 0.
	 qgz(k) = 0.
	 qsatwz(k)=0.
	 relhz(k)=0.
	 prof1(k)=0.
	 prof2(k)=0.
	 prof3(k)=0.
	 zero(k)=0.
	 do j=1,ny
	  do i=1,nx
           if(.not.dokruegermicro) then !bloss  
	   omn = omegan(tabs(i,j,k))
	   omp = omegap(tabs(i,j,k))
	   omg = omegag(tabs(i,j,k))
	   lstarn=fac_cond+(1.-omn)*fac_fus
	   lstarp=fac_cond+(1.-omp)*fac_fus
	   qcc=qn(i,j,k)*omn
	   qii=qn(i,j,k)*(1.-omn)
	   qrr=qp(i,j,k)*omp
	   qss=qp(i,j,k)*(1.-omp)*(1.-omg)
	   qgg=qp(i,j,k)*(1.-omp)*omg
           elseif(dokruegermicro) then !bloss
            qcc=qc(i,j,k)
            qii=qi(i,j,k)
            qrr=qr(i,j,k)
            qss=qs(i,j,k)
            qgg=qg(i,j,k)
            if (qcc+qii.gt.0.) then
             omn = qcc/(qcc+qii)
            else
             omn = max(0.,min(1.,(tabs(i,j,k)-ttfrz)*a_bgkru))
            end if
            lstarn=fac_cond+(1.-omn)*fac_fus
            lstarp=fac_cond+fac_fus*(qss+qgg)/(qss+qgg+qrr+eps)
           end if
           qrz(k)=qrz(k)+qrr
           qsz(k)=qsz(k)+qss
           qgz(k)=qgz(k)+qgg
           qcz(k)=qcz(k)+qcc
           qiz(k)=qiz(k)+qii
	   prof1(k)=prof1(k)+qcc+qii
	   prof2(k)=prof2(k)+qrr+qss+qgg
	   prof3(k)=prof3(k)+qcc+qii+qrr+qss+qgg
	   tmp(1)=tabs(i,j,k)*prespot(k)
           tpz(k)=tpz(k)+tmp(1)
	   tlz(k)=tlz(k)+tmp(1)*(1.-fac_cond*qn(i,j,k)/tabs(i,j,k))
           tvirt(i,j,k)=tmp(1)*(1.+0.61*q(i,j,k)-1.61*(qn(i,j,k))-qp(i,j,k))
	   tvz(k)=tvz(k)+tvirt(i,j,k)
           tez(k)=tez(k)+t(i,j,k)+lstarn*qn(i,j,k)+lstarp*qp(i,j,k) &
                     +fac_cond*(q(i,j,k)-qcc-qii)-fac_fus*(qii+qss+qgg)
	   qvz(k) =qvz(k)+q(i,j,k)-qn(i,j,k)
     	   dse(k)=dse(k)+t(i,j,k)+lstarn*qn(i,j,k)+lstarp*qp(i,j,k)
	   mse(k)=mse(k)+t(i,j,k)+lstarn*qn(i,j,k)+lstarp*qp(i,j,k) &
                	+fac_cond*(q(i,j,k)-qn(i,j,k))
	   qsat = omn*qsatw(tabs(i,j,k),pres(k))+ &
             (1.-omn)*qsati(tabs(i,j,k),pres(k)) 
	   sse(k)=sse(k)+t(i,j,k)+lstarn*qn(i,j,k)+lstarp*qp(i,j,k) &
                        +fac_cond*qsat
	   qsatwz(k) = qsatwz(k)+qsatw(tabs(i,j,k),pres(k))
	   relhz(k)=relhz(k)+(q(i,j,k)-qn(i,j,k))/qsatw(tabs(i,j,k),pres(k))
	  end do
	 end do
	end do	
	

	call hbuf_avg_put('TL',t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)
	call hbuf_avg_put('TABS',tabs,1,nx, 1,ny, nzm,1.)
	call hbuf_avg_put('U',u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,1.)
	call hbuf_avg_put('V',v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,1.)
	call hbuf_avg_put('QT',q,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.e3)
	call hbuf_avg_put('TK_XY',tk_xy,0,nxp1, (1-YES3D),nyp1, nzm,1.)
	call hbuf_avg_put('TK_Z',tk_z,0,nxp1, (1-YES3D),nyp1, nzm,1.)
	call hbuf_avg_put('TKH_XY',tkh_xy,0,nxp1, (1-YES3D),nyp1, nzm,1.)
	call hbuf_avg_put('TKH_Z',tkh_z,0,nxp1, (1-YES3D),nyp1, nzm,1.)
	call hbuf_put('TABSOBS',tg0,1.)
	call hbuf_put('QVOBS',qg0,1.e3)
	call hbuf_put('UOBS',ug0,1.)
	call hbuf_put('VOBS',vg0,1.)
        call hbuf_put('WOBS',wsub,1.)
	call hbuf_put('TTEND',ttend,86400.)
	call hbuf_put('QTEND',qtend,86400.*1.e3)
        call hbuf_put('WWVCONV',wwave_conv,1.)
        call hbuf_put('W_WAVE',w_wave,1.)
        call hbuf_put('H_WAVE',heating_wave,86400.)
        call hbuf_put('Q_WAVE',q_wave,1.)
        call hbuf_put('T_WAVE',t_wave,1.)
        call hbuf_put('H_WAVEBG',heating_wavebg,86400.)
        call hbuf_put('Q_WAVEBG',q_wavebg,1.)
        call hbuf_put('T_WAVEBG',t_wavebg,1.)
        call hbuf_put('P_WAVE',pres_wave,1.)
        call hbuf_put('T_TARGET',t_target,1.)
        call hbuf_put('TTENDWV',ttend_wave,86400.)
        call hbuf_put('QTENDWV',qtend_wave,86400.*1.e3)

	call hbuf_put('DSE',dse,factor_xy)
	call hbuf_put('MSE',mse,factor_xy)
	call hbuf_put('SSE',sse,factor_xy)
	call hbuf_put('THETA',tpz,factor_xy)
	call hbuf_put('THETAL',tlz,factor_xy)
	call hbuf_put('THETAV',tvz,factor_xy)
	call hbuf_put('THETAE',tez,factor_xy)

	call hbuf_put('PRES',pres,1.)
	call hbuf_put('RHO',rho,1.)
	call hbuf_put('QV',qvz,1.e3*factor_xy)
	call hbuf_put('QC',qcz,1.e3*factor_xy)
	call hbuf_put('QI',qiz,1.e3*factor_xy)
	call hbuf_put('QR',qrz,1.e3*factor_xy)
	call hbuf_put('QS',qsz,1.e3*factor_xy)
	call hbuf_put('QG',qgz,1.e3*factor_xy)
	call hbuf_put('QN',prof1,1.e3*factor_xy)
	call hbuf_put('QP',prof2,1.e3*factor_xy)
	call hbuf_put('QCOND',prof3,1.e3*factor_xy)
	call hbuf_put('QSAT',qsatwz,1.e3*factor_xy)
	call hbuf_put('RELH',relhz,100.*factor_xy)

!-------------------------------------------------------------
!	Fluxes:
!-------------------------------------------------------------

	do k=1,nzm
	  tmp(1) = dz/rhow(k)
	  tmp(2) = tmp(1) / dtn
	  uwsb(k) = uwsb(k) * tmp(1)
	  vwsb(k) = vwsb(k) * tmp(1)
	  twsb(k) = twsb(k) * tmp(1) * rhow(k) * cp
	  qwsb(k) = qwsb(k) * tmp(1) * rhow(k) * lcond
	  uwle(k) = uwle(k)*tmp(1) + uwsb(k)
	  vwle(k) = vwle(k)*tmp(1) + vwsb(k)
	  twle(k) = twle(k)*tmp(2)*rhow(k)*cp + twsb(k)
	  qwle(k) = qwle(k)*tmp(2)*rhow(k)*lcond + qwsb(k)
	  if(docloud.and.doprecip) then
	    qpwsb(k) = qpwsb(k) * tmp(1)*rhow(k)*lcond
     	    qpwle(k) = qpwle(k)*tmp(2)*rhow(k)*lcond + qpwsb(k)
	  else
	    qpwsb(k) = 0.
	    qpwle(k) = 0.
	    precflux(k) = 0.
	  endif
	end do
	uwle(nz) = 0.
	vwle(nz) = 0.
	uwsb(nz) = 0.
	vwsb(nz) = 0.

	call hbuf_put('UW',uwle,factor_xy)
	call hbuf_put('VW',vwle,factor_xy)
	call hbuf_put('UWSB',uwsb,factor_xy)
	call hbuf_put('VWSB',vwsb,factor_xy)
	call hbuf_put('TLFLUX',twle,factor_xy)
	call hbuf_put('QTFLUX',qwle,factor_xy)
	call hbuf_put('TLFLUXS',twsb,factor_xy)
	call hbuf_put('QTFLUXS',qwsb,factor_xy)
	call hbuf_put('QPFLUXS',qpwsb,factor_xy)
	call hbuf_put('QPFLUX',qpwle,factor_xy)
	call hbuf_put('PRECIP',precflux,factor_xy/dt*dz*86400./(nstatis+1.e-5))
	
	do j=1,ny
	 do i=1,nx
	  precsfc(i,j)=precsfc(i,j)*dz/dt*3600./(nstatis+1.e-5)
	 end do
	end do

	call stat_precip(precsfc, prec_conv, prec_strat)
	p_conv = p_conv + prec_conv
	p_strat = p_strat + prec_strat

	do k=1,nzm
	 tvz(k) = 0.
	 qcz(k) = 0.
	 qiz(k) = 0.
	 qsatwz(k) = 0.
	 prof1(k)=0.
	 prof2(k)=0.
         if(.not.dokruegermicro) then !bloss
	 do j=1,ny
	  do i=1,nx
	    omn = omegan(tabs(i,j,k))
	    tvz(k) = tvz(k) + tvirt(i,j,k)
	    qcz(k) = qcz(k) + qn(i,j,k)*omn
	    qiz(k) = qiz(k) + qn(i,j,k)*(1.-omn)
	    qsat = omn*qsatw(tabs(i,j,k),pres(k))+ &
                       (1.-omn)*qsati(tabs(i,j,k),pres(k)) 
	    qsatwz(k) = qsatwz(k)+qsat
	  end do
	 end do
         elseif(dokruegermicro) then !bloss
          do j=1,ny
	   do i=1,nx
	     tvz(k) = tvz(k) + tvirt(i,j,k)
	     qcz(k) = qcz(k) + qc(i,j,k)
	     qiz(k) = qiz(k) + qi(i,j,k)
             if (qc(i,j,k)+qi(i,j,k).gt.0.) then
              omn = qc(i,j,k)/(qc(i,j,k)+qi(i,j,k))
             else
              omn = max(0.,min(1.,(tabs(i,j,k)-ttfrz)*a_bgkru))
             end if
	     omn = omegan(tabs(i,j,k))
	     qsat = omn*qsatw(tabs(i,j,k),pres(k))+ &
                        (1.-omn)*qsati(tabs(i,j,k),pres(k)) 
	     qsatwz(k) = qsatwz(k)+qsat
 	   end do
	  end do
         end if
	 tvz(k) = tvz(k)*factor_xy
	 qcz(k) = qcz(k)*factor_xy
	 qiz(k) = qiz(k)*factor_xy	 
	 qsatwz(k) = qsatwz(k)*factor_xy
	end do
	if(dompi) then
	  coef1 = 1./float(nsubdomains)
	  do k=1,nzm
	    buffer(k,1) = tvz(k)
	    buffer(k,2) = qcz(k)
	    buffer(k,3) = qiz(k)
	    buffer(k,4) = qsatwz(k)
	  end do
	  call task_sum_real(buffer,buffer1,nzm*4)
	  do k=1,nzm
	    tvz(k) = buffer1(k,1) * coef1
	    qcz(k) = buffer1(k,2) * coef1
	    qiz(k) = buffer1(k,3) * coef1
	    qsatwz(k) = buffer1(k,4) * coef1
	  end do
	end if ! dompi

	tvwle(1) = 0.
	qccwle(1) = 0.
	qiiwle(1) = 0.
	do k=2,nzm
	 tvwle(k) = 0.
	 qccwle(k) = 0.
	 qiiwle(k) = 0.
         if(.not.dokruegermicro) then !bloss
	 do j=1,ny
	  do i=1,nx
	    omn = omegan(tabs(i,j,k))
	    omn1 = omegan(tabs(i,j,k-1))
	    tvwle(k) = tvwle(k) + 0.5*w(i,j,k)* &
		(tvirt(i,j,k-1)-tvz(k-1)+tvirt(i,j,k)-tvz(k))
	    qccwle(k) = qccwle(k) + 0.5*w(i,j,k)* &
		(qn(i,j,k-1)*omn1-qcz(k-1)+ qn(i,j,k)*omn-qcz(k))	  
	    qiiwle(k) = qiiwle(k) + 0.5*w(i,j,k)* &
                (qn(i,j,k-1)*(1.-omn1)-qiz(k-1)+qn(i,j,k)*(1.-omn)-qiz(k))
	    prof1(k)=prof1(k)+rho(k)*0.5* &
                (w(i,j,k)**2+w(i,j,k+1)**2)*(t(i,j,k)-t0(k))
          end do
	 end do
         elseif(dokruegermicro) then !bloss
          do j=1,ny
	   do i=1,nx
	     tvwle(k) = tvwle(k) + 0.5*w(i,j,k)* &
		 (tvirt(i,j,k-1)-tvz(k-1)+tvirt(i,j,k)-tvz(k))
	     qccwle(k) = qccwle(k) + 0.5*w(i,j,k)* &
		 (qc(i,j,k-1)-qcz(k-1)+qc(i,j,k)-qcz(k))	  
	     qiiwle(k) = qiiwle(k) + 0.5*w(i,j,k)* &
                 (qi(i,j,k-1)-qiz(k-1)+qi(i,j,k)-qiz(k))
	     prof1(k)=prof1(k)+rho(k)*0.5* &
                 (w(i,j,k)**2+w(i,j,k+1)**2)*(t(i,j,k)-t0(k))
           end do
	  end do
         end if
	 tvwle(k) = tvwle(k)*rhow(k)*cp
	 qccwle(k) = qccwle(k)*rhow(k)*lcond
	 qiiwle(k) = qiiwle(k)*rhow(k)*lcond
	end do	

	call hbuf_put('TVFLUX',tvwle,factor_xy)
	call hbuf_put('QCFLUX',qccwle,factor_xy)
	call hbuf_put('QIFLUX',qiiwle,factor_xy)

!---------------------------------------------------------
!	Mean turbulence related profiles:
!-----------------------------------------------------------


	do k=1,nzm

	 u2z(k) = 0.
	 v2z(k) = 0.
	 w2z(k) = 0.
	 w22(k) = 0.
	 w3z(k) = 0.
	 aup(k) = 0.
	 t2z(k) = 0.
	 tqz(k) = 0.
	 q2z(k) = 0.
	 s2z(k) = 0.
	 qc2z(k) = 0.
	 qi2z(k) = 0.
	 qs2z(k) = 0.
	 do j=1,ny
	  do i=1,nx
	    u2z(k) = u2z(k)+(u(i,j,k)-u0(k))**2	  
	    v2z(k) = v2z(k)+(v(i,j,k)-v0(k))**2
	    w2z(k) = w2z(k)+0.5*(w(i,j,k+1)**2+w(i,j,k)**2)
	    w22(k) = w22(k)+w(i,j,k)**2
	    w3z(k) = w3z(k)+0.5*(w(i,j,k+1)**3+w(i,j,k)**3)	  
	    t2z(k) = t2z(k)+(t(i,j,k)-t0(k))**2	  
	    tqz(k) = tqz(k)+(t(i,j,k)-t0(k))*(q(i,j,k)-q0(k))	  
	    q2z(k) = q2z(k)+(q(i,j,k)-q0(k))**2
	    s2z(k) = s2z(k)+(tke(i,j,k)-tke0(k))**2
	    if(w(i,j,k)+w(i,j,k+1).gt.0) aup(k) = aup(k) + 1
	  end do
	 end do
	 skw(k) = w3z(k)/(w2z(k)*factor_xy+1.e-5)**1.5
	 tkez(k)= 0.5*(u2z(k)+v2z(k)*YES3D+w2z(k))
	 tvwle(k) = tvwle(k) * bet(k) /(rho(k)*cp)
         if(.not.dokruegermicro) then !bloss
	 do j=1,ny
	  do i=1,nx
	    omn = omegan(tabs(i,j,k))
	    qc2z(k) = qc2z(k)+(qn(i,j,k)*omn-qcz(k))**2
	    qi2z(k) = qi2z(k)+(qn(i,j,k)*(1.-omn)-qiz(k))**2
	    qsat = omn*qsatw(tabs(i,j,k),pres(k))+ &
                    (1.-omn)*qsati(tabs(i,j,k),pres(k)) 
	    qs2z(k) = qs2z(k)+(qsat-qsatwz(k))**2
	  end do
	 end do
	 elseif(dokruegermicro) then !bloss
	  do j=1,ny
	   do i=1,nx
	     qc2z(k) = qc2z(k)+(qc(i,j,k)-qcz(k))**2
	     qi2z(k) = qi2z(k)+(qi(i,j,k)-qiz(k))**2
             if (qcc+qii.gt.0.) then
              omn = qcc/(qcc+qii)
             else
              omn = max(0.,min(1.,(tabs(i,j,k)-ttfrz)*a_bgkru))
             end if
	     qsat = omn*qsatw(tabs(i,j,k),pres(k))+ &
                     (1.-omn)*qsati(tabs(i,j,k),pres(k)) 
	     qs2z(k) = qs2z(k)+(qsat-qsatwz(k))**2
	   end do
	  end do
         end if
	end do

	call hbuf_put('U2',u2z,factor_xy)
	call hbuf_put('V2',v2z,factor_xy)
	call hbuf_put('W2',w2z,factor_xy)
	call hbuf_put('W3',w3z,factor_xy)
	call hbuf_put('WSKEW',skw,factor_xy)
	call hbuf_put('AUP',aup,factor_xy)

	call hbuf_put('TL2',t2z,factor_xy)
	call hbuf_put('TQ',tqz,factor_xy)
	call hbuf_put('QT2',q2z,1.e6*factor_xy)
	call hbuf_put('QC2',qc2z,1.e6*factor_xy)
	call hbuf_put('QI2',qi2z,1.e6*factor_xy)
	call hbuf_put('QS2',qs2z,1.e6*factor_xy)
	
	call hbuf_put('TKE',tkez,factor_xy)
	
	fadv(1)=0.
	fadv(nz)=0.

!-----------------------------------------------------------------
!  TKE balance:

	shear(1)=0.
	shear(nz)=0.
	do k=2,nzm
          shear(k)=-( (uwle(k)-uwsb(k))*(u0(k)-u0(k-1)) &
            +(vwle(k)-vwsb(k))*(v0(k)-v0(k-1))*YES3D )*factor_xy /(dz*adzw(k))
	end do	
	do k=1,nzm
	  shear(k)=0.5*(shear(k)+shear(k+1)) 
	  tkeleadv(k)=tkeleadv(k)-shear(k)
	  tkelediff(k)=tkelediff(k)-tkelediss(k)
	end do

	call hbuf_put('ADVTR',tkeleadv,1.)
	call hbuf_put('PRESSTR',tkelepress,1.)
	call hbuf_put('BUOYA',tkelebuoy,1.)
	call hbuf_put('SHEAR',shear,1.)
	call hbuf_put('DISSIP',tkelediss,1.)
	call hbuf_put('DIFTR',tkelediff,1.)

!-----------------------------------------------------------------
!  Momentum flux balance:

!  UW advection d(w'w'u')/dz:

	do k=2,nzm
	 fadv(k)=0.
	 do j=1,ny
	  do i=1,nx
	    fadv(k)=fadv(k)+w(i,j,k)**2*rhow(k)*0.5* &
                         ( u(i,j,k-1)-u0(k-1)+u(i,j,k)-u0(k))
	  end do
	 end do
	 fadv(k)=fadv(k)*factor_xy
	end do
	do k=1,nzm
	 coef=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	 shearx(k)=momleadv(k,1)-coef
	 momleadv(k,1)=coef
	end do	


!  VW advection d(w'w'v')/dz:

	do k=2,nzm
	 fadv(k)=0.
	 do j=1,ny
	  do i=1,nx
	    fadv(k)=fadv(k)+w(i,j,k)**2*rhow(k)*0.5* &	
                         ( v(i,j,k-1)-v0(k-1)+v(i,j,k)-v0(k))
	  end do
	 end do
	 fadv(k)=fadv(k)*factor_xy
	end do
	do k=1,nzm
	 coef=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	 sheary(k)=momleadv(k,2)-coef
	 momleadv(k,2)=coef
	end do	


!  UW advection d(p'u')/dz:

	do k=1,nz
	 fadv(k)=0.
	 if(k.eq.1) then
	   do j=1,ny
	    do i=1,nx
	      fadv(k)=fadv(k)+(1.5*(u(i,j,k)-u0(k))*p(i,j,k)*rho(k)- &
                     0.5*(u(i,j,k+1)-u0(k+1))*p(i,j,k+1)*rho(k+1))
	    end do
	   end do
	 else if(k.eq.nz) then
	   do j=1,ny
	    do i=1,nx
	     fadv(k)=fadv(k)+(1.5*(u(i,j,k-1)-u0(k-1))*p(i,j,k-1)*rho(k-1)- &
                  0.5*(u(i,j,k-2)-u0(k-2))*p(i,j,k-2)*rho(k-2))
	    end do
	   end do
	 else
	   do j=1,ny
	    do i=1,nx
	      fadv(k)=fadv(k)+0.5*((u(i,j,k)-u0(k))*p(i,j,k)*rho(k)+ &
                       (u(i,j,k-1)-u0(k-1))*p(i,j,k-1)*rho(k-1))
	    end do
	   end do
	 end if
	 fadv(k)=fadv(k)*factor_xy
	end do
	do k=1,nzm
	 presx(k)=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	end do	


!  VW advection d(p'v')/dz:

	do k=1,nz
	 fadv(k)=0.
	 if(k.eq.1) then
	   do j=1,ny
	    do i=1,nx
	      fadv(k)=fadv(k)+(1.5*(v(i,j,k)-v0(k))*p(i,j,k)*rho(k)- & 
                    0.5*(v(i,j,k+1)-v0(k+1))*p(i,j,k+1)*rho(k+1))
	    end do
	   end do
	 else if(k.eq.nz) then
	   do j=1,ny
	    do i=1,nx
             fadv(k)=fadv(k)+(1.5*(v(i,j,k-1)-v0(k-1))*p(i,j,k-1)*rho(k-1)- &
                     0.5*(v(i,j,k-2)-v0(k-2))*p(i,j,k-2)*rho(k-2))
	    end do
	   end do
	 else
	   do j=1,ny
	    do i=1,nx
	      fadv(k)=fadv(k)+0.5*((v(i,j,k)-v0(k))*p(i,j,k)*rho(k)+ &
                       (v(i,j,k-1)-v0(k-1))*p(i,j,k-1)*rho(k-1))
	    end do
	   end do
	 end if
	 fadv(k)=fadv(k)*factor_xy
	end do
	do k=1,nzm
	 presy(k)=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	end do	

	do k=1,nzm
	  momlepress(k,1)=momlepress(k,1)-presx(k)
	  momlepress(k,2)=momlepress(k,2)-presy(k)
	  momlepress(k,3)=momlepress(k,3)-tkelepress(k)
	end do

	call hbuf_put('WUADV',momleadv(1,1),1.)
	call hbuf_put('WUANIZ',momlepress(1,1),1.)
	call hbuf_put('WUBUOY',momlebuoy(1,1),1.)
	call hbuf_put('WUSHEAR',shearx,1.)
	call hbuf_put('WUPRES',presx,1.)
	call hbuf_put('WUDIFF',momlediff(1,1),1.)

	call hbuf_put('WVADV',momleadv(1,2),1.)
	call hbuf_put('WVANIZ',momlepress(1,2),1.)
	call hbuf_put('WVBUOY',momlebuoy(1,2),1.)
	call hbuf_put('WVSHEAR',sheary,1.)
	call hbuf_put('WVPRES',presy,1.)
	call hbuf_put('WVDIFF',momlediff(1,2),1.)

	call hbuf_put('W2BUOY',momlebuoy(1,3),2.)
	call hbuf_put('W2ADV',momleadv(1,3),2.)
	call hbuf_put('W2REDIS',momlepress(1,3),2.)
	call hbuf_put('W2PRES',tkelepress,2.)
	call hbuf_put('W2DIFF',momlediff(1,3),2.)

!-----------------------------------------------------------
! T2 and Q2 variance budget:


	do k=1,nzm
	  q2lediff(k)=q2lediff(k)-q2lediss(k)
	  t2lediff(k)=t2lediff(k)-t2lediss(k)
	end do	
	

	call hbuf_put('T2ADVTR',t2leadv,1.)
	call hbuf_put('T2GRAD',t2legrad,1.)
	call hbuf_put('T2DISSIP',t2lediss,1.)
	call hbuf_put('T2DIFTR',t2lediff,1.)
	call hbuf_put('T2PREC',t2leprec,1.)

	call hbuf_put('Q2ADVTR',q2leadv,1.)
	call hbuf_put('Q2GRAD',q2legrad,1.)
	call hbuf_put('Q2DISSIP',q2lediss,1.)
	call hbuf_put('Q2DIFTR',q2lediff,1.)
	call hbuf_put('Q2PREC',q2leprec,1.)

!------------------------------------------------------------------
! HW and QW budgets:


	fadv(1)=0.
	fadv(nz)=0.


!  HW advection d(w'w'h')/dz:

	do k=2,nzm
	 fadv(k)=0.
	 do j=1,ny
	  do i=1,nx
	    fadv(k)=fadv(k)+w(i,j,k)**2*rhow(k)*0.5* &
                 	      ( t(i,j,k-1)-t0(k-1)+t(i,j,k)-t0(k))
	  end do
	 end do
	end do
	do k=1,nzm
	 coef=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	 twgrad(k)=twleadv(k)-coef
	 twleadv(k)=coef
	end do	


!  QW advection d(w'w'h')/dz:

	do k=2,nzm
	 fadv(k)=0.
	 do j=1,ny
	  do i=1,nx
	    fadv(k)=fadv(k)+w(i,j,k)**2*rhow(k)*0.5* &
                	      ( q(i,j,k-1)-q0(k-1)+q(i,j,k)-q0(k))
	  end do
	 end do
	end do
	do k=1,nzm
	 coef=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	 qwgrad(k)=qwleadv(k)-coef
	 qwleadv(k)=coef
	end do	



	call hbuf_put('TWADV',twleadv,factor_xy)
	call hbuf_put('TWDIFF',twlediff,factor_xy)
	call hbuf_put('TWGRAD',twgrad,factor_xy)
	call hbuf_put('TWBUOY',twlebuoy,factor_xy)
	call hbuf_put('TWPRES',twlepres,factor_xy)
	call hbuf_put('TWPREC',twleprec,factor_xy)

	call hbuf_put('QWADV',qwleadv,factor_xy)
	call hbuf_put('QWDIFF',qwlediff,factor_xy)
	call hbuf_put('QWGRAD',qwgrad,factor_xy)
	call hbuf_put('QWBUOY',qwlebuoy,factor_xy)
	call hbuf_put('QWPRES',qwlepres,factor_xy)
	call hbuf_put('QWPREC',qwleprec,factor_xy)


!  SW advection d(w'w's')/dz:

	if(doscalar) then

	  do k=2,nzm
	   fadv(k)=0.
	   do j=1,ny
	    do i=1,nx
	      fadv(k)=fadv(k)+w(i,j,k)**2*rhow(k)*0.5* &
          	        ( tke(i,j,k-1)-tke0(k-1)+tke(i,j,k)-tke0(k))
	    end do
	   end do
	  end do
	  do k=1,nzm
	   coef=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	   swgrad(k)=swleadv(k)-coef
	   swleadv(k)=coef
	  end do	

	  call hbuf_put('SWADV',swleadv,factor_xy)
	  call hbuf_put('SWDIFF',swlediff,factor_xy)
	  call hbuf_put('SWGRAD',swgrad,factor_xy)
	  call hbuf_put('SWBUOY',swlebuoy,factor_xy)
	  call hbuf_put('SWPRES',swlepres,factor_xy)
	  call hbuf_put('S2ADVTR',s2leadv,1.)
	  call hbuf_put('S2GRAD',s2legrad,1.)
	  call hbuf_put('S2DISSIP',s2lediss,1.)
	  call hbuf_put('S2DIFTR',s2lediff,1.)

	  do k=1,nzm
	   tmp(1) = dz/rhow(k)
	   tmp(2) = tmp(1) / dtn
	   tkewsb(k) = tkewsb(k) * tmp(1) * rhow(k)
	   tkewle(k) = tkewle(k)*tmp(2)*rhow(k) + tkewsb(k)
	  end do
	  call hbuf_put('SW',tkewle,factor_xy)
	  call hbuf_put('SWSB',tkewsb,factor_xy)
	  call hbuf_put('S2',s2z,1.e6*factor_xy)
	  call hbuf_avg_put('S',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)

	 else

	  call hbuf_put('SWADV',zero,factor_xy)
	  call hbuf_put('SWDIFF',zero,factor_xy)
	  call hbuf_put('SWGRAD',zero,factor_xy)
	  call hbuf_put('SWBUOY',zero,factor_xy)
	  call hbuf_put('SWPRES',zero,factor_xy)
	  call hbuf_put('S2ADVTR',zero,1.)
	  call hbuf_put('S2GRAD',zero,1.)
	  call hbuf_put('S2DISSIP',zero,1.)
	  call hbuf_put('S2DIFTR',zero,1.)
	  call hbuf_put('SW',zero,factor_xy)
	  call hbuf_put('SWSB',zero,factor_xy)
	  call hbuf_put('S',zero,factor_xy)
	  call hbuf_put('S2',zero,factor_xy)

	end if


!---------------------------------------------------------
! SGS TKE Budget:

	if(.not.dosmagor) then
 	 call hbuf_put('ADVTRS',tkewle,factor_xy)
	 call hbuf_avg_put('TKES',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)
	else
	 if(doscalar) then
 	   call hbuf_put('TKES',zero,factor_xy)
	 else
	   call hbuf_avg_put('TKES',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)
	 end if
 	 call hbuf_put('ADVTRS',zero,factor_xy)
	end if
	 call hbuf_put('BUOYAS',tkesbbuoy,factor_xy)
	 call hbuf_put('SHEARS',tkesbshear,factor_xy)
	 call hbuf_put('DISSIPS',tkesbdiss,factor_xy)
	
	
!-------------------------------------------------------------
!	Cloud statistics:
!-------------------------------------------------------------

	do k=1,nzm
	 cld(k) = 0.
	 hydro(k) = 0.
	 wcl(k) = 0.
	 ucl(k) = 0.
	 vcl(k) = 0.
	 wacl(k) = 0.
	 tcl(k) = 0.
	 tacl(k) = 0.
	 tvcl(k)= 0.
	 tvcla(k)= 0.
	 qcl(k) = 0.
	 qccl(k)= 0.
	 qicl(k)= 0.
	 qpcl(k)= 0.
	 tvwcl(k)= 0.
	 twcl(k)= 0.
	 swcl(k)= 0.
	 qwcl(k)= 0.
	 qcwcl(k)= 0.
	 qiwcl(k)= 0.
	 scl(k) = 0.
	 prof1(k)=0.
	 prof2(k)=0.
	 prof3(k)=0.
	 prof4(k)=0.
	 dse(k)=0.
	 mse(k)=0.
	 sse(k)=0.
	 if(LES) then
	  coef=0.
	 else
	  coef=1.3e-6/rho(k)
!!$	  coef=min(1.e-5,0.01*qsatw(tabs0(k),pres(k)))
	 endif
	 do j=1,ny
	  do i=1,nx
	    if(qn(i,j,k).gt.coef) then
              if(.not.dokruegermicro) then ! bloss
	      omn = omegan(tabs(i,j,k))
	      omp = omegap(tabs(i,j,k))
	      lstarn=fac_cond+(1.-omn)*fac_fus
	      lstarp=fac_cond+(1.-omp)*fac_fus
	      qcc=qn(i,j,k)*omn
	      qii=qn(i,j,k)*(1.-omn)
              elseif(.not.dokruegermicro) then ! bloss
               qcc=qc(i,j,k)
               qii=qi(i,j,k)
               qrr=qr(i,j,k)
               qss=qs(i,j,k)
               qgg=qg(i,j,k)
               if (qcc+qii.gt.0.) then
                omn = qcc/(qcc+qii)
               else
                omn = max(0.,min(1.,(tabs(i,j,k)-ttfrz)*a_bgkru))
               end if
               lstarn=fac_cond+(1.-omn)*fac_fus
               lstarp=fac_cond+fac_fus*(qss+qgg)/(qss+qgg+qrr+eps)
              end if
	      cld(k)=cld(k) + 1
	      hydro(k) = hydro(k) + 1
	      tmp(1)=0.5*(w(i,j,k+1)+w(i,j,k))
	      wcl(k) = wcl(k) + tmp(1)
	      ucl(k) = ucl(k) + u(i,j,k) 
	      vcl(k) = vcl(k) + v(i,j,k) 
	      if(tmp(1).gt.0.) then
		prof1(k)=prof1(k)+rho(k)*tmp(1)
	      else
	        prof2(k)=prof2(k)+rho(k)*tmp(1)
	      endif	
	      tmp(1)=t(i,j,k)+lstarn*qn(i,j,k)+lstarp*qp(i,j,k)
	      dse(k)=dse(k)+tmp(1)	
	      mse(k)=mse(k)+tmp(1)+fac_cond*(q(i,j,k)-qn(i,j,k))	
	      tcl(k) = tcl(k) + t(i,j,k)
	      qcl(k) = qcl(k) + q(i,j,k)  
	      scl(k) = scl(k) + tke(i,j,k)
	      qccl(k) = qccl(k) + qcc
	      qicl(k) = qicl(k) + qii
	      qpcl(k) = qpcl(k) + qp(i,j,k)
	      tvcl(k) = tvcl(k) + tvirt(i,j,k)	 
	      tvcla(k) = tvcla(k) + tvirt(i,j,k) - tvz(k)	 
	      tacl(k) = tacl(k) + tabs(i,j,k)	 
	      twcl(k) = twcl(k) + t(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      swcl(k) = swcl(k) + tke(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      qwcl(k) = qwcl(k) + q(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      tvwcl(k) = tvwcl(k)+tvirt(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      qcwcl(k) = qcwcl(k) + qcc*0.5*(w(i,j,k+1)+w(i,j,k))
	      qiwcl(k) = qiwcl(k) + qii*0.5*(w(i,j,k+1)+w(i,j,k))
	    elseif(qp(i,j,k).gt.1.e-4) then
	      hydro(k) = hydro(k) + 1
	      if(w(i,j,k)+w(i,j,k+1).lt.0.) &
 	         prof3(k)=prof3(k)+rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))      
	    endif
	  end do
	 end do
	 if(cld(k).gt.0.) then
	   wacl(k) = wcl(k)
	   wcl(k) = wcl(k)/cld(k)
	   ucl(k) = ucl(k)/cld(k)
	   vcl(k) = vcl(k)/cld(k)
	   dse(k)=dse(k)/cld(k)
	   mse(k)=mse(k)/cld(k)
	   tcl(k) = tcl(k)/cld(k)
	   tvcl(k) = tvcl(k)/cld(k)
	   tvcla(k) = tvcla(k)/cld(k)
	   tacl(k) = tacl(k)/cld(k)
	   qcl(k) = qcl(k)/cld(k)
	   qccl(k) = qccl(k)/cld(k)
	   qicl(k) = qicl(k)/cld(k)
	   qpcl(k) = qpcl(k)/cld(k)
	   scl(k) = scl(k)/cld(k)
	   cloud_factor(k) = cloud_factor(k)+1
	 endif

	 prof4(k)=prof1(k)+prof2(k)+prof3(k)

	end do
		
	call hbuf_put('CLD',cld,factor_xy)
	call hbuf_put('HYDRO',hydro,factor_xy)
	call hbuf_put('WCLD',wcl,1.)
	call hbuf_put('UCLD',ucl,1.)
	call hbuf_put('VCLD',vcl,1.)
	call hbuf_put('DSECLD',dse,1.)
	call hbuf_put('MSECLD',mse,1.)
	call hbuf_put('TLCLD',tcl,1.)
	call hbuf_put('TVCLD',tvcl,1.)
	call hbuf_put('TVCLDA',tvcla,1.)
	call hbuf_put('TACLD',tacl,1.)
	call hbuf_put('QTCLD',qcl,1.e3)
	call hbuf_put('QCCLD',qccl,1.e3)
	call hbuf_put('QICLD',qicl,1.e3)
	call hbuf_put('QPCLD',qpcl,1.e3)
	call hbuf_put('WCLDA',wacl,factor_xy)
	call hbuf_put('TLWCLD',twcl,factor_xy)
	call hbuf_put('TVWCLD',tvwcl,factor_xy)
	call hbuf_put('QTWCLD',qwcl,factor_xy*1.e3)
	call hbuf_put('QCWCLD',qcwcl,factor_xy*1.e3)
	call hbuf_put('QIWCLD',qiwcl,factor_xy*1.e3)
	call hbuf_put('MCUP',prof1,factor_xy)
	call hbuf_put('MCDNS',prof2,factor_xy)
	call hbuf_put('MCDNU',prof3,factor_xy)
	call hbuf_put('MC',prof4,factor_xy)

	if(doscalar) then
	 call hbuf_put('SCLD',tcl,1.)
	 call hbuf_put('SWCLD',swcl,factor_xy)
	else
	 call hbuf_put('SCLD',zero,1.)
	 call hbuf_put('SWCLD',zero,factor_xy)
	end if

!-------------------------------------------------------------
!	Updraft Core statistics:
!-------------------------------------------------------------


	do k=1,nzm
	 cld(k) = 0.
	 cldd(k) = 0.
	 wcl(k) = 0.
	 ucl(k) = 0.
	 vcl(k) = 0.
	 wacl(k) = 0.
	 tcl(k) = 0.
	 scl(k) = 0.
	 tvcl(k)= 0.
	 tvcla(k)= 0.
	 tacl(k)= 0.
	 qcl(k) = 0.
	 qccl(k)= 0.
	 qicl(k)= 0.
	 qpcl(k)= 0.
	 tvwcl(k)= 0.
	 twcl(k)= 0.
	 swcl(k)= 0.
	 qwcl(k)= 0.
	 qcwcl(k)= 0.
	 qiwcl(k)= 0.
	 dse(k)=0.
	 mse(k)=0.
	 prof1(k)=0.
         if(LES) then
          coef=0.
         else
	  coef=1.3e-6/rho(k)
!!$          coef=min(1.e-5,0.01*qsatw(tabs0(k),pres(k)))
         endif
	 do j=1,ny
	  do i=1,nx
	    condition_cl = qn(i,j,k).gt.coef
     	    condition = tvirt(i,j,k).gt.tvz(k) 
	    if(CEM) condition=condition.and.w(i,j,k)+w(i,j,k+1).gt.2.
	    if(condition) then
              if(.not.dokruegermicro) then ! bloss
	      omn = omegan(tabs(i,j,k))
	      omp = omegap(tabs(i,j,k))
	      lstarn=fac_cond+(1.-omn)*fac_fus
	      lstarp=fac_cond+(1.-omp)*fac_fus
	      qcc=qn(i,j,k)*omn
	      qii=qn(i,j,k)*(1.-omn)
              elseif(.not.dokruegermicro) then ! bloss
               qcc=qc(i,j,k)
               qii=qi(i,j,k)
               qrr=qr(i,j,k)
               qss=qs(i,j,k)
               qgg=qg(i,j,k)
               if (qcc+qii.gt.0.) then
                omn = qcc/(qcc+qii)
               else
                omn = max(0.,min(1.,(tabs(i,j,k)-ttfrz)*a_bgkru))
               end if
               lstarn=fac_cond+(1.-omn)*fac_fus
               lstarp=fac_cond+fac_fus*(qss+qgg)/(qss+qgg+qrr+eps)
              end if
	      cld(k)=cld(k) + 1
	      wcl(k) = wcl(k) + 0.5*(w(i,j,k+1)+w(i,j,k))
	      ucl(k) = ucl(k) + u(i,j,k) 
	      vcl(k) = vcl(k) + v(i,j,k) 
	      tmp(1)=t(i,j,k)+lstarn*qn(i,j,k)+lstarp*qp(i,j,k)
	      dse(k)=dse(k)+tmp(1)	
	      mse(k)=mse(k)+tmp(1)+fac_cond*(q(i,j,k)-qn(i,j,k))	
	      tcl(k) = tcl(k) + t(i,j,k)
	      scl(k) = scl(k) + tke(i,j,k)
	      qcl(k) = qcl(k) + q(i,j,k)  
	      qccl(k) = qccl(k) + qcc
	      qicl(k) = qicl(k) + qii  
	      qpcl(k) = qpcl(k) + qp(i,j,k)  
	      tvcl(k) = tvcl(k) + tvirt(i,j,k)
	      tvcla(k) = tvcla(k) + tvirt(i,j,k)-tvz(k)
	      tacl(k) = tacl(k) + tabs(i,j,k)
	      twcl(k) = twcl(k) + t(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      swcl(k) = swcl(k) + tke(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      qwcl(k) = qwcl(k) + q(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      tvwcl(k) = tvwcl(k)+tvirt(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      qcwcl(k) = qcwcl(k) + qcc*0.5*(w(i,j,k+1)+w(i,j,k))
	      qiwcl(k) = qiwcl(k) + qii*0.5*(w(i,j,k+1)+w(i,j,k))
	      prof1(k)=prof1(k)+rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))      
	      if(condition_cl) then
	        cldd(k)=cldd(k)+1
	      end if
	    endif
	  end do
	 end do
	 if(cld(k).gt.0.) then
	   wacl(k) = wcl(k)
	   wcl(k) = wcl(k)/cld(k)
	   ucl(k) = ucl(k)/cld(k)
	   vcl(k) = vcl(k)/cld(k)
	   dse(k)=dse(k)/cld(k)
	   mse(k)=mse(k)/cld(k)
	   tcl(k) = tcl(k)/cld(k)
	   scl(k) = scl(k)/cld(k)
	   qcl(k) = qcl(k)/cld(k)
	   qccl(k) = qccl(k)/cld(k)
	   qicl(k) = qicl(k)/cld(k)
	   qpcl(k) = qpcl(k)/cld(k)
	   tvcl(k) = tvcl(k)/cld(k)
	   tvcla(k) = tvcla(k)/cld(k)
	   tacl(k) = tacl(k)/cld(k)
	   core_factor(k) = core_factor(k)+1
	 endif
	end do
		
	call hbuf_put('CORE',cld,factor_xy)
	call hbuf_put('CORECL',cldd,factor_xy)
	call hbuf_put('WCOR',wcl,1.)
	call hbuf_put('UCOR',ucl,1.)
	call hbuf_put('VCOR',vcl,1.)
	call hbuf_put('DSECOR',dse,1.)
	call hbuf_put('MSECOR',mse,1.)
	call hbuf_put('TLCOR',tcl,1.)
	call hbuf_put('TVCOR',tvcl,1.)
	call hbuf_put('TVCORA',tvcla,1.)
	call hbuf_put('TACOR',tacl,1.)
	call hbuf_put('QTCOR',qcl,1.e3)
	call hbuf_put('QCCOR',qccl,1.e3)
	call hbuf_put('QICOR',qicl,1.e3)
	call hbuf_put('QPCOR',qpcl,1.e3)
	call hbuf_put('WCORA',wacl,factor_xy)
	call hbuf_put('TLWCOR',twcl,factor_xy)
	call hbuf_put('TVWCOR',tvwcl,factor_xy)
	call hbuf_put('QTWCOR',qwcl,factor_xy*1.e3)
	call hbuf_put('QCWCOR',qcwcl,factor_xy*1.e3)
	call hbuf_put('QIWCOR',qiwcl,factor_xy*1.e3)

	if(doscalar) then
	 call hbuf_put('SCOR',tcl,1.)
	 call hbuf_put('SWCOR',swcl,factor_xy)
	else
	 call hbuf_put('SCOR',zero,1.)
	 call hbuf_put('SWCOR',zero,factor_xy)
	end if
!-------------------------------------------------------------
!	Cloud Downdraft Core statistics:
!-------------------------------------------------------------


	do k=1,nzm
	 cld(k) = 0.
	 cldd(k) = 0.
	 wcl(k) = 0.
	 ucl(k) = 0.
	 vcl(k) = 0.
	 wacl(k) = 0.
	 tcl(k) = 0.
	 scl(k) = 0.
	 tvcl(k)= 0.
	 tacl(k)= 0.
	 qcl(k) = 0.
	 qccl(k)= 0.
	 qicl(k)= 0.
	 qpcl(k)= 0.
	 tvwcl(k)= 0.
	 twcl(k)= 0.
	 swcl(k)= 0.
	 qwcl(k)= 0.
	 qcwcl(k)= 0.
	 qiwcl(k)= 0.
	 dse(k)=0.
	 mse(k)=0.
	 prof2(k)=0.
	 prof3(k)=0.
         if(LES) then
          coef=0.
         else
	  coef=1.3e-6/rho(k)
!!$          coef=min(1.e-5,0.01*qsatw(tabs0(k),pres(k)))
         endif
	 do j=1,ny
	  do i=1,nx
	    condition_cl = qn(i,j,k).gt.coef .or. qp(i,j,k).gt.1.e-4 
     	    condition = tvirt(i,j,k).lt.tvz(k) 
	    if(CEM) condition=condition.and.w(i,j,k)+w(i,j,k+1).lt.-2.
	    if(condition) then
              if(.not.dokruegermicro) then ! bloss
	      omn = omegan(tabs(i,j,k))
	      omp = omegap(tabs(i,j,k))
	      lstarn=fac_cond+(1.-omn)*fac_fus
	      lstarp=fac_cond+(1.-omp)*fac_fus
	      qcc=qn(i,j,k)*omn
	      qii=qn(i,j,k)*(1.-omn)
              elseif(.not.dokruegermicro) then ! bloss
               qcc=qc(i,j,k)
               qii=qi(i,j,k)
               qrr=qr(i,j,k)
               qss=qs(i,j,k)
               qgg=qg(i,j,k)
               if (qcc+qii.gt.0.) then
                omn = qcc/(qcc+qii)
               else
                omn = max(0.,min(1.,(tabs(i,j,k)-ttfrz)*a_bgkru))
               end if
               lstarn=fac_cond+(1.-omn)*fac_fus
               lstarp=fac_cond+fac_fus*(qss+qgg)/(qss+qgg+qrr+eps)
              end if
	      cld(k)=cld(k) + 1
	      wcl(k) = wcl(k) + 0.5*(w(i,j,k+1)+w(i,j,k))
	      ucl(k) = ucl(k) + u(i,j,k) 
	      vcl(k) = vcl(k) + v(i,j,k) 
	      tmp(1)=t(i,j,k)+lstarn*qn(i,j,k)+lstarp*qp(i,j,k)
	      dse(k)=dse(k)+tmp(1)	
	      mse(k)=mse(k)+tmp(1)+fac_cond*(q(i,j,k)-qn(i,j,k))	
	      tcl(k) = tcl(k) + t(i,j,k)
	      scl(k) = scl(k) + tke(i,j,k)
	      qcl(k) = qcl(k) + q(i,j,k)  
	      qccl(k) = qccl(k) + qcc
	      qicl(k) = qicl(k) + qii  
	      qpcl(k) = qpcl(k) + qp(i,j,k)  
	      tvcl(k) = tvcl(k) + tvirt(i,j,k)
	      tacl(k) = tacl(k) + tabs(i,j,k)
	      twcl(k) = twcl(k) + t(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      swcl(k) = swcl(k) + tke(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      qwcl(k) = qwcl(k) + q(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      tvwcl(k) = tvwcl(k)+tvirt(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      qcwcl(k) = qcwcl(k) + qcc*0.5*(w(i,j,k+1)+w(i,j,k))
	      qiwcl(k) = qiwcl(k) + qii*0.5*(w(i,j,k+1)+w(i,j,k))
	      if(condition_cl) then
	        prof2(k)=prof2(k)+rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))      
	        cldd(k)=cldd(k) + 1
	      else
	        prof3(k)=prof3(k)+rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))      
	      end if
	    endif
	  end do
	 end do
	 if(cld(k).gt.0.) then
	   wacl(k) = wcl(k)
	   wcl(k) = wcl(k)/cld(k)
	   ucl(k) = ucl(k)/cld(k)
	   vcl(k) = vcl(k)/cld(k)
	   dse(k)=dse(k)/cld(k)
	   mse(k)=mse(k)/cld(k)
	   tcl(k) = tcl(k)/cld(k)
	   scl(k) = scl(k)/cld(k)
	   qcl(k) = qcl(k)/cld(k)
	   qccl(k) = qccl(k)/cld(k)
	   qicl(k) = qicl(k)/cld(k)
	   qpcl(k) = qpcl(k)/cld(k)
	   tvcl(k) = tvcl(k)/cld(k)
	   tacl(k) = tacl(k)/cld(k)
	   coredn_factor(k) = coredn_factor(k)+1
	 endif

	 prof4(k)=prof1(k)+prof2(k)+prof3(k)

	end do
		
	call hbuf_put('COREDN',cld,factor_xy)
	call hbuf_put('COREDNCL',cldd,factor_xy)
	call hbuf_put('UCORDN',ucl,1.)
	call hbuf_put('VCORDN',vcl,1.)
	call hbuf_put('WCORDN',wcl,1.)
	call hbuf_put('DSECORDN',dse,1.)
	call hbuf_put('MSECORDN',mse,1.)
	call hbuf_put('TLCORDN',tcl,1.)
	call hbuf_put('TVCORDN',tvcl,1.)
	call hbuf_put('TACORDN',tacl,1.)
	call hbuf_put('QTCORDN',qcl,1.e3)
	call hbuf_put('QCCORDN',qccl,1.e3)
	call hbuf_put('QICORDN',qicl,1.e3)
	call hbuf_put('QPCORDN',qpcl,1.e3)
	call hbuf_put('WCORDNA',wacl,factor_xy)
	call hbuf_put('TLWCORDN',twcl,factor_xy)
	call hbuf_put('TVWCORDN',tvwcl,factor_xy)
	call hbuf_put('QTWCORDN',qwcl,factor_xy*1.e3)
	call hbuf_put('QCWCORDN',qcwcl,factor_xy*1.e3)
	call hbuf_put('QIWCORDN',qiwcl,factor_xy*1.e3)
	call hbuf_put('MCRUP',prof1,factor_xy)
	call hbuf_put('MCRDNS',prof2,factor_xy)
	call hbuf_put('MCRDNU',prof3,factor_xy)
	call hbuf_put('MCR',prof4,factor_xy)

	if(doscalar) then
	 call hbuf_put('SCORDN',tcl,1.)
	 call hbuf_put('SWCORDN',swcl,factor_xy)
	else
	 call hbuf_put('SCORDN',zero,1.)
	 call hbuf_put('SWCORDN',zero,factor_xy)
	end if
!---------------------------------------------------------
!  Radiation and other stuff

	do j=1,ny
	 do i=1,nx
	   cwp(i,j)=0.
	   cwpl(i,j)=0.
	   cwpm(i,j)=0.
	   cwph(i,j)=0.
	   topind(i,j)=1
	   tau_isccp(i,j)=0.
	 end do
	end do
	
	if(CEM) then
	  cwpmax=0.02
	else
	  cwpmax=0.0
	endif

	do k=nzm,1,-1
	 prof1(k)=(radqrlw(k)+radqrsw(k))*factor_xy
	 tmp(1)=rho(k)*adzw(k)*dz
	 do j=1,ny
	  do i=1,nx
	    cwp(i,j)=cwp(i,j)+tmp(1)*qn(i,j,k)
            if(pres(k).ge.700.) then
	      cwpl(i,j)=cwpl(i,j)+tmp(1)*qn(i,j,k)
            else if(pres(k).le.400.) then
	      cwph(i,j)=cwph(i,j)+tmp(1)*qn(i,j,k)
            else
	      cwpm(i,j)=cwpm(i,j)+tmp(1)*qn(i,j,k)
	    end if
	    if(cwp(i,j).gt.cwpmax.and.topind(i,j).eq.1)topind(i,j)=k
            if(.not.dokruegermicro) then ! bloss
             omn = omegan(tabs(i,j,k))
             tau_isccp(i,j) = tau_isccp(i,j)+ tmp(1)*qn(i,j,k)*1000.* &
               ( (2.817e-02 + 1.305/10.)*omn +  &
               (3.448e-03 + 2.431/25.)*(1.-omn))
            elseif(dokruegermicro) then ! bloss
             tau_isccp(i,j) = tau_isccp(i,j)+ tmp(1)*1000.* &
             ( (2.817e-02 + 1.305/10.)*qc(i,j,k) +  &
               (3.448e-03 + 2.431/25.)*qi(i,j,k))
            end if
	  end do
	 end do
	end do

	do j=1,ny
	 do i=1,nx
            !bloss: Use isccp approximation from SAM 6.2 for shaded cloud
            ! fraction so that both day and night data are available
            ! to compare (roughly) to the ISCCP OLR-based cloud
            ! fraction product
!bloss           if(cwp(i,j).gt.cwpmax) s_acld=s_acld+1.
	   if(tau_isccp(i,j).gt.0.3) s_acld=s_acld+1.
	   if(cwpl(i,j).gt.cwpmax) s_acldl=s_acldl+1.
	   if(cwpm(i,j).gt.cwpmax) s_acldm=s_acldm+1.
	   if(cwph(i,j).gt.cwpmax) s_acldh=s_acldh+1.
	   if(tabs(i,j,topind(i,j)).lt.245.) s_acldcold=s_acldcold+1
	 end do
	end do

	do k=1,nzm
	 prof2(k)=0.
	 prof3(k)=0.	 
	 n=0
	 if(dolongwave.or.doshortwave) then
	   do j=1,ny
	     do i=1,nx
	       if(cwp(i,j).gt.cwpmax) then
	         n=n+1
	         prof2(k)=prof2(k)+misc(i,j,k)
	       else
	         prof3(k)=prof3(k)+misc(i,j,k)
	       endif
	     end do
	   end do 
	 end if
	end do


	call hbuf_put('HLADV',tadv,factor_xy*86400./dtn)
	call hbuf_put('HLDIFF',tdiff,factor_xy*86400./dtn)
	call hbuf_put('HLLAT',tlat+tlatqi,factor_xy*86400./dtn)
	call hbuf_put('HLRAD',prof1,86400.)
	call hbuf_put('QTADV',qadv+qifall,factor_xy*86400000./dtn)
	call hbuf_put('QTDIFF',qdiff,factor_xy*86400000./dtn)
	call hbuf_put('QTSINK',qpsrc,-factor_xy*86400000./dtn)
	call hbuf_put('QTSRC',qpevp,-factor_xy*86400000./dtn)
        if (dokruegermicro) then
           qpfall = qrfall + qgfall + qsfall
           qpadv = qradv + qgadv + qsadv
           qpdiff = qrdiff + qgdiff + qsdiff
        end if
	call hbuf_put('QPADV',qpadv,factor_xy*86400000./dtn)
	call hbuf_put('QPDIFF',qpdiff,factor_xy*86400000./dtn)
	call hbuf_put('QPFALL',qpfall,factor_xy*86400000./dtn)
	call hbuf_put('QPSRC',qpsrc,factor_xy*86400000./dtn)
	call hbuf_put('QPEVP',qpevp,factor_xy*86400000.)


	call hbuf_put('RADLWUP',radlwup,factor_xy)
	call hbuf_put('RADLWDN',radlwdn,factor_xy)
	call hbuf_put('RADSWUP',radswup,factor_xy)
	call hbuf_put('RADSWDN',radswdn,factor_xy)
	call hbuf_put('RADQRLW',radqrlw,factor_xy*86400.) 	
	call hbuf_put('RADQRSW',radqrsw,factor_xy*86400.) 	
	call hbuf_put('RADQR',prof1,86400.) 	
	call hbuf_put('RADQRC',prof2,86400./(n+1.e-5)) 	
	call hbuf_put('RADQRS',prof3,86400./(nx*ny-n+1.e-5))

!---------------------------------------------------------
!  Apparent heat/moisture sources/sinks

	tmp(1)=1./dtn

	do k=1,nzm
	 prof2(k)=0.
	 prof3(k)=0.	 
	 n=0
	 do j=1,ny
	  do i=1,nx
	    prof2(k)=prof2(k)+(tabs(i,j,k)-t01(k))*tmp(1)-ttend(k)-ttend_wave(k)-prof1(k)
	    prof3(k)=prof3(k)-fac_cond*((q(i,j,k)-q01(k))*tmp(1)-qtend(k)-qtend_wave(k))
	  end do
	 end do 
	end do

	call hbuf_put('Q1C',prof2,factor_xy*86400.) 	
	call hbuf_put('Q2',prof3,factor_xy*86400.) 	

        call isccp_get()

!===================================================
! UW ADDITIONS
        
        !bloss: extra nudging/vertical large-scale advective tendency outputs
	call hbuf_put('TVTEND',tlsvadv,86400.)
        call hbuf_put('QVTEND',qlsvadv,86400.*1.e3)
        call hbuf_put('TNUDGE',tnudge,86400.)
        call hbuf_put('QNUDGE',qnudge,86400.*1.e3)

        !bloss: extra radiation outputs, clearsky sw/lw heating rates
	call hbuf_put('RADQRCLW',radqrclw,factor_xy*86400.) 	
	call hbuf_put('RADQRCSW',radqrcsw,factor_xy*86400.) 	

        !cw: uwmf from UW PBL scheme
	call hbuf_avg_put('UMF', umf, 1, nx, 1, ny, nzm, 1.)

        !cw: uwmf from UW PBL scheme
	call hbuf_avg_put('SLFLX', slflx_uwcu, 1, nx, 1, ny, nzm, 1.)

        !cw: uwmf from UW PBL scheme
	call hbuf_avg_put('QTFLX', qtflx_uwcu, 1, nx, 1, ny, nzm, 1.)

        !cw: tke_uw from UW PBL scheme
	call hbuf_avg_put('TKE_UW', tke_uw, 1, nx, 1, ny, nz, 1.)

        !cw: ttend from UW PBL scheme
	call hbuf_avg_put('TT_UWPBL', ttend_uwpbl, 1, nx, 1, ny, nz, 1.)
	call hbuf_avg_put('QT_UWPBL', qtend_uwpbl, 1, nx, 1, ny, nz, 1.)

        !cw: ttend from UW CU scheme
	call hbuf_avg_put('TTEND_UW', ttend_uwcu, 1, nx, 1, ny, nzm, 1.)

        !cw: qtend from UW CU scheme
	call hbuf_avg_put('QVTEN_UW', qvtend_uwcu, 1, nx, 1, ny, nzm, 1.)

        !cw: qtend from UW CU scheme
	call hbuf_avg_put('QLTEN_UW', qltend_uwcu, 1, nx, 1, ny, nzm, 1.)

        !cw: qtend from UW CU scheme
	call hbuf_avg_put('QITEN_UW', qitend_uwcu, 1, nx, 1, ny, nzm, 1.)

        !cw: cloud fraction from UW CU scheme
	call hbuf_avg_put('CFRC_UW', cldfrc, 1, nx, 1, ny, nzm, 1.)

	call hbuf_avg_put('TK_Z_UW', tk_z_uw, 1, nx, 1,ny, nzm,1.)
	call hbuf_avg_put('TKH_Z_UW', tkh_z_uw, 1, nx, 1,ny, nzm,1.)

        !kzm compute mean profile and variance of trx
	if(dotrx) then 
           do k=1,nzm
              tr0(k)=0.
              do j=1,ny
                 do i=1,nx
                    tr0(k)=tr0(k)+trx(i,j,k)
                 end do
              end do
              tr2(k)=0.
              do j=1,ny
                 do i=1,nx
                    tr2(k)=tr2(k)+(trx(i,j,k)-tr0(k))**2
                 end do
              end do
           end do
           call hbuf_put('TRX',tr0,factor_xy*1.e3)
           call hbuf_put('TRX2',tr2,factor_xy*1.e6)
	endif

        !kzm compute mean profile and variance of try
	if(dotry) then 
           do k=1,nzm
              tr0(k)=0.
              do j=1,ny
                 do i=1,nx
                    tr0(k)=tr0(k)+try(i,j,k)
                 end do
              end do
              tr2(k)=0.
              do j=1,ny
                 do i=1,nx
                    tr2(k)=tr2(k)+(try(i,j,k)-tr0(k))**2
                 end do
              end do
           end do
           call hbuf_put('TRY',tr0,factor_xy*1.e3)
           call hbuf_put('TRY2',tr2,factor_xy*1.e6)
	endif

        !kzm compute mean profile and variance of trz
	if(dotrz) then 
           do k=1,nzm
              tr0(k)=0.
              do j=1,ny
                 do i=1,nx
                    tr0(k)=tr0(k)+trz(i,j,k)
                 end do
              end do
              tr2(k)=0.
              do j=1,ny
                 do i=1,nx
                    tr2(k)=tr2(k)+(trz(i,j,k)-tr0(k))**2
                 end do
              end do
           end do
           call hbuf_put('TRZ',tr0,factor_xy*1.e3)
           call hbuf_put('TRZ2',tr2,factor_xy*1.e6)
	endif

        !kzm compute mean profile and variance of trzz
	if(dotrzz) then 
           do k=1,nzm
              tr0(k)=0.
              do j=1,ny
                 do i=1,nx
                    tr0(k)=tr0(k)+trzz(i,j,k)
                 end do
              end do
              tr2(k)=0.
              do j=1,ny
                 do i=1,nx
                    tr2(k)=tr2(k)+(trzz(i,j,k)-tr0(k))**2
                 end do
              end do
           end do
           call hbuf_put('TRT',tr0,factor_xy*1.e3)
           call hbuf_put('TRT2',tr2,factor_xy*1.e6)
	endif

        !kzm compute mean profile and variance of tro3
	if(dotro3) then 
           do k=1,nzm
              tr0(k)=0.
              do j=1,ny
                 do i=1,nx
                    tr0(k)=tr0(k)+tro3(i,j,k)
                 end do
              end do
              tr2(k)=0.
              do j=1,ny
                 do i=1,nx
                    tr2(k)=tr2(k)+(tro3(i,j,k)-tr0(k))**2
                 end do
              end do
           end do
           call hbuf_put('TRO3',tr0,factor_xy*1.e6)
           call hbuf_put('TRO32',tr2,factor_xy*1.e12)
	endif


        if (doreflectivity) then

!---------------------------------------------------------
!  Equivalent Radar reflectivity:
!     -- contoured frequency with altitude diagrams (CFAD's)
!     -- slices at 1km, 3km and 6 km for output with 2D slices.
        
        pii = 3.141592654
        fac_ice = 0.197/0.93 ! Reduction of equivalent reflectivity for ice.

        ! Coefficients which encapsulate effect of particle type
        ! on reflectivity. 

        ! Based on taking sixth moment of exponential drop size
        ! distribution for rain with fixed intercept parameter nzeror.
        ! Note that this coefficient will be multiplied by the density
        ! of rain (in grams per cubic meter of air) raised to the 1.75.
        tmpr = 7.2e20*nzeror/(1000*pii*rhor*nzeror)**(1.75)

        ! Based on taking sixth moment of exponential drop size
        ! distribution with fixed intercept parameter nzeror/nzerog.
        ! However, using the formulas from Heymsfield et al (2002)
        ! the sixth power of the equivalent melter diameter is
        ! used in place of the diameter itself.
        tmp_mu       = 0. ! Parameter for gamma distn -- zero for
        ! exponential distribution (as assumed here)
        tmp_exponent = (5.5+tmp_mu)/(3.2+tmp_mu)
        ! Note that the units of nzero are converted to cgs in the following:
        tmpg = 1.2e8*gammafff(5.5+tmp_mu)*(1.e-8*nzerog)**(1.-tmp_exponent) &
             /(5.7e3*gammafff(3.2+tmp_mu))**(tmp_exponent)
        tmps = 1.2e8*gammafff(5.5+tmp_mu)*(1.e-8*nzeros)**(1.-tmp_exponent) &
             /(5.7e3*gammafff(3.2+tmp_mu))**(tmp_exponent)

	do k=1,nzm
           kb = max(k-1,1)
           kc = min(k+1,nzm)
           cfadm10(k) = 0.
           cfadm05(k) = 0.
           cfad00(k) = 0.
           cfad05(k) = 0.
           cfad10(k) = 0.
           cfad15(k) = 0.
           cfad20(k) = 0.
           cfad25(k) = 0.
           cfad30(k) = 0.
           cfad35(k) = 0.
           cfad40(k) = 0.
           cfad45(k) = 0.
           cfad50(k) = 0.
           cfad55(k) = 0.
           cfad60(k) = 0.
           cfad65(k) = 0.
           do j=1,ny
              do i=1,nx
                 omp = omegap(tabs(i,j,k))
                 omg = omegag(tabs(i,j,k))
                 qrr=qp(i,j,k)*omp
                 qss=qp(i,j,k)*(1.-omp)*(1.-omg)
                 qgg=qp(i,j,k)*(1.-omp)*omg
                 dbZ(i,j) = 10.*log10(tmpr*(rho(k)*qrr*1.e3)**(1.75) &
                           + fac_ice*(tmpg*(rho(k)*qgg*1.e3)**(tmp_exponent) &
                                    + tmps*(rho(k)*qss*1.e3)**(tmp_exponent))&
                                    + 1.e-8)
                 select case (floor(dbZ(i,j)+2.5))
                 case (-10:-6)
                    cfadm10(k) = cfadm10(k) + 1.
                 case (-5:-1)
                    cfadm05(k) = cfadm05(k) + 1.
                 case (0:4)
                    cfad00(k) = cfad00(k) + 1.
                 case (5:9)
                    cfad05(k) = cfad05(k) + 1.
                 case (10:14)
                    cfad10(k) = cfad10(k) + 1.
                 case (15:19)
                    cfad15(k) = cfad15(k) + 1.
                 case (20:24)
                    cfad20(k) = cfad20(k) + 1.
                 case (25:29)
                    cfad25(k) = cfad25(k) + 1.
                 case (30:34)
                    cfad30(k) = cfad30(k) + 1.
                 case (35:39)
                    cfad35(k) = cfad35(k) + 1.
                 case (40:44)
                    cfad40(k) = cfad40(k) + 1.
                 case (45:49)
                    cfad45(k) = cfad45(k) + 1.
                 case (50:54)
                    cfad50(k) = cfad50(k) + 1.
                 case (55:59)
                    cfad55(k) = cfad55(k) + 1.
                 case (60:64)
                    cfad60(k) = cfad60(k) + 1.
                 case (65:69)
                    cfad65(k) = cfad65(k) + 1.
                 case default
                 end select
              end do
           end do

           ! Extract slices of reflectivity at 1km, 3km and 6km
           if ((z(k).lt.1000.).and.(z(kc).ge.1000.)) dbZe1km_xy = dbZ
           if ((z(k).lt.3000.).and.(z(kc).ge.3000.)) dbZe3km_xy = dbZ
           if ((z(k).lt.6000.).and.(z(kc).ge.6000.)) dbZe6km_xy = dbZ

	end do

	call hbuf_put('DBZEM10',cfadm10,factor_xy) 	
	call hbuf_put('DBZEM05',cfadm05,factor_xy) 	
	call hbuf_put('DBZE00',cfad00,factor_xy) 	
	call hbuf_put('DBZE05',cfad05,factor_xy) 	
	call hbuf_put('DBZE10',cfad10,factor_xy) 	
	call hbuf_put('DBZE15',cfad15,factor_xy) 	
	call hbuf_put('DBZE20',cfad20,factor_xy) 	
	call hbuf_put('DBZE25',cfad25,factor_xy) 	
	call hbuf_put('DBZE30',cfad30,factor_xy) 	
	call hbuf_put('DBZE35',cfad35,factor_xy) 	
	call hbuf_put('DBZE40',cfad40,factor_xy) 	
	call hbuf_put('DBZE45',cfad45,factor_xy) 	
	call hbuf_put('DBZE50',cfad50,factor_xy) 	
	call hbuf_put('DBZE55',cfad55,factor_xy) 	
	call hbuf_put('DBZE60',cfad60,factor_xy) 	
	call hbuf_put('DBZE65',cfad65,factor_xy) 	

        ! Compute hydrometeor mixture fractions based on strength of
        ! radar reflectivity at 6 km.  The notion here is that 
        ! convective regions with significant ice microphysics will
        ! have high reflectivity at 6 km, while stratiform regions
        ! will have weaker reflectivity.  We will use four bins:
        !
        ! dbZe \in [10,20), [20,30), [30,40) and [40,inf)
        prof4 = 0.
        do n = 1,4
           !initialize conditional statistics
           prof1 = 0.
           prof2 = 0.
           prof3 = 0.
           q2z = 0.
           dse=0.
           mse=0.
           sse=0.
           tpz = 0.
           tlz = 0.
           tvz = 0.
           taz = 0.
           tez = 0.
           qvz = 0.
           qcz = 0.
           qiz = 0.
           qrz = 0.
           qsz = 0.
           qgz = 0.
           qsatwz=0.
           relhz=0.

           wcl = 0.
           scl = 0.
           twcl = 0.
           qwcl = 0.
           tvwcl = 0.

           ! set up thresholds for each dbZe bin
           if (n.eq.1) then
              threshmin = 10.
              threshmax = 20.
           elseif (n.eq.2) then
              threshmin = 20.
              threshmax = 30.
           elseif (n.eq.3) then
              threshmin = 30.
              threshmax = 40.
           elseif (n.eq.4) then
              threshmin = 40.
              threshmax = 1.e6
           end if

           do j = 1,ny
              do i = 1,nx

                 ! accumulate statistics if threshmin < dbZe < threshmax
                 if ((dbZe6km_xy(i,j).gt.threshmin).and. &
                      (dbZe6km_xy(i,j).le.threshmax)) then

                    prof4(n) = prof4(n) + precsfc(i,j) ! binned precip

                    do k = 1,nzm
                       prof1(k) = prof1(k) + 1.
                       prof2(k) = prof2(k) + misc(i,j,k)
!!$                       q2z(k)   = q2z(k)   + misc2(i,j,k)

                       if(.not.dokruegermicro) then !bloss
                          omn = omegan(tabs(i,j,k))
                          omp = omegap(tabs(i,j,k))
                          omg = omegag(tabs(i,j,k))
                          lstarn=fac_cond+(1.-omn)*fac_fus
                          lstarp=fac_cond+(1.-omp)*fac_fus
                          qvv=q(i,j,k)-qn(i,j,k)
                          qcc=qn(i,j,k)*omn
                          qii=qn(i,j,k)*(1.-omn)
                          qrr=qp(i,j,k)*omp
                          qss=qp(i,j,k)*(1.-omp)*(1.-omg)
                          qgg=qp(i,j,k)*(1.-omp)*omg
                       else !bloss: dokruegermicro = .true.
                          qcc=qc(i,j,k)
                          qii=qi(i,j,k)
                          qrr=qr(i,j,k)
                          qss=qs(i,j,k)
                          qgg=qg(i,j,k)
                          if (qcc+qii.gt.0.) then
                             omn = qcc/(qcc+qii)
                          else
                             omn = max(0.,min(1.,(tabs(i,j,k)-ttfrz)*a_bgkru))
                          end if
                          lstarn=fac_cond+(1.-omn)*fac_fus
                          lstarp=fac_cond+fac_fus*(qss+qgg)/(qss+qgg+qrr+eps)
                       end if

                       qvz(k)=qvz(k)+qvv
                       qrz(k)=qrz(k)+qrr
                       qsz(k)=qsz(k)+qss
                       qgz(k)=qgz(k)+qgg
                       qcz(k)=qcz(k)+qcc
                       qiz(k)=qiz(k)+qii
                       tmp(1)=tabs(i,j,k)*prespot(k)
                       taz(k)=taz(k)+tabs(i,j,k)
                       tpz(k)=tpz(k)+tmp(1)
                       tlz(k)=tlz(k)+tmp(1)*(1.-fac_cond*(qcc+qii)/tabs(i,j,k))
                       tmp(3) = tmp(1)*(1.+0.61*qvv-(qcc+qii+qrr+qss+qgg))
                       tvz(k)=tvz(k) + tmp(3)
                       tez(k)=tez(k)+t(i,j,k) &
                            +fac_cond*(qvv+qcc+qii+qrr+qgg+qss)
                       dse(k)=dse(k) &
                            +t(i,j,k)+fac_cond*(qcc+qrr)+fac_sub*(qii+qgg+qss)
                       mse(k)=mse(k) &
                            +t(i,j,k)+fac_cond*(qcc+qrr)+fac_sub*(qii+qgg+qss)&
                            +fac_cond*qvv

                       qsat = omn*qsatw(tabs(i,j,k),pres(k))+ &
                            (1.-omn)*qsati(tabs(i,j,k),pres(k)) 
                       sse(k)=sse(k) &
                            +t(i,j,k)+fac_cond*(qcc+qrr)+fac_sub*(qii+qgg+qss)&
                            +fac_cond*qsat
                       qsatwz(k) = qsatwz(k)+qsatw(tabs(i,j,k),pres(k))
                       relhz(k)=relhz(k)+qvv/qsatw(tabs(i,j,k),pres(k))

                       tmp(2) = 0.5*(w(i,j,k+1)+w(i,j,k))
                       wcl(k) = wcl(k) + tmp(2)
                       scl(k) = scl(k) + tke(i,j,k)
                       twcl(k) = twcl(k) + tmp(2)*t(i,j,k)
                       qwcl(k) = qwcl(k) + tmp(2)*(qvv+qcc+qii)
                       tvwcl(k) = tvwcl(k) + tmp(2)*tmp(3)
                       prof3(k) = prof3(k) + tmp(2)*rho(k)
                    end do
                 end if
              end do
           end do
           
           if (n.eq.1) then
              call hbuf_put('FRAC10',prof1,factor_xy)
              call hbuf_put('QRAD10',prof2,factor_xy*86400.)
!!$              call hbuf_put('QLAT10',q2z,factor_xy*86400.)
              call hbuf_put('QV10',qvz,1.e3*factor_xy)
              call hbuf_put('QC10',qcz,1.e3*factor_xy)
              call hbuf_put('QI10',qiz,1.e3*factor_xy)
              call hbuf_put('QR10',qrz,1.e3*factor_xy)
              call hbuf_put('QS10',qsz,1.e3*factor_xy)
              call hbuf_put('QG10',qgz,1.e3*factor_xy)
              call hbuf_put('QSAT10',qsatwz,1.e3*factor_xy)
              call hbuf_put('RELH10',relhz,100.*factor_xy)

              call hbuf_put('DSE10',dse,factor_xy)
              call hbuf_put('MSE10',mse,factor_xy)
              call hbuf_put('SSE10',sse,factor_xy)
              call hbuf_put('TABS10',taz,factor_xy)
              call hbuf_put('THETA10',tpz,factor_xy)
              call hbuf_put('THETAL10',tlz,factor_xy)
              call hbuf_put('THETAV10',tvz,factor_xy)
              call hbuf_put('THETAE10',tez,factor_xy)

              call hbuf_put('W10',wcl,factor_xy)
              call hbuf_put('TKE10',scl,factor_xy)
              call hbuf_put('TLW10',twcl,factor_xy)
              call hbuf_put('TVW10',tvwcl,factor_xy)
              call hbuf_put('QTW10',qwcl,factor_xy*1.e3)
              call hbuf_put('MFLX10',prof3,factor_xy)
           elseif (n.eq.2) then
              call hbuf_put('FRAC20',prof1,factor_xy)
              call hbuf_put('QRAD20',prof2,factor_xy*86400.)
!!$              call hbuf_put('QLAT20',q2z,factor_xy*86400.)
              call hbuf_put('QV20',qvz,1.e3*factor_xy)
              call hbuf_put('QC20',qcz,1.e3*factor_xy)
              call hbuf_put('QI20',qiz,1.e3*factor_xy)
              call hbuf_put('QR20',qrz,1.e3*factor_xy)
              call hbuf_put('QS20',qsz,1.e3*factor_xy)
              call hbuf_put('QG20',qgz,1.e3*factor_xy)
              call hbuf_put('QSAT20',qsatwz,1.e3*factor_xy)
              call hbuf_put('RELH20',relhz,100.*factor_xy)

              call hbuf_put('DSE20',dse,factor_xy)
              call hbuf_put('MSE20',mse,factor_xy)
              call hbuf_put('SSE20',sse,factor_xy)
              call hbuf_put('TABS20',taz,factor_xy)
              call hbuf_put('THETA20',tpz,factor_xy)
              call hbuf_put('THETAL20',tlz,factor_xy)
              call hbuf_put('THETAV20',tvz,factor_xy)
              call hbuf_put('THETAE20',tez,factor_xy)

              call hbuf_put('W20',wcl,factor_xy)
              call hbuf_put('TKE20',scl,factor_xy)
              call hbuf_put('TLW20',twcl,factor_xy)
              call hbuf_put('TVW20',tvwcl,factor_xy)
              call hbuf_put('QTW20',qwcl,factor_xy*1.e3)
              call hbuf_put('MFLX20',prof3,factor_xy)
           elseif (n.eq.3) then
              call hbuf_put('FRAC30',prof1,factor_xy)
              call hbuf_put('QRAD30',prof2,factor_xy*86400.)
!!$              call hbuf_put('QLAT30',q2z,factor_xy*86400.)
              call hbuf_put('QV30',qvz,1.e3*factor_xy)
              call hbuf_put('QC30',qcz,1.e3*factor_xy)
              call hbuf_put('QI30',qiz,1.e3*factor_xy)
              call hbuf_put('QR30',qrz,1.e3*factor_xy)
              call hbuf_put('QS30',qsz,1.e3*factor_xy)
              call hbuf_put('QG30',qgz,1.e3*factor_xy)
              call hbuf_put('QSAT30',qsatwz,1.e3*factor_xy)
              call hbuf_put('RELH30',relhz,100.*factor_xy)

              call hbuf_put('DSE30',dse,factor_xy)
              call hbuf_put('MSE30',mse,factor_xy)
              call hbuf_put('SSE30',sse,factor_xy)
              call hbuf_put('TABS30',taz,factor_xy)
              call hbuf_put('THETA30',tpz,factor_xy)
              call hbuf_put('THETAL30',tlz,factor_xy)
              call hbuf_put('THETAV30',tvz,factor_xy)
              call hbuf_put('THETAE30',tez,factor_xy)

              call hbuf_put('W30',wcl,factor_xy)
              call hbuf_put('TKE30',scl,factor_xy)
              call hbuf_put('TLW30',twcl,factor_xy)
              call hbuf_put('TVW30',tvwcl,factor_xy)
              call hbuf_put('QTW30',qwcl,factor_xy*1.e3)
              call hbuf_put('MFLX30',prof3,factor_xy)
           elseif (n.eq.4) then
              call hbuf_put('FRAC40',prof1,factor_xy)
              call hbuf_put('QRAD40',prof2,factor_xy*86400.)
!!$              call hbuf_put('QLAT40',q2z,factor_xy*86400.)
              call hbuf_put('QV40',qvz,1.e3*factor_xy)
              call hbuf_put('QC40',qcz,1.e3*factor_xy)
              call hbuf_put('QI40',qiz,1.e3*factor_xy)
              call hbuf_put('QR40',qrz,1.e3*factor_xy)
              call hbuf_put('QS40',qsz,1.e3*factor_xy)
              call hbuf_put('QG40',qgz,1.e3*factor_xy)
              call hbuf_put('QSAT40',qsatwz,1.e3*factor_xy)
              call hbuf_put('RELH40',relhz,100.*factor_xy)

              call hbuf_put('DSE40',dse,factor_xy)
              call hbuf_put('MSE40',mse,factor_xy)
              call hbuf_put('SSE40',sse,factor_xy)
              call hbuf_put('TABS40',taz,factor_xy)
              call hbuf_put('THETA40',tpz,factor_xy)
              call hbuf_put('THETAL40',tlz,factor_xy)
              call hbuf_put('THETAV40',tvz,factor_xy)
              call hbuf_put('THETAE40',tez,factor_xy)

              call hbuf_put('W40',wcl,factor_xy)
              call hbuf_put('TKE40',scl,factor_xy)
              call hbuf_put('TLW40',twcl,factor_xy)
              call hbuf_put('TVW40',tvwcl,factor_xy)
              call hbuf_put('QTW40',qwcl,factor_xy*1.e3)
              call hbuf_put('MFLX40',prof3,factor_xy)
           end if
        end do
        call hbuf_put('BINPCP',prof4,factor_xy*24.)

        end if

        if (dokruegermicro) then
           !bloss output additional statistics for lin et al microphysics.
           call hbuf_put('RFRAC',rainfrac,factor_xy)
           call hbuf_put('SFRAC',snowfrac,factor_xy)
           call hbuf_put('GFRAC',graufrac,factor_xy)
           call hbuf_put('LFRAC',clwfrac,factor_xy)
           call hbuf_put('IFRAC',clifrac,factor_xy)

           call hbuf_put('VTR',vtrz,factor_xy)
           call hbuf_put('VTS',vtsz,factor_xy)
           call hbuf_put('VTG',vtgz,factor_xy)
           call hbuf_put('VTI',vtiz,factor_xy)

           call hbuf_put('PCLW',pclwz,factor_xy*86400.*1000.)
           call hbuf_put('PVAPOR',pvaporz,factor_xy*86400.*1000.)
           call hbuf_put('PCLI',pcliz,factor_xy*86400.*1000.)
           call hbuf_put('PIMLT',pimltz,factor_xy*86400.*1000.)
           call hbuf_put('PIHOM',pihomz,factor_xy*86400.*1000.)

           call hbuf_put('PIDW',pidwz,factor_xy*86400.*1000.)
           call hbuf_put('PRAIN',prainz,factor_xy*86400.*1000.)
           call hbuf_put('PRAUT',prautz,factor_xy*86400.*1000.)
           call hbuf_put('PRACW',pracwz,factor_xy*86400.*1000.)
           call hbuf_put('PREVP',prevpz,factor_xy*86400.*1000.)

           call hbuf_put('PSNOW',psnowz,factor_xy*86400.*1000.)
           call hbuf_put('PSAUT',psautz,factor_xy*86400.*1000.)
           call hbuf_put('PSFW',psfwz,factor_xy*86400.*1000.)
           call hbuf_put('PSFI',psfiz,factor_xy*86400.*1000.)
           call hbuf_put('PRACI',praciz,factor_xy*86400.*1000.)

           call hbuf_put('PIACR',piacrz,factor_xy*86400.*1000.)
           call hbuf_put('PSACI',psaciz,factor_xy*86400.*1000.)
           call hbuf_put('PSACW',psacwz,factor_xy*86400.*1000.)
           call hbuf_put('PSDEP',psdepz,factor_xy*86400.*1000.)
           call hbuf_put('PSSUB',pssubz,factor_xy*86400.*1000.)

           call hbuf_put('PRACS',pracsz,factor_xy*86400.*1000.)
           call hbuf_put('PSACR',psacrz,factor_xy*86400.*1000.)
           call hbuf_put('PSMLT',psmltz,factor_xy*86400.*1000.)
           call hbuf_put('PSMLTEVP',psmltevpz,factor_xy*86400.*1000.)
           call hbuf_put('PLADJ',pladjz,factor_xy*86400.*1000.)

           call hbuf_put('PIADJ',piadjz,factor_xy*86400.*1000.)
           call hbuf_put('PGRAUPEL',pgraupelz,factor_xy*86400.*1000.)
           call hbuf_put('PGAUT',pgautz,factor_xy*86400.*1000.)
           call hbuf_put('PGFR',pgfrz,factor_xy*86400.*1000.)
           call hbuf_put('PGACW',pgacwz,factor_xy*86400.*1000.)

           call hbuf_put('PGACI',pgaciz,factor_xy*86400.*1000.)
           call hbuf_put('PGACR',pgacrz,factor_xy*86400.*1000.)
           call hbuf_put('PGACS',pgacsz,factor_xy*86400.*1000.)
           call hbuf_put('PGACIP',pgacipz,factor_xy*86400.*1000.)
           call hbuf_put('PGACRP',pgacrpz,factor_xy*86400.*1000.)

           call hbuf_put('PGACSP',pgacspz,factor_xy*86400.*1000.)
           call hbuf_put('PGWET',pgwetz,factor_xy*86400.*1000.)
           call hbuf_put('PDRY',pgwordz-pgwetz,factor_xy*86400.*1000.)
           call hbuf_put('PGSUB',pgsubz,factor_xy*86400.*1000.)
           call hbuf_put('PGDEP',pgdepz,factor_xy*86400.*1000.)

           call hbuf_put('PGMLT',pgmltz,factor_xy*86400.*1000.)
           call hbuf_put('PGMLTEVP',pgmltevpz,factor_xy*86400.*1000.)
        end if

        if (dorceintercomparison) then

        do k = 1,nzm
           sfz(k) = 0.
           cfz(k) = 0.
           ufz(k) = 0.
           dfz(k) = 0.
           umfz(k) = 0.
           dmfz(k) = 0.
           uhz(k) = 0.
           dhz(k) = 0.
           do j=1,ny
              do i=1,nx
                 if(.not.dokruegermicro) then ! bloss
                    omn = omegan(tabs(i,j,k))
                    omp = omegap(tabs(i,j,k))
                    omg = omegag(tabs(i,j,k))
                    lstarn=fac_cond+(1.-omn)*fac_fus
                    lstarp=fac_cond+(1.-omp)*fac_fus
                    qcc=qn(i,j,k)*omn
                    qii=qn(i,j,k)*(1.-omn)
                    qrr=qp(i,j,k)*omp
                    qss=qp(i,j,k)*(1.-omp)*(1.-omg)
                    qgg=qp(i,j,k)*(1.-omp)*omg
                 elseif(.not.dokruegermicro) then ! bloss
                    qcc=qc(i,j,k)
                    qii=qi(i,j,k)
                    qrr=qr(i,j,k)
                    qss=qs(i,j,k)
                    qgg=qg(i,j,k)
                    if (qcc+qii.gt.0.) then
                       omn = qcc/(qcc+qii)
                    else
                       omn = max(0.,min(1.,(tabs(i,j,k)-ttfrz)*a_bgkru))
                    end if
                    lstarn=fac_cond+(1.-omn)*fac_fus
                    lstarp=fac_cond+fac_fus*(qss+qgg)/(qss+qgg+qrr+eps)
                 end if
                 ! Compute moist static energy and vert. vel.
                 tmpmse = cp*(tabs(i,j,k)+gamaz(k))+lcond*(q(i,j,k)-qn(i,j,k))
                 tmpw = 0.5*(w(i,j,k+1)+w(i,j,k))
                 if ((qcc+qii).gt.5.e-6) cfz(k) = cfz(k) + 1.
                 if ((qcc+qii).gt.min(1.e-6,0.01*qsatw(tabs(i,j,k),pres(k)))) &
                      sfz(k) = sfz(k) + 1.
                 if (tmpw.gt.1.) then
                    ufz(k) = ufz(k) + 1.
                    umfz(k) = umfz(k) + rho(k)*tmpw
                    uhz(k) = uhz(k) + tmpmse
                 elseif (tmpw.lt.-1.) then
                    dfz(k) = dfz(k) + 1.
                    dmfz(k) = dmfz(k) + rho(k)*tmpw
                    dhz(k) = dhz(k) + tmpmse
                 end if
              end do
           end do
        end do

        ! Create a vector which holds time series of scalar values for rce
        ! intercomparison
        rcei = 0.
        do j = 1,ny
           do i = 1,nx
              rcei(1) = rcei(1) + prec_inst(i,j)*86400.*1000./rhor
              rcei(2) = rcei(2) + fluxbq(i,j)*rhow(1)*86400.*1000.&
                   &/rhor
              rcei(3) = rcei(3) + swntxy(i,j)
              rcei(4) = rcei(4) + lwntxy(i,j)
              rcei(5) = rcei(5) + swnsxy(i,j)
              rcei(6) = rcei(6) + lwnsxy(i,j)
           end do
        end do

        ! Precipitable water and surface-100 mbar mass-weighted dry
        ! static energy
        do k = 1,nzm
           coef=max(0.,min(presi(k)-100.,presi(k)-presi(k+1)))/(presi(1)-100.)
           do j = 1,ny
              do i = 1,nx
                 rcei(7) = rcei(7) + dz*adz(k)*rho(k)*q(i,j,k)
                 rcei(8) = rcei(8) + coef*cp*(tabs(i,j,k)+gamaz(k))
              end do
           end do
        end do
        
        rcei(9) = rcei(9) + float(nx*ny)*presi(1) ! surface pressure

        ! lowest-level temperature, moisture, wind speed
        do j = 1,ny
           do i = 1,nx
              rcei(10) = rcei(10) + tabs(i,j,1)
              rcei(11) = rcei(11) + q(i,j,1) - qn(i,j,1)
              rcei(12) = rcei(12) + sqrt(u(i,j,1)**2 + v(i,j,1)**2)
           end do
        end do
	
        ! satellite-detectable cloud fraction
        do j = 1,ny
           do i = 1,nx
              coef = 0.
              coef1 = 0.
              tmpw = 0.
              do k = 1,nzm
                 if(.not.dokruegermicro) then ! bloss
                    omn = omegan(tabs(i,j,k))
                    omp = omegap(tabs(i,j,k))
                    omg = omegag(tabs(i,j,k))
                    qcc=qn(i,j,k)*omn
                    qii=qn(i,j,k)*(1.-omn)
                    qss=qp(i,j,k)*(1.-omp)*(1.-omg)
                    qgg=qp(i,j,k)*(1.-omp)*omg
                 elseif(.not.dokruegermicro) then ! bloss
                    qcc=qc(i,j,k)
                    qii=qi(i,j,k)
                    qrr=qr(i,j,k)
                    qss=qs(i,j,k)
                    qgg=qg(i,j,k)
                    if (qcc+qii.gt.0.) then
                       omn = qcc/(qcc+qii)
                    else
                       omn = max(0.,min(1.,(tabs(i,j,k)-ttfrz)*a_bgkru))
                    end if
                    lstarn=fac_cond+(1.-omn)*fac_fus
                    lstarp=fac_cond+fac_fus*(qss+qgg)/(qss+qgg+qrr+eps)
                 end if
                 coef1 = coef1 + dz*adz(k)*rho(k)*(qcc+qii)
                 if (qcc+qii.gt.5.e-6) coef = 1.
                 tmpw = tmpw + dz*adz(k)*rho(k)*(0.15*qcc+ (0.005 &
                      +1./(30.-20.*max(0.,min(1.,(pres(k)/pres0-.4)/.4))))*qii)
              end do
              if (coef1.gt.6.5e-4) rcei(13) = rcei(13) + 1.
              if (tmpw.gt.1.e-4) rcei(15) = rcei(15) + 1.
              rcei(14) = rcei(14) + coef
           end do
        end do

	call hbuf_put('SF',sfz,factor_xy) 	
	call hbuf_put('CF',cfz,factor_xy) 	
	call hbuf_put('UF',ufz,factor_xy) 	
	call hbuf_put('DF',dfz,factor_xy) 	
	call hbuf_put('UMF',umfz,factor_xy) 	
	call hbuf_put('DMF',dmfz,factor_xy) 	
	call hbuf_put('UH',uhz,factor_xy) 	
	call hbuf_put('DH',dhz,factor_xy) 	
	call hbuf_put('HFLX',twle+qwle+qpwle-twsb-qwsb-qpwsb,factor_xy)
	call hbuf_put('HFLXS',twsb+qwsb+qpwsb,factor_xy) 	
	call hbuf_put('RCEI',rcei,factor_xy) 	

        end if

        if(docasefive) then

           thflux = 0.
           thflxs = 0.
           qvflux = 0.
           qvflxs = 0.

           thflxs(1) = SUM(fluxbt)
           qvflxs(1) = SUM(fluxbq)

           do k = 2,nzm
              do j = 1,ny
                 do i = 1,nx
                    thp = tabs(i,j,k)*prespot(k)
                    thm = tabs(i,j,k-1)*prespot(k-1)
                    qvp = q(i,j,k) - qn(i,j,k)
                    qvm = q(i,j,k-1) - qn(i,j,k-1)
                    tkhm = 0.5*(tkh_z(i,j,k)+tkh_z(i,j,k-1))
                    if(LES) tkhm = sqrt(tkh_z(i,j,k)*tkh_z(i,j,k-1))

                    thflux(k) = thflux(k) + w(i,j,k)*0.5*(thp+thm)
                    thflxs(k) = thflxs(k) - tkhm*(thp-thm)/dz/adzw(k)
                    qvflux(k) = qvflux(k) + w(i,j,k)*0.5*(qvp+qvm)
                    qvflxs(k) = qvflxs(k) - tkhm*(qvp-qvm)/dz/adzw(k)
                 end do
              end do
           end do
           
           do k = 1,nzm
              dthadv(k) = -(rhow(k+1)*thflux(k+1)-rhow(k)*thflux(k)) &
                   /rho(k)/dz/adz(k)
              dthdif(k) = -(rhow(k+1)*thflxs(k+1)-rhow(k)*thflxs(k)) &
                   /rho(k)/dz/adz(k)
              dqvadv(k) = -(rhow(k+1)*qvflux(k+1)-rhow(k)*qvflux(k)) &
                   /rho(k)/dz/adz(k)
              dqvdif(k) = -(rhow(k+1)*qvflxs(k+1)-rhow(k)*qvflxs(k)) &
                   /rho(k)/dz/adz(k)
           end do

           call hbuf_put('THFLUX',thflux,factor_xy)
           call hbuf_put('THFLXS',thflxs,factor_xy)
           call hbuf_put('QVFLUX',qvflux,factor_xy*1000.)
           call hbuf_put('QVFLXS',qvflxs,factor_xy*1000.)
           call hbuf_put('THADV',dthadv,factor_xy*86400.)
           call hbuf_put('THDIFF',dthdif,factor_xy*86400.)
           call hbuf_put('QVADV',dqvadv,factor_xy*86400000.)
           call hbuf_put('QVDIFF',dqvdif,factor_xy*86400000.)
           call hbuf_put('QVCOND',qvcond,factor_xy*86400000./dtn)
           call hbuf_put('TLCOND',tlcond,factor_xy*86400./dtn)
           call hbuf_put('TLEVP',tlevp,factor_xy*86400.)




           cloudxy = 0.
           cttemp = sstxy(1:nx,1:ny)
           topind = -1

           clfracz = 0.
           hmfracz = 0.
           cufracz = 0. ! cloudy updrafts area frac
           cumflxz = 0. ! cloudy updrafts mass flux
           wufracz = 0. ! updraft core area frac
           wumflxz = 0. ! updraft core mass flux
           bcufracz = 0. ! buoyant cloudy updraft area frac
           bcumflxz = 0. ! buoyant cloudy updraft mass flux
           do k = 1,nzm
              do j = 1,ny
                 do i = 1,nx
                    tmpw = 0.5*(w(i,j,k)+w(i,j,k+1))
                    if(qn(i,j,k).gt.1.e-5) then 
                       clfracz(k) = clfracz(k) + 1. ! cloudy
                       cloudxy(i,j) = 1.
                       cttemp(i,j) = tabs(i,j,k)
                       topind(i,j) = k

                       if(tmpw.gt.0.) then 
                          cufracz(k) = cufracz(k) + 1. ! cloudy updraft
                          cumflxz(k) = cumflxz(k) + rho(k)*tmpw

                          if(tvirt(i,j,k).gt.tvz(k)) then 
                             bcufracz(k) = bcufracz(k) + 1. ! buoy cldy updraft
                             bcumflxz(k) = bcumflxz(k) + rho(k)*tmpw
                          end if
                       end if
                    end if

                    if(tmpw.gt.5.) then
                       wufracz(k) = wufracz(k) + 1. ! updraft core
                       wumflxz(k) = wumflxz(k) + rho(k)*tmpw
                    end if

                    if(qn(i,j,k)+qp(i,j,k).gt.1.e-5) then
                       hmfracz(k) = hmfracz(k) + 1.
                    end if
                 end do
              end do
           end do

           call hbuf_put('CLFRAC',clfracz,factor_xy) 	
           call hbuf_put('HMFRAC',hmfracz,factor_xy) 	
           call hbuf_put('CUFRAC',cufracz,factor_xy) 	
           call hbuf_put('CUMFLX',cumflxz,factor_xy) 	
           call hbuf_put('WUFRAC',wufracz,factor_xy) 	
           call hbuf_put('WUMFLX',wumflxz,factor_xy) 	
           call hbuf_put('BCUFRAC',bcufracz,factor_xy) 	
           call hbuf_put('BCUMFLX',bcumflxz,factor_xy) 	

           ! TIMESERIES
           c5t = 0.
           c5sst =0.
           c5dse =0.
           c5qv = 0.
           c5mse =0.
           c5u = 0.
           c5v = 0.
           c5shf = 0.
           c5lhf = 0.
           c5tau = 0.
           c5tav = 0.

           c5swds =0.
           c5lwds =0.
           c5swus =0.
           c5lwus =0.
           c5swdt =0.
           c5swut =0.
           c5lwut =0.

           c5cld =0.
           c5pcp =0.

           c5wcf = 0.
           c5wpcp = 0.
           c5ifc = 0.
           c5ifpcp = 0.
           c5hcf = 0.
           c5hpcp = 0.
           c5lcf = 0.
           c5lpcp = 0.
           c5mcf = 0.
           c5mpcp = 0.

           c5pw = 0.
           c5lwp = 0.
           c5iwp = 0.
           c5rp = 0.
           c5sp = 0.
           c5gp = 0.

           nts = mod(nstep,nstat)/nstatis
           if (nts.eq.0) nts = nstatfrq
           c5t(nts) = (time/3600.)*float(nx*ny)
           c5sst(nts) = SUM(sstxy(1:nx,1:ny))
           c5dse(nts) = SUM(t(1:nx,1:ny,1))*cp/1000. ! kJ/kg
           c5qv(nts) = SUM(q(1:nx,1:ny,1))*1000. ! g/kg
           c5mse(nts) = SUM(cp*t(1:nx,1:ny,1)+lcond*q(1:nx,1:ny,1))/1000.
           c5u(nts) = SUM(u(1:nx,1:ny,1)) 
           c5v(nts) = SUM(v(1:nx,1:ny,1)) 
           c5shf(nts) = SUM(rhow(1)*cp*fluxbt)
           c5lhf(nts) = SUM(rhow(1)*lcond*fluxbq)
           c5tau(nts) = SUM(fluxbu)
           c5tav(nts) = SUM(fluxbv)

           c5swds(nts) = SUM(swdsxy)
           c5lwds(nts) = SUM(lwdsxy)
           c5swus(nts) = SUM(swdsxy-swnsxy)
           c5lwus(nts) = SUM(lwnsxy+lwdsxy)
           c5swdt(nts) = SUM(solinxy)
           c5swut(nts) = SUM(solinxy-swntxy)
           c5lwut(nts) = SUM(lwutxy)
           
           c5cld(nts) = SUM(cloudxy)
           c5pcp(nts) = SUM(precsfc)

           if(dokruegermicro) then
              do k = 1,nzm
                 fac = dz*adz(k)*rho(k)
                 c5pw(nts) = c5pw(nts) + fac*SUM(q(1:nx,1:ny,k))
                 c5lwp(nts)= c5lwp(nts)+ fac*SUM(qc(1:nx,1:ny,k))
                 c5iwp(nts)= c5iwp(nts)+ fac*SUM(qi(1:nx,1:ny,k))
                 c5rp(nts) = c5rp(nts) + fac*SUM(qr(1:nx,1:ny,k))
                 c5sp(nts) = c5sp(nts) + fac*SUM(qs(1:nx,1:ny,k))
                 c5gp(nts) = c5gp(nts) + fac*SUM(qg(1:nx,1:ny,k))
              end do
           else
              do k = 1,nzm
                 fac = dz*adz(k)*rho(k)
                 c5pw(nts) = c5pw(nts) + fac*SUM(q(1:nx,1:ny,k))
                 do j = 1,ny
                    do i = 1,nx
                       omn = omegan(tabs(i,j,k))
                       omp = omegap(tabs(i,j,k))
                       omg = omegag(tabs(i,j,k))
                       c5lwp(nts)= c5lwp(nts)+ fac*qn(i,j,k)*omn
                       c5iwp(nts)= c5iwp(nts)+ fac*qn(i,j,k)*(1.-omn)
                       c5rp(nts) = c5rp(nts) + fac*qp(i,j,k)*omp
                       c5sp(nts) = c5sp(nts) + fac*qp(i,j,k)*(1.-omp)*(1.-omg)
                       c5gp(nts) = c5gp(nts) + fac*qp(i,j,k)*(1.-omp)*omg
                    end do
                 end do
              end do
           end if
           
           do j = 1,ny
              do i = 1,nx
                 if (cloudxy(i,j).eq.1.) then
                    if(cttemp(i,j).gt.273.15) then ! WARM CLOUD
                       c5wcf(nts) = c5wcf(nts) + 1.
                       c5wpcp(nts) = c5wpcp(nts) + precsfc(i,j)
                    end if
                    if(c5iwp(nts)+c5sp(nts) &
                         +c5gp(nts).lt.0.02) then ! ICE-FREE CLOUD
                       c5ifc(nts) = c5ifc(nts) + 1.
                       c5ifpcp(nts) = c5ifpcp(nts) + precsfc(i,j)
                    end if
                    if(z(topind(i,j)).gt.9000.) then ! DEEP CLOUD
                       c5hcf(nts) = c5hcf(nts) + 1.
                       c5hpcp(nts) = c5hpcp(nts) + precsfc(i,j)
                    elseif(z(topind(i,j)).lt.4000.) then ! SHALLOW CLOUD
                       c5lcf(nts) = c5lcf(nts) + 1.
                       c5lpcp(nts) = c5lpcp(nts) + precsfc(i,j)
                    else ! MIDDLE CLOUD
                       c5mcf(nts) = c5mcf(nts) + 1.
                       c5mpcp(nts) = c5mpcp(nts) + precsfc(i,j)
                    end if
                 end if
              end do
           end do

           fac = factor_xy*float(nstatfrq)
           call hbuf_put('TSTIME',c5t,fac)
           call hbuf_put('TSSST',c5sst,fac)
           call hbuf_put('TSDSE',c5dse,fac)
           call hbuf_put('TSQV',c5qv,fac)
           call hbuf_put('TSMSE',c5mse,fac)
           call hbuf_put('TSU',c5u,fac)
           call hbuf_put('TSV',c5v,fac)
           call hbuf_put('TSSHF',c5shf,fac)
           call hbuf_put('TSLHF',c5lhf,fac)
           call hbuf_put('TSTAUX',c5tau,fac)
           call hbuf_put('TSTAUY',c5tav,fac)

           call hbuf_put('TSSWDS',c5swds,fac)
           call hbuf_put('TSLWDS',c5lwds,fac)
           call hbuf_put('TSSWUS',c5swus,fac)
           call hbuf_put('TSLWUS',c5lwus,fac)
           call hbuf_put('TSSOLIN',c5swdt,fac)
           call hbuf_put('TSSWUT',c5swut,fac)
           call hbuf_put('TSLWUT',c5lwut,fac)

           call hbuf_put('TSCLD',c5cld,fac)
           call hbuf_put('TSPCP',c5pcp,fac*24.)

           call hbuf_put('TSWCLF',c5wcf,fac)
           call hbuf_put('TSWPCP',c5wpcp,fac*24.)
           call hbuf_put('TSIFCLF',c5ifc,fac)
           call hbuf_put('TSIFPCP',c5ifpcp,fac*24.)
           call hbuf_put('TSHCLF',c5hcf,fac)
           call hbuf_put('TSHPCP',c5hpcp,fac*24.)
           call hbuf_put('TSLCLF',c5lcf,fac)
           call hbuf_put('TSLPCP',c5lpcp,fac*24.)
           call hbuf_put('TSMCLF',c5mcf,fac)
           call hbuf_put('TSMPCP',c5mpcp,fac*24.)

           call hbuf_put('TSPW',c5pw,fac)
           call hbuf_put('TSLWP',c5lwp,fac)
           call hbuf_put('TSIWP',c5iwp,fac)
           call hbuf_put('TSRWP',c5rp,fac)
           call hbuf_put('TSSWP',c5sp,fac)
           call hbuf_put('TSGWP',c5gp,fac)


        end if


! END UW ADDITIONS
!===================================================

end
	
