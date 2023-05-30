
! peters.  utility subroutines for calling surface scheme
!
!  takes location i,j to call surface scheme at and gridflag designating
!  the grid to call at

subroutine surface_scheme(i,j,gridflag,zz00)

  use grid, only : nxp1, nyp1
  use domain, only : YES3D
  use simple_land, only : landmask_frac
  implicit none

  integer i, j       ! input: current location to compute at
  integer gridflag   ! input: 1=u grid, 2=v grid, 3=scalar grid
  real zz00(0:nxp1,1-YES3D:nyp1) ! input: surface roughness

  ! local
  real tmplandmask, fluxu, fluxv, fluxt, fluxq, tfluxu, tfluxv, tfluxt, tfluxq
  integer ib,ic,jb,jc


  !------------------------------------------
  ib=i-1
  ic=i+1
  jc=j+YES3D
  jb=j-YES3D

  ! make the landmask for this grid
  select case(gridflag)
    case(1)
      tmplandmask = 0.5*(landmask_frac(i,j)+landmask_frac(ib,j))
    case(2)
      tmplandmask = 0.5*(landmask_frac(i,j)+landmask_frac(i,jb))
    case(3)
      tmplandmask = landmask_frac(i,j)
  end select
  !------------------------------------------

  !------------------------------------------
  !  call appropriate land/ocean scheme
  if(tmplandmask.gt.0.9999) then       ! all land

    call call_landflx(gridflag,i,j,ib,ic,jb,jc,fluxu,fluxv,fluxt,fluxq,zz00)
    call assign_fluxes(i,j,gridflag,fluxu,fluxv,fluxt,fluxq)

  elseif(tmplandmask.lt.1e-4) then     ! all ocean

    call call_oceflx(gridflag,i,j,ib,ic,jb,jc,fluxu,fluxv,fluxt,fluxq)
    call assign_fluxes(i,j,gridflag,fluxu,fluxv,fluxt,fluxq)

  else

    call call_landflx(gridflag,i,j,ib,ic,jb,jc,tfluxu,tfluxv,tfluxt,tfluxq,zz00)
    call call_oceflx(gridflag,i,j,ib,ic,jb,jc,fluxu,fluxv,fluxt,fluxq)
 
    ! take weighted average of land/ocean fluxes
    fluxu = (1.-tmplandmask)*fluxu+tmplandmask*tfluxu
    fluxv = (1.-tmplandmask)*fluxv+tmplandmask*tfluxv
    fluxt = (1.-tmplandmask)*fluxt+tmplandmask*tfluxt
    fluxq = (1.-tmplandmask)*fluxq+tmplandmask*tfluxq

    call assign_fluxes(i,j,gridflag,fluxu,fluxv,fluxt,fluxq)

  end if

end subroutine surface_scheme

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine call_landflx(gridflag,i,j,ib,ic,jb,jc,fluxu,fluxv,fluxt,fluxq,zz00)

  use domain, only : YES3D
  use grid, only : nxp1, nyp1
  use vars, only : u, v, sstxy, t, q, gamaz
  use simple_land, only : xysoil_wet, xysurf_rough
  use grid, only : z, pres, presi
  use params
  implicit none

  integer gridflag,i,j,ib,ic,jb,jc              ! input
  real fluxu,fluxv,fluxt,fluxq         ! output
  real zz00(0:nxp1,1-YES3D:nyp1)       ! input

  ! local
  real t_s, q_s, tb, tmpu, tmpv, tmpq, tmpz0, xlmo
  real coef, coef1
  real qsatw
  external qsatw


  ! make the input to landflx
  coef = (1000./pres0)**(rgas/cp)
  coef1 = (1000./pres(1))**(rgas/cp)

  select case(gridflag)
    case(1)
      tmpu = u(i,j,1) + ug
      tmpv = 0.25*(v(i,j,1) +    v(i,jc,1) + &
                   v(ib,jc,1) + v(ib,j,1))+vg
      t_s = 0.5*(sstxy(i,j)*coef + sstxy(ib,j)*coef)
      q_s = 0.5*(xysoil_wet(i,j)*qsatw(sstxy(i,j),pres(1)) + &
                 xysoil_wet(ib,j)*qsatw(sstxy(ib,j),pres(1)))
      tb = 0.5*( t(i,j,1)-gamaz(1) + &
                 t(ib,j,1)-gamaz(1) )
      tmpq = 0.5*(q(i,j,1)+q(ib,j,1))
      tmpz0 = 0.5*(zz00(i,j)+zz00(ib,j))

    case(2)
      tmpu = 0.25*(u(i,j,1) + u(ic,j,1) &
               +  u(i,jb,1)+ u(ic,jb,1))+ug
      tmpv = v(i,j,1)+vg
      t_s = 0.5*(sstxy(i,j)*coef + sstxy(i,jb)*coef)
      q_s = 0.5*(xysoil_wet(i,j)*qsatw(sstxy(i,j),pres(1)) + &
                 xysoil_wet(i,jb)*qsatw(sstxy(i,jb),pres(1)))
      tb = 0.5*( t(i,j,1)-gamaz(1) + &
                 t(i,jb,1)-gamaz(1) )
      tmpq = 0.5*(q(i,j,1)+q(i,jb,1))
      tmpz0 = 0.5*(zz00(i,j)+zz00(i,jb))

    case(3)
      tmpu = 0.5*(u(ic,j,1)+u(i,j,1))+ug
      tmpv = 0.5*(v(i,jc,1)+v(i,j,1))+vg
      t_s = sstxy(i,j)*coef
      q_s = xysoil_wet(i,j)*qsatw(sstxy(i,j),pres(1))
      tb = t(i,j,1)-gamaz(1)
      tmpq = q(i,j,1)
      tmpz0 = zz00(i,j)
  end select

  call landflx(pres(1),tb*coef1, t_s, tmpq, q_s, tmpu, tmpv, z(1), tmpz0,&
               fluxt, fluxq, fluxu, fluxv, xlmo)

end subroutine call_landflx

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine call_oceflx(gridflag,i,j,ib,ic,jb,jc,fluxu,fluxv,fluxt,fluxq)

  use vars, only : u, v, sstxy, t, q, gamaz
  use grid, only : z, pres, sfcum
  use params
  implicit none

  integer gridflag,i,j,ib,ic,jb,jc              ! input
  real fluxu,fluxv,fluxt,fluxq         ! output

  ! local
  real tb, tmpu, tmpv, tmpq, q_s, t_s

  select case(gridflag)
    case(1)
      tmpu = u(i,j,1) + ug
      tmpv = 0.25*(v(i,j,1)   +    v(i,jc,1) + &
                   v(ib,jc,1) +    v(ib,j,1))+vg
      tb = 0.5*(t(i,j,1)-gamaz(1) + t(ib,j,1)-gamaz(1) )
      tmpq = 0.5*(q(i,j,1)+q(ib,j,1))
      t_s = 0.5*(sstxy(i,j)+sstxy(ib,j))

    case(2)
      tmpu = 0.25*(u(i,j,1) + u(ic,j,1) &
               +  u(i,jb,1)+ u(ic,jb,1))+ug
      tmpv = v(i,j,1)+vg
      tb = 0.5*(t(i,j,1)-gamaz(1) + t(i,jb,1)-gamaz(1) )
      tmpq = 0.5*(q(i,j,1)+q(i,jb,1))
      t_s = 0.5*(sstxy(i,j)+sstxy(i,jb))

    case(3)
      tmpu = 0.5*(u(ic,j,1)+u(i,j,1))+ug
      tmpv = 0.5*(v(i,jc,1)+v(i,j,1))+vg
      tb = t(i,j,1)-gamaz(1)
      tmpq = q(i,j,1)
      t_s = sstxy(i,j)
  end select

  call oceflx(pres(1), tmpu, tmpv, tb, tmpq, t(i,j,1), z(1), t_s, &
              fluxt, fluxq, fluxu, fluxv, q_s, sfcum)

end subroutine call_oceflx

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine assign_fluxes(i,j,gridflag,fluxu,fluxv,fluxt,fluxq)

  use vars, only : fluxbu, fluxbv, fluxbq, fluxbt

  implicit none

  integer i,j,gridflag
  real fluxu,fluxv,fluxt,fluxq

  select case(gridflag)
    case(1)
      fluxbu(i,j) = fluxu
    case(2)
      fluxbv(i,j) = fluxv
    case(3)
      fluxbq(i,j) = fluxq
      fluxbt(i,j) = fluxt
  end select

end subroutine assign_fluxes

