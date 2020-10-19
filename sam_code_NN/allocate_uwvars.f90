subroutine allocate_uwvars()
  use vars
  implicit none

  integer ierr

  ! ALLOCATE MICROPHYSICAL VARIABLES FOR KRUEGER (LIN/LORD) SCHEME
  if(.not.ALLOCATED(qv).and.dokruegermicro) then
     ALLOCATE(qv(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm),STAT=ierr)
     if (ierr.ne.0) call failed_allocation('qv')       
  end if

  if(.not.ALLOCATED(qc).and.dokruegermicro) then
     ALLOCATE(qc(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm),STAT=ierr)
     if (ierr.ne.0) call failed_allocation('qc')       
  end if

  if(.not.ALLOCATED(qi).and.dokruegermicro) then
     ALLOCATE(qi(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm),STAT=ierr)
     if (ierr.ne.0) call failed_allocation('qi')       
  end if

  if(.not.ALLOCATED(qr).and.dokruegermicro) then
     ALLOCATE(qr(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm),STAT=ierr)
     if (ierr.ne.0) call failed_allocation('qr')       
  end if

  if(.not.ALLOCATED(qs).and.dokruegermicro) then
     ALLOCATE(qs(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm),STAT=ierr)
     if (ierr.ne.0) call failed_allocation('qs')       
  end if

  if(.not.ALLOCATED(qg).and.dokruegermicro) then
     ALLOCATE(qg(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm),STAT=ierr)
     if (ierr.ne.0) call failed_allocation('qg')       
  end if

  if(.not.ALLOCATED(qcwle).and.dokruegermicro) then
     ALLOCATE(qcwle(nz),qcwsb(nz),qcadv(nz),qcdiff(nz), & 
            qcsrc(nzm),qcfall(nz),qcevp(nzm),&
          qiwle(nz),qiwsb(nz),qiadv(nz),qidiff(nz), &
            qisrc(nzm),qievp(nzm),&
          qrwle(nz),qrwsb(nz),qradv(nz),qrdiff(nz), &
            qrsrc(nzm),qrfall(nz),qrevp(nzm),&
          qgwle(nz),qgwsb(nz),qgadv(nz),qgdiff(nz), &
            qgsrc(nzm),qgfall(nz),qgevp(nzm),&
          qswle(nz),qswsb(nz),qsadv(nz),qsdiff(nz), &
            qssrc(nzm),qsfall(nz),qsevp(nzm),&
          rainfrac(nzm), snowfrac(nzm), graufrac(nzm),  &
            clifrac(nzm), clwfrac(nzm), &
          vtrz(nzm), vtsz(nzm), vtgz(nzm), vtiz(nzm), &
          rainflux(nz), snowflux(nz), grauflux(nz), &
          pclwz(nzm), pvaporz(nzm), pcliz(nzm), pimltz(nzm), pihomz(nzm), &
          pidwz(nzm), prainz(nzm), prautz(nzm), pracwz(nzm), prevpz(nzm), &
          psnowz(nzm), psautz(nzm), psfwz(nzm), psfiz(nzm), praciz(nzm), &
          piacrz(nzm), psaciz(nzm), psacwz(nzm), psdepz(nzm), pssubz(nzm), &
          pracsz(nzm), psacrz(nzm), psmltz(nzm), psmltevpz(nzm), pladjz(nzm), &
          piadjz(nzm), pgraupelz(nzm), pgautz(nzm),pgfrz(nzm), pgacwz(nzm), &
          pgaciz(nzm), pgacrz(nzm), pgacsz(nzm), pgacipz(nzm), pgacrpz(nzm), &
          pgacspz(nzm), pgwetz(nzm), pdryz(nzm), pgsubz(nzm), pgdepz(nzm), &
          pgmltz(nzm), pgmltevpz(nzm), pgwordz(nzm),STAT=ierr)
     if (ierr.ne.0) call failed_allocation('lin/lord microphysical stats')
     qcadv=0.; qiadv=0.; qradv=0.; qsadv=0.; qgadv=0. ! zero out q*adv
     qcdiff=0.; qidiff=0.; qrdiff=0.; qsdiff=0.; qgdiff=0. ! zero out q*diff
     qcwle=0.; qiwle=0.; qrwle=0.; qswle=0.; qgwle=0. ! zero out q*wle
     qcwsb=0.; qiwsb=0.; qrwsb=0.; qswsb=0.; qgwsb=0. ! zero out q*adv
     qrfall=0.; qsfall=0.; qgfall=0. ! zero out q*fall
     qcsrc=0.; qisrc=0.; qrsrc=0.; qssrc=0.; qgsrc=0. ! zero out q*src
     qcevp=0.; qievp=0.; qrevp=0.; qsevp=0.; qgevp=0. ! zero out q*evp
  end if


!!$  ! ALLOCATE MICROPHYSICAL VARIABLES FOR ISOTOPES
!!$  if(.not.ALLOCATED(iso_qv).and.doisotopes) then
!!$     ALLOCATE(iso_qv(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm),STAT=ierr)
!!$     if (ierr.ne.0) call failed_allocation('iso_qv')       
!!$  end if
!!$
!!$  if(.not.ALLOCATED(iso_qc).and.doisotopes) then
!!$     ALLOCATE(iso_qc(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm),STAT=ierr)
!!$     if (ierr.ne.0) call failed_allocation('iso_qc')       
!!$  end if
!!$
!!$  if(.not.ALLOCATED(iso_qi).and.doisotopes) then
!!$     ALLOCATE(iso_qi(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm),STAT=ierr)
!!$     if (ierr.ne.0) call failed_allocation('iso_qi')       
!!$  end if
!!$
!!$  if(.not.ALLOCATED(iso_qr).and.doisotopes) then
!!$     ALLOCATE(iso_qr(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm),STAT=ierr)
!!$     if (ierr.ne.0) call failed_allocation('iso_qr')       
!!$  end if
!!$
!!$  if(.not.ALLOCATED(iso_qs).and.doisotopes) then
!!$     ALLOCATE(iso_qs(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm),STAT=ierr)
!!$     if (ierr.ne.0) call failed_allocation('iso_qs')       
!!$  end if
!!$
!!$  if(.not.ALLOCATED(iso_qg).and.doisotopes) then
!!$     ALLOCATE(iso_qg(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm),STAT=ierr)
!!$     if (ierr.ne.0) call failed_allocation('iso_qg')       
!!$  end if

  !tracers added by kzm

  ! allocate variables related to tracer of x grid
  if(.not.ALLOCATED(trx).and.(dotrx.or.dooldrestart)) then
     ALLOCATE(trx(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm), &
          fluxbtrx(nx,ny),trxadv(nz),trxdiff(nz),trxwsb(nz), &
          trxwle(nz),trxwleadv(nz),trxwlediff(nz),trx2leadv(nz), &
          trx2lediff(nz),trx2legrad(nz),trx2lediss(nz),STAT=ierr) 
     if (ierr.ne.0) call failed_allocation('x grid tracer')
     fluxbtrx = 0.
  end if

  ! allocate variables related to tracer of y grid
  if(.not.ALLOCATED(try).and.(dotry.or.dooldrestart)) then
     ALLOCATE(try(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm), &
          fluxbtry(nx,ny),tryadv(nz),trydiff(nz),trywsb(nz), &
          trywle(nz),trywleadv(nz),trywlediff(nz),try2leadv(nz), &
          try2lediff(nz),try2legrad(nz),try2lediss(nz),STAT=ierr) 
     if (ierr.ne.0) call failed_allocation('y grid tracer')
     fluxbtry = 0.
  end if

  ! allocate variables related to tracer of z grid
  if(.not.ALLOCATED(trz).and.(dotrz.or.dooldrestart)) then
     ALLOCATE(trz(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm), &
          fluxbtrz(nx,ny),trzadv(nz),trzdiff(nz),trzwsb(nz), &
          trzwle(nz),trzwleadv(nz),trzwlediff(nz),trz2leadv(nz), &
          trz2lediff(nz),trz2legrad(nz),trz2lediss(nz),STAT=ierr) 
     if (ierr.ne.0) call failed_allocation('z grid tracer')
     fluxbtrz = 0.
  end if

  ! allocate variables related to tracer of zz grid
  if(.not.ALLOCATED(trzz).and.(dotrzz.or.dooldrestart)) then
     ALLOCATE(trzz(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm), &
          fluxbtrzz(nx,ny),trzzadv(nz),trzzdiff(nz),trzzwsb(nz), &
          trzzwle(nz),trzzwleadv(nz),trzzwlediff(nz),trzz2leadv(nz), &
          trzz2lediff(nz),trzz2legrad(nz),trzz2lediss(nz),STAT=ierr) 
     if (ierr.ne.0) call failed_allocation('zz grid tracer')
     fluxbtrzz = 0.
  end if

  ! allocate variables related to tracer of o3 grid
  if(.not.ALLOCATED(tro3).and.(dotro3.or.dooldrestart)) then
     ALLOCATE(tro3(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm), &
          fluxbtro3(nx,ny),tro3adv(nz),tro3diff(nz),tro3wsb(nz), &
          tro3wle(nz),tro3wleadv(nz),tro3wlediff(nz),tro32leadv(nz), &
          tro32lediff(nz),tro32legrad(nz),tro32lediss(nz),STAT=ierr) 
     if (ierr.ne.0) call failed_allocation('o3 grid tracer')
     fluxbtro3 = 0.
  end if

end subroutine allocate_uwvars

subroutine deallocate_uwvars()
  use vars
  implicit none

  integer ierr

  if(ALLOCATED(qv)) then
     DEALLOCATE(qv,STAT=ierr)
     if (ierr.ne.0) call failed_allocation('qv')       
  end if

  if(ALLOCATED(qc)) then
     DEALLOCATE(qc,STAT=ierr)
     if (ierr.ne.0) call failed_allocation('qc')   
  end if

  if(ALLOCATED(qi)) then
     DEALLOCATE(qi,STAT=ierr)
     if (ierr.ne.0) call failed_allocation('qi')     
  end if

  if(ALLOCATED(qr)) then
     DEALLOCATE(qr,STAT=ierr)
     if (ierr.ne.0) call failed_allocation('qr')      
  end if

  if(ALLOCATED(qs)) then
     DEALLOCATE(qs,STAT=ierr)
     if (ierr.ne.0) call failed_allocation('qs')      
  end if

  if(ALLOCATED(qg)) then
     DEALLOCATE(qg,STAT=ierr)
     if (ierr.ne.0) call failed_allocation('qg')       
  end if

  if(ALLOCATED(qcwle)) then
     DEALLOCATE( &
          qcwle,qcwsb,qcadv,qcdiff, qcsrc,qcfall,qcevp,&
          qiwle,qiwsb,qiadv,qidiff, qisrc,qievp,&
          qrwle,qrwsb,qradv,qrdiff, qrsrc,qrfall,qrevp,&
          qgwle,qgwsb,qgadv,qgdiff, qgsrc,qgfall,qgevp,&
          qswle,qswsb,qsadv,qsdiff, qssrc,qsfall,qsevp,&
          rainfrac, snowfrac, graufrac, clifrac, clwfrac, &
          vtrz, vtsz, vtgz, vtiz, rainflux, snowflux, grauflux, &
          pclwz, pvaporz, pcliz, pimltz, pihomz, &
          pidwz, prainz, prautz, pracwz, prevpz, &
          psnowz, psautz, psfwz, psfiz, praciz, &
          piacrz, psaciz, psacwz, psdepz, pssubz, &
          pracsz, psacrz, psmltz, psmltevpz, pladjz, &
          piadjz, pgraupelz, pgautz,pgfrz, pgacwz, &
          pgaciz, pgacrz, pgacsz, pgacipz, pgacrpz, &
          pgacspz, pgwetz, pdryz, pgsubz, pgdepz, &
          pgmltz, pgmltevpz, pgwordz,STAT=ierr)
     if (ierr.ne.0) call failed_allocation('lin/lord microphysical stats')
  end if

  !ISOTOPES
  if(ALLOCATED(iso_qv)) then
     DEALLOCATE(iso_qv,STAT=ierr)
     if (ierr.ne.0) call failed_allocation('iso_qv')  
  end if

  if(ALLOCATED(iso_qc)) then
     DEALLOCATE(iso_qc,STAT=ierr)
     if (ierr.ne.0) call failed_allocation('iso_qc')  
  end if

  if(ALLOCATED(iso_qi)) then
     DEALLOCATE(iso_qi,STAT=ierr)
     if (ierr.ne.0) call failed_allocation('iso_qi')  
  end if

  if(ALLOCATED(iso_qr)) then
     DEALLOCATE(iso_qr,STAT=ierr)
     if (ierr.ne.0) call failed_allocation('iso_qr') 
  end if

  if(ALLOCATED(iso_qs)) then
     DEALLOCATE(iso_qs,STAT=ierr)
     if (ierr.ne.0) call failed_allocation('iso_qs')  
  end if

  if(ALLOCATED(iso_qg)) then
     DEALLOCATE(iso_qg,STAT=ierr)
     if (ierr.ne.0) call failed_allocation('iso_qg')  
  end if

  !kzm: tracers
  if(ALLOCATED(trx)) then
     DEALLOCATE(&
          trx,fluxbtrx,trxadv,trxdiff,trxwsb,trxwle,trxwleadv, &
          trxwlediff,trx2leadv,trx2lediff,trx2legrad,trx2lediss,STAT=ierr) 
     if (ierr.ne.0) call failed_allocation('x grid tracer')
  end if

  if(ALLOCATED(try)) then
     DEALLOCATE(&
          try,fluxbtry,tryadv,trydiff,trywsb,trywle,trywleadv, &
          trywlediff,try2leadv,try2lediff,try2legrad,try2lediss,STAT=ierr) 
     if (ierr.ne.0) call failed_allocation('y grid tracer')
  end if

  if(ALLOCATED(trz)) then
     DEALLOCATE(&
          trz,fluxbtrz,trzadv,trzdiff,trzwsb,trzwle,trzwleadv, &
          trzwlediff,trz2leadv,trz2lediff,trz2legrad,trz2lediss,STAT=ierr) 
     if (ierr.ne.0) call failed_allocation('z grid tracer')
  end if

  if(ALLOCATED(trzz)) then
     DEALLOCATE(&
          trzz,fluxbtrzz,trzzadv,trzzdiff,trzzwsb,trzzwle,trzzwleadv, &
          trzzwlediff,trzz2leadv,trzz2lediff,trzz2legrad,trzz2lediss,STAT=ierr)
     if (ierr.ne.0) call failed_allocation('zz grid tracer')
  end if

  if(ALLOCATED(tro3)) then
     DEALLOCATE(&
          tro3,fluxbtro3,tro3adv,tro3diff,tro3wsb,tro3wle,tro3wleadv, &
          tro3wlediff,tro32leadv,tro32lediff,tro32legrad,tro32lediss,STAT=ierr)
     if (ierr.ne.0) call failed_allocation('o3 grid tracer')
  end if

end subroutine deallocate_uwvars

subroutine failed_allocation(output_string)
  use grid, only: rank
  implicit none

  character(LEN=*), intent(in) :: output_string

  write(*,*) 'Failed to (de)allocate ', output_string, &
       ' variables on processor ', rank
  call task_abort()

end subroutine failed_allocation

