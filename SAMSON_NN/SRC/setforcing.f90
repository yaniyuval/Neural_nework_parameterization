subroutine setforcing()
	
use vars
use params
implicit none
	
integer ndmax,i,n,iz
parameter (ndmax = 1000)
real zz(ndmax),tt(ndmax),qq(ndmax),pp(ndmax)
real uu(ndmax),vv(ndmax),ww(ndmax),rad(ndmax)
real ratio1, ratio2, ratio_t1, ratio_t2
logical zgrid

!-------------------------------------------------
!	Read sounding arrays:


if(masterproc) print*,'Initializeing forcing data...'

open(77,file='./'//trim(case)//'/snd', status='old',form='formatted') 
read(77,*)
nsnd=1
do while(.true.)
  read(77,err=44,end=44,fmt=*) daysnd(nsnd),n
  do i=1,n
      read(77,*) zz(i),pp(i),tt(i),qq(i),uu(i),vv(i)
  end do
  if(zz(2).gt.zz(1)) then
   zgrid=.true.
  else if(pp(2).lt.pp(1)) then
   zgrid=.false.
  else
   if(masterproc) print*,'error in grid in snd'
   stop
  end if
  do iz = 1,nzm
   if(zgrid) then	
     do i = 2,n
      if(z(iz).le.zz(i)) then
        tsnd(iz,nsnd)=tt(i-1)+(tt(i)-tt(i-1))/(zz(i)-zz(i-1))*(z(iz)-zz(i-1))	
        tsnd(iz,nsnd) = tsnd(iz,nsnd) /prespot(iz)
        qsnd(iz,nsnd)=qq(i-1)+(qq(i)-qq(i-1))/(zz(i)-zz(i-1))*(z(iz)-zz(i-1))
        usnd(iz,nsnd)=uu(i-1)+(uu(i)-uu(i-1))/(zz(i)-zz(i-1))*(z(iz)-zz(i-1))
        vsnd(iz,nsnd)=vv(i-1)+(vv(i)-vv(i-1))/(zz(i)-zz(i-1))*(z(iz)-zz(i-1))
        goto 11
      endif
     end do
   else
     do i = 2,n
      if(pres(iz).ge.pp(i)) then
        tsnd(iz,nsnd)=tt(i-1)+(tt(i)-tt(i-1))/(pp(i)-pp(i-1))*(pres(iz)-pp(i-1))
        tsnd(iz,nsnd) = tsnd(iz,nsnd) /prespot(iz)
        qsnd(iz,nsnd)=qq(i-1)+(qq(i)-qq(i-1))/(pp(i)-pp(i-1))*(pres(iz)-pp(i-1))
        usnd(iz,nsnd)=uu(i-1)+(uu(i)-uu(i-1))/(pp(i)-pp(i-1))*(pres(iz)-pp(i-1))
        vsnd(iz,nsnd)=vv(i-1)+(vv(i)-vv(i-1))/(pp(i)-pp(i-1))*(pres(iz)-pp(i-1))
	goto 11
      endif
     end do
   end if

   call atmosphere(z(iz-1)/1000.,ratio1,ratio2,ratio_t1)
   call atmosphere(z(iz)/1000.,ratio1,ratio2,ratio_t2)

   tsnd(iz,nsnd)=ratio_t2/ratio_t1*tsnd(iz-1,nsnd)
   qsnd(iz,nsnd)=max(0.,2.*qsnd(iz-1,nsnd)-qsnd(iz-2,nsnd))
   usnd(iz,nsnd)=usnd(iz-1,nsnd)
   vsnd(iz,nsnd)=vsnd(iz-1,nsnd)

11 continue

  end do ! iz

  nsnd=nsnd+1
  if(nsnd.gt.nmaxsnd.and.masterproc) then
    print*,'Maximum dimension of sounding arrays is smaller'// &
		   'than the number of sounding time-points'
    call task_abort()
  endif
end do ! while...
44 continue
nsnd=nsnd-1  

if(nsnd.eq.1.and.masterproc) then
  print*,'Error: minimum two sounding profiles are needed.'
  call task_abort()
endif

if(masterproc)print*,'Observed sounding interval (days):',daysnd(1),daysnd(nsnd)

!-------------------------------------------------
!	Read Large-scale forcing arrays:



if(dolargescale.or.dosubsidence) then

open(77,file='./'//trim(case)//'/lsf',status='old',form='formatted') 
read(77,*)
nlsf=1
do while(.true.)
  read(77,err=55,end=55,fmt=*) dayls(nlsf),n, pres0ls(nlsf)
  do i=1,n
      read(77,*) zz(i),pp(i),tt(i),qq(i),uu(i),vv(i),ww(i)
  end do
!kzm accelerate the tendencies 
  tt=tt*sqrt(betafactor)
  qq=qq*sqrt(betafactor)
  ww=ww*sqrt(betafactor)
 

  if(zz(2).gt.zz(1)) then
   zgrid=.true.
  else if(pp(2).lt.pp(1)) then
   zgrid=.false.
  else
   if(masterproc) print*,'error in grid in lsf'
   stop
  end if
  do iz = 1,nzm
   if(zgrid) then
    do i = 2,n
     if(z(iz).le.zz(i)) then
       dtls(iz,nlsf)=tt(i-1)+(tt(i)-tt(i-1))/(zz(i)-zz(i-1))*(z(iz)-zz(i-1))	
       dqls(iz,nlsf)=qq(i-1)+(qq(i)-qq(i-1))/(zz(i)-zz(i-1))*(z(iz)-zz(i-1))
       ugls(iz,nlsf)=uu(i-1)+(uu(i)-uu(i-1))/(zz(i)-zz(i-1))*(z(iz)-zz(i-1)) 	
       vgls(iz,nlsf)=vv(i-1)+(vv(i)-vv(i-1))/(zz(i)-zz(i-1))*(z(iz)-zz(i-1)) 	
       wgls(iz,nlsf)=ww(i-1)+(ww(i)-ww(i-1))/(zz(i)-zz(i-1))*(z(iz)-zz(i-1)) 	
       dosubsidence = dosubsidence .or. wgls(iz,nlsf).ne.0. 
       goto 12
     endif
    end do
   else
    do i = 2,n
     if(pres(iz).ge.pp(i)) then
       dtls(iz,nlsf)=tt(i-1)+(tt(i)-tt(i-1))/(pp(i)-pp(i-1))*(pres(iz)-pp(i-1))	
       dqls(iz,nlsf)=qq(i-1)+(qq(i)-qq(i-1))/(pp(i)-pp(i-1))*(pres(iz)-pp(i-1))
       ugls(iz,nlsf)=uu(i-1)+(uu(i)-uu(i-1))/(pp(i)-pp(i-1))*(pres(iz)-pp(i-1)) 
       vgls(iz,nlsf)=vv(i-1)+(vv(i)-vv(i-1))/(pp(i)-pp(i-1))*(pres(iz)-pp(i-1)) 
       wgls(iz,nlsf)=ww(i-1)+(ww(i)-ww(i-1))/(pp(i)-pp(i-1))*(pres(iz)-pp(i-1))	
       dosubsidence = dosubsidence .or. wgls(iz,nlsf).ne.0. 
       goto 12
     endif
    end do
   end if
   dtls(iz,nlsf)=0.
   dqls(iz,nlsf)=0.
   ugls(iz,nlsf)=ugls(iz-1,nlsf)
   vgls(iz,nlsf)=vgls(iz-1,nlsf)
   wgls(iz,nlsf)=0.
12 continue

  end do
  
  nlsf=nlsf+1
  if(nlsf.gt.nmaxlsf.and.masterproc) then
    print*,'Maximum dimension of forcing arrays is smaller'// &
		   'than the number of l.s. forcing time-points'
    call task_abort()
  endif
end do
55 continue
if(nlsf.gt.2) then
   dayls(nlsf)=dayls(nlsf-1)+dayls(nlsf-1)-dayls(nlsf-2)
else
   dayls(nlsf)=dayls(nlsf-1)+100.
end if
pres0ls(nlsf) = pres0ls(nlsf-1)
do iz = 1,nzm
   dtls(iz,nlsf) = dtls(iz,nlsf-1)
   dqls(iz,nlsf) = dqls(iz,nlsf-1)  
   ugls(iz,nlsf) = ugls(iz,nlsf-1)
   vgls(iz,nlsf) = vgls(iz,nlsf-1)
   wgls(iz,nlsf) = wgls(iz,nlsf-1) 
end do

if(masterproc)print*,'Large-Scale Forcing interval (days):',dayls(1),dayls(nlsf)

endif

!-------------------------------------------------
!	Read Radiation forcing arrays:



if(doradforcing) then

open(77,file='./'//trim(case)//'/rad',status='old',form='formatted') 

read(77,*)
nrfc=1
do while(.true.)
  read(77,err=66,end=66,fmt=*) dayrfc(nrfc),n
  do i=1,n
      read(77,*) pp(i),rad(i)
  end do

  !kzm accelerate rad
  rad=rad*sqrt(betafactor)

  if(pp(2).gt.pp(1)) then
   zgrid=.true.
  else if(pp(2).lt.pp(1)) then
   zgrid=.false.
  else
   if(masterproc) print*,'error in grid in rad'
   stop
  end if
  do iz = 1,nzm
   if(zgrid) then
    do i = 2,n
     if(z(iz).le.pp(i)) then
      dtrfc(iz,nrfc)=rad(i-1)+(rad(i)-rad(i-1))/(pp(i)-pp(i-1))*(z(iz)-pp(i-1))	
      goto 13
     endif
    end do
   else
    do i = 2,n
     if(pres(iz).ge.pp(i)) then
      dtrfc(iz,nrfc)=rad(i-1)+(rad(i)-rad(i-1))/(pp(i)-pp(i-1)) &
                                                         *(pres(iz)-pp(i-1))	
      goto 13
     endif
    end do
   end if
   dtrfc(iz,nrfc)=0.
13 continue
  end do
  nrfc=nrfc+1
  if(nrfc.gt.nmaxrfc.and.masterproc) then
    print*,'Maximum dimension of forcing arrays is smaller'//  &
                       'than the number of rad.forcing time-points'
	    call task_abort()
  endif
end do
66 continue

if(nrfc.gt.2) then
   dayrfc(nrfc)=dayrfc(nrfc-1)+dayrfc(nrfc-1)-dayrfc(nrfc-2)
else
   dayrfc(nrfc)=dayrfc(nrfc-1)+100.
end if
do iz = 1,nzm
   dtrfc(iz,nrfc) =  dtrfc(iz,nrfc-1)
end do
if(masterproc)print*,'Radiative Forcing interval (days):',dayrfc(1),dayrfc(nrfc)

endif ! doradforcing

!-------------------------------------------------------
! Read Surface Forcing Arrays:
!

if(dosfcforcing) then

open(77,file='./'//trim(case)//'/sfc',status='old',form='formatted') 
read(77,*)

nsfc=1
do while(.true.)
  read(77,err=77,end=77,fmt=*) daysfc(nsfc),  &
                               sstsfc(nsfc),hsfc(nsfc),lesfc(nsfc),tausfc(nsfc)
  nsfc=nsfc+1
  if(nsfc.gt.nmaxsfc.and.masterproc) then
    print*,'Maximum dimension of forcing arrays is smaller'// &
 		   'than the number of surface forcing time-points'
    call task_abort()
  endif
end do
77 continue
close(77)

if(nsfc.gt.2) then
   daysfc(nsfc)=daysfc(nsfc-1)+daysfc(nsfc-1)-daysfc(nsfc-2)
else
   daysfc(nsfc)=daysfc(nsfc-1)+100.
end if
do iz = 1,nzm
   sstsfc(nsfc) = sstsfc(nsfc-1)
   hsfc(nsfc) = hsfc(nsfc-1)
   lesfc(nsfc) = lesfc(nsfc-1)
   tausfc(nsfc) = tausfc(nsfc-1)  
end do


if(masterproc)print*,'Surface Forcing interval (days):',daysfc(1),daysfc(nsfc)

end if

!kzm
daysnd=daysnd/sqrt(betafactor)
dayls=dayls/sqrt(betafactor)
dayrfc=dayrfc/sqrt(betafactor)
daysfc=daysfc/sqrt(betafactor)

end

