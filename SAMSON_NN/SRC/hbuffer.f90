module hbuffer

use grid
use params

implicit none

	integer hbuf_length, HBUF_MAX_LENGTH
	parameter(HBUF_MAX_LENGTH = 500)
	
	real hbuf(HBUF_MAX_LENGTH*nzm+30)	
	character *8 namelist(HBUF_MAX_LENGTH)	
	character *80 deflist(HBUF_MAX_LENGTH)	
	character *10 unitlist(HBUF_MAX_LENGTH)	
	integer status(HBUF_MAX_LENGTH)
	integer average_type(HBUF_MAX_LENGTH)

CONTAINS

!------------------------------------------------------------

subroutine hbuf_average(n)

!       Average the profiles in buffer

use vars        

implicit none
integer l, n 
real coef

coef=1./float(n)
do l = 1,hbuf_length*nzm
 hbuf(l) = hbuf(l) *coef
end do  

end subroutine hbuf_average

!------------------------------------------------------------


subroutine hbuf_avg_put(name, f, dimx1,dimx2,dimy1,dimy2,dimz,factor)

!       Write to buffer an averaged 3D field (not to file yet)

implicit none

integer dimx1, dimx2, dimy1, dimy2, dimz
real f(dimx1:dimx2, dimy1:dimy2, dimz), fz(nzm), factor
character *(*) name
integer l, begin, k
logical flag

flag=.false.
do  l = 1,hbuf_length
  if(.not.(lgt(name,namelist(l)).or.llt(name,namelist(l)))) then
     flag=.true.
     if(status(l).gt.0) then
         status(l) = 999
         call averageXY(f,dimx1,dimx2,dimy1,dimy2,dimz,fz)
         begin = (l-1)*nzm
         do k = 1,nzm
           hbuf(begin+k) = hbuf(begin+k) + fz(k) * factor
         end do
     endif
  endif
end do
if(.not.flag.and.masterproc) print*,name

end subroutine hbuf_avg_put

!------------------------------------------------------------

subroutine hbuf_flush

!       Flush the buffer


use vars
implicit none
integer l,k

do l=1,hbuf_length*nzm
  hbuf(l) = 0.
end do
do k=1,nzm
 cloud_factor(k) = 0.
 core_factor(k) = 0.
 coredn_factor(k) = 0.
end do

s_acld=0.
s_acldl=0.
s_acldm=0.
s_acldh=0.
s_acldcold=0.
s_acldisccp=0.
s_acldlisccp=0.
s_acldmisccp=0.
s_acldhisccp=0.
s_ar=0.
w_max=0.
p_conv=0.
p_strat=0.
s_flns=0.
s_flnt=0.
s_flnsc=0.
s_flntc=0.
s_flds=0.
s_fsns=0.
s_fsnt=0.
s_fsnsc=0.
s_fsntc=0.
s_fsds=0.
s_solin=0.
lhobs=0.
shobs=0.

end subroutine hbuf_flush

!------------------------------------------------------------
subroutine hbuf_init

!       Read list of vertical profile names to store in file

use grid, only: case, masterproc

implicit none
character *8 nm
character *80 def
character *10 un
integer stat,count,type,i
integer lenstr
external lenstr
character*3 filestatus

open(66,file='./'//case(1:lenstr(case))//'/lst',&
                          status='old',form='formatted')

! first determine how many entries are gonna be:

hbuf_length = 0
111    continue
read(66,err=111,end=222,fmt=*) stat,type,nm
if(stat.gt.0) hbuf_length = hbuf_length+1
goto 111
222    continue

if(hbuf_length.gt.HBUF_MAX_LENGTH) then
   print *,'Error: hbuf_length > HBUF_MAX_LENGTH.'
   call task_abort()
endif

! fill the buffers:

rewind(66)
count = 0
333    continue
read(66,err=333,end=444,fmt=*) stat,type,nm,def,un
if(stat.gt.0) then
   count = count + 1
   namelist(count) = nm
   deflist(count) = def
   unitlist(count) = un
   status(count) = stat
   average_type(count) = type
endif
goto 333
444    continue
if(masterproc) then
         print *,'Number of statistics profiles:', hbuf_length
         print *,'Statistics profiles to save:'
         print *,(namelist(i),i=1,hbuf_length)

! make sure that the stat file doesn't exist if a new run to prevent
! accidental overwrite

  filestatus='old'
  if(nrestart.eq.0.or.nrestart.eq.2) then
    filestatus='new'
  end if


  open (55,file='./'//case(1:lenstr(case))//'/'// &
                  case(1:lenstr(case))//'_'// &
                  caseid(1:lenstr(caseid))//'.stat', &
                  status=filestatus,form='unformatted')

  close(55)

end if

close (66)

call hbuf_flush()

end subroutine hbuf_init


!-----------------------------------------------------------------

subroutine hbuf_put(name, f, factor)

!       Write to buffer (not to file yet)

use grid, only: masterproc
implicit none

real f(nzm), factor
character *(*) name
integer l, begin, k
logical flag

flag=.false.
do  l = 1,hbuf_length
   if(.not.(lgt(name,namelist(l)).or.llt(name,namelist(l)))) then
       flag=.true.
       if(status(l).gt.0) then
          status(l) = 999
          begin = (l-1)*nzm
          do k = 1,nzm
            hbuf(begin+k) = hbuf(begin+k) + f(k)*factor
          end do
       endif
   endif
end do
if(.not.flag.and.masterproc) print*,'name ', name,' is missing in "lst" file.'

end subroutine hbuf_put

!----------------------------------------------------------------

subroutine hbuf_write(n)

!       Write buffer to file

use vars
implicit none
integer l, k, n, ntape, length
real coef,aver,factor,tmp(nzm),tmp1(nzm)
real cloud_f(nzm),core_f(nzm),coredn_f(nzm)
real hbuf1(HBUF_MAX_LENGTH*nzm+30)
integer nbuf

data ntape/55/

aver=1./float(n)
factor=1./float(nx*ny)

if(dompi) then
  coef = 1./float(nsubdomains)
  do k=1,nzm
     cloud_f(k)=cloud_factor(k)/(cloud_factor(k)+1.e-10)
     core_f(k)=core_factor(k)/(core_factor(k)+1.e-10)
     coredn_f(k)=coredn_factor(k)/(coredn_factor(k)+1.e-10)
  end do
  do k=1,nzm
     tmp(k)=cloud_f(k)
  end do
  call task_sum_real(tmp,cloud_f,nzm)
  do k=1,nzm
     tmp(k)=core_f(k)
  end do
  call task_sum_real(tmp,core_f,nzm)
  do k=1,nzm
     tmp(k)=coredn_f(k)
  end do
  call task_sum_real(tmp,coredn_f,nzm)
  do k=1,nzm
     cloud_factor(k)=cloud_factor(k)*cloud_f(k)*coef
     core_factor(k)=core_factor(k)*core_f(k)*coef
     coredn_factor(k)=coredn_factor(k)*coredn_f(k)*coef
  end do
end if

do k=1,nzm
  cloud_factor(k) = float(n)/(cloud_factor(k)+1.e-10)
  core_factor(k) = float(n)/(core_factor(k)+1.e-10)
  coredn_factor(k) = float(n)/(coredn_factor(k)+1.e-10)
end do

length = 0
do l = 1,hbuf_length
  if(status(l).eq. 999) then
    length = length+1
    if(average_type(l).eq.1) then
      do k=1,nzm
        hbuf((l-1)*nzm+k) = hbuf((l-1)*nzm+k)*cloud_factor(k)
      end do
    endif
    if(average_type(l).eq.2) then
      do k=1,nzm
        hbuf((l-1)*nzm+k) = hbuf((l-1)*nzm+k)*core_factor(k)
      end do
    endif
    if(average_type(l).eq.3) then
      do k=1,nzm
        hbuf((l-1)*nzm+k) = hbuf((l-1)*nzm+k)*coredn_factor(k)
      end do
    endif
  endif
end do


!  Get statistics buffer from different processes, add them together
!  and average

if(dompi) then

   coef = 1./float(nsubdomains)

   tmp1(1)=w_max
   call task_max_real(tmp1,tmp,1)
   w_max=tmp(1)
   k=hbuf_length*nzm
   hbuf(k+1)=s_acld
   hbuf(k+2)=s_acldl
   hbuf(k+3)=s_acldm
   hbuf(k+4)=s_acldh
   hbuf(k+5)=s_acldcold
   hbuf(k+6)=s_ar
   hbuf(k+7)=s_flns
   hbuf(k+8)=s_flnt
   hbuf(k+9)=s_flnsc
   hbuf(k+10)=s_flntc
   hbuf(k+11)=s_flds
   hbuf(k+12)=s_fsns
   hbuf(k+13)=s_fsnt
   hbuf(k+14)=s_fsnsc
   hbuf(k+15)=s_fsntc
   hbuf(k+16)=s_fsds
   hbuf(k+17)=s_solin
   hbuf(k+18)=s_acldisccp
   hbuf(k+19)=s_acldlisccp
   hbuf(k+20)=s_acldmisccp
   hbuf(k+21)=s_acldhisccp
   nbuf = k+21
   call task_sum_real(hbuf,hbuf1,nbuf)
   hbuf(1:nbuf) = hbuf1(1:nbuf)*coef
   s_acld=hbuf(k+1)
   s_acldl=hbuf(k+2)
   s_acldm=hbuf(k+3)
   s_acldh=hbuf(k+4)
   s_acldcold=hbuf(k+5)
   s_ar=hbuf(k+6)
   s_flns=hbuf(k+7)
   s_flnt=hbuf(k+8)
   s_flnsc=hbuf(k+9)
   s_flntc=hbuf(k+10)
   s_flds=hbuf(k+11)
   s_fsns=hbuf(k+12)
   s_fsnt=hbuf(k+13)
   s_fsnsc=hbuf(k+14)
   s_fsntc=hbuf(k+15)
   s_fsds=hbuf(k+16)
   s_solin=hbuf(k+17)
   s_acldisccp=hbuf(k+18)
   s_acldlisccp=hbuf(k+19)
   s_acldmisccp=hbuf(k+20)
   s_acldhisccp=hbuf(k+21)

endif
if(masterproc) then

  open (ntape,file='./'//case(1:lenstr(case))//'/'// &
                  case(1:lenstr(case))//'_'// &
                  caseid(1:lenstr(caseid))//'.stat', &
                  status='unknown',form='unformatted')
  if(nstep.ne.nstat) then
    do while(.true.)
       read (ntape,end=222)
    end do
  endif
222  continue
  backspace(ntape)

  print *,'Writting history file ',caseid
  write(ntape) caseid
  write(ntape)  day,dt,nstep,nx,ny,nz,nzm,dx,dy,dz,adz,z,pres,tabs_s,pres0, &
                s_acld*aver*factor,s_ar*factor/(nstat+1.e-5), &
                s_acldcold*aver*factor,w_max,p_conv*aver*24.,p_strat*aver*24.,&
                s_flns*aver*factor,s_flnt*aver*factor,s_flnsc*aver*factor, &
                s_flntc*aver*factor,s_flds*aver*factor,s_fsns*aver*factor, &
                s_fsnt*aver*factor,s_fsnsc*aver*factor,s_fsntc*aver*factor, &
                s_fsds*aver*factor,s_solin*aver*factor, &
                sstobs,lhobs*aver,shobs*aver, &
                s_acldl*aver*factor,s_acldm*aver*factor,s_acldh*aver*factor, &
                s_acldisccp, &
                s_acldlisccp,s_acldmisccp,s_acldhisccp

  write(ntape) length
  do l = 1,hbuf_length
     if(status(l).eq. 999) then
         write(ntape) namelist(l)
         write(ntape) deflist(l)
         write(ntape) unitlist(l)
         write(ntape) (hbuf((l-1)*nzm+k),k=1,nzm)
     end if
  end do
  close (ntape)

end if

end subroutine hbuf_write



end module hbuffer
