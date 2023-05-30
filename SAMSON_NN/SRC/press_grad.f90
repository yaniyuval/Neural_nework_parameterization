
subroutine press_grad
	
!       pressure term of the momentum equations

use vars
implicit none
	
real *8 rdx,rdy,rdz
integer i,j,k,kb,jb,ib
real du(nx,ny,nz,3)
real wrk

wrk=1.
if(dorave) then
   wrk=ravefactor*ravefactor
endif

rdx=1./dx
rdy=1./dy

if(dostatis) then
        
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        du(i,j,k,1)=dudt(i,j,k,na)
        du(i,j,k,2)=dvdt(i,j,k,na)
        du(i,j,k,3)=dwdt(i,j,k,na)
      end do
    end do 
  end do  
          
endif
	
	 
do k=1,nzm
 kb=max(1,k-1)
 rdz = 1./(dz*adzw(k))
 do j=1,ny
  jb=j-YES3D
  do i=1,nx
   ib=i-1 
   dudt(i,j,k,na)=dudt(i,j,k,na)-(p(i,j,k)-p(ib,j,k))*rdx	
   dvdt(i,j,k,na)=dvdt(i,j,k,na)-(p(i,j,k)-p(i,jb,k))*rdy	
   dwdt(i,j,k,na)=dwdt(i,j,k,na)-(p(i,j,k)-p(i,j,kb))*rdz/wrk
  end do ! i
 end do ! j	
end do ! k
	
do k=1,nzm
 do j=1,ny
  do i=1,nx
    p(i,j,k)=p(i,j,k)*rho(k)  ! convert p'/rho to p'
  end do
 end do 
end do  

if(dowallx.and.mod(rank,nsubdomains_x).eq.0) then

    do k=1,nzm
     do j=1,ny
      dudt(1,j,k,na) = 0.
     end do
    end do

end if

if(dowally.and.RUN3D.and.rank.lt.nsubdomains_x) then

    do k=1,nzm
     do i=1,nx
      dvdt(i,1,k,na) = 0.
     end do
    end do

end if

if(dompi) then
   call task_bound_duvdt()
else
   call bound_duvdt()	   
endif


if(dostatis) then
        
 do k=1,nzm
  do j=1,ny
   do i=1,nx
     du(i,j,k,1)=dudt(i,j,k,na)-du(i,j,k,1)
     du(i,j,k,2)=dvdt(i,j,k,na)-du(i,j,k,2)
     du(i,j,k,3)=dwdt(i,j,k,na)-du(i,j,k,3)
   end do
  end do 
 end do  
          
 call stat_tke(du,tkelepress)
 call stat_mom(du,momlepress)
 call setvalue(twlepres,nzm,0.)
 call setvalue(qwlepres,nzm,0.)
 call setvalue(swlepres,nzm,0.)
 call stat_sw1(du,twlepres,qwlepres,swlepres)          

endif
	
call task_barrier()
	
end subroutine press_grad



