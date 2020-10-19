
subroutine coriolis

use vars

implicit none
	
real u_av, v_av, w_av
integer i,j,k,ib,ic,jb,jc,kc
	
if(RUN3D) then

do k=1,nzm
 kc=k+1
 do j=1,ny
  jb=j-1
  jc=j+1
  do i=1,nx
   ib=i-1
   ic=i+1
    v_av=0.25*(v(i,j,k)+v(i,jc,k)+v(ib,j,k)+v(ib,jc,k))
    w_av=0.25*(w(i,j,kc)+w(ib,j,kc)+w(i,j,k)+w(ib,j,k))
!    dudt(i,j,k,na)=dudt(i,j,k,na)+0.5*(fcory(j) + fcory(jc))*(v_av-vg0(k))-fcorzy(j)*w_av !JY YANI 
    dudt(i,j,k,na)=dudt(i,j,k,na)+(fcory(j))*(v_av-vg0(k))-fcorzy(j)*w_av
    u_av=0.25*(u(i,j,k)+u(ic,j,k)+u(i,jb,k)+u(ic,jb,k))
!    dvdt(i,j,k,na)=dvdt(i,j,k,na)-fcory(j)*(u_av-ug0(k))
 !   dvdt(i,j,k,na)=dvdt(i,j,k,na)-0.5*(fcory(j)+fcory(jb))*(u_av-ug0(k)) !Yani modified
    !had to put weighted averages on the cory parameter because it was calculated in a wrong way in the high resolution simulatuibs
    dvdt(i,j,k,na)=dvdt(i,j,k,na)-(0.5625*fcory(j)+0.4375*fcory(jb))*(u_av-ug0(k)) !Yani modified
!    print *,'cor parameter:', 0.625*fcory(j)+0.375*fcory(jb), 0.5*fcory(j)+0.5*fcory(jb), fcory(j)
!    if (i.eq.1.and.k.eq.1) then
!    print *,'cor parameter:', fcory(j), fcory(jc), j 
!    end if
  end do ! i
 end do ! j
end do ! k

else

do k=1,nzm
 kc=k+1
 do j=1,ny
  do i=1,nx
   ib=i-1
   ic=i+1
   w_av=0.25*(w(i,j,kc)+w(ib,j,kc)+w(i,j,k)+w(ib,j,k))
   dudt(i,j,k,na)=dudt(i,j,k,na)+fcory(j)*(v(i,j,k)-vg0(k))-fcorzy(j)*w_av
   dvdt(i,j,k,na)=dvdt(i,j,k,na)-fcory(j)*(u(i,j,k)-ug0(k))
  end do ! i
 end do ! i
end do ! k

endif
	
end subroutine coriolis

