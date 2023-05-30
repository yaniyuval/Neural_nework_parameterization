!this program renews parts of the tracerfield based on the tracer option
	subroutine renewtracer()
	
	use vars
	use params
	implicit none
	integer i, j, k, kk, it, jt
	
!destruct o3 in the subcloud layer (z<600m)	  	
	if(dotro3) then 
		k=1
	  	do while(z(k).lt.600)
	  	 do j=1,ny
	   	  do i=1,nx
			tro3(i,j,k)=o3g0(k)	
	    	  end do
	  	 end do
		 k=k+1
	  	end do
	endif
!put tracer trz in the subcloud layer (z<600m)	  	
	if(tracer_option .eq. 1) then 
		k=1
	  	do while(z(k).lt.600)
	  	 do j=1,ny
	   	  do i=1,nx
			trz(i,j,k)=1
	    	  end do
	  	 end do
		 k=k+1
	  	end do
	endif
!put tracer trz in the subcloud layer and trzz in the bulk troposphere (z<=10km)
	if(tracer_option .eq. 2) then 
		k=1
	  	do while(z(k).lt.600)
	  	 do j=1,ny
	   	  do i=1,nx
			trz(i,j,k)=1
	    	  end do
	  	 end do
		 k=k+1
	  	end do
		k=1
	  	do while(z(k).lt.10000.)
	  	 do j=1,ny
	   	  do i=1,nx
			trzz(i,j,k)=1
	    	  end do
	  	 end do
		 k=k+1
	  	end do
	endif
!put tracer trz and trzz in the subcloud layer (z<600m)	  	
	if(tracer_option .eq. 3) then 
		k=1
	  	do while(z(k).lt.600)
	  	 do j=1,ny
	   	  do i=1,nx
			trz(i,j,k)=1
			trzz(i,j,k)=1
	    	  end do
	  	 end do
		 k=k+1
	  	end do
	endif
	
!put tracer trz, try and trzz in the  bulk troposphere (z<=10km)	  	
	if(tracer_option .eq. 4) then 
		k=1
	  	do while(z(k).lt.600)
	  	 do j=1,ny
	   	  do i=1,nx
			trz(i,j,k)=1
			try(i,j,k)=1
			trzz(i,j,k)=1
	    	  end do
	  	 end do
		 k=k+1
	  	end do
	endif
	
        if(dotrz) then
           call setvalue(trzwleadv,nzm,0.)
           call setvalue(trzwlediff,nzm,0.)
        end if
        if(dotrx) then
           call setvalue(trxwleadv,nzm,0.)
           call setvalue(trxwlediff,nzm,0.)
        end if
        if(dotry) then
           call setvalue(trywleadv,nzm,0.)
           call setvalue(trywlediff,nzm,0.)
        end if
        if(dotrzz) then
           call setvalue(trzzwleadv,nzm,0.)
           call setvalue(trzzwlediff,nzm,0.)
        end if
        if(dotro3) then
           call setvalue(tro3wleadv,nzm,0.)
           call setvalue(tro3wlediff,nzm,0.)
        end if

	return
	end
