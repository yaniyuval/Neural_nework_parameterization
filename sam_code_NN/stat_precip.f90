	

subroutine stat_precip(precsfc, prec_conv, prec_strat)

use grid
implicit none

real precsfc(nx,ny)
real prec_conv, prec_strat, tmp(2)

real precip(nx_gl,ny_gl), buf(nx,ny,nsubdomains)
integer i,j,m,count,it,jt,reqs_in(nsubdomains),rf,tag
logical flag(nsubdomains)

real prec_bgnd
integer nrain,nradius_x,nradius_y,ii,jj,iii,jjj
logical convective(nx_gl, ny_gl)

	if(dompi) then

	  if(rank.ne.0) then

  	    prec_conv=-1.
	    prec_strat=-1.

            call task_bsend_float(0,precsfc,nx*ny,66)

	  else

	    do m = 1,nsubdomains-1 
              call task_receive_float(buf(1,1,m),nx*ny,reqs_in(m))
              flag(m) = .false.
            end do ! m         
                     
            count = nsubdomains-1
            do while (count .gt. 0)
              do m = 1,nsubdomains-1
               if(.not.flag(m)) then
                call task_test(reqs_in(m), flag(m), rf, tag) 
                if(flag(m)) then
                   count=count-1
                   call task_rank_to_index(rf,it,jt)
                    do j = 1,ny
                     do i = 1,nx
                       precip(i+it,j+jt) = buf(i,j,m)
                     end do
                    end do
                endif
               endif
              end do    
            end do


            call task_rank_to_index(0,it,jt)
            do j = 1,ny
              do i = 1,nx
               precip(i+it,j+jt) = precsfc(i,j)
              end do
            end do
           

          endif ! rank.ne.0


	else

	
          do j = 1,ny
            do i = 1,nx
               precip(i,j) = precsfc(i,j)
            end do
          end do


	end if ! dompi


	if(masterproc) then
	  
	  nradius_x = min(nx_gl/2.,11000./dx)
	  nradius_y = min(ny_gl/2.,11000./dy)
	  convective = .false.

	  do j=1,ny_gl
	   do i=1,nx_gl
	     if(precip(i,j).lt.0.01) CYCLE
	     if(precip(i,j).ge.10.) then
	       convective(i,j) = .true.
	       CYCLE 
	     endif
	     prec_bgnd=0.
	     nrain=0
	     do jj=j-nradius_y,j+nradius_y
	      jjj=mod(jj-1+3*ny_gl,ny_gl)+1
	      do ii=i-nradius_x,i+nradius_x
	        iii=mod(ii-1+3*nx_gl,nx_gl)+1
		if(((ii-i)*dx)**2+((jj-j)*dy)**2.le.11000.**2 &
			.and. precip(iii,jjj).ge.0.01) then
		   nrain=nrain+1
		   prec_bgnd=prec_bgnd+precip(iii,jjj)
	        endif 
	      end do
	     end do
	     prec_bgnd=prec_bgnd/float(nrain)
	     if(precip(i,j).gt.2.*prec_bgnd) convective(i,j)=.true.	
	   end do
	  end do

	  do j=1,ny_gl
	   do i=1,nx_gl

	     if(convective(i,j)) then

	     do jj=j-nradius_y,j+nradius_y
	      jjj=mod(jj-1+3*ny_gl,ny_gl)+1
	      do ii=i-nradius_x,i+nradius_x
	        iii=mod(ii-1+3*nx_gl,nx_gl)+1
	        
	       if(precip(iii,jjj).ge.1.2.and.precip(iii,jjj).lt.2.5 &
		 .and.((ii-i)*dx)**2+((jj-j)*dy)**2.le.2000.**2) then
		      convective(iii,jjj) = .true.
	       elseif(precip(iii,jjj).ge.2.5.and.precip(iii,jjj).lt.5.0 &
		 .and.((ii-i)*dx)**2+((jj-j)*dy)**2.le.3000.**2) then
		      convective(iii,jjj) = .true.
	       elseif(precip(iii,jjj).ge.5.0.and.precip(iii,jjj).lt.10.  &
		 .and.((ii-i)*dx)**2+((jj-j)*dy)**2.le.4000.**2) then
		      convective(iii,jjj) = .true.
	       endif

	      end do
	     end do

	     endif

	   end do
	  end do

	  prec_conv=0.
	  prec_strat=0.
	  do j=1,ny_gl
           do i=1,nx_gl
	    if(convective(i,j)) then
	     prec_conv = prec_conv + precip(i,j)	     
	    else
	     prec_strat = prec_strat + precip(i,j)	     
	    endif
	   end do
	  end do

	  prec_conv = prec_conv/float(nx_gl*ny_gl)
	  prec_strat = prec_strat/float(nx_gl*ny_gl)

	endif



	if(dompi) then

	  tmp(1) = prec_conv
	  call task_max_real(tmp(1),tmp(2),1)
	  prec_conv = tmp(2)
	  tmp(1) = prec_strat
	  call task_max_real(tmp(1),tmp(2),1)
	  prec_strat = tmp(2)

	endif

end



