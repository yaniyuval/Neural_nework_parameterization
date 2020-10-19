!this program sets the tracers based on the tracer option
	subroutine settracer()
	
	use vars
	use params
        use rad, only: o3
	implicit none
	integer i, j, k, kk, it, jt

        ! Initialize tro3 and o3g0 sounding here.
        if(dotro3) then
           call tracesini() ! get o3 sounding
           o3=o3*o3factor ! scale if desired
           do k=1,nzm
              o3g0(k)=o3(nz-k) !kzm save background o3 
           end do
        end if

        !bloss: all tracers are now saved in restart file.
        if((dorestart_tracer).or.(firststep.and.(nrestart.eq.0))) then
           !Below resets the values
           if(ALLOCATED(trx)) then
              trx=0.
              fluxbtrx = 0.
           end if

           if(ALLOCATED(try)) then
              try=0.
              fluxbtry = 0.
           end if

           if(ALLOCATED(trz)) then
              trz=0.
              fluxbtrz = 0.
           end if

           if(ALLOCATED(trzz)) then
              trzz=0.
              fluxbtrzz = 0.
           end if

           if(ALLOCATED(tro3)) then
              do k=1,nzm ! interactive ozone initialization.
                 do j=1,ny
                    do i=1,nx
                       tro3(i,j,k)=o3(nz-k) ! note that o3 order is reversed
                    end do
                 end do
              end do
              fluxbtro3 = 0.
           end if
        end if

        ! Initialize surface flux for tracer_option=6
        if(tracer_option.eq.6) fluxbtrx = rhow(1)*itautrx

        if(tracer_option .eq. 99) then 
           call task_rank_to_index(rank,it,jt)
           do k=1,nzm
              do j=1,ny
                 do i=1,nx
                    trx(i,j,k)=i+it-1
                    try(i,j,k)=j+jt-1
                    trz(i,j,k)=k
                    trzz(i,j,k)=float(k)*k
                 end do
              end do
           end do
        endif

	if(tracer_option .eq. 5) then 
           k=1
           do while(z(k).lt.2000)
              trx(:,:,k)=1
              k=k+1
           end do
           k=1
           do while(z(k).ge.2000.and.z(k).lt.10000)
              try(:,:,k)=1
              k=k+1
           end do
           k=1
           do while(z(k).ge.10000.and.z(k).lt.15000)
              trz(:,:,k)=1
              k=k+1
           end do
	endif
	
	return
	end
	
