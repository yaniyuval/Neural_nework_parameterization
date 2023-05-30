
subroutine init_linear_wave

!
! parameterized large-scale linear wave dynamics
!added the following global variables
!alfafactor: for the meridional scale difference of convection and wave
!wsub_conv(nzm): vertical velocity on ITCZ scale
!nstartlinearwave, nsteplinearwave,nsteplinearwavebg
use vars
use params
implicit none

!   doupperbound=.false.


   t_wavebg_hist=0.
   q_wavebg_hist=0.
   heating_wavebg_hist=0.
   heating_anom_hist=0.
   t_wave_local=0.
   q_wave_local=0.
   h_wave_local=0.
   t_wave=0.
   q_wave=0.
   w_wave=0.
   wwave_conv=0.
   heating_wave=0.
   t_wavebg=0.
   q_wavebg=0.
   heating_wavebg=0.
   wavecounter=1
   hanomcounter=1
!!$   if(.not.doevolvewavebg) then
!!$      t_wave_local=tabs0
!!$      q_wave_local=q0
!!$      t_wavebg=tabs0
!!$      q_wavebg=q0
!!$      t_wave=tabs0
!!$      q_wave=q0
!!$   endif

end subroutine init_linear_wave
