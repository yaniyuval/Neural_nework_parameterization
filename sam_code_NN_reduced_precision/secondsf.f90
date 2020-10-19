subroutine secondsf(time)
  implicit none
  double precision, intent(inout) :: time
  integer, dimension(2) :: timeofday

  call gtodfff(timeofday)
  time = timeofday(1) + 1.d-6*timeofday(2)

end subroutine secondsf
