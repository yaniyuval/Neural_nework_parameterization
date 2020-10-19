! Initialize some arrays from vars module:

subroutine init()

  use vars

  implicit none

  radlwup = 0.
  radlwdn = 0.
  radswup = 0.
  radswdn = 0.
  radqrlw = 0.
  radqrsw = 0.
 
  tlat = 0. 
  tlatqi = 0.
  tadv = 0.
  tdiff = 0.
  qifall = 0.
  qadv = 0.
  qdiff = 0.
  qpadv = 0.
  qpdiff = 0.
  qpsrc = 0.
  qpfall = 0.
  qpevp = 0.

  usnd = 0.
  vsnd = 0.
  qsnd = 0.
  tsnd = 0.
  daysnd = 0.
  nsnd = 0

  dqls = 0.
  dtls = 0.
  ugls = 0.
  vgls = 0.
  wgls = 0.
  pres0ls = 0.
  dayls = 0.
  nlsf = 0

  dayrfc = 0.
  dtrfc = 0.
  nrfc = 0

  daysfc = 0.
  sstsfc = 0.
  hsfc = 0.
  lesfc = 0.
  tausfc = 0.
  nsfc = 0

  gamt0 = 0.  
  gamq0 = 0.
  ttend_wave = 0.
  qtend_wave = 0.

end subroutine init



