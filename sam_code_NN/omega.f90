real function omegan(tabs)
use params
implicit none
real tabs
omegan = max(0.,min(1.,(tabs-tbgmin)*a_bg))
return
end

real function omegap(tabs)
use params
implicit none
real tabs
omegap = max(0.,min(1.,(tabs-tprmin)*a_pr))
return
end

real function omegag(tabs)
use params
implicit none
real tabs
omegag = max(0.,min(1.,(tabs-tgrmin)*a_gr))
return
end

