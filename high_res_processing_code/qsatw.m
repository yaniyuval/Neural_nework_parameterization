function q = qsatw(t,p)
% t is temperature (K)
% p is pressure    (mb)
% from sat.f90

esat = esatw(t);

q = 0.622 * esat./max(esat,p-esat);
        
