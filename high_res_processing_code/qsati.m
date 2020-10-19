function q = qsati(t,p)
% t is temperature (K)
% p is pressure    (mb)
% from sat.f90

esat=esati(t);

q = 0.622 * esat./max(esat,p-esat);
        
        
      
