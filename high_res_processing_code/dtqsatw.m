function dq = dtqsatw(t,p)
% t is temperature (K)
% p is pressure    (mb)
% from sat.f90

dq=0.622*dtesatw(t)/p;
        
