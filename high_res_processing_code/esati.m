function e = esati(t)
% t is temperature (K)
% e is in hPa
% from sat.f90


a0 = 6.11147274;
a1 = 0.503160820;
a2 = 0.188439774e-1;
a3 = 0.420895665e-3;
a4 = 0.615021634e-5;
a5 = 0.602588177e-7;
a6 = 0.385852041e-9;
a7 = 0.146898966e-11;
a8 = 0.252751365e-14;

dt = max(-80.,t-273.16);

e = a0 + dt.*(a1+dt.*(a2+dt.*(a3+dt.*(a4+dt.*(a5+dt.*(a6+dt.*(a7+a8.*dt)))))));

%disp('mod single precision')
%
%a0 = single(6.11147274);
%a1 = single(0.503160820);
%a2 = single(0.188439774e-1);
%a3 = single(0.420895665e-3);
%a4 = single(0.615021634e-5);
%a5 = single(0.602588177e-7);
%a6 = single(0.385852041e-9);
%a7 = single(0.146898966e-11);
%a8 = single(0.252751365e-14);
%
%t = single(t);
%dt = max(-80.,t-273.16);
%
%e = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))));
        
