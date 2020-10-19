function e = esatw(t)
% from sat.f90
% t is temperature (K)
% e is in hPa

a0 = 6.105851;
a1 = 0.4440316;
a2 = 0.1430341e-1;
a3 = 0.2641412e-3;
a4 = 0.2995057e-5;
a5 = 0.2031998e-7;
a6 = 0.6936113e-10;
a7 = 0.2564861e-13;
a8 = -0.3704404e-15;

dt = max(-80.,t-273.16);

e = a0 + dt.*(a1+dt.*(a2+dt.*(a3+dt.*(a4+dt.*(a5+dt.*(a6+dt.*(a7+a8.*dt)))))));
        
