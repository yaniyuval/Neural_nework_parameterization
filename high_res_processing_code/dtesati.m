function de = dtesati(t)
% t is temperature (K)
% from sat.f90

a0 = 0.503223089;
a1 = 0.377174432e-1;
a2 = 0.126710138e-2;
a3 = 0.249065913e-4;
a4 = 0.312668753e-6;
a5 = 0.255653718e-8;
a6 = 0.132073448e-10;
a7 = 0.390204672e-13;
a8 = 0.497275778e-16;

dt = max(-80.,t-273.16);

de = a0 + dt.*(a1+dt.*(a2+dt.*(a3+dt.*(a4+dt.*(a5+dt.*(a6+dt.*(a7+a8.*dt)))))));

