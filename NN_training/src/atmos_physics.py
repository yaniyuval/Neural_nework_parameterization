import numpy as np

# physical constants
L = 2.5104e6  # latent heat of condensation [J/kg]
g = 9.81  # gravitational acceleration [m/s^2]
cp = 1004.0  # specific heat capacity dry air [J/kg/K]
Lf = 0.3336e+06 # Latent heat of fusion [J/kg]
tprmin = 268.16 # Minimum temperature for rain, K
tprmax = 283.16 # Maximum temperature for snow+graupel, K
tbgmin = 253.16 # Minimum temperature for cloud water, K
tbgmax = 273.16 # Maximum temperature for cloud ice, K



# scripts to calculate atmospheric variables

def vertical_diff(rho, z):

    # follow vertical differencing from setgrid.f90 in SAM
    # changed indexing from starting at 1 to 0
    nzm = z.size

    adz = np.zeros(nzm)
    
    dz = 0.5*(z[0]+z[1]) 
    adz[0] = 1.

    for k in range(1,nzm-1): # range doesn't include stopping number
       adz[k] = 0.5*(z[k+1]-z[k-1])/dz

    adz[nzm-1] = (z[nzm-1]-z[nzm-2])/dz

    rho_dz = adz*dz*rho

    return rho_dz

def vertical_integral(data, rho, z):
# vertical integral with respect to sigma

    rho_dz = vertical_diff(rho, z)
    if rho_dz[None,:].shape[1] != data.shape[1]:
        print('The vertical integral is not performed on the whole column, using only ' + str(data.shape[1]) + ' levels')
    rho_dz_int = rho_dz[None,:data.shape[1]]
    # int_data = np.sum(data * rho_dz[None,:], axis=1)
    int_data = np.sum(data * rho_dz_int, axis=1)


    return int_data

def energy_tendency_residual(tabs, t_tend, q_tend, rho, z, output_vert_vars, o_dict):
# error in tendency of column-averaged liquid/ice static energy tendency in K/s 

# t corresponds to hL but excluding precipitating condensate:
# cp*t = cp*Tabs -Lc*qc-Ls*qi
# Equation A3 of Kharioutdinov et al 2003 gives energy conservation
# for hL that includes precipitating condensate
# We are working in limit where qp immediately falls out as precipitation
# and thus \partial qp/ partial time is zero such that \partial hL/\partial time
# is just cp*\partial t/partial time and we will refer to the column integral
# of this divided by cp as LHS. 
# Furthermore, equation 3 gives that
# the column integrated change in hL is L_sfc*P_tot_sfc which we will
# call RHS after dividing by cp 
# (noting that SAM actually uses L*P_tot where L = Lc+(1-omp)*Lf rather than what is written in the precip divergence term in equation 3)

   lhs_tend = vertical_integral(t_tend,rho,z)

   a_pr = 1.0/(tprmax-tprmin)
   omp = np.maximum(0.0,np.minimum(1.0,(tabs-tprmin)*a_pr))
   fac = (L + Lf*(1.0-omp))/cp
   # the following is correct since the sgs flux divergence in q_tend integrates to zero
   if 'qpout' in output_vert_vars:
        # precip = 0 # No need to correct for precipitation since it is in t_tend
        rhs_tend = -fac[:, 0] * vertical_integral(q_tend + o_dict['qpout'], rho, z) #I think it is a problem that I use the same fac here for both q's, but I have no other option since I want the cancelation of the mic tend.
   else:
        precip = calc_precip(q_tend,rho,z, output_vert_vars, o_dict)
        rhs_tend = fac[:,1]*precip

   vert_coordinate_norm = np.zeros([1,q_tend.shape[1]]) + 1.0
   tend_normalized = (lhs_tend-rhs_tend)/vertical_integral(vert_coordinate_norm, rho, z)

   return tend_normalized

def calc_precip(q, rho, z, output_vert_vars,o_dict):
# surface precipitation rate given tendency of specific humidity
    if 'qpout' in output_vert_vars:
        precip = -vertical_integral(q + o_dict['qpout'], rho, z)
    else:
        precip = -vertical_integral(q, rho, z)
    rho_dz = vertical_diff(rho, z)
    if 'qsout' in output_vert_vars:
        # precip=precip+np.squeeze(o_dict['qsout'])*rho_dz[0]
        precip = precip + np.squeeze(o_dict['qsout']) * rho_dz[0]
    return precip


def sam_qsat(T,p):
# formulation used in SAM

    esatw = sam_esatw(T)
    esati = sam_esati(T)

    qw = 0.622 * esatw/np.maximum(esatw,p-esatw);
    qi = 0.622 * esati/np.maximum(esati,p-esati);

    a_bg = 1.0/(tbgmax-tbgmin)
    omn = np.maximum(0.0,np.minimum(1.0,(T-tbgmin)*a_bg))
    qsat = omn*qw+(1-omn)*qi

    return qsat


def sam_esatw(T):
# SAM saturation vapor pressure for liquid water
# from sat.f90
# T is in K
# e is converted here to Pa

    a0 = 6.105851
    a1 = 0.4440316
    a2 = 0.1430341e-1
    a3 = 0.2641412e-3
    a4 = 0.2995057e-5
    a5 = 0.2031998e-7
    a6 = 0.6936113e-10
    a7 = 0.2564861e-13
    a8 = -0.3704404e-15

    dt = np.maximum(-80.,T-273.16)

    e = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))

    return e*100


def sam_esati(T):
# SAM saturation vapor pressure for ice
# from sat.f90
# T is in K
# e is converted here to Pa

    a0 = 6.11147274
    a1 = 0.503160820
    a2 = 0.188439774e-1
    a3 = 0.420895665e-3
    a4 = 0.615021634e-5
    a5 = 0.602588177e-7
    a6 = 0.385852041e-9
    a7 = 0.146898966e-11
    a8 = 0.252751365e-14

    dt = np.maximum(-80.,T-273.16)

    e = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))

    return e*100

def sam_qsatw(t,p):
# % t is temperature (K)
# % p is pressure    (mb)
# % from sat.f90

    esat = sam_esatw(t) /100 #(the 100 factor is to be consistant with the matlab)
    q = 0.622 * esat/np.maximum(esat,p-esat)

    return q
