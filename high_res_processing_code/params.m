
% from params.f90

%   Constants:

cp = 1004.;             % Specific heat of air, J/kg/K
ggr = 9.81;             % Gravity acceleration, m/s2
lcond = 2.5104e+06;     % Latent heat of condensation, J/kg
lfus = 0.3336e+06;      % Latent heat of fusion, J/kg
lsub = 2.8440e+06;      % Latent heat of sublimation, J/kg
rv = 461.;              % Gas constant for water vapor, J/kg/K
rgas = 287.;            % Gas constant for dry air, J/kg/K
diffelq = 2.21e-05;     % Diffusivity of water vapor, m2/s
therco = 2.40e-02;      % Thermal conductivity of air, J/m/s/K
muelq = 1.717e-05;      % Dynamic viscosity of air

fac_cond = lcond/cp;
fac_fus = lfus/cp;
fac_sub = lsub/cp;

%  Microphysics stuff:

% Densities of hydrometeors

rhor = 1000.; % Density of water, kg/m3
rhos = 100.;  % Density of snow, kg/m3
rhog = 400.;  % Density of graupel, kg/m3

% Temperatures limits for various hydrometeors

tbgmin = 253.16;    % Minimum temperature for cloud water., K
tbgmax = 273.16;    % Maximum temperature for cloud ice, K
tprmin = 268.16;    % Minimum temperature for rain, K
tprmax = 283.16;    % Maximum temperature for snow+graupel, K
tgrmin = 223.16;    % Minimum temperature for snow, K
tgrmax = 283.16;    % Maximum temperature for graupel, K

% Terminal velocity coefficients

a_rain = 842.; % Coeff.for rain term vel 
b_rain = 0.8;  % Fall speed exponent for rain
a_snow = 4.84; % Coeff.for snow term vel
b_snow = 0.25; % Fall speed exponent for snow
a_grau = 94.5; % Lin (1983) (rhog=400)
b_grau = 0.5;  % Fall speed exponent for graupel


% Autoconversion

qcw0 = 1.e-3;      % Threshold for water autoconversion, g/g  
qci0 = 1.e-4;      % Threshold for ice autoconversion, g/g
alphaelq = 1.e-3;  % autoconversion of cloud water rate coef
betaelq = 1.e-3;   % autoconversion of cloud ice rate coef

% Accretion

erccoef = 1.0;   % Rain/Cloud water collection efficiency
esccoef = 1.0;   % Snow/Cloud water collection efficiency
esicoef = 0.1;   % Snow/cloud ice collection efficiency
egccoef = 1.0;   % Graupel/Cloud water collection efficiency
egicoef = 0.1;   % Graupel/Cloud ice collection efficiency

% Interseption parameters for exponential size spectra

nzeror = 8.e6;   % Intercept coeff. for rain  
nzeros = 3.e6;   % Intersept coeff. for snow
nzerog = 4.e6;   % Intersept coeff. for graupel

% Cloud droplet distribution properties, used in sedimentation scheme.
Nc0 = 65.; % cloud droplet concentration [cm^{-3}]

qp_threshold = 1.e-8; % minimal rain/snow water content

% from setparm.f90
a_bg = 1./(tbgmax-tbgmin);
a_pr = 1./(tprmax-tprmin);
a_gr = 1./(tgrmax-tgrmin);

