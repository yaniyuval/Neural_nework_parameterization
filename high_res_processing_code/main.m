%This code takes the high resolution fields from Fortran, and propagate it
%a time step, calculate tendencies and fluxes, coarse grain it, and
%calculates tendencies and fluxes also for the coarse grain fields

%%%%%%%%%%
% Step 01:
% step 01-01: definitions and flags
% step 01-02: path to data
% step 01-03: read base z data from files
% step 01-04: import params and init_precip
% step 01-05: Initialize path to files/ coarse sst x y

% Step 02: Main loop - for the high resolution data:
% Step 02:01 define path of file
% Step 02:02 announce variables
% Step 02:02b read input variables
% Step 02:03 calculate derived quantities
% Step 02:03b get coarse grained variables before they change
% Step 02:04 Set grid (setgrid.f90) for calculations
% Step 02:05:Calculate shear_prod3D.f90 for the tke/tkz calculation
% Step 02:06: tke/tkz
% Step 02:07: Advect variables
% Step 02:07a: Vertiacl momentum advection
% Step 02:07b precip_fall calculation
% Step 02:08: surface fluxes
% Step 02:09: vertical diffusion
% Step 02:10: horizontal diffusion
% Step 02:11: do cloud (for recalculation of qn, and sedimentation)
% Step 02:12: calculate microphysical tendency of qp

% Step 03: Main loop over coarse grained quanteties:
% Step 03:03 calculate coarse quantities
% Step 03:05:Calculate shear_prod3D.f90 for the tke/tkz calculation
% Step 03:06: tke/tkz
% Step 03:07: Advect variables
% Step 03:07a: Vertiacl momentum advection
% Step 03:07b precip_fall calculation
% Step 03:08: surface fluxes
% Step 03:09: vertical diffusion
% Step 03:10: horizontal diffusion
% Step 03:11: do cloud (for recalculation of qn, and sedimentation)
% Step 03:12: calculate microphysical tendency of qp
 % Step 03:13: calculate rh in the end of the time step - for
            % coarse grained tendencies and the resolved tendencies:

% Step 04 output to netcdf
%Step 04:01 calculate resituals
%Step 04:02 permute fields
%Step 04:03 write to netcdf

%%%%%%%%%

% Step 1

clear all
close all
clc
is_cheyenne = 1;
is_test = 0; %flag when I want to test the script  - can compare to the results of the test code
plot_plots  = 0;
plot_coarse_vs_residuals = 0;
do_show_times = 1;
dimfactor_as_in_samson = 0; %1 if we want the calculations to follow samson - 0 if we want the coarse and high to have the same dimfactor.
ssty_from_txt = 1; %This is a new thing - to read the SST and not interpolate it.

sprintf('good luck')
if is_test
    compare_advection_t_flux = 1; %Good acuracy
    compare_advection_t_flux_z = 1; %Good acuracy
    compare_advection_t_flux_x = 1; %Good acuracy
    compare_advection_t_x_tend = 1; %Good acuracy
    compare_advection_t_flux_y = 1; %Good acuracy
    compare_advection_t_y_tend = 1; %Good acuracy
    dqp_fall_tendencies = 1;
    t_fall_tendencies = 1;
    do_compare_tkz_vs_approx = 1;
    do_approx_tke_vs_real = 1;
    def2_compare = 1;
    t_diff_plot = 1; %I have large errors in the  vertical diffusion. I am not sure if we could really use it - to think!
    t_z_diff_tend_plot = 1;
    compare_x_diffusion=1;
    compare_y_diffusion =1;
    compare_x_diffusion_tend = 1;
    compare_y_diffusion_tend = 1;
    plot_dqp_mic_tend = 1;
    
    %From full advection
    compare_qt_total_advection_tendency =1;
    compare_qp_total_advection_tendency =1;
    
    compare_advection_qp_flux_y =1;
    compare_advection_qp_flux_x =1;
    compare_advection_qp_flux_z = 1;
    compare_advection_qt_flux_y =1;
    compare_advection_qt_flux_x =1;
    compare_advection_qt_flux_z =1;
end

do_sedimentation = 1;
calc_advection_tend = 0;
calc_diffusive_tend = 0;

% timestep in seconds (needed in microphysics routines, but outputs are tendencies)
% step 01-01: definitions and flags
dtn = 24;
do_fall_tend_tfull = 1;
add_cloud_tfull_tend = 1;
do_horizontal_diffusion_at_all = 1;
do_fall_tend_qp = 1;
advect_fields = 1; % If I want to actually advect the fields
diffuse_fields =1;

if is_test
    start_time = 1008450;
    end_time = start_time; %Yani- try
else
    if is_cheyenne
       start_time = 810000 +(1-1)*31*450;
       end_time = min(start_time + 30*450, 2520000);  
    else
    start_time =360000;
    end_time = start_time + 1000.*450; % 2520000;
    end
end

interval = 450;
ravefactor = 4.0; % used in diffusion

if is_test
    resolutions = 16;
    fac_redefine_dx_dy = 16;
else
    resolutions = [4,8,16,32];
    fac_redefine_dx_dy = 1;
end

dx =   12000.0.*fac_redefine_dx_dy;
dy =   12000.0.*fac_redefine_dx_dy;

exper = {'qobs'};
filename_profile = 'sounding_z_rho_rhow.txt';
delimiterIn = ' ';
headerlinesIn = 1;
sounding = importdata(filename_profile,delimiterIn,headerlinesIn);
rhow_input = flip(sounding.data(:,3));
rho = flip(sounding.data(:,2));
zin = flip(sounding.data(:,1));

filename_profile = 'sounding_z_pres0_tabs0.txt';
sounding = importdata(filename_profile,delimiterIn,headerlinesIn);
pres0 = flip(sounding.data(:,2));
tabs0 = flip(sounding.data(:,3));
num_z = length(rho);
rhow(1:num_z,1) = rhow_input;
rhow(num_z+1,1) = 2*rhow_input(num_z)-rhow_input(num_z-1);

% get SST distribution for use in diffusive fsslux calculation
filename_profile = 'ssty.txt';
profile = importdata(filename_profile,delimiterIn,headerlinesIn);
ssty = profile.data(:,1);

if is_test
    exper_path = ['/net/aimsir/archive1/janniy/sam_tests/DATA3D/'];
elseif is_cheyenne
    exper_path = ['/glade/scratch/janniy/bill_crm_data/'];
else
    exper_path = ['/net/aimsir/archive1/pog/bill_crm_data/'];
end

% step 01-03: read base z data from files
% get rho, rhow, pres0 and tabs0 profiles
% pres0 is the inital pres and tabs0 is the initial tabs0; both are used in precip_init

% step 01-04: import params and init_precip
% parameters
params

% initialize microphysics
precip_init

% step 01-05: Initialize path to files/ coarse sst x y
%% initialize variables in loop
for i_exper = 1:length(exper)
    time_index = start_time;
    cycle_init = 1;
    do_loop_init = 1;
    while do_loop_init
        if is_test
            filename_base_init = [exper_path, exper{i_exper}, '_32x90x48_12_', num2str(time_index, '%010d'), '_000', num2str(cycle_init)];
            filename_init = [filename_base_init, '.nc'];
        else
            filename_base_init = [exper_path, exper{i_exper}, 'km12x576/', exper{i_exper}, 'km12x576_576x1440x48_ctl_288_', num2str(time_index, '%010d'), '_000', num2str(cycle_init)];
            filename_init = [filename_base_init, '.nc4'];
        end
        if exist(filename_init,'file')
            do_loop_init = 0;
        else
            cycle_init = cycle_init+1;
            if cycle_init>8
                disp(filename_init)
                error('no filename')
            end
        end
    end
    
    disp(filename_init)
    if i_exper==1 & time_index==start_time
        
        x     = my_ncread(filename_init, 'x'); % m
        y     = my_ncread(filename_init, 'y'); % m
        z     = my_ncread(filename_init, 'z'); % m
        pres  = my_ncread(filename_init, 'p'); % hPa
        
        if max(abs(zin-z))>0.1
            error('inconsistent sounding file')
        end
        
        num_x = length(x);
        num_y = length(y);
        num_z = length(z);
        
        
        sstxy = repmat(ssty, [1, num_x]);
        if is_test
            multiple_space = 1;
            num_blocks_x = num_x/multiple_space;
            num_blocks_y = num_y/multiple_space;
            x_coarse = zeros(num_blocks_x,1);
            y_coarse = zeros(num_blocks_y,1);
            filename_profile = sprintf('ssty_16.txt');
            profile = importdata(filename_profile);
            ssty_coarse = profile(:);
            sstxy_coarse = repmat(ssty_coarse, [1, length(x_coarse)]); %TO CHECK
            sstxy = sstxy_coarse;
        end
    end
end


%% Main loop
% Step 02: Main loop - for the high resolution data
% Step 02:01 define path of file
times_1 = start_time:interval:end_time;
sprintf('Main loop')
for i_exper = 1:length(exper)
%     parfor dummy_time = 1:length(times_1)
    for dummy_time = 1:length(times_1)
        sprintf('starting measuring time')
        tStart = tic;
        tic
        time_index = times_1(dummy_time);
        cycle = 1;
        do_loop = 1;
        while do_loop
            if is_test
                filename_base = [exper_path, exper{i_exper}, '_32x90x48_12_', num2str(time_index, '%010d'), '_000', num2str(cycle_init)];
                filename = [filename_base, '.nc'];
            else
                filename_base = [exper_path, exper{i_exper}, 'km12x576/', exper{i_exper}, 'km12x576_576x1440x48_ctl_288_', num2str(time_index, '%010d'), '_000', num2str(cycle)];
                filename = [filename_base, '.nc4'];
            end
            if exist(filename,'file')
                do_loop = 0;
            else
                cycle = cycle+1;
                if cycle>8
                    disp(filename)
                    error('no filename')
                end
            end
        end
        disp(filename)
        if is_cheyenne
            outfilename_janni =['/glade/scratch/janniy/ML_convection_data/', exper{i_exper}, '/', exper{i_exper}, 'km12x576_576x1440x48_ctl_288_', num2str(time_index, '%010d'), '_000', num2str(cycle), '_diff_coarse_space_corrected_tkz'];
        else
            outfilename_janni =['/net/aimsir/archive1/janniy/ML_convection_data_cheyenne/', exper{i_exper}, 'km12x576_576x1440x48_ctl_288_', num2str(time_index, '%010d'), '_000', num2str(cycle), '_diff_coarse_space_corrected_tkz'];
        end
        sum_files = 0;
        for res_i = resolutions
            outfilename_janni_tmp = [outfilename_janni, num2str(res_i), '.nc4' ];
            if exist(outfilename_janni_tmp, 'file') == 2
                sum_files = sum_files + 1; %
                sprintf('the file exist:')
                disp(outfilename_janni_tmp)
            end
        end
        if sum_files == length(resolutions)
            sprintf('Already calculated all the fields - skipping calc for this time step')
            continue
        end
        sprintf('the dummy index is %i and time is %i',dummy_time,time_index)
        if do_show_times, sprintf('step 02:01: %f', toc), end
        tic
        %%
        % Step 02:02 announce variables        
        tkz_out = zeros(num_z,num_y,num_x);
        
        %dqp tendencies
        dqp = zeros(num_z,num_y,num_x);
        
        %advection fluxes and tendencies:
        tfull_flux_x = zeros(num_z,num_y,num_x);
        tfull_flux_y = zeros(num_z,num_y,num_x);
        tfull_flux_z = zeros(num_z,num_y,num_x);
        t_flux_x = zeros(num_z,num_y,num_x);
        t_flux_y = zeros(num_z,num_y,num_x);
        t_flux_z = zeros(num_z,num_y,num_x);
        if calc_advection_tend
            tfull_flux_x_tend = zeros(num_z,num_y,num_x);
            tfull_flux_y_tend = zeros(num_z,num_y,num_x);
            tfull_flux_z_tend = zeros(num_z,num_y,num_x);
            t_flux_x_tend = zeros(num_z,num_y,num_x);
            t_flux_y_tend = zeros(num_z,num_y,num_x);
            t_flux_z_tend = zeros(num_z,num_y,num_x);
        end
        
        %diffusion fluxes and tendencies
        t_diff_flx_x = zeros(num_z,num_y,num_x);
        t_diff_flx_y = zeros(num_z,num_y,num_x);
        t_diff_flx_z = zeros(num_z,num_y,num_x);
        
        tfull_diff_flx_x = zeros(num_z,num_y,num_x);
        tfull_diff_flx_y = zeros(num_z,num_y,num_x);
        tfull_diff_flx_z = zeros(num_z,num_y,num_x);
        
        qt_diff_flx_x = zeros(num_z,num_y,num_x);
        qt_diff_flx_y = zeros(num_z,num_y,num_x);
        qt_diff_flx_z = zeros(num_z,num_y,num_x);
        
        qp_diff_flx_x = zeros(num_z,num_y,num_x);
        qp_diff_flx_y = zeros(num_z,num_y,num_x);
        qp_diff_flx_z = zeros(num_z,num_y,num_x);
        
        if calc_diffusive_tend
            t_diff_x_tend = zeros(num_z,num_y,num_x);
            t_diff_y_tend = zeros(num_z,num_y,num_x);
            t_diff_z_tend = zeros(num_z,num_y,num_x);
            
            tfull_diff_x_tend = zeros(num_z,num_y,num_x);
            tfull_diff_y_tend = zeros(num_z,num_y,num_x);
            tfull_diff_z_tend = zeros(num_z,num_y,num_x);
            
            qt_diff_x_tend = zeros(num_z,num_y,num_x);
            qt_diff_y_tend = zeros(num_z,num_y,num_x);
            qt_diff_z_tend = zeros(num_z,num_y,num_x);
            
            qp_diff_x_tend = zeros(num_z,num_y,num_x);
            qp_diff_y_tend = zeros(num_z,num_y,num_x);
            qp_diff_z_tend = zeros(num_z,num_y,num_x);
        end
        
        %cloud
        cloud_lat_heat = zeros(num_z,num_y,num_x);
        cloud_qt_tend = zeros(num_z,num_y,num_x);
        
        %         tke_approx = zeros(num_z,num_y,num_x); %since I calculate the tke
        %         with iterations I don't use it (can't coarse grain it!)
        %         tke_approx_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
        fz = zeros(num_z,num_y,num_x);
        fzt = zeros(num_z,num_y,num_x);
        Pr1 =  zeros(num_z,num_y,num_x);
        
        windspeed_high =  zeros(num_y,num_x);
        fluxbu = zeros(num_y,num_x);
        fluxbv = zeros(num_y,num_x);
        fluxbu_samson = zeros(num_y,num_x);
        fluxbv_samson = zeros(num_y,num_x);
        if do_show_times, sprintf('step 02:02 : %f', toc), end
        tic
        %%
        % Step 02:02b read input variables
        u = my_ncread(filename, 'U'); % m/s
        v = my_ncread(filename, 'V'); % m/s
        w = my_ncread(filename, 'W'); % m/s
        qv = my_ncread(filename, 'Q')/1000.0; % water vapor (kg/kg)
        qn = my_ncread(filename, 'QN')/1000.0; % non precip cond (water+ice) (kg/kg)
        qp = my_ncread(filename, 'QP')/1000.0; % precip (kg/kg)
        Qrad = my_ncread(filename, 'QRAD'); % rad heating rate (K/day)
        tabs = my_ncread(filename, 'TABS'); % absolute temperature (K)
        tkz = my_ncread(filename, 'tk_z'); % vertical diffusivity (m^2/s)
        if is_test
            u_i = my_ncread(filename, 'u_i'); % m/s
            v_i = my_ncread(filename, 'v_i'); % m/s
            w_i = my_ncread(filename, 'w_i'); % m/s
            qv_i = my_ncread(filename, 'q_i')/1000.0; % water vapor (kg/kg)
            qn_i = my_ncread(filename, 'qn_i')/1000.0; % non precip cond (water+ice) (kg/kg)
            qp_i = my_ncread(filename, 'qp_i')/1000.0; % precip (kg/kg)
            tabs_i = my_ncread(filename, 'tabs_i'); % absolute temperature (K)
            tkz_i = my_ncread(filename, 'tk_z_i'); % vertical diffusivity (m^2/s)
            tkhz_i = my_ncread(filename, 'tkh_z_i'); % vertical diffusivity (m^2/s)
            t_i = my_ncread(filename, 't_i'); % absolute temperature (K)
            tke_i = my_ncread(filename, 'tke_i'); % absolute temperature (K
            [u_i, u] = switch_vars(u_i,u); %Assigning the beginning of step values to the end of step and vice versa
            [v_i, v] = switch_vars(v_i,v);
            [w_i, e] = switch_vars(w_i,w);
            [qv_i, qv] = switch_vars(qv_i,qv);
            [qn_i, qn] = switch_vars(qn_i,qn);
            [qp_i, qp] = switch_vars(qp_i,qp);
            [tabs_i, tabs] = switch_vars(tabs_i,tabs);
            [tkz_i, tkz] = switch_vars(tkz_i,tkz);
            
            t_diff_tend = my_ncread(filename, 't_v_diff'); %compare to Python tendency
            t_diff_flux = my_ncread(filename, 't_d_flux'); %compare to Python tendency
            q_diff_flux = my_ncread(filename, 'q_v_diff'); %compare to Python tendency
            t_bottom_flux = my_ncread(filename, 't_bot_di');
            dummy11 = my_ncread(filename, 'dummy1');
            dummy22 = my_ncread(filename, 'dummy2');
            dummy33 = my_ncread(filename, 'dummy3');
            dummy44 = my_ncread(filename, 'dummy4');
            dqp_test = my_ncread(filename, 'dqp_test');
            dummy55 = my_ncread(filename, 'dummy5');
            
            dummy100 = my_ncread(filename, 'dummy10');
            dummy111 = my_ncread(filename, 'dummy11');
            dummy122 = my_ncread(filename, 'dummy12');
            dummy133 = my_ncread(filename, 'dummy13');
            dummy144 = my_ncread(filename, 'dummy14');
            dummy155 = my_ncread(filename, 'dummy15');
            dummy166 = my_ncread(filename, 'dummy16');
            dummy177 = my_ncread(filename, 'dummy17');
            dummy188 = my_ncread(filename, 'dummy18');
            dummy199 = my_ncread(filename, 'dummy19');
            dummy200 = my_ncread(filename, 'dummy20');
            dummy211 = my_ncread(filename, 'dummy21');
            dummy222 = my_ncread(filename, 'dummy22');
            
            dummy233 = my_ncread(filename, 'dummy23');
            dummy244 = my_ncread(filename, 'dummy24');
            
            y_flux_adv = my_ncread(filename, 'y_f_adv');
            y_tend_adv = my_ncread(filename, 'y_f_tend');
            
            qp_x_flux_adv = my_ncread(filename, 'qp_x_f_a');
            qp_x_tend_adv = my_ncread(filename, 'qp_x_t_a');
            qp_y_flux_adv = my_ncread(filename, 'qp_y_f_a');
            qp_y_tend_adv = my_ncread(filename, 'qp_y_t_a');
            qp_z_flux_adv = my_ncread(filename, 'qp_z_f_a');
            qp_z_tend_adv = my_ncread(filename, 'qp_z_t_a');
            
            qt_x_flux_adv = my_ncread(filename, 'qt_x_f_a');
            qt_x_tend_adv = my_ncread(filename, 'qt_x_t_a');
            qt_y_flux_adv = my_ncread(filename, 'qt_y_f_a');
            qt_y_tend_adv = my_ncread(filename, 'qt_y_t_a');
            qt_z_flux_adv = my_ncread(filename, 'qt_z_f_a');
            qt_z_tend_adv = my_ncread(filename, 'qt_z_t_a');
        end
        
        if do_show_times, sprintf('step 02:02b: %f', toc), end
        tic
        %%
        % Step 02:03 calculate derived quantities:
        
        %         liquid static energy and tendency of qp
        qt = qv+qn; % total non-precipitating water (referred to as simply q in SAM)
        
        % calculate liquid/ice static energy
        gamaz = zeros(size(tabs));
        for k=1:num_z
            gamaz(k,:,:)=ggr/cp*z(k);
        end
        
        % liquid/ice water static energy h_L divided by cp
        % but here set qp to zero in SAM with rf
        omn_3d  = max(0.,min(1.,(tabs-tbgmin)*a_bg));
        t = tabs + gamaz - (fac_cond+(1.-omn_3d).*fac_fus).*qn;
        
        % full t including precipitating condensates
        omp_3d  = max(0.,min(1.,(tabs-tprmin)*a_pr));
        tfull = tabs + gamaz - (fac_cond+(1.-omn_3d).*fac_fus).*qn - ...
            (fac_cond+(1.-omp_3d).*fac_fus).*qp;
        if do_show_times, sprintf('step 02:03: %f', toc), end
        tic
        %%
        % Step 02:04 Set grid (setgrid.f90) for calculations
        % from setgrid.f90
        
        dz = 0.5*(z(1)+z(2));
        adzw = zeros(length(num_z));% Yani added
        for k=2:num_z
            adzw(k) = (z(k)-z(k-1))/dz;
        end
        adzw(1) = 1.;
        adzw(num_z+1) = adzw(num_z);
        adz = zeros(length(num_z));% Yani added
        
        for k=2:num_z-1
            adz(k) = 0.5*(z(k+1)-z(k-1))/dz;
        end
        adz(1) = 1.;
        adz(num_z) = adzw(num_z);
        
        rdz2=1./(dz*dz);
        rdz=1./dz;
        
        if do_show_times, sprintf('step 02:04: %f', toc), end
        tic
        %%
        % Step 02:05:Calculate shear_prod3D.f90 for the tke/tkz calculation
        def2 = shear_prod3D(u,v,w,dx,dy,dz,adz,adzw,num_x,num_y,num_z,ravefactor);
        
        %%
        if do_show_times, sprintf('step 02:05: %f', toc), end
        tic
        % Step 02:06: tke/tkz
        %tke_full.f90 - correct its location as in the fortran code!
        
        dxdzfactor  = 1.0;
        dimfactor=max(1.,log10(sqrt(dx*dy)/dz/dxdzfactor/ravefactor));
        Ck=0.1;
        Cs = 0.1944;
        Cs1 = 0.14;
        Pr = 3.0;
        Ce=Ck^3/Cs^4;
        Ces=Ce/0.7*3.0;
        
        for k=1:num_z
            if k>1
                kb=k-1;
            else
                kb=k;
            end
            kc=min(k+1,num_z); 
            grd=dz*adz(k)*dimfactor;
            rhoi = rhow(kc)/adzw(kc);
            rdz5=0.5*rdz2;
            % following from setdata.f90
            bet = ggr/tabs0(k);
            % following from tke_full.f90
            betdz=bet/dz/(adzw(kc)+adzw(k));
            Ce1=Ce/0.7*0.19;
            Ce2=Ce/0.7*0.51;
            
            if k==1
                betdz=bet/dz/adzw(kc);
                Ce1=Ces/0.7*0.19;
                Ce2=Ces/0.7*0.51;
            end
            
            if k==num_z-1
                betdz=bet/dz/adzw(k);
                Ce1=Ces/0.7*0.19;
                Ce2=Ces/0.7*0.51;
            end
            
            for i=1:num_x
                for j=1:num_y
                    
                    omn = max(0.,min(1.,(tabs(k,j,i)-tbgmin)*a_bg));
                    omp = max(0.,min(1.,(tabs(k,j,i)-tprmin)*a_pr));
                    
                    lstarn = fac_cond+(1.-omn)*fac_fus;
                    lstarp = fac_cond+(1.-omp)*fac_fus;
                    
                    
                    if qn(k,j,i)>0.0
                        
                        dqsat = omn*dtqsatw(tabs(k,j,i),pres(k))+ ...
                            (1.-omn)*dtqsati(tabs(k,j,i),pres(k));
                        qsat = omn*qsatw(tabs(k,j,i),pres(k))+(1.-omn)*qsati(tabs(k,j,i),pres(k));
                        bbb = 1. + 0.61*qsat-qn(k,j,i) -qp(k,j,i)+1.61*tabs(k,j,i)*dqsat;
                        bbb = bbb / (1.+lstarn*dqsat);
                        % note replaced t with tfull in following line
                        buoy_sgs=betdz*(bbb*(tfull(kc,j,i)-tfull(kb,j,i)) ...
                            +(bbb*lstarn - (1.+lstarn*dqsat)*tabs(k,j,i))*(qt(kc,j,i)-qt(kb,j,i)) ...
                            +(bbb*lstarp - (1.+lstarp*dqsat)*tabs(k,j,i))*(qp(kc,j,i)-qp(kb,j,i)) );
                    else
                        
                        bbb = 1.+0.61*qt(k,j,i)-qp(k,j,i);
                        % note replaced t with tfull in following line
                        buoy_sgs=betdz*( bbb*(tfull(kc,j,i)-tfull(kb,j,i)) ...
                            +0.61*tabs(k,j,i)*(qt(kc,j,i)-qt(kb,j,i)) ...
                            +(bbb*lstarp-tabs(k,j,i))*(qp(kc,j,i)-qp(kb,j,i)) );
                    end
                 
                    if buoy_sgs<=0.0 % unstable/neutral
                        smix=grd;
                    else % stable
                        smix = sqrt(0.76*tkz(k,j,i)/sqrt(buoy_sgs)/Ck); 
                        smix = min(grd,max(0.1*grd,smix));
                    end
                    ratio=smix/grd;
                    Pr1(k,j,i)=1.+2.*ratio;
                    Cee=Ce1+Ce2*ratio;
                    tkz(k,j,i)=sqrt(Ck^3/Cee*max(0.,def2(k,j,i)-Pr1(k,j,i)*buoy_sgs))*smix^2;
                    %%
                    if Pr1(k,j,i)>3 || Pr1(k,j,i)<1.2
                        error('Pr1 out of range')
                    end
                    
                end
            end
        end
        
        
        if do_show_times, sprintf('step 02:06: %f', toc), end
        tic
        %%
        
        % Step 02:07: Advect variables
        %simple version of scalar advection Advection
        
        % vertical fluxes
        % note w is on half levels but only first nzm are output by SAM
        % w is set to zero at top level n_z+1
        % w is first multiplied by rhow*dtn/dz in adams.f90 (note w and rhow
        % are grouped together in the differencing)
        % from advect_scalar3D.f90:
        % kb=max(1,k-1)
        % www(i,j,k)=max(0.,w(i,j,k))*f(i,j,kb)+min(0.,w(i,j,k))*f(i,j,k)
        % assuming dt->0 such that not iterative
        % leaving out non-oscillatory option,
        % see smolarkiewicz 2006 for general discussion
        
        
        % from setgrid.f90 - I think I need it again since some of the
        % vars were changed
        dz = 0.5*(z(1)+z(2));
        adzw = zeros(length(num_z));% Yani added
        for k=2:num_z
            adzw(k) = (z(k)-z(k-1))/dz;
        end
        adzw(1) = 1.;
        adzw(num_z+1) = adzw(num_z);
        adz = zeros(length(num_z));% Yani added
        
        for k=2:num_z-1
            adz(k) = 0.5*(z(k+1)-z(k-1))/dz;
        end
        adz(1) = 1.;
        adz(num_z) = adzw(num_z);
        
        
        rdz2=1./(dz*dz);
        rdz=1./dz;
        
        
        for k=1:num_z
            kb = max(1,k-1);
            w_rhow = w(k,:,:)*rhow(k);
            t_flux_z(k,:,:) = 0.5*w_rhow.*(t(k,:,:) + t(kb,:,:));%Yani simplified
            tfull_flux_z(k,:,:)  = 0.5*w_rhow.*(tfull(k,:,:) + tfull(kb,:,:));
        end
        rhow_num_x = repmat(rho,[1,num_y]);
        for i=1:num_x
            ib = i-1;
            if ib == 0
                ib = num_x; %Try to impose periodic domain - Yani ?
            end
            u_rhow = u(:,:,i).*rhow_num_x;
            t_flux_x(:,:,i) = 0.5*u_rhow.*(t(:,:,ib) + t(:,:,i));
            tfull_flux_x(:,:,i) = 0.5*u_rhow.*(tfull(:,:,ib) + tfull(:,:,i));
        end
        rhow_num_y = repmat(rho,[1,num_x]);
        for j=1:num_y
            jb = max(1,j-1);
            v_rhow = squeeze(v(:,j,:)).*rhow_num_y;
            t_flux_y(:,j,:) = 0.5*v_rhow.*squeeze(t(:,jb,:) + t(:,j,:));
            tfull_flux_y(:,j,:) = 0.5*v_rhow.*squeeze(tfull(:,jb,:) + tfull(:,j,:));
        end
        
        
        %advect fields without calculation of tendencies
        if (calc_advection_tend == 0 && advect_fields) %If I want to advect things without calculating them
            for k = 1:num_z
                if k < num_z
                    t(k,:,:) =t(k,:,:) - dtn.*(t_flux_z(k+1,:,:) - t_flux_z(k,:,:))./(adz(k).*dz.*rho(k));
                    tfull(k,:,:) =tfull(k,:,:) - dtn.*(tfull_flux_z(k+1,:,:) - tfull_flux_z(k,:,:))./(adz(k).*dz.*rho(k));
                elseif k == num_z
                    t(k,:,:) =t(k,:,:) - dtn.*(0.0 - t_flux_z(k,:,:))./(adz(k).*dz.*rho(k));
                    tfull(k,:,:) =tfull(k,:,:) - dtn.*(0.0 - tfull_flux_z(k,:,:))./(adz(k).*dz.*rho(k));
                else
                    sprintf('there is some error in the z indices advection')
                end
            end
            for i = 1:num_x
                if i < num_x
                    t(:,:,i) =t(:,:,i) - dtn.*(t_flux_x(:,:,i+1) - t_flux_x(:,:,i))./(dx.*rhow_num_x);
                    tfull(:,:,i) =tfull(:,:,i) - dtn.*(tfull_flux_x(:,:,i+1) - tfull_flux_x(:,:,i))./(dx.*rhow_num_x);
                elseif i == num_x
                    t(:,:,i) =t(:,:,i) - dtn.*(t_flux_x(:,:,1) - t_flux_x(:,:,i))./(dx.*rhow_num_x);
                    tfull(:,:,i) =tfull(:,:,i) - dtn.*(tfull_flux_x(:,:,1) - tfull_flux_x(:,:,i))./(dx.*rhow_num_x);
                else
                    sprintf('there is some error in the x indices advection')
                end
                
            end
            
            for j = 1:num_y
                if j < num_y
                    t(:,j,:) =squeeze(t(:,j,:)) - dtn.*squeeze(t_flux_y(:,j+1,:) - t_flux_y(:,j,:))./(dy.*rhow_num_y);
                    tfull(:,j,:) =squeeze(tfull(:,j,:)) - dtn.*squeeze(tfull_flux_y(:,j+1,:) - tfull_flux_y(:,j,:))./(dy.*rhow_num_y);
                elseif j == num_y
                    t(:,j,:) =squeeze(t(:,j,:)) - dtn.*(0.0 - squeeze(t_flux_y(:,num_y,:)))./(dy.*rhow_num_y);
                    tfull(:,j,:) =squeeze(tfull(:,j,:)) - dtn.*(0.0 - squeeze(tfull_flux_y(:,num_y,:)))./(dy.*rhow_num_y);
                else
                    sprintf('there is some error in the y indices advection')
                    
                end
            end
            
            
            
        end
        
        
        
        
        if calc_advection_tend %If I realy want to calculate the tendency itself
            % {
            for k = 1:num_z-1
                t_flux_z_tend(k,:,:) = - (t_flux_z(k+1,:,:) - t_flux_z(k,:,:))./(adz(k).*dz.*rho(k));
                tfull_flux_z_tend(k,:,:) = - (tfull_flux_z(k+1,:,:) - tfull_flux_z(k,:,:))./(adz(k).*dz.*rho(k));
            end
            k = num_z;
            t_flux_z_tend(k,:,:) = - (0.0 - t_flux_z(k,:,:))./(adz(k).*dz.*rho(k));
            tfull_flux_z_tend(k,:,:) = - (0.0 - tfull_flux_z(k,:,:))./(adz(k).*dz.*rho(k));
            
            for i = 1:num_x-1
                t_flux_x_tend(:,:,i) = - (t_flux_x(:,:,i+1) - t_flux_x(:,:,i))./(dx.*rhow_num_x);
                tfull_flux_x_tend(:,:,i) =  - (tfull_flux_x(:,:,i+1) - tfull_flux_x(:,:,i))./(dx.*rhow_num_x);
            end
            i = num_x;
            t_flux_x_tend(:,:,i) = - (t_flux_x(:,:,1) - t_flux_x(:,:,i))./(dx.*rhow_num_x);
            tfull_flux_x_tend(:,:,i) = - (tfull_flux_x(:,:,1) - tfull_flux_x(:,:,i))./(dx.*rhow_num_x);
            
            for j = 1:num_y-1
                t_flux_y_tend(:,j,:) = - squeeze(t_flux_y(:,j+1,:) - t_flux_y(:,j,:))./(dy.*rhow_num_y);
                tfull_flux_y_tend(:,j,:) = - squeeze(tfull_flux_y(:,j+1,:) - tfull_flux_y(:,j,:))./(dy.*rhow_num_y);
            end
            j = num_y;
            %             qt_flux_y_tend(:,j,:) = -(0.0 - squeeze(qt_flux_y(:,num_y,:)))./(dy.*rhow_num_y);
            t_flux_y_tend(:,j,:) = -(0.0 - squeeze(t_flux_y(:,num_y,:)))./(dy.*rhow_num_y);
            tfull_flux_y_tend(:,j,:) =-(0.0 - squeeze(tfull_flux_y(:,num_y,:)))./(dy.*rhow_num_y);
            
            if advect_fields
                t = t + (t_flux_x_tend + t_flux_y_tend + t_flux_z_tend)*dtn;
                tfull = tfull + (tfull_flux_x_tend + tfull_flux_y_tend + tfull_flux_z_tend)*dtn;
            end
            
        end
        
        
        
        %%02:07-try
        sprintf('try full advection for qp')
        
        %taken from fortran advect_scalar3D.f90 - when approximation for advection
        
        
        sprintf('full adv red')
        tic
        [f_out,qt_flux_x,qt_flux_y,qt_flux_z,qt_adv_tend] = advect_scalar3D_func_try_reduce_complexity(qt,u,v,w,rho,rhow,dx,dy,dz,dtn,num_x,num_y, num_z,adz);
        
        if advect_fields
            qt = f_out;
        end
        if do_show_times, toc, end
        sprintf('full adv ')
        tic
        [f_out,qp_flux_x,qp_flux_y,qp_flux_z,qp_adv_tend] = advect_scalar3D_func_try_reduce_complexity(qp,u,v,w,rho,rhow,dx,dy,dz,dtn,num_x,num_y, num_z,adz);
        if advect_fields
            qp = f_out;
        end
        if do_show_times, toc, end
        sprintf('the minimum qt, qp are:')
        min(min(min(qt)))
        min(min(min(qp)))
        
        qp = max(0,qp); 
        qt = max(0,qt);
        
        % 02:07a: vertical momentum advection
        [dudt_advect,dvdt_advect,dwdt_advect] = advect2_mom_z(num_x,num_y,num_z,dz,rho,rhow,adz,adzw,u,v,w);
        
        %%
        
        if do_show_times, sprintf('step 02:07: %f', toc), end
        tic
        %%
        %Step 02:07b precip_fall calculation
        
        % calculate precipitation flux
        %%precip_fall
        % from precip.f90
        
        
        sprintf('the precip fall func')
        tic
        [dqp_fall,t_fall_tend,precip]= ...
            precip_fall(qp,tabs,rho,rhow,rhofac,num_x,num_y,num_z,dz,adz,dtn,tprmin,a_pr,fac_fus,...
            fac_cond,crain,vrain,tgrmin,a_gr,qp_threshold,vgrau,cgrau,vsnow,csnow);%,do_fall_tend_qp,do_fall_tend_tfull);
        
        precip_energy = precip.*(fac_cond+(1.-omp_3d).*fac_fus); %I need to consider to recalc omp_3d due to changes in tabs
        
        if do_fall_tend_tfull
            tfull = tfull + t_fall_tend; % I found that it is better not to change
            qp = qp + dqp_fall;
        end
        dqp_fall = dqp_fall./dtn; %(all the tendencies are devided by dtn - and multiplied in the RF in SAMSON).
        t_fall_tend = t_fall_tend./dtn;
        
        if do_show_times, sprintf('step 02:07b: %f', toc), end
        tic
        %%
        % Step 02:08: surface fluxes
        % First find surface fluxes of t, tfull and qt following surface.f90
        umin = 1.0;
        cd=1.1e-3;
        wrk=(log(10/1.e-4)/log(z(1)/1.e-4))^2;
        fluxbt = zeros(num_y,num_x); %Yani added for parfor
        fluxbtfull = zeros(num_y,num_x);
        fluxbqt = zeros(num_y,num_x);
        for i=1:num_x
            for j=1:num_y
                
                if i<num_x
                    ic=i+1;
                else
                    ic=1;
                end
                
                if j<num_y
                    jc=j+1;
                else
                    jc=j;
                end
                
                if i == 1
                    ib = num_x;
                else
                    ib = i - 1;
                end
                
                if j == 1
                    jb = 1;
                else
                    jb = j -1;
                end
                
                ubot=0.5*(u(1,j,ic)+u(1,j,i));
                vbot=0.5*(v(1,jc,i)+v(1,j,i));
                windspeed=sqrt(ubot^2+vbot^2+umin^2);
                delt     = t(1,j,i)-gamaz(1) - sstxy(j,i);
                deltfull = tfull(1,j,i)-gamaz(1) - sstxy(j,i);
                ssq = qsatw(sstxy(j,i),pres(1));
                delqt   = qt(1,j,i)  - ssq;
                fluxbt(j,i) = -cd*windspeed*delt*wrk;
                fluxbtfull(j,i) = -cd*windspeed*deltfull*wrk;
                fluxbqt(j,i) = -cd*windspeed*delqt*wrk;
                windspeed_high(j,i) = windspeed;
                fluxbu_samson(j,i) = -rho(1)*(u(1,j,i))*cd*windspeed*wrk*rdz*rhow(1);
                fluxbv_samson(j,i) = -rho(1)*(v(1,j,i))*cd*windspeed*wrk*rdz*rhow(1);
                ubot=u(1,j,i);
                vbot=0.25*(v(1,j,i) + v(1,jc,i) + v(1,jc,ib) + v(1,j,ib));
                windspeed=sqrt(ubot^2+vbot^2+umin^2);
                fluxbu(j,i) = -rho(1)*(u(1,j,i))*cd*windspeed*wrk*rdz*rhow(1); %Note that we need to  multiply by things to get tendency (in fortran by - 1./(rho(k)*adz(k))) 
                vbot=v(1,j,i);
                ubot=0.25*(u(1,j,i) + u(1,jb,i) + u(1,jb,ic) + u(1,j,ic));
                windspeed=sqrt(ubot^2+vbot^2+umin^2);
                fluxbv(j,i) = -rho(1)*(v(1,j,i))*cd*windspeed*wrk*rdz*rhow(1);
            end
        end
        if do_show_times, sprintf('step 02:08: %f', toc), end
        tic
        % Step 02:09: vertical diffusion
        
        for i=1:num_x
            for j=1:num_y
                t_diff_flx_z(1,j,i)=fluxbt(j,i)*rdz*rhow(1);
                tfull_diff_flx_z(1,j,i)=fluxbtfull(j,i)*rdz*rhow(1);
                qt_diff_flx_z(1,j,i)=fluxbqt(j,i)*rdz*rhow(1);
                qp_diff_flx_z(1,j,i)= 0.0;
            end
        end
        
        sprintf("vert diff calc mat")
        tic
        rhoi = rhow(2:end-1)./adzw(2:end-1)';
        rhoi_rep = repmat(rhoi,[1,num_y,num_x]);
        tkh_z=rdz5*(tkz(1:num_z-1,:,:).*Pr1(1:num_z-1,:,:)+tkz(2:num_z,:,:).*Pr1(2:num_z,:,:));
        t_diff_flx_z(2:num_z,:,:)=-tkh_z.*(t(2:num_z,:,:)-t(1:num_z-1,:,:)).*rhoi/ravefactor;
        tfull_diff_flx_z(2:num_z,:,:)=-tkh_z.*(tfull(2:num_z,:,:)-tfull(1:num_z-1,:,:)).*rhoi/ravefactor;
        qt_diff_flx_z(2:num_z,:,:)=-tkh_z.*(qt(2:num_z,:,:)-qt(1:num_z-1,:,:)).*rhoi/ravefactor;
        qp_diff_flx_z(2:num_z,:,:)=-tkh_z.*(qp(2:num_z,:,:)-qp(1:num_z-1,:,:)).*rhoi/ravefactor;
        
        tkz_out(1:num_z-1,:,:) = tkz(1:num_z-1,:,:); 
        
        
        if do_show_times, toc, end
        
        t_diff_flx_z = t_diff_flx_z*dz;
        tfull_diff_flx_z = tfull_diff_flx_z*dz;
        qt_diff_flx_z = qt_diff_flx_z*dz;
        qp_diff_flx_z = qp_diff_flx_z*dz;
        
        if do_show_times, sprintf('step 02:09 (before vertical diffusive tendencies: %f', toc), end
        irhoadz = zeros([num_z,1]);
        for k=1:num_z
            irhoadz(k) = 1./(rho(k)*adz(k));
        end
        %Calculate  diffusive vertical tendencies:
        if calc_diffusive_tend
            for k=1:num_z
                irhoadz_dz = irhoadz(k)./dz;
                for j=1:num_y
                    for i=1:num_x
                        if k == num_z
                            t_diff_z_tend(k,j,i)=-(0-t_diff_flx_z(k,j,i)).*irhoadz_dz;
                            tfull_diff_z_tend(k,j,i)=-(0-tfull_diff_flx_z(k,j,i)).*irhoadz_dz;
                            qt_diff_z_tend(k,j,i)=-(0-qt_diff_flx_z(k,j,i)).*irhoadz_dz;
                            qp_diff_z_tend(k,j,i)=-(0-qp_diff_flx_z(k,j,i)).*irhoadz_dz;
                        else
                            kt=k + 1;
                            t_diff_z_tend(k,j,i)=-(t_diff_flx_z(kt,j,i)-t_diff_flx_z(k,j,i)).*irhoadz_dz;
                            tfull_diff_z_tend(k,j,i)=-(tfull_diff_flx_z(kt,j,i)-tfull_diff_flx_z(k,j,i)).*irhoadz_dz;
                            qt_diff_z_tend(k,j,i)=-(qt_diff_flx_z(kt,j,i)-qt_diff_flx_z(k,j,i)).*irhoadz_dz;
                            qp_diff_z_tend(k,j,i)=-(qp_diff_flx_z(kt,j,i)-qp_diff_flx_z(k,j,i)).*irhoadz_dz;
                        end
                    end
                end
            end
        end
        
        
        if do_show_times, sprintf('step 02:09: %f', toc), end
        tic
        %%
        % Step 02:10: horizontal diffusion
        
        
        if do_horizontal_diffusion_at_all
            %Horizontal diffusion
            rdx2=1./(dx*dx);
            rdy2=1./(dy*dy);
            tkh_xy = tkz;
            rdx5=0.5*rdx2;
            rdy5=0.5*rdy2;
            
            % x diffusion
            for k=1:num_z
                for j=1:num_y
                    for i=1:num_x
                        ic = i + 1;
                        if i == num_x
                            ic = 1;
                        end
                        tkx=rdx5*(tkh_xy(k,j,i)*Pr1(k,j,i)+tkh_xy(k,j,ic)*Pr1(k,j,ic));
                        t_diff_flx_x(k,j,i)=-tkx*(t(k,j,ic)-t(k,j,i))*ravefactor;
                        tfull_diff_flx_x(k,j,i)=-tkx*(tfull(k,j,ic)-tfull(k,j,i))*ravefactor;
                        qt_diff_flx_x(k,j,i)=-tkx*(qt(k,j,ic)-qt(k,j,i))*ravefactor;
                        qp_diff_flx_x(k,j,i)=-tkx*(qp(k,j,ic)-qp(k,j,i))*ravefactor;
                    end
                end
                
            end
            
            if calc_diffusive_tend
                for k=1:num_z
                    for j=1:num_y
                        for i=1:num_x
                            if i == 1
                                ib = num_x;
                            else
                                ib=i-1;
                            end
                            t_diff_x_tend(k,j,i)=-(t_diff_flx_x(k,j,i)-t_diff_flx_x(k,j,ib));
                            tfull_diff_x_tend(k,j,i)=-(tfull_diff_flx_x(k,j,i)-tfull_diff_flx_x(k,j,ib));
                            qt_diff_x_tend(k,j,i)=-(qt_diff_flx_x(k,j,i)-qt_diff_flx_x(k,j,ib));
                            qp_diff_x_tend(k,j,i)=-(qp_diff_flx_x(k,j,i)-qp_diff_flx_x(k,j,ib));
                        end
                    end
                end
            end
            
            % y diffusion
            for k=1:num_z
                for j=1:num_y
                    for i=1:num_x
                        jc = j + 1;
                        if j == num_y
                            t_diff_flx_y(k,j,i)=0;
                            tfull_diff_flx_y(k,j,i)=0;
                            qt_diff_flx_y(k,j,i)=0;
                            qp_diff_flx_y(k,j,i)=0;
                        else
                            tkx=rdy5*(tkh_xy(k,j,i)*Pr1(k,j,i)+tkh_xy(k,jc,i)*Pr1(k,jc,i));
                            t_diff_flx_y(k,j,i)=-tkx*(t(k,jc,i)-t(k,j,i))*ravefactor;
                            tfull_diff_flx_y(k,j,i)=-tkx*(tfull(k,jc,i)-tfull(k,j,i))*ravefactor;
                            qt_diff_flx_y(k,j,i)=-tkx*(qt(k,jc,i)-qt(k,j,i))*ravefactor;
                            qp_diff_flx_y(k,j,i)=-tkx*(qp(k,jc,i)-qp(k,j,i))*ravefactor;
                        end
                        
                    end
                end
                
            end
            
            if do_show_times, sprintf('step 02:10 (before horizontal diffusive tendencies): %f', toc),end
            
            if calc_diffusive_tend
                
                for k=1:num_z
                    for j=1:num_y
                        for i=1:num_x
                            if j == 1
                                t_diff_y_tend(k,j,i)=-(t_diff_flx_y(k,j,i)-0);
                                tfull_diff_y_tend(k,j,i)=-(tfull_diff_flx_y(k,j,i)-0);
                                qt_diff_y_tend(k,j,i)=-(qt_diff_flx_y(k,j,i)-0);
                                qp_diff_y_tend(k,j,i)=-(qp_diff_flx_y(k,j,i)-0);
                            else
                                jb=j-1;
                                t_diff_y_tend(k,j,i)=-(t_diff_flx_y(k,j,i)-t_diff_flx_y(k,jb,i));
                                tfull_diff_y_tend(k,j,i)=-(tfull_diff_flx_y(k,j,i)-tfull_diff_flx_y(k,jb,i));
                                qt_diff_y_tend(k,j,i)=-(qt_diff_flx_y(k,j,i)-qt_diff_flx_y(k,jb,i));
                                qp_diff_y_tend(k,j,i)=-(qp_diff_flx_y(k,j,i)-qp_diff_flx_y(k,jb,i));
                            end
                        end
                    end
                end
            end
        end
        
        if (diffuse_fields==1 && calc_diffusive_tend==1)
            t = t + (t_diff_x_tend + t_diff_y_tend + t_diff_z_tend).*dtn;
            tfull = tfull + (tfull_diff_x_tend + tfull_diff_y_tend + tfull_diff_z_tend).*dtn;
            qt = qt + (qt_diff_x_tend + qt_diff_y_tend + qt_diff_z_tend).*dtn;
            qp = qp + (qp_diff_x_tend + qp_diff_y_tend + qp_diff_z_tend).*dtn;
        end
        
        if (diffuse_fields==1 && calc_diffusive_tend==0) %Diffuse fields without calculating seperately the tendencies
            
            for k=1:num_z
                irhoadz_dz = irhoadz(k)./dz;
                for j=1:num_y
                    for i=1:num_x
                        if k == num_z
                            t(k,j,i)=t(k,j,i) - dtn.*(0-t_diff_flx_z(k,j,i)).*irhoadz_dz;
                            tfull(k,j,i)=tfull(k,j,i) - dtn.*(0-tfull_diff_flx_z(k,j,i)).*irhoadz_dz;
                            qt(k,j,i)=qt(k,j,i) - dtn.*(0-qt_diff_flx_z(k,j,i)).*irhoadz_dz;
                            qp(k,j,i)=qp(k,j,i) - dtn.*(0-qp_diff_flx_z(k,j,i)).*irhoadz_dz;
                        else
                            kt=k + 1;
                            t(k,j,i)=t(k,j,i) - dtn.*(t_diff_flx_z(kt,j,i)-t_diff_flx_z(k,j,i)).*irhoadz_dz;
                            tfull(k,j,i)=tfull(k,j,i) - dtn.*(tfull_diff_flx_z(kt,j,i)-tfull_diff_flx_z(k,j,i)).*irhoadz_dz;
                            qt(k,j,i)=qt(k,j,i) - dtn.*(qt_diff_flx_z(kt,j,i)-qt_diff_flx_z(k,j,i)).*irhoadz_dz;
                            qp(k,j,i)=qp(k,j,i) - dtn.*(qp_diff_flx_z(kt,j,i)-qp_diff_flx_z(k,j,i)).*irhoadz_dz;
                        end
                    end
                end
            end
            
            if do_horizontal_diffusion_at_all
                for k=1:num_z
                    for j=1:num_y
                        for i=1:num_x
                            if i == 1
                                ib = num_x;
                            else
                                ib=i-1;
                            end
                            t(k,j,i)=t(k,j,i) - dtn.*(t_diff_flx_x(k,j,i)-t_diff_flx_x(k,j,ib));
                            tfull(k,j,i)=tfull(k,j,i) - dtn.*(tfull_diff_flx_x(k,j,i)-tfull_diff_flx_x(k,j,ib));
                            qt(k,j,i)=qt(k,j,i) - dtn.*(qt_diff_flx_x(k,j,i)-qt_diff_flx_x(k,j,ib));
                            qp(k,j,i)=qp(k,j,i) - dtn.*(qp_diff_flx_x(k,j,i)-qp_diff_flx_x(k,j,ib));
                        end
                    end
                end
                
                
                
                for k=1:num_z
                    for j=1:num_y
                        for i=1:num_x
                            if j == 1
                                t(k,j,i)=t(k,j,i) - dtn.*(t_diff_flx_y(k,j,i)-0);
                                tfull(k,j,i)=tfull(k,j,i) - dtn.*(tfull_diff_flx_y(k,j,i)-0);
                                qt(k,j,i)=qt(k,j,i) - dtn.*(qt_diff_flx_y(k,j,i)-0);
                                qp(k,j,i)=qp(k,j,i) - dtn.*(qp_diff_flx_y(k,j,i)-0);
                            else
                                jb=j-1;
                                t(k,j,i)=t(k,j,i) - dtn.*(t_diff_flx_y(k,j,i)-t_diff_flx_y(k,jb,i));
                                tfull(k,j,i)=tfull(k,j,i) - dtn.*(tfull_diff_flx_y(k,j,i)-tfull_diff_flx_y(k,jb,i));
                                qt(k,j,i)=qt(k,j,i) - dtn.*(qt_diff_flx_y(k,j,i)-qt_diff_flx_y(k,jb,i));
                                qp(k,j,i)=qp(k,j,i) - dtn.*(qp_diff_flx_y(k,j,i)-qp_diff_flx_y(k,jb,i));
                            end
                        end
                    end
                end
            end
            
        end
        t_diff_flx_x = t_diff_flx_x.*dx;
        tfull_diff_flx_x = tfull_diff_flx_x.*dx;
        qt_diff_flx_x = qt_diff_flx_x.*dx;
        qp_diff_flx_x = qp_diff_flx_x.*dx;
        
        t_diff_flx_y = t_diff_flx_y.*dy;
        tfull_diff_flx_y = tfull_diff_flx_y.*dy;
        qt_diff_flx_y = qt_diff_flx_y.*dy;
        qp_diff_flx_y = qp_diff_flx_y.*dy;
        
        if do_show_times, sprintf('step 02:10: %f', toc), end
        tic
        
        %%
        % Step 02:11: do cloud (for recalculation of qn, and sedimentation)
        %%Do cloud before microphysics! cloud.f90
        an = 1./(tbgmax-tbgmin) ;
        bn = tbgmin * an;
        ap = 1./(tprmax-tprmin) ;
        bp = tprmin * ap;
        fac1 = fac_cond+(1+bp)*fac_fus;
        fac2 = fac_fus*ap;
        ag = 1./(tgrmax-tgrmin);
        
        kmax=0;
        kmin=num_z+1;
        
        for i=1:num_x
            for j=1:num_y
                for k=1:num_z
                    qn0 = qn(k,j,i);
                    qt(k,j,i)=max(0.,qt(k,j,i));
                    
                    %             ! Initail guess for temperature assuming no cloud water/ice:
                    tabs(k,j,i) = tfull(k,j,i)-gamaz(k); % Yani - modified to tfull.
                    tabs1=(tabs(k,j,i)+fac1*qp(k,j,i))/(1.+fac2*qp(k,j,i));
                    
                    %             ! Warm cloud:
                    
                    if(tabs1 >= tbgmax)
                        
                        tabs1=tabs(k,j,i)+fac_cond*qp(k,j,i);
                        qsat = qsatw(tabs1,pres(k));
                        
                        %                 ! Ice cloud:
                        
                    elseif(tabs1 <= tbgmin)
                        
                        tabs1=tabs(k,j,i)+fac_sub*qp(k,j,i);
                        qsat = qsati(tabs1,pres(k));
                        
                        %                 ! Mixed-phase cloud:
                        
                    else
                        
                        om = an*tabs1-bn;
                        qsat = om*qsatw(tabs1,pres(k))+(1.-om)*qsati(tabs1,pres(k));
                        
                    end
                    
                    
                    if(qt(k,j,i) > qsat)
                        
                        niter=0;
                        dtabs = 100.;
                        while(abs(dtabs)>0.01 && niter < 10)
                            if(tabs1>=tbgmax)
                                om=1.;
                                lstarn=fac_cond;
                                dlstarn=0.;
                                qsat=qsatw(tabs1,pres(k));
                                dqsat=dtqsatw(tabs1,pres(k));
                            elseif(tabs1<=tbgmin)
                                om=0.;
                                lstarn=fac_sub;
                                dlstarn=0.;
                                qsat=qsati(tabs1,pres(k));
                                dqsat=dtqsati(tabs1,pres(k));
                            else
                                om=an*tabs1-bn;
                                lstarn=fac_cond+(1.-om)*fac_fus;
                                dlstarn=an;
                                qsat=om*qsatw(tabs1,pres(k))+(1.-om)*qsati(tabs1,pres(k));
                                dqsat=om*dtqsatw(tabs1,pres(k))+(1.-om)*dtqsati(tabs1,pres(k));
                            end
                            if(tabs1>=tprmax)
                                omp=1.;
                                lstarp=fac_cond;
                                dlstarp=0.;
                            elseif(tabs1<=tprmin)
                                omp=0.;
                                lstarp=fac_sub;
                                dlstarp=0.;
                            else
                                omp=ap*tabs1-bp;
                                lstarp=fac_cond+(1.-omp)*fac_fus;
                                dlstarp=ap;
                            end
                            fff = tabs(k,j,i)-tabs1+lstarn*(qt(k,j,i)-qsat)+lstarp*qp(k,j,i);
                            dfff=dlstarn*(qt(k,j,i)-qsat)+dlstarp*qp(k,j,i)-lstarn*dqsat-1.;
                            dtabs=-fff/dfff;
                            niter=niter+1;
                            tabs1=tabs1+dtabs;
                        end
                        qsat = qsat + dqsat * dtabs;
                        qn(k,j,i) = max(0.,qt(k,j,i)-qsat);
                    else
                        qn(k,j,i) = 0.;
                    end
                    tabs(k,j,i) = tabs1;
                    qp(k,j,i) = max(0.,qp(k,j,i)); %! just in case
                    
                    if(qn(k,j,i)>qp_threshold)
                        kmin = min(kmin,k);
                        kmax = max(kmax,k);
                    end
                end
            end
        end
        
        if do_show_times, sprintf('step 02:11 before sedimetation): %f', toc), end
        
        if do_sedimentation
            % Sedimentation of ice and water:
            qifall = zeros(num_z,1);
            tlatqi = zeros(num_z,1);
            for k = 1:num_z
                qifall(k) = 0.;
                tlatqi(k) = 0.;
            end
            
            coef_cl = 1.19e8*(3./(4.*3.1415*1000.*Nc0*1.e6))^(2./3.)*exp(5.*log(1.5)^2);
            
            for k = max(1,kmin-1):kmax
                %    ! Set up indices for x-y planes above and below current plane.
                kc = min(num_z,k+1);
                kb = max(1,k-1);
                for j = 1:num_y
                    for i = 1:num_x
                        coef = dtn/(0.5*(adz(kb)+adz(k))*dz);
                        %
                        %          ! Compute cloud ice density in this cell and the ones above/below.
                        %          ! Since cloud ice is falling, the above cell is u (upwind),
                        %          ! this cell is c (center) and the one below is d (downwind).
                        omnu = max(0.,min(1.,(tabs(kc,j,i)-tbgmin)*a_bg));
                        omnc = max(0.,min(1.,(tabs(k,j,i) -tbgmin)*a_bg));
                        omnd = max(0.,min(1.,(tabs(kb,j,i)-tbgmin)*a_bg));
                        
                        qiu = rho(kc)*qn(kc,j,i)*(1.-omnu);
                        qic = rho(k) *qn(k,j,i) *(1.-omnc);
                        qid = rho(kb)*qn(kb,j,i)*(1.-omnd);
                        
                        %          ! Ice sedimentation velocity depends on ice content. The fiting is
                        %          ! based on the data by Heymsfield (JAS,2003). -Marat
                        %          ! 0.1 m/s low bound was suggested by Chris Bretherton
                        vt_ice = max(0.1,0.5*log10(qic+1.e-12)+3.);
                        
                        %          ! Use MC flux limiter in computation of flux correction.
                        %          ! (MC = monotonized centered difference).
                        if (qic==qid)
                            tmp_phi = 0.;
                        else
                            tmp_theta = (qiu-qic)/(qic-qid);
                            tmp_phi = max(0.,min([0.5*(1.+tmp_theta),2.,2.*tmp_theta]));
                        end
                        
                        %          ! Compute limited flux.
                        %          ! Since falling cloud ice is a 1D advection problem, this
                        %          ! flux-limited advection scheme is monotonic.
                        fluxi = -vt_ice*(qic - 0.5*(1.-coef*vt_ice)*tmp_phi*(qic-qid));
                        
                        doclouddropsed = 0;%This is our config- Yani
                        if(doclouddropsed)
                            %             ! Compute cloud water density in this cell and the ones above/below
                            %             ! Since cloud water is falling, the above cell is u (upwind),
                            %             ! this cell is c (center) and the one below is d (downwind).
                            qiu = rho(kc)*qn(kc,j,i)*omnu;
                            qic = rho(k) *qn(k,j,i) *omnc;
                            qid = rho(kb)*qn(kb,j,i)*omnd;
                            
                            vt_cl = coef_cl*(qic+1.e-12)^(2./3.);
                            
                            %             ! Use MC flux limiter in computation of flux correction.
                            %             ! (MC = monotonized centered difference).
                            if (qic==qid)
                                tmp_phi = 0.;
                            else
                                tmp_theta = (qiu-qic)/(qic-qid);
                                tmp_phi = max(0.,min(0.5*(1.+tmp_theta),2.,2.*tmp_theta));
                            end
                            
                            %             ! Compute limited flux.
                            %             ! Since falling cloud water is a 1D advection problem, this
                            %             ! flux-limited advection scheme is monotonic.
                            fluxc = -vt_cl*(qic - 0.5*(1.-coef*vt_cl)*tmp_phi*(qic-qid));
                        else
                            fluxc = 0.;
                        end
                        fz(k,j,i) = fluxi + fluxc;
                        fzt(k,j,i) = -(fac_cond+fac_fus)*fluxi - fac_cond*fluxc;
                    end
                end
            end
            
            for k=max(1,kmin-2):kmax
                for j=1:num_y
                    for i=1:num_x
                        coef=dtn/(dz*adz(k)*rho(k));
                        
                        %          ! The cloud ice increment is the difference of the fluxes.
                        dqi=coef*(fz(k,j,i)-fz(k+1,j,i));
                        %          ! Add this increment to both non-precipitating and total water.
                        qn(k,j,i) = qn(k,j,i) + dqi;
                        qt(k,j,i)  = qt(k,j,i)  + dqi;
                        
                        cloud_qt_tend(k,j,i) =  dqi;
                        %          ! Include this effect in the total moisture budget.
                        %          ! The latent heat flux induced by the falling cloud ice enters
                        %          ! the liquid-ice static energy budget in the same way as the
                        %          ! precipitation.  Note: use latent heat of sublimation.
                        lat_heat  = coef*(fzt(k,j,i)-fzt(k+1,j,i));
                        cloud_lat_heat(k,j,i) = lat_heat;
                        %          ! Add divergence of latent heat flux to liquid-ice static energy.
                        if add_cloud_tfull_tend
                            tfull(k,j,i)  = tfull(k,j,i)  + lat_heat; % Yani Need to think when need to change also t not full
                        end
                    end
                end
            end
            
            cloud_qt_tend = cloud_qt_tend./dtn; %all tendencies are divided by dtn (multiplied in RF of fortran).
            cloud_lat_heat = cloud_lat_heat./dtn;
            
        end
        if do_show_times, sprintf('step 02:11: %f', toc), end
        tic
        
        %%
        % Step 02:12: calculate microphysical tendency of qp
        % calculate microphysical tendency of qp
        %%precip_proc
        
        % from precip_proc.f90
        
        powr1 = (3 + b_rain) / 4.;
        powr2 = (5 + b_rain) / 8.;
        pows1 = (3 + b_snow) / 4.;
        pows2 = (5 + b_snow) / 8.;
        powg1 = (3 + b_grau) / 4.;
        powg2 = (5 + b_grau) / 8.;
        
        for i=1:num_x
            for j=1:num_y
                for k=1:num_z
                    
                    if (qn(k,j,i)+qp(k,j,i)>0.)
                        
                        omn = max(0.,min(1.,(tabs(k,j,i)-tbgmin)*a_bg));
                        omp = max(0.,min(1.,(tabs(k,j,i)-tprmin)*a_pr));
                        omg = max(0.,min(1.,(tabs(k,j,i)-tgrmin)*a_gr));
                        qrr = qp(k,j,i) * omp;
                        qss = qp(k,j,i) * (1.-omp)*(1.-omg); %Eq. A11 in the paper of MARAT. ;
                        qgg = qp(k,j,i) * (1.-omp)*omg;
                        
                        if (qn(k,j,i)>0.) %-------     Autoconversion/accretion
                            
                            qcc = qn(k,j,i) * omn;
                            qii = qn(k,j,i) * (1.-omn);
                            
                            if (qcc > qcw0)
                                autor = alphaelq;
                            else
                                autor = 0.;
                            end
                            
                            if (qii > qci0)
                                autos = betaelq*coefice(k);
                            else
                                autos = 0.;
                            end
                            accrr = accrrc(k) * qrr^powr1;
                            tmp = qss^pows1;
                            accrcs = accrsc(k) * tmp;
                            accris = accrsi(k) * tmp;
                            tmp = qgg^powg1;
                            accrcg = accrgc(k) * tmp;
                            accrig = accrgi(k) * tmp;
                            qcc = (qcc+dtn*autor*qcw0)/(1.+dtn*(accrr+accrcs+accrcg+autor));
                            qii = (qii+dtn*autos*qci0)/(1.+dtn*(accris+accrig+autos));
                            
                            
                            dqp(k,j,i) = dtn *(accrr*qcc + autor*(qcc-qcw0)+ ...
                                (accris+accrig)*qii + (accrcs+accrcg)*qcc + autos*(qii-qci0));
                            
                            dqp(k,j,i) = min(dqp(k,j,i),qn(k,j,i));
                            
                        elseif(qp(k,j,i)>qp_threshold && qn(k,j,i)==0.)  % evaporation
                            % I think that there is a missing condition!
                            % if(tabs(i,j,k).gt.tmin_evap) then !kzm limit evaporation to temperatures
                            % above tmin_evap tmin_evap = 0 I think in our case
                            
                            qsat = omn*qsatw(tabs(k,j,i),pres(k))+ ...
                                (1.-omn)*qsati(tabs(k,j,i),pres(k));
                            dqp(k,j,i) = dtn * (evapr1(k) * sqrt(qrr) + ...
                                evapr2(k) * qrr^powr2 + ...
                                evaps1(k) * sqrt(qss) + ...
                                evaps2(k) * qss^pows2 + ...
                                evapg1(k) * sqrt(qgg) + ...
                                evapg2(k) * qgg^powg2)* ...
                                (qt(k,j,i) /qsat-1.);
                            dqp(k,j,i) = max(-0.5*qp(k,j,i),dqp(k,j,i));
                            
                        else % complete evaporation
                            
                            dqp(k,j,i) = -qp(k,j,i);
                            
                        end
                        
                    end
                    qp(k,j,i) = qp(k,j,i)  + dqp(k,j,i);
                    qp(k,j,i) = max(0,qp(k,j,i));% Avoid negative values for
                    %             precip - should use only if I really modify qp here... which
                    %             probably I should! For some reason the advection of qp give
                    %             negative qp values - shouldn't happen - to check !
                    
                end  % for
            end % for
        end % for
        % get tendency in kg/kg/s
        dqp = dqp/dtn;
        omp_3d  = max(0.,min(1.,(tabs-tprmin)*a_pr));
        fac_dqp_tfull = (fac_cond+(1.-omp_3d)*fac_fus); 
        dqp_t_tend = dqp .*fac_dqp_tfull;
        
        if do_show_times, sprintf('step 02:12: %f', toc), end
        tic
        
        %% Step 03: Main loop over coarse grained quanteties:
        
        
        %03:01: define resolutions for calculations
        for res = resolutions
            if do_show_times, sprintf('step 03:01, with resolution %f):',res), end
            
            if is_test
                multiple_space = 1;
            else
                multiple_space = res; % average over how many grid cells
            end
            
            %
            % if is_test
            coarse_test = max(multiple_space,fac_redefine_dx_dy);
            ssty_coarse = zeros(length(ssty)./coarse_test,1);
            if ssty_from_txt 
                filename_profile = sprintf('ssty_%i.txt',res);
                profile = importdata(filename_profile);
                ssty_coarse = profile(:);
            else
                for aa = 1:length(ssty)./coarse_test
                    ssty_coarse(aa) = (mean(ssty(((aa-1).*coarse_test + 1):coarse_test.*aa),1)) ;
                end
            end
            num_blocks_x = num_x/multiple_space;
            num_blocks_y = num_y/multiple_space;
            
            %03:02: Allocate coarse grain values and resolved
            
            x_coarse = zeros(num_blocks_x,1);
            y_coarse = zeros(num_blocks_y,1);
            sstxy_coarse = repmat(ssty_coarse, [1, length(x_coarse)]); %TO CHECK
            
            
            u_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            v_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            w_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            tabs_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            t_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            tfull_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            qv_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            qn_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            qsat_resolved_init = zeros(num_z,num_blocks_y,num_blocks_x);
            qp_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            Qrad_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            dqp_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            dqp_t_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            dqp_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            tfull_flux_x_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            tfull_flux_y_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            tfull_flux_z_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            t_flux_x_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            t_flux_y_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            t_flux_z_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            qt_flux_x_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            qt_flux_y_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            qt_flux_z_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            qp_flux_x_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            qp_flux_y_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            qp_flux_z_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            
            if calc_advection_tend
                tfull_flux_x_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                tfull_flux_y_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                tfull_flux_z_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                t_flux_x_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                t_flux_y_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                t_flux_z_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            end
            tfull_flux_x_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            tfull_flux_y_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            tfull_flux_z_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            t_flux_x_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            t_flux_y_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            t_flux_z_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            
            if calc_advection_tend
                tfull_flux_x_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                tfull_flux_y_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                tfull_flux_z_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                t_flux_x_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                t_flux_y_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                t_flux_z_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            end
            
            
            t_diff_flx_x_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            t_diff_flx_y_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            t_diff_flx_z_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            
            tfull_diff_flx_x_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            tfull_diff_flx_y_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            tfull_diff_flx_z_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            
            qt_diff_flx_x_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            qt_diff_flx_y_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            qt_diff_flx_z_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            
            qp_diff_flx_x_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            qp_diff_flx_y_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            qp_diff_flx_z_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            
            if calc_diffusive_tend
                t_diff_x_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                t_diff_y_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                t_diff_z_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                
                tfull_diff_x_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                tfull_diff_y_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                tfull_diff_z_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                
                
                qt_diff_x_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                qt_diff_y_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                qt_diff_z_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                
                qp_diff_x_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                qp_diff_y_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
                qp_diff_z_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            end
            
            t_diff_flx_x_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            t_diff_flx_y_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            t_diff_flx_z_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            
            tfull_diff_flx_x_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            tfull_diff_flx_y_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            tfull_diff_flx_z_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            
            qt_diff_flx_x_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            qt_diff_flx_y_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            qt_diff_flx_z_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            
            qp_diff_flx_x_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            qp_diff_flx_y_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            qp_diff_flx_z_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            
            if calc_diffusive_tend
                t_diff_x_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                t_diff_y_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                t_diff_z_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                
                tfull_diff_x_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                tfull_diff_y_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                tfull_diff_z_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                
                qt_diff_x_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                qt_diff_y_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                qt_diff_z_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                
                qp_diff_x_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                qp_diff_y_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                qp_diff_z_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            end
            
            
            dqp_fall_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            t_fall_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            
            precip_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            precip_energy_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            
            cloud_lat_heat_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            cloud_qt_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            
            cloud_lat_heat_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            cloud_qt_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            
            fzt_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            fz_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            
            fzt_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            fz_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            
            Pr1_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            Pr1_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            
            tkz_out_coarse= zeros(num_z,num_blocks_y,num_blocks_x);
            tkz_out_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            
            windspeed_coarse= zeros(num_blocks_y,num_blocks_x);
            windspeed_resolved= zeros(num_blocks_y,num_blocks_x);
            fluxbu_coarse= zeros(num_blocks_y,num_blocks_x);
            fluxbv_coarse= zeros(num_blocks_y,num_blocks_x);
            fluxbu_resolved= zeros(num_blocks_y,num_blocks_x);
            fluxbv_resolved= zeros(num_blocks_y,num_blocks_x);
            
            fluxbu_samson_coarse= zeros(num_blocks_y,num_blocks_x);
            fluxbv_samson_coarse= zeros(num_blocks_y,num_blocks_x);
            fluxbu_samson_resolved= zeros(num_blocks_y,num_blocks_x);
            fluxbv_samson_resolved= zeros(num_blocks_y,num_blocks_x);
            
            
            dudt_advect_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            dvdt_advect_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            dwdt_advect_coarse = zeros(num_z,num_blocks_y,num_blocks_x);

            dudt_advect_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            dvdt_advect_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
            dwdt_advect_resolved = zeros(num_z,num_blocks_y,num_blocks_x);

%             tabs_coarse_end = zeros(num_z,num_blocks_y,num_blocks_x);
            
            
            % Step 03:02b get coarse grained variables before they change
            % check block size
            %read again vars that changed so I could coarse grain them (Yani: why not coarse in the beginning ? )            
            qv = my_ncread(filename, 'Q')/1000.0; % water vapor (kg/kg)
            qn = my_ncread(filename, 'QN')/1000.0; % non precip cond (water+ice) (kg/kg)
            qp = my_ncread(filename, 'QP')/1000.0; % precip (kg/kg)
            Qrad = my_ncread(filename, 'QRAD'); % rad heating rate (K/day)
            tabs = my_ncread(filename, 'TABS'); % absolute temperature (K) - Note that this is not an 
            qt = qv+qn; % total non-precipitating water (referred to as simply q in SAM)
            
            % calculate liquid/ice static energy
            gamaz = zeros(size(tabs));
            for k=1:num_z
                gamaz(k,:,:)=ggr/cp*z(k);
            end
            
            % liquid/ice water static energy h_L divided by cp
            % but here set qp to zero in SAM with rf
            omn_3d  = max(0.,min(1.,(tabs-tbgmin)*a_bg));
            t = tabs + gamaz - (fac_cond+(1.-omn_3d).*fac_fus).*qn;
            
            % full t including precipitating condensates
            omp_3d  = max(0.,min(1.,(tabs-tprmin)*a_pr));
            tfull = tabs + gamaz - (fac_cond+(1.-omn_3d).*fac_fus).*qn - ...
                (fac_cond+(1.-omp_3d).*fac_fus).*qp;
            
            if mod(num_x, multiple_space)~=0
                error('number of x values not multiple of spatial block size')
            end
            if mod(num_y, multiple_space)~=0
                error('number of y values not multiple of spatial block size')
            end
            
            % change to mean over spatial blocks
            multiple_space_weight1 = (multiple_space-1)./multiple_space ;
            multiple_space_weight2 = (1)./multiple_space ;
            for i=1:num_blocks_x
                i_indices = [(i-1)*multiple_space+1:i*multiple_space];
                x_coarse(i) = mean(x(i_indices));
                        %winds are defined on different grids than
                        %thermodynamic vars (C-grid), so need to coarse
                        %grain them differnetly. 
                if i==1
                    i_indices_x = [1:multiple_space/2,size(u,3)-multiple_space./2+2:size(u,3)]; 
                    i_indices_x2 = [multiple_space/2+1,size(u,3)-multiple_space./2 + 1];
                else
                    i_indices_x = [(i-1)*multiple_space+1 - multiple_space/2 + 1:(i*multiple_space-multiple_space/2)]; 
                    i_indices_x2 = [(i-1)*multiple_space+1 - multiple_space/2,i*multiple_space - multiple_space/2 + 1];
                end

                for j=1:num_blocks_y
                    j_indices = [(j-1)*multiple_space+1:j*multiple_space];
                    y_coarse(j) = mean(y(j_indices));
                    
                    if j==1
                        j_indices_y = 1;
                        j_indices_y2 =1; 
                    else
                        j_indices_y = [(j-1)*multiple_space+1 - multiple_space/2 + 1:(j*multiple_space - multiple_space/2)]; 
                        j_indices_y2 = [(j-1)*multiple_space+1 - multiple_space/2,j*multiple_space - multiple_space/2 + 1];
                    end
                    
                    if is_test ==1
                        %Surface wind fluxes (to correct momentum).
                        select = fluxbu(j_indices,i_indices);
                        fluxbu_coarse(j,i) = mean(select(:));
                        
                        select = fluxbv(j_indices,i_indices);
                        fluxbv_coarse(j,i) = mean(select(:));
                        
                        select = fluxbu_samson(j_indices,i_indices);
                        fluxbu_samson_coarse(j,i) = mean(select(:));
                        
                        select = fluxbv_samson(j_indices,i_indices);
                        fluxbv_samson_coarse(j,i) = mean(select(:));
                    else
                        select = fluxbu(j_indices,i_indices_x); %Modified x indices by shifting it
                        select2 =  fluxbu(j_indices,i_indices_x2);
                        fluxbu_coarse(j,i) = multiple_space_weight1*mean(select(:)) + multiple_space_weight2*mean(select2(:));
                        
                        select = fluxbv(j_indices_y,i_indices);
                        select2 =  fluxbv(j_indices_y2,i_indices);
                        fluxbv_coarse(j,i) = multiple_space_weight1*mean(select(:)) + multiple_space_weight2*mean(select2(:));
                        
                        select = fluxbu_samson(j_indices,i_indices_x); %Modified x indices by shifting it
                        select2 =  fluxbu_samson(j_indices,i_indices_x2);
                        fluxbu_samson_coarse(j,i) = multiple_space_weight1*mean(select(:)) + multiple_space_weight2*mean(select2(:));
                        
                        select = fluxbv_samson(j_indices_y,i_indices);
                        select2 =  fluxbv_samson(j_indices_y2,i_indices);
                        fluxbv_samson_coarse(j,i) = multiple_space_weight1*mean(select(:)) + multiple_space_weight2*mean(select2(:));
                    end

                    for k=1:num_z
                        if is_test ==1
                        select = u(k,j_indices,i_indices);
                        u_coarse(k,j,i) = mean(select(:));
                        
                        select = v(k,j_indices,i_indices);
                        v_coarse(k,j,i) = mean(select(:));
                        
                        
                        %Momentum vertical advection 
                        select = dudt_advect(k,j_indices,i_indices);
                        dudt_advect_coarse(k,j,i) = mean(select(:));
                        
                        select = dvdt_advect(k,j_indices,i_indices);
                        dvdt_advect_coarse(k,j,i) = mean(select(:));
                        

                        else
                            select = u(k,j_indices,i_indices_x); %Modified x indices by shifting it 
                            select2 =  u(k,j_indices,i_indices_x2);
                            u_coarse(k,j,i) = multiple_space_weight1*mean(select(:)) + multiple_space_weight2*mean(select2(:));

                            select = v(k,j_indices_y,i_indices);
                            select2 =  v(k,j_indices_y2,i_indices);
                            v_coarse(k,j,i) = multiple_space_weight1*mean(select(:)) + multiple_space_weight2*mean(select2(:));
                            
                            select = dudt_advect(k,j_indices,i_indices_x); %Modified x indices by shifting it 
                            select2 =  dudt_advect(k,j_indices,i_indices_x2);
                            dudt_advect_coarse(k,j,i) = multiple_space_weight1*mean(select(:)) + multiple_space_weight2*mean(select2(:));

                            select = dvdt_advect(k,j_indices_y,i_indices);
                            select2 =  dvdt_advect(k,j_indices_y2,i_indices);
                            dvdt_advect_coarse(k,j,i) = multiple_space_weight1*mean(select(:)) + multiple_space_weight2*mean(select2(:));
                        end
                        
                        select = w(k,j_indices,i_indices);
                        w_coarse(k,j,i) = mean(select(:));
                        
                        select = tabs(k,j_indices,i_indices); % This is a non prognostic variable - NOT SUPPOSE TO USE, BUT RECALCULATE ! 
                        tabs_coarse(k,j,i) = mean(select(:));
                        
                        select = t(k,j_indices,i_indices);
                        t_coarse(k,j,i) = mean(select(:));
                        
                        select = tfull(k,j_indices,i_indices);
                        tfull_coarse(k,j,i) = mean(select(:));
                        
                        select = qv(k,j_indices,i_indices);
                        qv_coarse(k,j,i) = mean(select(:));
                        
                        select = qn(k,j_indices,i_indices); 
                        qn_coarse(k,j,i) = mean(select(:));
                        
                        select = qp(k,j_indices,i_indices);
                        qp_coarse(k,j,i) = mean(select(:));
                    end
                end
            end
            
            tfull_coarse_init = tfull_coarse; %This is a prognostic variable
            t_coarse_init = t_coarse; %This is a prognostic variable when qp is not included. 
%         end
%     end
% end
% {
            
            %03:03
            % check block size
            if mod(num_x, multiple_space)~=0
                error('number of x values not multiple of spatial block size')
            end
            if mod(num_y, multiple_space)~=0
                error('number of y values not multiple of spatial block size')
            end
            
            % change to mean over spatial blocks
            for i=1:num_blocks_x
                i_indices = [(i-1)*multiple_space+1:i*multiple_space];
                x_coarse(i) = mean(x(i_indices));
                for j=1:num_blocks_y
                    j_indices = [(j-1)*multiple_space+1:j*multiple_space];
                    y_coarse(j) = mean(y(j_indices));
                     
                    select = windspeed_high(j_indices,i_indices);
                    windspeed_coarse(j,i) = mean(select(:));

                    
                    for k=1:num_z
                        
                        select = Qrad(k,j_indices,i_indices); 
                        Qrad_coarse(k,j,i) = mean(select(:));
                        
                        %Advection fluxes
                        
                        select = tfull_flux_x(k,j_indices,i_indices);
                        tfull_flux_x_coarse(k,j,i) = mean(select(:));
                        
                        select = tfull_flux_y(k,j_indices,i_indices);
                        tfull_flux_y_coarse(k,j,i) = mean(select(:));
                        
                        select = tfull_flux_z(k,j_indices,i_indices);
                        tfull_flux_z_coarse(k,j,i) = mean(select(:));
                        
                        select = t_flux_x(k,j_indices,i_indices);
                        t_flux_x_coarse(k,j,i) = mean(select(:));
                        
                        select = t_flux_y(k,j_indices,i_indices);
                        t_flux_y_coarse(k,j,i) = mean(select(:));
                        
                        select = t_flux_z(k,j_indices,i_indices);
                        t_flux_z_coarse(k,j,i) = mean(select(:));
                        
                        select = qt_flux_x(k,j_indices,i_indices);
                        qt_flux_x_coarse(k,j,i) = mean(select(:));
                        
                        select = qt_flux_y(k,j_indices,i_indices);
                        qt_flux_y_coarse(k,j,i) = mean(select(:));
                        
                        select = qt_flux_z(k,j_indices,i_indices);
                        qt_flux_z_coarse(k,j,i) = mean(select(:));
                        
                        select = qp_flux_x(k,j_indices,i_indices);
                        qp_flux_x_coarse(k,j,i) = mean(select(:));
                        
                        select = qp_flux_y(k,j_indices,i_indices);
                        qp_flux_y_coarse(k,j,i) = mean(select(:));
                        
                        select = qp_flux_z(k,j_indices,i_indices);
                        qp_flux_z_coarse(k,j,i) = mean(select(:));
                        
                        %Advection tendencies
                        if calc_advection_tend
                            select = tfull_flux_x_tend(k,j_indices,i_indices);
                            tfull_flux_x_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = tfull_flux_y_tend(k,j_indices,i_indices);
                            tfull_flux_y_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = tfull_flux_z_tend(k,j_indices,i_indices);
                            tfull_flux_z_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = t_flux_x_tend(k,j_indices,i_indices);
                            t_flux_x_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = t_flux_y_tend(k,j_indices,i_indices);
                            t_flux_y_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = t_flux_z_tend(k,j_indices,i_indices);
                            t_flux_z_tend_coarse(k,j,i) = mean(select(:));
                            
                        end
                        %diffusion fluxes
                        select = tfull_diff_flx_x(k,j_indices,i_indices);
                        tfull_diff_flx_x_coarse(k,j,i) = mean(select(:));
                        
                        select = tfull_diff_flx_y(k,j_indices,i_indices);
                        tfull_diff_flx_y_coarse(k,j,i) = mean(select(:));
                        
                        select = tfull_diff_flx_z(k,j_indices,i_indices);
                        tfull_diff_flx_z_coarse(k,j,i) = mean(select(:));
                        
                        select = t_diff_flx_x(k,j_indices,i_indices);
                        t_diff_flx_x_coarse(k,j,i) = mean(select(:));
                        
                        select = t_diff_flx_y(k,j_indices,i_indices);
                        t_diff_flx_y_coarse(k,j,i) = mean(select(:));
                        
                        select = t_diff_flx_z(k,j_indices,i_indices);
                        t_diff_flx_z_coarse(k,j,i) = mean(select(:));
                        
                        select = qt_diff_flx_x(k,j_indices,i_indices);
                        qt_diff_flx_x_coarse(k,j,i) = mean(select(:));
                        
                        select = qt_diff_flx_y(k,j_indices,i_indices);
                        qt_diff_flx_y_coarse(k,j,i) = mean(select(:));
                        
                        select = qt_diff_flx_z(k,j_indices,i_indices);
                        qt_diff_flx_z_coarse(k,j,i) = mean(select(:));
                        
                        select = qp_diff_flx_x(k,j_indices,i_indices);
                        qp_diff_flx_x_coarse(k,j,i) = mean(select(:));
                        
                        select = qp_diff_flx_y(k,j_indices,i_indices);
                        qp_diff_flx_y_coarse(k,j,i) = mean(select(:));
                        
                        select = qp_diff_flx_z(k,j_indices,i_indices);
                        qp_diff_flx_z_coarse(k,j,i) = mean(select(:));
                        
                        %diffusion tendencies
                        if calc_diffusive_tend
                            select = tfull_diff_x_tend(k,j_indices,i_indices);
                            tfull_diff_x_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = tfull_diff_y_tend(k,j_indices,i_indices);
                            tfull_diff_y_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = tfull_diff_z_tend(k,j_indices,i_indices);
                            tfull_diff_z_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = t_diff_x_tend(k,j_indices,i_indices);
                            t_diff_x_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = t_diff_y_tend(k,j_indices,i_indices);
                            t_diff_y_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = t_diff_z_tend(k,j_indices,i_indices);
                            t_diff_z_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = qt_diff_x_tend(k,j_indices,i_indices);
                            qt_diff_x_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = qt_diff_y_tend(k,j_indices,i_indices);
                            qt_diff_y_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = qt_diff_z_tend(k,j_indices,i_indices);
                            qt_diff_z_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = qp_diff_x_tend(k,j_indices,i_indices);
                            qp_diff_x_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = qp_diff_y_tend(k,j_indices,i_indices);
                            qp_diff_y_tend_coarse(k,j,i) = mean(select(:));
                            
                            select = qp_diff_z_tend(k,j_indices,i_indices);
                            qp_diff_z_tend_coarse(k,j,i) = mean(select(:));
                        end
                        %Precip_fall tendencies
                        select = dqp_fall(k,j_indices,i_indices);
                        dqp_fall_coarse(k,j,i) = mean(select(:));
                        
                        select = t_fall_tend(k,j_indices,i_indices);
                        t_fall_tend_coarse(k,j,i) = mean(select(:));
                        
                        %cloud
                        select = cloud_lat_heat(k,j_indices,i_indices);
                        cloud_lat_heat_coarse(k,j,i) = mean(select(:));
                        
                        select = cloud_qt_tend(k,j_indices,i_indices);
                        cloud_qt_tend_coarse(k,j,i) = mean(select(:));
                        
                        %qp micro tendency
                        select = dqp(k,j_indices,i_indices);
                        dqp_coarse(k,j,i) = mean(select(:));
                        
                        select = dqp_t_tend(k,j_indices,i_indices);
                        dqp_t_tend_coarse(k,j,i) = mean(select(:));
                        select = precip(k,j_indices,i_indices);
                        precip_coarse(k,j,i) = mean(select(:));
                        
                        select = precip_energy(k,j_indices,i_indices);
                        precip_energy_coarse(k,j,i) = mean(select(:));
                        
                        %The cloud fluxes:
                        select = fz(k,j_indices,i_indices);
                        fz_coarse(k,j,i) = mean(select(:));
                        
                        select = fzt(k,j_indices,i_indices);
                        fzt_coarse(k,j,i) = mean(select(:));
                        
                        %diffusivity.
                        select = tkz_out(k,j_indices,i_indices);
                        tkz_out_coarse(k,j,i) = mean(select(:));
                        
                        select = Pr1(k,j_indices,i_indices);
                        Pr1_coarse(k,j,i) = mean(select(:));
                        
                        %Momentum vertical advection 
                        select = dwdt_advect(k,j_indices,i_indices);
                        dwdt_advect_coarse(k,j,i) = mean(select(:));
                    end
                end
            end
            
            

            
            
            qt_coarse = qv_coarse+qn_coarse; % qt_coarse is the prognostic variable that in any case we can access
            qv_coarse_init = qv_coarse;
            qn_coarse_init = qn_coarse;
            qp_coarse_init = qp_coarse;


[qn_resolved, tabs_resolved, qt_coarse, qp_coarse,~,~] = update_qn_tabs(tbgmax, tbgmin, tprmax, tprmin, fac_cond, fac_fus, fac_sub, tgrmax, tgrmin, qp_threshold, num_z, ...
    num_blocks_x, num_blocks_y, qn_coarse, qt_coarse, qp_coarse, tfull_coarse, gamaz, pres);    

            %qn_coarse and tabs_coarse to be the resolved quanteties.             
            qn_coarse = qn_resolved;
            tabs_coarse = tabs_resolved;
            
            %Calcualting relative humidity at initial time step
            omp_3d_coarse  = max(0.,min(1.,(tabs_resolved-tprmin)*a_pr));
            omn_3d_coarse  = max(0.,min(1.,(tabs_resolved-tbgmin)*a_bg));
            for k = 1:num_z
                qsat_resolved_init(k,:,:) = omn_3d_coarse(k,:,:).*qsatw(tabs_resolved(k,:,:),pres(k))+ ...
                    (1.-omn_3d_coarse(k,:,:)).*qsati(tabs_resolved(k,:,:),pres(k));
            end
            rh_init = qt_coarse./qsat_resolved_init;

            
            
            %%
            % Step 03:04 Set grid (setgrid.f90) for calculations
            
            % liquid/ice water static energy h_L divided by cp
            % but here set qp to zero in SAM with rf
%             omn_3d_coarse  = max(0.,min(1.,(tabs_coarse-tbgmin)*a_bg));
            dx_coarse= dx.*multiple_space;
            dy_coarse= dy.*multiple_space;
            
            
            
            % from setgrid.f90
            dz = 0.5*(z(1)+z(2));
            adzw = zeros(length(num_z));% Yani added
            for k=2:num_z
                adzw(k) = (z(k)-z(k-1))/dz;
            end
            adzw(1) = 1.;
            adzw(num_z+1) = adzw(num_z);
            adz = zeros(length(num_z));% Yani added
            
            for k=2:num_z-1
                adz(k) = 0.5*(z(k+1)-z(k-1))/dz;
            end
            adz(1) = 1.;
            adz(num_z) = adzw(num_z);
            
            
            rdz2=1./(dz*dz);
            rdz=1./dz;
            
            %%
            % Step 03:05:Calculate shear_prod3D.f90 for the tke/tkz calculation
            def2_coarse = shear_prod3D(u_coarse,v_coarse,w_coarse,dx_coarse,dy_coarse,dz,adz,adzw,num_blocks_x,num_blocks_y,num_z,ravefactor);
            
            
            %%
            
            %%
            % Step 03:06: tke/tkz
            %tke_full.f90
            
            dxdzfactor  = 1.0;
            if dimfactor_as_in_samson
                dimfactor=max(1.,log10(sqrt(dx_coarse*dy_coarse)/dz/dxdzfactor/ravefactor));
            else
                dimfactor=max(1.,log10(sqrt(dx*dy)/dz/dxdzfactor/ravefactor));
            end
            Ck=0.1;
            
            
            Cs = 0.1944;
            Cs1 = 0.14;
            Pr = 3.0;
            % Ck=0.1;
            Ce=Ck^3/Cs^4;
            Ces=Ce/0.7*3.0;
            
            
            for k=1:num_z
                
                if k>1
                    kb=k-1;
                else
                    kb=k;
                end
                kc=min(k+1,num_z); %Yani added and changed the loop to reach till num_z
                
                grd=dz*adz(k)*dimfactor;
                
                rhoi = rhow(kc)/adzw(kc);
                rdz5=0.5*rdz2;
                
                % following from setdata.f90
                bet = ggr/tabs0(k);
                
                % following from tke_full.f90
                
                
                
                betdz=bet/dz/(adzw(kc)+adzw(k));
                Ce1=Ce/0.7*0.19;
                Ce2=Ce/0.7*0.51;
                
                if k==1
                    betdz=bet/dz/adzw(kc);
                    Ce1=Ces/0.7*0.19;
                    Ce2=Ces/0.7*0.51;
                end
                
                if k==num_z-1
                    betdz=bet/dz/adzw(k);
                    Ce1=Ces/0.7*0.19;
                    Ce2=Ces/0.7*0.51;
                end
                
                for i=1:num_blocks_x
                    for j=1:num_blocks_y
                        
                        omn = max(0.,min(1.,(tabs_coarse(k,j,i)-tbgmin)*a_bg));
                        omp = max(0.,min(1.,(tabs_coarse(k,j,i)-tprmin)*a_pr));
                        
                        lstarn = fac_cond+(1.-omn)*fac_fus;
                        lstarp = fac_cond+(1.-omp)*fac_fus;
                        
                        
                        if qn_coarse(k,j,i)>0.0
                            
                            dqsat = omn*dtqsatw(tabs_coarse(k,j,i),pres(k))+ ...
                                (1.-omn)*dtqsati(tabs_coarse(k,j,i),pres(k));
                            qsat = omn*qsatw(tabs_coarse(k,j,i),pres(k))+(1.-omn)*qsati(tabs_coarse(k,j,i),pres(k));
                            bbb = 1. + 0.61*qsat-qn_coarse(k,j,i) -qp_coarse(k,j,i)+1.61*tabs_coarse(k,j,i)*dqsat;
                            bbb = bbb / (1.+lstarn*dqsat);
                            % note replaced t with tfull in following line
                            buoy_sgs=betdz*(bbb*(tfull_coarse(kc,j,i)-tfull_coarse(kb,j,i)) ...
                                +(bbb*lstarn - (1.+lstarn*dqsat)*tabs_coarse(k,j,i))*(qt_coarse(kc,j,i)-qt_coarse(kb,j,i)) ...
                                +(bbb*lstarp - (1.+lstarp*dqsat)*tabs_coarse(k,j,i))*(qp_coarse(kc,j,i)-qp_coarse(kb,j,i)) );
                        else
                            
                            bbb = 1.+0.61*qt_coarse(k,j,i)-qp_coarse(k,j,i);
                            % note replaced t with tfull in following line
                            buoy_sgs=betdz*( bbb*(tfull_coarse(kc,j,i)-tfull_coarse(kb,j,i)) ...
                                +0.61*tabs_coarse(k,j,i)*(qt_coarse(kc,j,i)-qt_coarse(kb,j,i)) ...
                                +(bbb*lstarp-tabs_coarse(k,j,i))*(qp_coarse(kc,j,i)-qp_coarse(kb,j,i)) );
                        end
                        
                        
                        %try iterative solution.
                        smix_resolved = grd;
                        hhh = 1;
                        while hhh<100
                            smix_prev = smix_resolved;
                            if buoy_sgs<=0.0 % unstable/neutral
                                smix_resolved=grd;
                            else % stable
                                if hhh>1
                                    smix_resolved = sqrt(0.76*tkz_out_resolved(k,j,i)/sqrt(buoy_sgs)/Ck); 
                                end
                                smix_resolved = min(grd,max(0.1*grd,smix_resolved));
                            end
                            
                            if hhh > 2
                                smix_resolved = smix_prev + (smix_resolved - smix_prev )./10;
                            end
                            ratio_i=smix_resolved/grd;
                            Pr1_resolved(k,j,i)=1.+2.*ratio_i;
                            Cee_i=Ce1+Ce2*ratio_i;
                            tkz_out_resolved(k,j,i)=sqrt(Ck^3/Cee_i*max(0.,def2_coarse(k,j,i)-Pr1_resolved(k,j,i)*buoy_sgs))*smix_resolved^2;
                            if (buoy_sgs<0)
                                break
                            elseif (hhh>1 && abs(smix_resolved-smix_prev)./(abs(smix_prev + smix_resolved))<0.00001)
                                break
                            end
                            hhh = hhh + 1;
                        end
                    end
                end
            end
            
            
            %%
            
            % Step 03:07: Advect variables
            %simple version of scalar advection Advection
            
            % vertical fluxes
            % note w is on half levels but only first nzm are output by SAM
            % w is set to zero at top level n_z+1
            % w is first multiplied by rhow*dtn/dz in adams.f90 (note w and rhow
            % are grouped together in the differencing)
            % from advect_scalar3D.f90:
            % kb=max(1,k-1)
            % www(i,j,k)=max(0.,w(i,j,k))*f(i,j,kb)+min(0.,w(i,j,k))*f(i,j,k)
            % assuming dt->0 such that not iterative
            % leaving out non-oscillatory option,
            % see smolarkiewicz 2006 for general discussion
            
            
            % from setgrid.f90 - I think I need it again since some of the
            % vars were changed
            dz = 0.5*(z(1)+z(2));
            adzw = zeros(length(num_z));% Yani added
            for k=2:num_z
                adzw(k) = (z(k)-z(k-1))/dz;
            end
            adzw(1) = 1.;
            adzw(num_z+1) = adzw(num_z);
            adz = zeros(length(num_z));% Yani added
            
            for k=2:num_z-1
                adz(k) = 0.5*(z(k+1)-z(k-1))/dz;
            end
            adz(1) = 1.;
            adz(num_z) = adzw(num_z);
            
            
            rdz2=1./(dz*dz);
            rdz=1./dz;
            
            
            %Yani: need to understand why he defines the flux with the 0.5 factor + what is rhow ? (Could it be related to the vert coordinate, but then why isn't it in the denominator?).
            % Paul mentioned that this is a contiuum form of a simplification of some scheme. There might be a reference in the MARAT paper to the scheme.
            for k=1:num_z
                kb = max(1,k-1);
                w_rhow = w_coarse(k,:,:)*rhow(k);
                t_flux_z_resolved(k,:,:) = 0.5*w_rhow.*(t_coarse(k,:,:) + t_coarse(kb,:,:));%Yani simplified
                tfull_flux_z_resolved(k,:,:)  = 0.5*w_rhow.*(tfull_coarse(k,:,:) + tfull_coarse(kb,:,:));
            end
            
            rhow_num_x = repmat(rho,[1,num_blocks_y]);
            for i=1:num_blocks_x
                ib = i-1;
                if ib == 0
                    ib = num_blocks_x; %Try to impose periodic domain - Yani ?
                end
                u_rhow = u_coarse(:,:,i).*rhow_num_x;
                t_flux_x_resolved(:,:,i) = 0.5*u_rhow.*(t_coarse(:,:,ib) + t_coarse(:,:,i));
                tfull_flux_x_resolved(:,:,i) = 0.5*u_rhow.*(tfull_coarse(:,:,ib) + tfull_coarse(:,:,i));
            end
            rhow_num_y = repmat(rho,[1,num_blocks_x]);
            for j=1:num_blocks_y
                jb = max(1,j-1);
                v_rhow = squeeze(v_coarse(:,j,:)).*rhow_num_y;
                t_flux_y_resolved(:,j,:) = 0.5*v_rhow.*squeeze(t_coarse(:,jb,:) + t_coarse(:,j,:));
                tfull_flux_y_resolved(:,j,:) = 0.5*v_rhow.*squeeze(tfull_coarse(:,jb,:) + tfull_coarse(:,j,:));
            end
            
            
            
            
            
            
            if calc_advection_tend %If I realy want to calculate the tendency itself
                % {
                for k = 1:num_z-1
                    t_flux_z_tend_resolved(k,:,:) = - (t_flux_z_resolved(k+1,:,:) - t_flux_z_resolved(k,:,:))./(adz(k).*dz.*rho(k));
                    tfull_flux_z_tend_resolved(k,:,:) = - (tfull_flux_z_resolved(k+1,:,:) - tfull_flux_z_resolved(k,:,:))./(adz(k).*dz.*rho(k));
                end
                k = num_z;
                t_flux_z_tend_resolved(k,:,:) = - (0.0 - t_flux_z_resolved(k,:,:))./(adz(k).*dz.*rho(k));
                tfull_flux_z_tend_resolved(k,:,:) = - (0.0 - tfull_flux_z_resolved(k,:,:))./(adz(k).*dz.*rho(k));
                
                for i = 1:num_blocks_x-1
                    t_flux_x_tend_resolved(:,:,i) = - (t_flux_x_resolved(:,:,i+1) - t_flux_x_resolved(:,:,i))./(dx_coarse.*rhow_num_x);
                    tfull_flux_x_tend_resolved(:,:,i) =  - (tfull_flux_x_resolved(:,:,i+1) - tfull_flux_x_resolved(:,:,i))./(dx_coarse.*rhow_num_x);
                end
                i = num_blocks_x;
                t_flux_x_tend_resolved(:,:,i) = - (t_flux_x_resolved(:,:,1) - t_flux_x_resolved(:,:,i))./(dx_coarse.*rhow_num_x);
                tfull_flux_x_tend_resolved(:,:,i) = - (tfull_flux_x_resolved(:,:,1) - tfull_flux_x_resolved(:,:,i))./(dx_coarse.*rhow_num_x);
                
                for j = 1:num_blocks_y-1
                    t_flux_y_tend_resolved(:,j,:) = - squeeze(t_flux_y_resolved(:,j+1,:) - t_flux_y_resolved(:,j,:))./(dy_coarse.*rhow_num_y);
                    tfull_flux_y_tend_resolved(:,j,:) = - squeeze(tfull_flux_y_resolved(:,j+1,:) - tfull_flux_y_resolved(:,j,:))./(dy_coarse.*rhow_num_y);
                end
                j = num_blocks_y;
                t_flux_y_tend_resolved(:,j,:) = -(0.0 - squeeze(t_flux_y_resolved(:,num_blocks_y,:)))./(dy_coarse.*rhow_num_y);
                tfull_flux_y_tend_resolved(:,j,:) =-(0.0 - squeeze(tfull_flux_y_resolved(:,num_blocks_y,:)))./(dy_coarse.*rhow_num_y);
                
                
                if advect_fields
                    t_coarse = t_coarse + (t_flux_x_tend_resolved + t_flux_y_tend_resolved + t_flux_z_tend_resolved)*dtn;
                    tfull_coarse = tfull_coarse + (tfull_flux_x_tend_resolved + tfull_flux_y_tend_resolved + tfull_flux_z_tend_resolved)*dtn;
                end
            end
            
            
            [f_out_resolved,qt_flux_x_resolved,qt_flux_y_resolved,qt_flux_z_resolved,qt_adv_tend_resolved] = advect_scalar3D_func_try_reduce_complexity(qt_coarse,u_coarse,v_coarse,w_coarse,rho,rhow,dx_coarse,dy_coarse,dz,dtn,num_blocks_x,num_blocks_y, num_z,adz);
            if advect_fields
                qt_coarse = f_out_resolved;
            end
            
            [f_out_resolved,qp_flux_x_resolved,qp_flux_y_resolved,qp_flux_z_resolved,qp_adv_tend_resolved] = advect_scalar3D_func_try_reduce_complexity(qp_coarse,u_coarse,v_coarse,w_coarse,rho,rhow,dx_coarse,dy_coarse,dz,dtn,num_blocks_x,num_blocks_y, num_z,adz);
            if advect_fields
                qp_coarse = f_out_resolved;
            end
            
            
            
            qp_coarse = max(0,qp_coarse); % It is possible that due to the approximate advection scheme I use I have errors that lead to negative qp values
            qt_coarse = max(0,qt_coarse);
            
             %% 03:07a: vertical momentum advection
            [dudt_advect_resolved,dvdt_advect_resolved,dwdt_advect_resolved] = advect2_mom_z(num_blocks_x,num_blocks_y,num_z,dz,rho,rhow,adz,adzw,u_coarse,v_coarse,w_coarse);
       
            
            %%
            %Step 03:07b precip_fall calculation
            
            % calculate precipitation flux
            %%precip_fall
            % from precip.f90
            
            % only kept parts needed to calculate precipitation flux
            % and neglected non-oscillatory option for speed
            
            sprintf('the precip fall func coarse')
            
            [dqp_fall_resolved,t_fall_tend_resolved,precip_resolved]= ...
                precip_fall(qp_coarse,tabs_coarse,rho,rhow,rhofac,num_blocks_x,num_blocks_y,num_z,dz,adz,dtn,tprmin,a_pr,fac_fus,...
                fac_cond,crain,vrain,tgrmin,a_gr,qp_threshold,vgrau,cgrau,vsnow,csnow);
            
            precip_energy_resolved = precip_resolved.*(fac_cond+(1.-omp_3d_coarse).*fac_fus); 
            
            if do_fall_tend_tfull
                tfull_coarse = tfull_coarse + t_fall_tend_resolved; 
                qp_coarse = qp_coarse + dqp_fall_resolved;
            end
            dqp_fall_resolved = dqp_fall_resolved./dtn; 
            t_fall_tend_resolved = t_fall_tend_resolved./dtn;
            
            
            
            
            %%
            % Step 03:08: surface fluxes
            % First find surface fluxes of t, tfull and qt following surface.f90
            umin = 1.0;
            cd=1.1e-3;
            wrk=(log(10/1.e-4)/log(z(1)/1.e-4))^2;
            fluxbt = zeros(num_blocks_y,num_blocks_x); %Yani added for parfor
            fluxbtfull = zeros(num_blocks_y,num_blocks_x);
            fluxbqt = zeros(num_blocks_y,num_blocks_x);
            
            resc_fact_res = 1/res;
            resc_fact_res_2 = resc_fact_res/2;
            resc_fact_res_2_1min = 1- resc_fact_res_2;
            resc_fact_res_4 = resc_fact_res/4;
            fact1_uv = 0.25 + resc_fact_res_4;
            fact1_uv_min = 0.25 - resc_fact_res_4;
%             resc_fact_res_1min_2 = resc_fact_res_1min/2;
            
            for i=1:num_blocks_x
                for j=1:num_blocks_y
                    
                    if i<num_blocks_x
                        ic=i+1;
                    else
                        ic=1;
                    end
                    
                    if j<num_blocks_y
                        jc=j+1;
                    else
                        jc=j;
                    end
                    
                    if i == 1
                        ib = num_blocks_x;
                    else
                        ib = i - 1;
                    end
                    
                    if j == 1
                        jb = 1;
                    else
                        jb = j -1;
                    end
                    
                    ubot=0.5*(u_coarse(1,j,ic)+u_coarse(1,j,i));
                    vbot=0.5*(v_coarse(1,jc,i)+v_coarse(1,j,i));
                    windspeed=sqrt(ubot^2+vbot^2+umin^2);
                    delt     = t_coarse(1,j,i)-gamaz(1) - sstxy_coarse(j,i);
                    deltfull = tfull_coarse(1,j,i)-gamaz(1) - sstxy_coarse(j,i);
                    ssq = qsatw(sstxy_coarse(j,i),pres(1));
                    delqt   = qt_coarse(1,j,i)  - ssq;
                    fluxbt(j,i) = -cd*windspeed*delt*wrk;
                    fluxbtfull(j,i) = -cd*windspeed*deltfull*wrk;
                    fluxbqt(j,i) = -cd*windspeed*delqt*wrk;
                    windspeed_resolved(j,i) = windspeed;
                    ubot=resc_fact_res_2_1min .* u_coarse(1,j,i) + resc_fact_res_2 * u_coarse(1,j,ic);
                    vbot = (fact1_uv .*v_coarse(1,j,i) + fact1_uv.*v_coarse(1,jc,i) + fact1_uv_min.*v_coarse(1,jc,ib) + fact1_uv_min.*v_coarse(1,j,ib));
                    windspeed=sqrt(ubot^2+vbot^2+umin^2);
                    fluxbu_samson_resolved(j,i) = -rho(1)*(u_coarse(1,j,i))*cd*windspeed*wrk*rdz*rhow(1); %Note that we need to  multiply by things to get tendency

                    vbot= resc_fact_res_2_1min .*v_coarse(1,j,i) + resc_fact_res_2 * v_coarse(1,jc,i);
                    ubot=(fact1_uv .*u_coarse(1,j,i) +fact1_uv_min.*u_coarse(1,jb,i) + fact1_uv_min.*u_coarse(1,jb,ic) + fact1_uv .*u_coarse(1,j,ic));
                    windspeed=sqrt(ubot^2+vbot^2+umin^2);
                    fluxbv_samson_resolved(j,i) = -rho(1)*(v_coarse(1,j,i))*cd*windspeed*wrk*rdz*rhow(1);
                    
                    ubot=u_coarse(1,j,i);
                    vbot=0.25*(v_coarse(1,j,i) + v_coarse(1,jc,i) + v_coarse(1,jc,ib) + v_coarse(1,j,ib));
                    windspeed=sqrt(ubot^2+vbot^2+umin^2);
                    fluxbu_resolved(j,i) = -rho(1)*(u_coarse(1,j,i))*cd*windspeed*wrk*rdz*rhow(1); %Note that we need to  multiply by things to get tendency
                    vbot=v_coarse(1,j,i);
                    ubot=0.25*(u_coarse(1,j,i) + u_coarse(1,jb,i) + u_coarse(1,jb,ic) + u_coarse(1,j,ic));
                    windspeed=sqrt(ubot^2+vbot^2+umin^2);
                    fluxbv_resolved(j,i) = -rho(1)*(v_coarse(1,j,i))*cd*windspeed*wrk*rdz*rhow(1);
                    
                end
            end
            
            % Step 03:09: vertical diffusion
            
            for i=1:num_blocks_x
                for j=1:num_blocks_y
                    t_diff_flx_z_resolved(1,j,i)=fluxbt(j,i)*rdz*rhow(1);
                    tfull_diff_flx_z_resolved(1,j,i)=fluxbtfull(j,i)*rdz*rhow(1);
                    qt_diff_flx_z_resolved(1,j,i)=fluxbqt(j,i)*rdz*rhow(1);
                    qp_diff_flx_z_resolved(1,j,i)= 0.0;
                end
            end
            
            for k=1:num_z-1
                kc = k + 1;
                rhoi = rhow(kc)/adzw(kc);
                for i=1:num_blocks_x
                    for j=1:num_blocks_y
                        tkh_z=rdz5*(tkz_out_resolved(k,j,i)*Pr1_resolved(k,j,i)+tkz_out_resolved(kc,j,i)*Pr1_resolved(kc,j,i)); 
                        t_diff_flx_z_resolved(kc,j,i)=-tkh_z*(t_coarse(kc,j,i)-t_coarse(k,j,i))*rhoi/ravefactor;
                        tfull_diff_flx_z_resolved(kc,j,i)=-tkh_z*(tfull_coarse(kc,j,i)-tfull_coarse(k,j,i))*rhoi/ravefactor;
                        qt_diff_flx_z_resolved(kc,j,i)=-tkh_z*(qt_coarse(kc,j,i)-qt_coarse(k,j,i))*rhoi/ravefactor;
                        qp_diff_flx_z_resolved(kc,j,i)=-tkh_z*(qp_coarse(kc,j,i)-qp_coarse(k,j,i))*rhoi/ravefactor;
                        
                        
                    end
                end
            end
            
            % the above includes an additional 1/dz (in rdz and rdz5), so multiply by dz to get the actual fluxes
            t_diff_flx_z_resolved = t_diff_flx_z_resolved*dz;
            tfull_diff_flx_z_resolved = tfull_diff_flx_z_resolved*dz;
            qt_diff_flx_z_resolved = qt_diff_flx_z_resolved*dz;
            qp_diff_flx_z_resolved = qp_diff_flx_z_resolved*dz;
            
            
            
            
            %Calculate  diffusive vertical tendencies:
            if calc_diffusive_tend
                for k=1:num_z
                    irhoadz_dz = irhoadz(k)./dz;
                    for j=1:num_blocks_y
                        for i=1:num_blocks_x
                            if k == num_z
                                t_diff_z_tend_resolved(k,j,i)=-(0-t_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                                tfull_diff_z_tend_resolved(k,j,i)=-(0-tfull_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                                qt_diff_z_tend_resolved(k,j,i)=-(0-qt_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                                qp_diff_z_tend_resolved(k,j,i)=-(0-qp_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                            else
                                kt=k + 1;
                                t_diff_z_tend_resolved(k,j,i)=-(t_diff_flx_z_resolved(kt,j,i)-t_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                                tfull_diff_z_tend_resolved(k,j,i)=-(tfull_diff_flx_z_resolved(kt,j,i)-tfull_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                                qt_diff_z_tend_resolved(k,j,i)=-(qt_diff_flx_z_resolved(kt,j,i)-qt_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                                qp_diff_z_tend_resolved(k,j,i)=-(qp_diff_flx_z_resolved(kt,j,i)-qp_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                            end
                        end
                    end
                end
            end
            
            %%
            % Step 03:10: horizontal diffusion
            
            
            if do_horizontal_diffusion_at_all
                %Horizontal diffusion
                rdx2=1./(dx_coarse*dx_coarse);
                rdy2=1./(dy_coarse*dy_coarse);
                tkh_xy = tkz_out_resolved; 
                rdx5=0.5*rdx2;
                rdy5=0.5*rdy2;
                
                % x diffusion
                for k=1:num_z
                    for j=1:num_blocks_y
                        for i=1:num_blocks_x
                            ic = i + 1;
                            if i == num_blocks_x
                                ic = 1;
                            end
                            tkx=rdx5*(tkh_xy(k,j,i)*Pr1_resolved(k,j,i)+tkh_xy(k,j,ic)*Pr1_resolved(k,j,ic));
                            t_diff_flx_x_resolved(k,j,i)=-tkx*(t_coarse(k,j,ic)-t_coarse(k,j,i))*ravefactor;
                            tfull_diff_flx_x_resolved(k,j,i)=-tkx*(tfull_coarse(k,j,ic)-tfull_coarse(k,j,i))*ravefactor;
                            qt_diff_flx_x_resolved(k,j,i)=-tkx*(qt_coarse(k,j,ic)-qt_coarse(k,j,i))*ravefactor;
                            qp_diff_flx_x_resolved(k,j,i)=-tkx*(qp_coarse(k,j,ic)-qp_coarse(k,j,i))*ravefactor;
                        end
                    end
                    
                end
                if calc_diffusive_tend
                    for k=1:num_z
                        for j=1:num_blocks_y
                            for i=1:num_blocks_x
                                if i == 1
                                    ib = num_blocks_x;
                                else
                                    ib=i-1;
                                end
                                t_diff_x_tend_resolved(k,j,i)=-(t_diff_flx_x_resolved(k,j,i)-t_diff_flx_x_resolved(k,j,ib));
                                tfull_diff_x_tend_resolved(k,j,i)=-(tfull_diff_flx_x_resolved(k,j,i)-tfull_diff_flx_x_resolved(k,j,ib));
                                qt_diff_x_tend_resolved(k,j,i)=-(qt_diff_flx_x_resolved(k,j,i)-qt_diff_flx_x_resolved(k,j,ib));
                                qp_diff_x_tend_resolved(k,j,i)=-(qp_diff_flx_x_resolved(k,j,i)-qp_diff_flx_x_resolved(k,j,ib));
                            end
                        end
                    end
                end
                
                % y diffusion
                for k=1:num_z
                    for j=1:num_blocks_y
                        for i=1:num_blocks_x
                            jc = j + 1;
                            if j == num_blocks_y
                                %             tkx=rdy5*(tkh_xy(k,j,i)*Pr1(k,j,i)+0);
                                t_diff_flx_y_resolved(k,j,i)=0;%-tkx*(0-t(k,j,i))*ravefactor;
                                tfull_diff_flx_y_resolved(k,j,i)=0;%-tkx*(0-tfull(k,j,i))*ravefactor;
                                qt_diff_flx_y_resolved(k,j,i)=0;%-tkx*(0-qt(k,j,i))*ravefactor;
                                qp_diff_flx_y_resolved(k,j,i)=0;%-tkx*(0-qp(k,j,i))*ravefactor;
                            else
                                tkx=rdy5*(tkh_xy(k,j,i)*Pr1_resolved(k,j,i)+tkh_xy(k,jc,i)*Pr1_resolved(k,jc,i));
                                t_diff_flx_y_resolved(k,j,i)=-tkx*(t_coarse(k,jc,i)-t_coarse(k,j,i))*ravefactor;
                                tfull_diff_flx_y_resolved(k,j,i)=-tkx*(tfull_coarse(k,jc,i)-tfull_coarse(k,j,i))*ravefactor;
                                qt_diff_flx_y_resolved(k,j,i)=-tkx*(qt_coarse(k,jc,i)-qt_coarse(k,j,i))*ravefactor;
                                qp_diff_flx_y_resolved(k,j,i)=-tkx*(qp_coarse(k,jc,i)-qp_coarse(k,j,i))*ravefactor;
                            end
                            
                        end
                    end
                    
                end
                
                if calc_diffusive_tend
                    
                    for k=1:num_z
                        for j=1:num_blocks_y
                            for i=1:num_blocks_x
                                if j == 1
                                    t_diff_y_tend_resolved(k,j,i)=-(t_diff_flx_y_resolved(k,j,i)-0);
                                    tfull_diff_y_tend_resolved(k,j,i)=-(tfull_diff_flx_y_resolved(k,j,i)-0);
                                    qt_diff_y_tend_resolved(k,j,i)=-(qt_diff_flx_y_resolved(k,j,i)-0);
                                    qp_diff_y_tend_resolved(k,j,i)=-(qp_diff_flx_y_resolved(k,j,i)-0);
                                else
                                    jb=j-1;
                                    t_diff_y_tend_resolved(k,j,i)=-(t_diff_flx_y_resolved(k,j,i)-t_diff_flx_y_resolved(k,jb,i));
                                    tfull_diff_y_tend_resolved(k,j,i)=-(tfull_diff_flx_y_resolved(k,j,i)-tfull_diff_flx_y_resolved(k,jb,i));
                                    qt_diff_y_tend_resolved(k,j,i)=-(qt_diff_flx_y_resolved(k,j,i)-qt_diff_flx_y_resolved(k,jb,i));
                                    qp_diff_y_tend_resolved(k,j,i)=-(qp_diff_flx_y_resolved(k,j,i)-qp_diff_flx_y_resolved(k,jb,i));
                                end
                            end
                        end
                    end
                end
            end
            
            if (diffuse_fields==1 && calc_diffusive_tend==1)
                t_coarse = t_coarse + (t_diff_x_tend_resolved + t_diff_y_tend_resolved + t_diff_z_tend_resolved).*dtn;
                tfull_coarse = tfull_coarse + (tfull_diff_x_tend_resolved + tfull_diff_y_tend_resolved + tfull_diff_z_tend_resolved).*dtn;
                qt_coarse = qt_coarse + (qt_diff_x_tend_resolved + qt_diff_y_tend_resolved + qt_diff_z_tend_resolved).*dtn;
                qp_coarse = qp_coarse + (qp_diff_x_tend_resolved + qp_diff_y_tend_resolved + qp_diff_z_tend_resolved).*dtn;
            end
            
            
            if (diffuse_fields==1 && calc_diffusive_tend==0) %Diffuse fields without calculating seperately the tendencies
                
                for k=1:num_z
                    irhoadz_dz = irhoadz(k)./dz;
                    for j=1:num_blocks_y
                        for i=1:num_blocks_x
                            if k == num_z
                                t_coarse(k,j,i)=t_coarse(k,j,i) - dtn.*(0-t_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                                tfull_coarse(k,j,i)=tfull_coarse(k,j,i) - dtn.*(0-tfull_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                                qt_coarse(k,j,i)=qt_coarse(k,j,i) - dtn.*(0-qt_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                                qp_coarse(k,j,i)=qp_coarse(k,j,i) - dtn.*(0-qp_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                            else
                                kt=k + 1;
                                t_coarse(k,j,i)=t_coarse(k,j,i) - dtn.*(t_diff_flx_z_resolved(kt,j,i)-t_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                                tfull_coarse(k,j,i)=tfull_coarse(k,j,i) - dtn.*(tfull_diff_flx_z_resolved(kt,j,i)-tfull_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                                qt_coarse(k,j,i)=qt_coarse(k,j,i) - dtn.*(qt_diff_flx_z_resolved(kt,j,i)-qt_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                                qp_coarse(k,j,i)=qp_coarse(k,j,i) - dtn.*(qp_diff_flx_z_resolved(kt,j,i)-qp_diff_flx_z_resolved(k,j,i)).*irhoadz_dz;
                            end
                        end
                    end
                end
                
                if do_horizontal_diffusion_at_all
                    for k=1:num_z
                        for j=1:num_blocks_y
                            for i=1:num_blocks_x
                                if i == 1
                                    ib = num_blocks_x;
                                else
                                    ib=i-1;
                                end
                                t_coarse(k,j,i)=t_coarse(k,j,i) - dtn.*(t_diff_flx_x_resolved(k,j,i)-t_diff_flx_x_resolved(k,j,ib));
                                tfull_coarse(k,j,i)=tfull_coarse(k,j,i) - dtn.*(tfull_diff_flx_x_resolved(k,j,i)-tfull_diff_flx_x_resolved(k,j,ib));
                                qt_coarse(k,j,i)=qt_coarse(k,j,i) - dtn.*(qt_diff_flx_x_resolved(k,j,i)-qt_diff_flx_x_resolved(k,j,ib));
                                qp_coarse(k,j,i)=qp_coarse(k,j,i) - dtn.*(qp_diff_flx_x_resolved(k,j,i)-qp_diff_flx_x_resolved(k,j,ib));
                            end
                        end
                    end
                    
                    
                    
                    for k=1:num_z
                        for j=1:num_blocks_y
                            for i=1:num_blocks_x
                                if j == 1
                                    t_coarse(k,j,i)=t_coarse(k,j,i) - dtn.*(t_diff_flx_y_resolved(k,j,i)-0);
                                    tfull_coarse(k,j,i)=tfull_coarse(k,j,i) - dtn.*(tfull_diff_flx_y_resolved(k,j,i)-0);
                                    qt_coarse(k,j,i)=qt_coarse(k,j,i) - dtn.*(qt_diff_flx_y_resolved(k,j,i)-0);
                                    qp_coarse(k,j,i)=qp_coarse(k,j,i) - dtn.*(qp_diff_flx_y_resolved(k,j,i)-0);
                                else
                                    jb=j-1;
                                    t_coarse(k,j,i)=t_coarse(k,j,i) - dtn.*(t_diff_flx_y_resolved(k,j,i)-t_diff_flx_y_resolved(k,jb,i));
                                    tfull_coarse(k,j,i)=tfull_coarse(k,j,i) - dtn.*(tfull_diff_flx_y_resolved(k,j,i)-tfull_diff_flx_y_resolved(k,jb,i));
                                    qt_coarse(k,j,i)=qt_coarse(k,j,i) - dtn.*(qt_diff_flx_y_resolved(k,j,i)-qt_diff_flx_y_resolved(k,jb,i));
                                    qp_coarse(k,j,i)=qp_coarse(k,j,i) - dtn.*(qp_diff_flx_y_resolved(k,j,i)-qp_diff_flx_y_resolved(k,jb,i));
                                end
                            end
                        end
                    end
                end
                
            end
            
            t_diff_flx_x_resolved = t_diff_flx_x_resolved.*dx_coarse;
            tfull_diff_flx_x_resolved = tfull_diff_flx_x_resolved.*dx_coarse;
            qt_diff_flx_x_resolved = qt_diff_flx_x_resolved.*dx_coarse;
            qp_diff_flx_x_resolved = qp_diff_flx_x_resolved.*dx_coarse;
            
            t_diff_flx_y_resolved = t_diff_flx_y_resolved.*dy_coarse;
            tfull_diff_flx_y_resolved = tfull_diff_flx_y_resolved.*dy_coarse;
            qt_diff_flx_y_resolved = qt_diff_flx_y_resolved.*dy_coarse;
            qp_diff_flx_y_resolved = qp_diff_flx_y_resolved.*dy_coarse;
            
            
            %%
            % Step 03:11: do cloud (for recalculation of qn_coarse, and sedimentation)
            %%Do cloud before microphysics! cloud.f90
            [qn_coarse, tabs_coarse, qt_coarse, qp_coarse,kmin,kmax] = update_qn_tabs(tbgmax, tbgmin, tprmax, tprmin, fac_cond, fac_fus, fac_sub, tgrmax, tgrmin, qp_threshold, num_z, ...
    num_blocks_x, num_blocks_y, qn_coarse, qt_coarse, qp_coarse, tfull_coarse, gamaz, pres);
    
        
            
            if do_sedimentation
                % Sedimentation of ice and water:
                qifall = zeros(num_z,1);
                tlatqi = zeros(num_z,1);
                for k = 1:num_z
                    qifall(k) = 0.;
                    tlatqi(k) = 0.;
                end
                
                coef_cl = 1.19e8*(3./(4.*3.1415*1000.*Nc0*1.e6))^(2./3.)*exp(5.*log(1.5)^2);
                
                for k = max(1,kmin-1):kmax
                    %    ! Set up indices for x-y planes above and below current plane.
                    kc = min(num_z,k+1);
                    kb = max(1,k-1);
                    for j = 1:num_blocks_y
                        for i = 1:num_blocks_x
                            coef = dtn/(0.5*(adz(kb)+adz(k))*dz);
                            %
                            %          ! Compute cloud ice density in this cell and the ones above/below.
                            %          ! Since cloud ice is falling, the above cell is u (upwind),
                            %          ! this cell is c (center) and the one below is d (downwind).
                            omnu = max(0.,min(1.,(tabs_coarse(kc,j,i)-tbgmin)*a_bg));
                            omnc = max(0.,min(1.,(tabs_coarse(k,j,i) -tbgmin)*a_bg));
                            omnd = max(0.,min(1.,(tabs_coarse(kb,j,i)-tbgmin)*a_bg));
                            
                            qiu = rho(kc)*qn_coarse(kc,j,i)*(1.-omnu);
                            qic = rho(k) *qn_coarse(k,j,i) *(1.-omnc);
                            qid = rho(kb)*qn_coarse(kb,j,i)*(1.-omnd);
                            
                            %          ! Ice sedimentation velocity depends on ice content. The fiting is
                            %          ! based on the data by Heymsfield (JAS,2003). -Marat
                            %          ! 0.1 m/s low bound was suggested by Chris Bretherton
                            vt_ice = max(0.1,0.5*log10(qic+1.e-12)+3.);
                            
                            %          ! Use MC flux limiter in computation of flux correction.
                            %          ! (MC = monotonized centered difference).
                            if (qic==qid)
                                tmp_phi = 0.;
                            else
                                tmp_theta = (qiu-qic)/(qic-qid);
                                tmp_phi = max(0.,min([0.5*(1.+tmp_theta),2.,2.*tmp_theta]));
                            end
                            
                            %          ! Compute limited flux.
                            %          ! Since falling cloud ice is a 1D advection problem, this
                            %          ! flux-limited advection scheme is monotonic.
                            fluxi = -vt_ice*(qic - 0.5*(1.-coef*vt_ice)*tmp_phi*(qic-qid));
                            
                            doclouddropsed = 0;%This is our config- Yani
                            if(doclouddropsed)
                                %             ! Compute cloud water density in this cell and the ones above/below
                                %             ! Since cloud water is falling, the above cell is u (upwind),
                                %             ! this cell is c (center) and the one below is d (downwind).
                                qiu = rho(kc)*qn_coarse(kc,j,i)*omnu;
                                qic = rho(k) *qn_coarse(k,j,i) *omnc;
                                qid = rho(kb)*qn_coarse(kb,j,i)*omnd;
                                
                                vt_cl = coef_cl*(qic+1.e-12)^(2./3.);
                                
                                %             ! Use MC flux limiter in computation of flux correction.
                                %             ! (MC = monotonized centered difference).
                                if (qic==qid)
                                    tmp_phi = 0.;
                                else
                                    tmp_theta = (qiu-qic)/(qic-qid);
                                    tmp_phi = max(0.,min(0.5*(1.+tmp_theta),2.,2.*tmp_theta));
                                end
                                
                                %             ! Compute limited flux.
                                %             ! Since falling cloud water is a 1D advection problem, this
                                %             ! flux-limited advection scheme is monotonic.
                                fluxc = -vt_cl*(qic - 0.5*(1.-coef*vt_cl)*tmp_phi*(qic-qid));
                            else
                                fluxc = 0.;
                            end
                            fz_resolved(k,j,i) = fluxi + fluxc;
                            fzt_resolved(k,j,i) = -(fac_cond+fac_fus)*fluxi - fac_cond*fluxc;
                        end
                    end
                end
                
                
                % { It seems that these modifications does not make dqp  more accurate -
                % need to check why.
                for k=max(1,kmin-2):kmax
                    % !   coef=dtn/(dz*adz(k)*rho(k))
                    for j=1:num_blocks_y
                        for i=1:num_blocks_x
                            coef=dtn/(dz*adz(k)*rho(k));
                            
                            %          ! The cloud ice increment is the difference of the fluxes.
                            dqi=coef*(fz_resolved(k,j,i)-fz_resolved(k+1,j,i));
                            %          ! Add this increment to both non-precipitating and total water.
                            qn_coarse(k,j,i) = qn_coarse(k,j,i) + dqi;
                            qt_coarse(k,j,i)  = qt_coarse(k,j,i)  + dqi;
                            
                            cloud_qt_tend_resolved(k,j,i) =  dqi;
                            %          ! Include this effect in the total moisture budget.
                            %          qifall(k) = qifall(k) + dqi
                            
                            %          ! The latent heat flux induced by the falling cloud ice enters
                            %          ! the liquid-ice static energy budget in the same way as the
                            %          ! precipitation.  Note: use latent heat of sublimation.
                            lat_heat  = coef*(fzt_resolved(k,j,i)-fzt_resolved(k+1,j,i));
                            cloud_lat_heat_resolved(k,j,i) = lat_heat;
                            %          ! Add divergence of latent heat flux to liquid-ice static energy.
                            if add_cloud_tfull_tend
                                tfull_coarse(k,j,i)  = tfull_coarse(k,j,i)  + lat_heat; % Yani Need to think when need to change also t not full
                            end
                            %          ! Add divergence to liquid-ice static energy budget.
                            %          tlatqi(k) = tlatqi(k) + lat_heat
                        end
                    end
                end
                
                cloud_qt_tend_resolved = cloud_qt_tend_resolved./dtn; %all tendencies are divided by dtn (multiplied in RF of fortran).
                cloud_lat_heat_resolved = cloud_lat_heat_resolved./dtn;
            end
            %%
            % Step 03:12: calculate microphysical tendency of qp
            % calculate microphysical tendency of qp
            %%precip_proc
            
            % from precip_proc.f90
            
            powr1 = (3 + b_rain) / 4.;
            powr2 = (5 + b_rain) / 8.;
            pows1 = (3 + b_snow) / 4.;
            pows2 = (5 + b_snow) / 8.;
            powg1 = (3 + b_grau) / 4.;
            powg2 = (5 + b_grau) / 8.;
            
            for i=1:num_blocks_x
                for j=1:num_blocks_y
                    for k=1:num_z
                        
                        if (qn_coarse(k,j,i)+qp_coarse(k,j,i)>0.)
                            
                            omn = max(0.,min(1.,(tabs_coarse(k,j,i)-tbgmin)*a_bg));
                            omp = max(0.,min(1.,(tabs_coarse(k,j,i)-tprmin)*a_pr));
                            omg = max(0.,min(1.,(tabs_coarse(k,j,i)-tgrmin)*a_gr));
                            qrr = qp_coarse(k,j,i) * omp;
                            qss = qp_coarse(k,j,i) * (1.-omp)*(1.-omg); %Eq. A11 in the paper of MARAT. ;
                            qgg = qp_coarse(k,j,i) * (1.-omp)*omg;
                            
                            if (qn_coarse(k,j,i)>0.) %-------     Autoconversion/accretion
                                
                                qcc = qn_coarse(k,j,i) * omn;
                                qii = qn_coarse(k,j,i) * (1.-omn);
                                
                                if (qcc > qcw0)
                                    autor = alphaelq;
                                else
                                    autor = 0.;
                                end
                                
                                if (qii > qci0)
                                    autos = betaelq*coefice(k);
                                else
                                    autos = 0.;
                                end
                                accrr = accrrc(k) * qrr^powr1;
                                tmp = qss^pows1;
                                accrcs = accrsc(k) * tmp;
                                accris = accrsi(k) * tmp;
                                tmp = qgg^powg1;
                                accrcg = accrgc(k) * tmp;
                                accrig = accrgi(k) * tmp;
                                qcc = (qcc+dtn*autor*qcw0)/(1.+dtn*(accrr+accrcs+accrcg+autor));
                                qii = (qii+dtn*autos*qci0)/(1.+dtn*(accris+accrig+autos));
                                
                                
                                dqp_resolved(k,j,i) = dtn *(accrr*qcc + autor*(qcc-qcw0)+ ...
                                    (accris+accrig)*qii + (accrcs+accrcg)*qcc + autos*(qii-qci0));
                                
                                dqp_resolved(k,j,i) = min(dqp_resolved(k,j,i),qn_coarse(k,j,i));
                                
                            elseif(qp_coarse(k,j,i)>qp_threshold && qn_coarse(k,j,i)==0.)  % evaporation
                                % I think that there is a missing condition!
                                % if(tabs_coarse(i,j,k).gt.tmin_evap) then !kzm limit evaporation to temperatures
                                % above tmin_evap tmin_evap = 0 I think in our case
                                
                                qsat = omn*qsatw(tabs_coarse(k,j,i),pres(k))+ ...
                                    (1.-omn)*qsati(tabs_coarse(k,j,i),pres(k));
                                dqp_resolved(k,j,i) = dtn * (evapr1(k) * sqrt(qrr) + ...
                                    evapr2(k) * qrr^powr2 + ...
                                    evaps1(k) * sqrt(qss) + ...
                                    evaps2(k) * qss^pows2 + ...
                                    evapg1(k) * sqrt(qgg) + ...
                                    evapg2(k) * qgg^powg2)* ...
                                    (qt_coarse(k,j,i) /qsat-1.);
                                dqp_resolved(k,j,i) = max(-0.5*qp_coarse(k,j,i),dqp_resolved(k,j,i));
                                
                            else % complete evaporation
                                
                                dqp_resolved(k,j,i) = -qp_coarse(k,j,i);
                                
                            end
                            
                        end
                        qp_coarse(k,j,i) = qp_coarse(k,j,i)  + dqp_resolved(k,j,i);
                        qp_coarse(k,j,i) = max(0,qp_coarse(k,j,i));% Avoid negative values for
                        %             precip - should use only if I really modify qp here... which
                        %             probably I should! For some reason the advection of qp give
                        %             negative qp values - shouldn't happen - to check !
                        
                    end  % for
                end % for
            end % for
            % get tendency in kg/kg/s
            dqp_resolved = dqp_resolved/dtn;
            omp_3d_res1  = max(0.,min(1.,(tabs_coarse-tprmin)*a_pr));
            fac_dqp_tfull_res1 = (fac_cond+(1.-omp_3d_res1)*fac_fus); %This is approximated because I assume that the temperature 
            dqp_t_tend_resolved = dqp_resolved .*fac_dqp_tfull_res1;

            
            if do_show_times, sprintf('who;e step 03 (coarse grained calc: %f', toc), end
            tic
            
            
            
            % Step 04 output to netcdf
            if do_show_times, sprintf('writing vars to netcdf'), end
            if do_show_times, toc, end
            
            %%
            % Step 03:13: calculate rh in the end of the time step - for
            % coarse grained tendencies and the resolved tendencies:
            
            qt_vert_adv_resolved_tend = zeros(size(qt_flux_z_resolved));
            qt_vert_adv_coarse_tend = zeros(size(qt_flux_z_resolved));
            
            if calc_advection_tend ==0
                t_flux_z_tend_resolved = zeros(num_z,num_blocks_y,num_blocks_x);
                t_flux_z_tend_coarse = zeros(num_z,num_blocks_y,num_blocks_x);
            end
                
            for k=1:num_z - 1
                irhoadzdz = 1./(adz(k).*dz.*rho(k));
                qt_vert_adv_resolved_tend(k,:,:) = -(qt_flux_z_resolved(k+1,:,:) - qt_flux_z_resolved(k,:,:))*irhoadzdz;
                qt_vert_adv_coarse_tend(k,:,:) = -(qt_flux_z_coarse(k+1,:,:) - qt_flux_z_coarse(k,:,:))*irhoadzdz;

                t_flux_z_tend_coarse(k,:,:) = - (t_flux_z_coarse(k+1,:,:) - t_flux_z_coarse(k,:,:))./(adz(k).*dz.*rho(k));
                t_flux_z_tend_resolved(k,:,:) = - (tfull_flux_z_resolved(k+1,:,:) - tfull_flux_z_resolved(k,:,:))./(adz(k).*dz.*rho(k));
            end
            
            k = num_z;
            irhoadzdz = 1./(adz(k).*dz.*rho(k));
            qt_vert_adv_resolved_tend(k,:,:) = -(0.0 - qt_flux_z_resolved(k,:,:)).*irhoadzdz;
            qt_vert_adv_coarse_tend(k,:,:) = -(0.0 - qt_flux_z_coarse(k,:,:)).*irhoadzdz;
            t_flux_z_tend_coarse(k,:,:) = - (0.0 - t_flux_z_coarse(k,:,:)).*irhoadzdz;
            t_flux_z_tend_resolved(k,:,:) = - (0.0 - t_flux_z_resolved(k,:,:)).*irhoadzdz;
            
            

            t_coarse_init_plus_coarse_tend = t_coarse_init + (t_flux_z_tend_coarse + cloud_lat_heat_coarse + Qrad_coarse/86400+ dqp_t_tend_coarse)*dtn; 
            t_coarse_init_plus_resolved_tend = t_coarse_init + (t_flux_z_tend_resolved + cloud_lat_heat_resolved)*dtn;
            
            %Note that I don't consider the effect of QRAD on q
            qt_coarse_init = qv_coarse_init+qn_coarse_init;
            qt_coarse_init_init_plus_coarse_tend = qt_coarse_init + (cloud_qt_tend_coarse - dqp_coarse + qt_vert_adv_coarse_tend)*dtn;
            qt_coarse_init_init_plus_resolved_tend = qt_coarse_init + (cloud_qt_tend_resolved + qt_vert_adv_resolved_tend)*dtn;
            
            %calculate tabs for the case we propagated with coarse grained
            %high res tendencies
            [~, tabs_resolved_end_coarse, qt_coarse_init_init_plus_coarse_tend, ~,~,~] = update_qn_tabs(tbgmax, tbgmin, tprmax, tprmin, fac_cond, fac_fus, fac_sub, tgrmax, tgrmin, qp_threshold, num_z, ...
    num_blocks_x, num_blocks_y, qn_coarse, qt_coarse_init_init_plus_coarse_tend, zeros(size(qp_coarse)), t_coarse_init_plus_coarse_tend, gamaz, pres);
            
%calculate tabs for the case we propagated with coarse grained
            %high res tendencies
            [~, tabs_resolved_end_resolved, qt_coarse_init_init_plus_resolved_tend, ~,~,~] = update_qn_tabs(tbgmax, tbgmin, tprmax, tprmin, fac_cond, fac_fus, fac_sub, tgrmax, tgrmin, qp_threshold, num_z, ...
    num_blocks_x, num_blocks_y, qn_coarse, qt_coarse_init_init_plus_resolved_tend, zeros(size(qp_coarse)), t_coarse_init_plus_resolved_tend, gamaz, pres);
            
            %Calcualting relative humidity at end of time step - for coarse
            %and resolved tendencies
            omn_3d_coarse  = max(0.,min(1.,(tabs_resolved_end_coarse-tbgmin)*a_bg));
            qsat_resolved_end_coarse = zeros(size(omn_3d_coarse));
            for k = 1:num_z
                qsat_resolved_end_coarse(k,:,:) = omn_3d_coarse(k,:,:).*qsatw(tabs_resolved_end_coarse(k,:,:),pres(k))+ ...
                    (1.-omn_3d_coarse(k,:,:)).*qsati(tabs_resolved_end_coarse(k,:,:),pres(k));
            end
            rh_end_coarse  = qt_coarse_init_init_plus_coarse_tend./qsat_resolved_end_coarse;
            
            omn_3d_coarse  = max(0.,min(1.,(tabs_resolved_end_resolved-tbgmin)*a_bg));
            qsat_resolved_end_resolved = zeros(size(omn_3d_coarse));
            for k = 1:num_z
                qsat_resolved_end_resolved(k,:,:) = omn_3d_coarse(k,:,:).*qsatw(tabs_resolved_end_resolved(k,:,:),pres(k))+ ...
                    (1.-omn_3d_coarse(k,:,:)).*qsati(tabs_resolved_end_resolved(k,:,:),pres(k));
            end
            rh_end_resolved  = qt_coarse_init_init_plus_resolved_tend./qsat_resolved_end_resolved;

            
            %%
            %Step 04:01 calculate resituals
            
            dqp_resolved = dqp_coarse - dqp_resolved; %getting the residuals...
            
            %advection
            tfull_flux_x_resolved = tfull_flux_x_coarse - tfull_flux_x_resolved;
            tfull_flux_y_resolved = tfull_flux_y_coarse-tfull_flux_y_resolved;
            tfull_flux_z_resolved = tfull_flux_z_coarse-tfull_flux_z_resolved;
            
            t_flux_x_resolved = t_flux_x_coarse - t_flux_x_resolved;
            t_flux_y_resolved = t_flux_y_coarse-t_flux_y_resolved;
            t_flux_z_resolved = t_flux_z_coarse-t_flux_z_resolved;
            
            qt_flux_x_resolved = qt_flux_x_coarse - qt_flux_x_resolved;
            qt_flux_y_resolved = qt_flux_y_coarse-qt_flux_y_resolved;
            qt_flux_z_resolved = qt_flux_z_coarse-qt_flux_z_resolved;
            
            qp_flux_x_resolved = qp_flux_x_coarse - qp_flux_x_resolved;
            qp_flux_y_resolved = qp_flux_y_coarse-qp_flux_y_resolved;
            qp_flux_z_resolved = qp_flux_z_coarse-qp_flux_z_resolved;
            if calc_advection_tend
                tfull_flux_x_tend_resolved = tfull_flux_x_tend_coarse - tfull_flux_x_tend_resolved;
                tfull_flux_y_tend_resolved = tfull_flux_y_tend_coarse - tfull_flux_y_tend_resolved;
                tfull_flux_z_tend_resolved = tfull_flux_z_tend_coarse - tfull_flux_z_tend_resolved;
                
                
            end
            %diffusion fluxes and tendencies
            
            t_diff_flx_x_resolved = t_diff_flx_x_coarse - t_diff_flx_x_resolved;
            t_diff_flx_y_resolved = t_diff_flx_y_coarse - t_diff_flx_y_resolved;
            t_diff_flx_z_resolved = t_diff_flx_z_coarse - t_diff_flx_z_resolved;
            
            tfull_diff_flx_x_resolved = tfull_diff_flx_x_coarse - tfull_diff_flx_x_resolved;
            tfull_diff_flx_y_resolved = tfull_diff_flx_y_coarse - tfull_diff_flx_y_resolved;
            tfull_diff_flx_z_resolved = tfull_diff_flx_z_coarse - tfull_diff_flx_z_resolved;
            
            qt_diff_flx_x_resolved = qt_diff_flx_x_coarse - qt_diff_flx_x_resolved;
            qt_diff_flx_y_resolved = qt_diff_flx_y_coarse - qt_diff_flx_y_resolved;
            qt_diff_flx_z_resolved = qt_diff_flx_z_coarse - qt_diff_flx_z_resolved;
            
            qp_diff_flx_x_resolved = qp_diff_flx_x_coarse - qp_diff_flx_x_resolved;
            qp_diff_flx_y_resolved = qp_diff_flx_y_coarse - qp_diff_flx_y_resolved;
            qp_diff_flx_z_resolved = qp_diff_flx_z_coarse - qp_diff_flx_z_resolved;
            if calc_diffusive_tend
                t_diff_x_tend_resolved = t_diff_x_tend_coarse - t_diff_x_tend_resolved;
                t_diff_y_tend_resolved = t_diff_y_tend_coarse - t_diff_y_tend_resolved;
                t_diff_z_tend_resolved = t_diff_z_tend_coarse - t_diff_z_tend_resolved;
                
                tfull_diff_x_tend_resolved = tfull_diff_x_tend_coarse - tfull_diff_x_tend_resolved;
                tfull_diff_y_tend_resolved = tfull_diff_y_tend_coarse - tfull_diff_y_tend_resolved;
                tfull_diff_z_tend_resolved = tfull_diff_z_tend_coarse - tfull_diff_z_tend_resolved;
                
                qt_diff_x_tend_resolved = qt_diff_x_tend_coarse - qt_diff_x_tend_resolved;
                qt_diff_y_tend_resolved = qt_diff_y_tend_coarse - qt_diff_y_tend_resolved;
                qt_diff_z_tend_resolved = qt_diff_z_tend_coarse - qt_diff_z_tend_resolved;
                
                qp_diff_x_tend_resolved = qp_diff_x_tend_coarse - qp_diff_x_tend_resolved;
                qp_diff_y_tend_resolved = qp_diff_y_tend_coarse - qp_diff_y_tend_resolved;
                qp_diff_z_tend_resolved = qp_diff_z_tend_coarse - qp_diff_z_tend_resolved;
            end
            dqp_fall_resolved = dqp_fall_coarse - dqp_fall_resolved;
            t_fall_tend_resolved = t_fall_tend_coarse - t_fall_tend_resolved;
            
            %precip
            %         precip_resolved = precip_coarse - precip_resolved;
            %         precip_energy_resolved = precip_energy_coarse - precip_energy_resolved;
            
            %cloud
            cloud_lat_heat_resolved = cloud_lat_heat_coarse - cloud_lat_heat_resolved;
            cloud_qt_tend_resolved = cloud_qt_tend_coarse - cloud_qt_tend_resolved;
            
            %cloud fluxes
            fzt_resolved = fzt_coarse -fzt_resolved;
            fz_resolved = fz_coarse -fz_resolved;
            
            %diffusivity
            tkz_out_resolved = tkz_out_coarse -tkz_out_resolved;
            Pr1_resolved = Pr1_coarse -Pr1_resolved;
            
            windspeed_resolved = windspeed_coarse - windspeed_resolved;
            
            %momentum vertical advection
            dudt_advect_resolved= dudt_advect_coarse - dudt_advect_resolved;
            dvdt_advect_resolved= dvdt_advect_coarse - dvdt_advect_resolved;
            dwdt_advect_resolved= dwdt_advect_coarse - dwdt_advect_resolved;
            
            %Momentum Surface wind fluxes
            fluxbu_resolved = fluxbu_coarse - fluxbu_resolved;
            fluxbv_resolved = fluxbv_coarse - fluxbv_resolved;
            
            fluxbu_samson_resolved = fluxbu_samson_coarse - fluxbu_samson_resolved;
            fluxbv_samson_resolved = fluxbv_samson_coarse - fluxbv_samson_resolved;
            
            
            dqp_t_tend_resolved = dqp_t_tend_coarse - dqp_t_tend_resolved;
            
            
            %%
            %Step 04:02 permute fields
            % reverse order so that follows convention when output
            
            perm = [3 2 1];
            u_coarse_flip = permute(u_coarse, perm);
            v_coarse_flip = permute(v_coarse, perm);
            w_coarse_flip = permute(w_coarse, perm);
            tabs_coarse_flip = permute(tabs_coarse, perm);
            t_coarse_flip = permute(t_coarse, perm);
            tfull_coarse_flip = permute(tfull_coarse, perm);
            tfull_coarse_init_flip = permute(tfull_coarse_init, perm);
            qv_coarse_flip = permute(qv_coarse, perm);
            qn_coarse_flip = permute(qn_coarse, perm);
            qp_coarse_flip = permute(qp_coarse, perm);
            Qrad_coarse_flip = permute(Qrad_coarse, perm);
            
            dqp_coarse_flip = permute(dqp_coarse, perm);
            dqp_resolved_flip = permute(dqp_resolved, perm);
            precip_coarse_flip = permute(precip_coarse, perm);
            precip_resolved_flip = permute(precip_resolved, perm);
            precip_energy_coarse_flip = permute(precip_energy_coarse, perm);
            precip_energy_resolved_flip = permute(precip_energy_resolved, perm);
            
            
            %advection
            tfull_flux_x_resolved_flip= permute(tfull_flux_x_resolved, perm);
            tfull_flux_y_resolved_flip = permute(tfull_flux_y_resolved, perm);
            tfull_flux_z_resolved_flip = permute(tfull_flux_z_resolved, perm);
            
            t_flux_x_resolved_flip = permute(t_flux_x_resolved, perm);
            t_flux_y_resolved_flip = permute(t_flux_y_resolved, perm);
            t_flux_z_resolved_flip = permute(t_flux_z_resolved, perm);
            
            qt_flux_x_resolved_flip = permute(qt_flux_x_resolved, perm);
            qt_flux_y_resolved_flip = permute(qt_flux_y_resolved, perm);
            qt_flux_z_resolved_flip = permute(qt_flux_z_resolved, perm);
            
            qp_flux_x_resolved_flip = permute(qp_flux_x_resolved, perm);
            qp_flux_y_resolved_flip = permute(qp_flux_y_resolved, perm);
            qp_flux_z_resolved_flip = permute(qp_flux_z_resolved, perm);
            
            %Also using the coarse advection - (necessary for one of the tendency terms)
            tfull_flux_x_coarse_flip= permute(tfull_flux_x_coarse, perm);
            tfull_flux_y_coarse_flip = permute(tfull_flux_y_coarse, perm);
            tfull_flux_z_coarse_flip = permute(tfull_flux_z_coarse, perm);
            
            t_flux_x_coarse_flip = permute(t_flux_x_coarse, perm);
            t_flux_y_coarse_flip = permute(t_flux_y_coarse, perm);
            t_flux_z_coarse_flip = permute(t_flux_z_coarse, perm);
            
            qt_flux_x_coarse_flip = permute(qt_flux_x_coarse, perm);
            qt_flux_y_coarse_flip = permute(qt_flux_y_coarse, perm);
            qt_flux_z_coarse_flip = permute(qt_flux_z_coarse, perm);
            
            qp_flux_x_coarse_flip = permute(qp_flux_x_coarse, perm);
            qp_flux_y_coarse_flip = permute(qp_flux_y_coarse, perm);
            qp_flux_z_coarse_flip = permute(qp_flux_z_coarse, perm);
            
            %advection tendency
            if calc_advection_tend
                tfull_flux_x_tend_resolved_flip = permute(tfull_flux_x_tend_resolved, perm);
                tfull_flux_y_tend_resolved_flip = permute(tfull_flux_y_tend_resolved, perm);
                tfull_flux_z_tend_resolved_flip = permute(tfull_flux_z_tend_resolved, perm);
                
                t_flux_x_tend_resolved_flip = permute(t_flux_x_tend_resolved, perm);
                t_flux_y_tend_resolved_flip = permute(t_flux_y_tend_resolved, perm);
                t_flux_z_tend_resolved_flip = permute(t_flux_z_tend_resolved, perm);
            end
            %diffusion fluxes and tendencies
            
            t_diff_flx_x_resolved_flip = permute(t_diff_flx_x_resolved, perm);
            t_diff_flx_y_resolved_flip = permute(t_diff_flx_y_resolved, perm);
            t_diff_flx_z_resolved_flip = permute(t_diff_flx_z_resolved, perm);
            
            tfull_diff_flx_x_resolved_flip = permute(tfull_diff_flx_x_resolved, perm);
            tfull_diff_flx_y_resolved_flip = permute(tfull_diff_flx_y_resolved, perm);
            tfull_diff_flx_z_resolved_flip = permute(tfull_diff_flx_z_resolved, perm);
            
            qt_diff_flx_x_resolved_flip = permute(qt_diff_flx_x_resolved, perm);
            qt_diff_flx_y_resolved_flip = permute(qt_diff_flx_y_resolved, perm);
            qt_diff_flx_z_resolved_flip = permute(qt_diff_flx_z_resolved, perm);
            
            qp_diff_flx_x_resolved_flip = permute(qp_diff_flx_x_resolved, perm);
            qp_diff_flx_y_resolved_flip = permute(qp_diff_flx_y_resolved, perm);
            qp_diff_flx_z_resolved_flip = permute(qp_diff_flx_z_resolved, perm);
            
            %Also using the coarse diffusion - (Undecided if it makes more sense than the coarse-resolved)
            t_diff_flx_x_coarse_flip = permute(t_diff_flx_x_coarse, perm);
            t_diff_flx_y_coarse_flip = permute(t_diff_flx_y_coarse, perm);
            t_diff_flx_z_coarse_flip = permute(t_diff_flx_z_coarse, perm);
            
            tfull_diff_flx_x_coarse_flip = permute(tfull_diff_flx_x_coarse, perm);
            tfull_diff_flx_y_coarse_flip = permute(tfull_diff_flx_y_coarse, perm);
            tfull_diff_flx_z_coarse_flip = permute(tfull_diff_flx_z_coarse, perm);
            
            qt_diff_flx_x_coarse_flip = permute(qt_diff_flx_x_coarse, perm);
            qt_diff_flx_y_coarse_flip = permute(qt_diff_flx_y_coarse, perm);
            qt_diff_flx_z_coarse_flip = permute(qt_diff_flx_z_coarse, perm);
            
            qp_diff_flx_x_coarse_flip = permute(qp_diff_flx_x_coarse, perm);
            qp_diff_flx_y_coarse_flip = permute(qp_diff_flx_y_coarse, perm);
            qp_diff_flx_z_coarse_flip = permute(qp_diff_flx_z_coarse, perm);
            
            
            if calc_diffusive_tend
                t_diff_x_tend_resolved_flip = permute(t_diff_x_tend_resolved, perm);
                t_diff_y_tend_resolved_flip = permute(t_diff_y_tend_resolved, perm);
                t_diff_z_tend_resolved_flip = permute(t_diff_z_tend_resolved, perm);
                
                tfull_diff_x_tend_resolved_flip = permute(tfull_diff_x_tend_resolved, perm);
                tfull_diff_y_tend_resolved_flip = permute(tfull_diff_y_tend_resolved, perm);
                tfull_diff_z_tend_resolved_flip = permute(tfull_diff_z_tend_resolved, perm);
                
                qt_diff_x_tend_resolved_flip = permute(qt_diff_x_tend_resolved, perm);
                qt_diff_y_tend_resolved_flip = permute(qt_diff_y_tend_resolved, perm);
                qt_diff_z_tend_resolved_flip = permute(qt_diff_z_tend_resolved, perm);
                
                qp_diff_x_tend_resolved_flip = permute(qp_diff_x_tend_resolved, perm);
                qp_diff_y_tend_resolved_flip = permute(qp_diff_y_tend_resolved, perm);
                qp_diff_z_tend_resolved_flip = permute(qp_diff_z_tend_resolved, perm);
            end
            % precip fall tedndencies
            dqp_fall_resolved_flip = permute(dqp_fall_resolved, perm);
            t_fall_tend_resolved_flip = permute(t_fall_tend_resolved, perm);
            
            dqp_fall_coarse_flip = permute(dqp_fall_coarse, perm);
            t_fall_tend_coarse_flip = permute(t_fall_tend_coarse, perm);
            
            %cloud
            cloud_lat_heat_resolved_flip = permute(cloud_lat_heat_resolved, perm);
            cloud_qt_tend_resolved_flip = permute(cloud_qt_tend_resolved, perm);
            
            cloud_lat_heat_coarse_flip = permute(cloud_lat_heat_coarse, perm);
            cloud_qt_tend_coarse_flip = permute(cloud_qt_tend_coarse, perm);
            
            fzt_resolved_flip = permute(fzt_resolved, perm);
            fz_resolved_flip = permute(fz_resolved, perm);
            
            fzt_coarse_flip = permute(fzt_coarse, perm);
            fz_coarse_flip = permute(fz_coarse, perm);
            
            %difusivity
            tkz_out_coarse_flip = permute(tkz_out_coarse, perm);
            tkz_out_resolved_flip = permute(tkz_out_resolved, perm);
            
            Pr1_coarse_flip = permute(Pr1_coarse, perm);
            Pr1_resolved_flip = permute(Pr1_resolved, perm);
            
            
            qp_coarse_init_flip  =  permute(qp_coarse_init, perm);
            qn_coarse_init_flip  =  permute(qn_coarse_init, perm);
            qv_coarse_init_flip  =  permute(qv_coarse_init, perm);
            tabs_resolved_init_flip  =  permute(tabs_resolved, perm);
            
            perm_2d = [2,1];
            windspeed_coarse_flip = permute(windspeed_coarse, perm_2d);
            windspeed_resolved_flip = permute(windspeed_resolved, perm_2d);
            
            

            %momentum vertical advection:
            dudt_advect_coarse_flip  =  permute(dudt_advect_coarse, perm);
            dvdt_advect_coarse_flip  =  permute(dvdt_advect_coarse, perm);
            dwdt_advect_coarse_flip  =  permute(dwdt_advect_coarse, perm);
            
            dudt_advect_resolved_flip  =  permute(dudt_advect_resolved, perm);
            dvdt_advect_resolved_flip  =  permute(dvdt_advect_resolved, perm);
            dwdt_advect_resolved_flip  =  permute(dwdt_advect_resolved, perm);
            %Momentum Surface wind fluxes
            fluxbu_coarse_flip = permute(fluxbu_coarse, perm_2d);
            fluxbv_coarse_flip = permute(fluxbv_coarse, perm_2d);
            fluxbu_resolved_flip = permute(fluxbu_resolved, perm_2d);
            fluxbv_resolved_flip = permute(fluxbv_resolved, perm_2d);
            
            fluxbu_samson_coarse_flip = permute(fluxbu_samson_coarse, perm_2d);
            fluxbv_samson_coarse_flip = permute(fluxbv_samson_coarse, perm_2d);
            fluxbu_samson_resolved_flip = permute(fluxbu_samson_resolved, perm_2d);
            fluxbv_samson_resolved_flip = permute(fluxbv_samson_resolved, perm_2d);
            % change place? 
            dqp_t_tend_resolved_flip =  permute(dqp_t_tend_resolved, perm);
            dqp_t_tend_coarse_flip =  permute(dqp_t_tend_coarse, perm);
            % Relative humedity and generalization tendencies:
            rh_init_flip =  permute(rh_init, perm);
            rh_end_resolved_flip =  permute(rh_end_resolved, perm);
            rh_end_coarse_flip =  permute(rh_end_coarse, perm);
            tabs_resolved_end_resolved_flip =  permute(tabs_resolved_end_resolved, perm);
            tabs_resolved_end_coarse_flip =  permute(tabs_resolved_end_coarse, perm);
            %%
            %Step 04:03 write to netcdf
            
            
            
                    if is_cheyenne
                        outfilename_janni_final =['/glade/scratch/janniy/ML_convection_data/', exper{i_exper}, '/', exper{i_exper},'km12x576_576x1440x48_ctl_288_', num2str(time_index, '%010d'), '_000', num2str(cycle), '_diff_coarse_space_corrected_tkz', num2str(multiple_space), '.nc4'];
                    else
                        outfilename_janni_final =['/net/aimsir/archive1/janniy/ML_convection_data/', exper{i_exper}, 'km12x576_576x1440x48_ctl_288_', num2str(time_index, '%010d'), '_000', num2str(cycle), '_diff_coarse_space_corrected', num2str(multiple_space), '.nc4'];
                    end
            
            disp(outfilename_janni_final)
            % create netcdf
            ncid = netcdf.create(outfilename_janni_final,'NETCDF4');
            % define dimensions
            xdimid = netcdf.defDim(ncid,'x',length(x_coarse));
            ydimid = netcdf.defDim(ncid,'y',length(y_coarse));
            zdimid = netcdf.defDim(ncid,'z',length(z));
            
            % define variables
            x_varid = netcdf.defVar(ncid,'x','NC_FLOAT',xdimid);
            y_varid = netcdf.defVar(ncid,'y','NC_FLOAT',ydimid);
            z_varid = netcdf.defVar(ncid,'z','NC_FLOAT',zdimid);
            p_varid = netcdf.defVar(ncid,'p','NC_FLOAT',zdimid);
            rho_varid = netcdf.defVar(ncid,'rho','NC_FLOAT',zdimid);
            
            u_varid = netcdf.defVar(ncid,'U','NC_FLOAT',[xdimid ydimid zdimid]);
            v_varid = netcdf.defVar(ncid,'V','NC_FLOAT',[xdimid ydimid zdimid]);
            w_varid = netcdf.defVar(ncid,'W','NC_FLOAT',[xdimid ydimid zdimid]);
            tabs_varid = netcdf.defVar(ncid,'TABS','NC_FLOAT',[xdimid ydimid zdimid]);
            

            t_varid = netcdf.defVar(ncid,'T','NC_FLOAT',[xdimid ydimid zdimid]);
            tfull_varid = netcdf.defVar(ncid,'TFULL','NC_FLOAT',[xdimid ydimid zdimid]);
            tfull_init_varid = netcdf.defVar(ncid,'TFULL_INIT','NC_FLOAT',[xdimid ydimid zdimid]);
            qv_varid = netcdf.defVar(ncid,'Q','NC_FLOAT',[xdimid ydimid zdimid]);
            qn_varid = netcdf.defVar(ncid,'QN','NC_FLOAT',[xdimid ydimid zdimid]);
            qp_varid = netcdf.defVar(ncid,'QP','NC_FLOAT',[xdimid ydimid zdimid]);
            Qrad_varid = netcdf.defVar(ncid,'QRAD','NC_FLOAT',[xdimid ydimid zdimid]);
            
            tfull_flux_x_resolved_varid = netcdf.defVar(ncid,'TFULL_FLUX_X','NC_FLOAT',[xdimid ydimid zdimid]);
            tfull_flux_y_resolved_varid = netcdf.defVar(ncid,'TFULL_FLUX_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            tfull_flux_z_resolved_varid = netcdf.defVar(ncid,'TFULL_FLUX_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            t_flux_x_resolved_varid = netcdf.defVar(ncid,'T_FLUX_X','NC_FLOAT',[xdimid ydimid zdimid]);
            t_flux_y_resolved_varid = netcdf.defVar(ncid,'T_FLUX_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            t_flux_z_resolved_varid = netcdf.defVar(ncid,'T_FLUX_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            qt_flux_x_resolved_varid = netcdf.defVar(ncid,'QT_FLUX_X','NC_FLOAT',[xdimid ydimid zdimid]);
            qt_flux_y_resolved_varid = netcdf.defVar(ncid,'QT_FLUX_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            qt_flux_z_resolved_varid = netcdf.defVar(ncid,'QT_FLUX_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            qp_flux_x_resolved_varid = netcdf.defVar(ncid,'QP_FLUX_X','NC_FLOAT',[xdimid ydimid zdimid]);
            qp_flux_y_resolved_varid = netcdf.defVar(ncid,'QP_FLUX_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            qp_flux_z_resolved_varid = netcdf.defVar(ncid,'QP_FLUX_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            
            %Also using the coarse advection - (necessary for one of the tendency terms)
            tfull_flux_x_coarse_varid = netcdf.defVar(ncid,'TFULL_FLUX_COARSE_X','NC_FLOAT',[xdimid ydimid zdimid]);
            tfull_flux_y_coarse_varid = netcdf.defVar(ncid,'TFULL_FLUX_COARSE_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            tfull_flux_z_coarse_varid = netcdf.defVar(ncid,'TFULL_FLUX_COARSE_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            t_flux_x_coarse_varid = netcdf.defVar(ncid,'T_FLUX_COARSE_X','NC_FLOAT',[xdimid ydimid zdimid]);
            t_flux_y_coarse_varid = netcdf.defVar(ncid,'T_FLUX_COARSE_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            t_flux_z_coarse_varid = netcdf.defVar(ncid,'T_FLUX_COARSE_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            qt_flux_x_coarse_varid = netcdf.defVar(ncid,'QT_FLUX_COARSE_X','NC_FLOAT',[xdimid ydimid zdimid]);
            qt_flux_y_coarse_varid = netcdf.defVar(ncid,'QT_FLUX_COARSE_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            qt_flux_z_coarse_varid = netcdf.defVar(ncid,'QT_FLUX_COARSE_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            qp_flux_x_coarse_varid = netcdf.defVar(ncid,'QP_FLUX_COARSE_X','NC_FLOAT',[xdimid ydimid zdimid]);
            qp_flux_y_coarse_varid = netcdf.defVar(ncid,'QP_FLUX_COARSE_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            qp_flux_z_coarse_varid = netcdf.defVar(ncid,'QP_FLUX_COARSE_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            
            if calc_advection_tend
                tfull_flux_x_tend_resolved_varid = netcdf.defVar(ncid,'TFULL_ADV_TEND_X','NC_FLOAT',[xdimid ydimid zdimid]);
                tfull_flux_y_tend_resolved_varid = netcdf.defVar(ncid,'TFULL_ADV_TEND_Y','NC_FLOAT',[xdimid ydimid zdimid]);
                tfull_flux_z_tend_resolved_varid = netcdf.defVar(ncid,'TFULL_ADV_TEND_Z','NC_FLOAT',[xdimid ydimid zdimid]);
                
                t_flux_x_tend_resolved_varid = netcdf.defVar(ncid,'T_ADV_TEND_X','NC_FLOAT',[xdimid ydimid zdimid]);
                t_flux_y_tend_resolved_varid = netcdf.defVar(ncid,'T_ADV_TEND_Y','NC_FLOAT',[xdimid ydimid zdimid]);
                t_flux_z_tend_resolved_varid = netcdf.defVar(ncid,'T_ADV_TEND_Z','NC_FLOAT',[xdimid ydimid zdimid]);
                
            end
            t_diff_flx_x_resolved_varid = netcdf.defVar(ncid,'T_DIFF_FLUX_X','NC_FLOAT',[xdimid ydimid zdimid]);
            t_diff_flx_y_resolved_varid = netcdf.defVar(ncid,'T_DIFF_FLUX_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            t_diff_flx_z_resolved_varid = netcdf.defVar(ncid,'T_DIFF_FLUX_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            tfull_diff_flx_x_resolved_varid = netcdf.defVar(ncid,'TFULL_DIFF_FLUX_X','NC_FLOAT',[xdimid ydimid zdimid]);
            tfull_diff_flx_y_resolved_varid = netcdf.defVar(ncid,'TFULL_DIFF_FLUX_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            tfull_diff_flx_z_resolved_varid = netcdf.defVar(ncid,'TFULL_DIFF_FLUX_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            qt_diff_flx_x_resolved_varid = netcdf.defVar(ncid,'QT_DIFF_FLUX_X','NC_FLOAT',[xdimid ydimid zdimid]);
            qt_diff_flx_y_resolved_varid = netcdf.defVar(ncid,'QT_DIFF_FLUX_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            qt_diff_flx_z_resolved_varid = netcdf.defVar(ncid,'QT_DIFF_FLUX_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            qp_diff_flx_x_resolved_varid = netcdf.defVar(ncid,'QP_DIFF_FLUX_X','NC_FLOAT',[xdimid ydimid zdimid]);
            qp_diff_flx_y_resolved_varid = netcdf.defVar(ncid,'QP_DIFF_FLUX_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            qp_diff_flx_z_resolved_varid = netcdf.defVar(ncid,'QP_DIFF_FLUX_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            
            
            %Also using the coarse diffusion - (Undecided if it makes more sense than the coarse-resolved)
            t_diff_flx_x_coarse_varid = netcdf.defVar(ncid,'T_DIFF_F_COARSE_X','NC_FLOAT',[xdimid ydimid zdimid]);
            t_diff_flx_y_coarse_varid = netcdf.defVar(ncid,'T_DIFF_F_COARSE_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            t_diff_flx_z_coarse_varid = netcdf.defVar(ncid,'T_DIFF_F_COARSE_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            tfull_diff_flx_x_coarse_varid = netcdf.defVar(ncid,'TFULL_DIFF_F_COARSE_X','NC_FLOAT',[xdimid ydimid zdimid]);
            tfull_diff_flx_y_coarse_varid = netcdf.defVar(ncid,'TFULL_DIFF_F_COARSE_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            tfull_diff_flx_z_coarse_varid = netcdf.defVar(ncid,'TFULL_DIFF_F_COARSE_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            qt_diff_flx_x_coarse_varid = netcdf.defVar(ncid,'QT_DIFF_F_COARSE_X','NC_FLOAT',[xdimid ydimid zdimid]);
            qt_diff_flx_y_coarse_varid = netcdf.defVar(ncid,'QT_DIFF_F_COARSE_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            qt_diff_flx_z_coarse_varid = netcdf.defVar(ncid,'QT_DIFF_F_COARSE_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            qp_diff_flx_x_coarse_varid = netcdf.defVar(ncid,'QP_DIFF_F_COARSE_X','NC_FLOAT',[xdimid ydimid zdimid]);
            qp_diff_flx_y_coarse_varid = netcdf.defVar(ncid,'QP_DIFF_F_COARSE_Y','NC_FLOAT',[xdimid ydimid zdimid]);
            qp_diff_flx_z_coarse_varid = netcdf.defVar(ncid,'QP_DIFF_F_COARSE_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            
            if calc_diffusive_tend
                tfull_diff_x_tend_resolved_varid = netcdf.defVar(ncid,'TFULL_DIFF_TEND_X','NC_FLOAT',[xdimid ydimid zdimid]);
                tfull_diff_y_tend_resolved_varid = netcdf.defVar(ncid,'TFULL_DIFF_TEND_Y','NC_FLOAT',[xdimid ydimid zdimid]);
                tfull_diff_z_tend_resolved_varid = netcdf.defVar(ncid,'TFULL_DIFF_TEND_Z','NC_FLOAT',[xdimid ydimid zdimid]);
                
                t_diff_x_tend_resolved_varid = netcdf.defVar(ncid,'T_DIFF_TEND_X','NC_FLOAT',[xdimid ydimid zdimid]);
                t_diff_y_tend_resolved_varid = netcdf.defVar(ncid,'T_DIFF_TEND_Y','NC_FLOAT',[xdimid ydimid zdimid]);
                t_diff_z_tend_resolved_varid = netcdf.defVar(ncid,'T_DIFF_TEND_Z','NC_FLOAT',[xdimid ydimid zdimid]);
                
                qt_diff_x_tend_resolved_varid = netcdf.defVar(ncid,'QT_DIFF_TEND_X','NC_FLOAT',[xdimid ydimid zdimid]);
                qt_diff_y_tend_resolved_varid = netcdf.defVar(ncid,'QT_DIFF_TEND_Y','NC_FLOAT',[xdimid ydimid zdimid]);
                qt_diff_z_tend_resolved_varid = netcdf.defVar(ncid,'QT_DIFF_TEND_Z','NC_FLOAT',[xdimid ydimid zdimid]);
                
                qp_diff_x_tend_resolved_varid = netcdf.defVar(ncid,'QP_DIFF_TEND_X','NC_FLOAT',[xdimid ydimid zdimid]);
                qp_diff_y_tend_resolved_varid = netcdf.defVar(ncid,'QP_DIFF_TEND_Y','NC_FLOAT',[xdimid ydimid zdimid]);
                qp_diff_z_tend_resolved_varid = netcdf.defVar(ncid,'QP_DIFF_TEND_Z','NC_FLOAT',[xdimid ydimid zdimid]);
            end
            % precip fall tedndencies
            dqp_fall_resolved_varid  = netcdf.defVar(ncid,'DQP_FALL_RES','NC_FLOAT',[xdimid ydimid zdimid]);
            dqp_fall_coarse_varid  = netcdf.defVar(ncid,'DQP_FALL_COARSE','NC_FLOAT',[xdimid ydimid zdimid]);
            
            t_fall_tend_resolved_varid  = netcdf.defVar(ncid,'T_FALL_RES','NC_FLOAT',[xdimid ydimid zdimid]);
            t_fall_tend_coarse_varid  = netcdf.defVar(ncid,'T_FALL_COARSE','NC_FLOAT',[xdimid ydimid zdimid]);
            
            %cloud
            cloud_lat_heat_resolved_varid  = netcdf.defVar(ncid,'LAT_HEAT_CLOUD_RES','NC_FLOAT',[xdimid ydimid zdimid]);
            cloud_lat_heat_coarse_varid  = netcdf.defVar(ncid,'LAT_HEAT_CLOUD_COARSE','NC_FLOAT',[xdimid ydimid zdimid]);
            
            cloud_qt_tend_resolved_varid  =netcdf.defVar(ncid,'QT_TEND_CLOUD_RES','NC_FLOAT',[xdimid ydimid zdimid]);
            cloud_qt_tend_coarse_varid  =netcdf.defVar(ncid,'QT_TEND_CLOUD_COARSE','NC_FLOAT',[xdimid ydimid zdimid]);
            
            %cloud flux
            fzt_resolved_varid  = netcdf.defVar(ncid,'CLOUD_FZT_FLUX_RES','NC_FLOAT',[xdimid ydimid zdimid]);
            fz_resolved_varid  = netcdf.defVar(ncid,'CLOUD_FZ_FLUX_RES','NC_FLOAT',[xdimid ydimid zdimid]);
            
            fzt_coarse_varid  = netcdf.defVar(ncid,'CLOUD_FZT_FLUX_COARSE','NC_FLOAT',[xdimid ydimid zdimid]);
            fz_coarse_varid  = netcdf.defVar(ncid,'CLOUD_FZ_FLUX_COARSE','NC_FLOAT',[xdimid ydimid zdimid]);
            
            %diffusivity
            tkz_out_resolved_varid  = netcdf.defVar(ncid,'TKZ_RES','NC_FLOAT',[xdimid ydimid zdimid]);
            tkz_out_coarse_varid  = netcdf.defVar(ncid,'TKZ_COARSE','NC_FLOAT',[xdimid ydimid zdimid]);
            
            Pr1_resolved_varid  = netcdf.defVar(ncid,'PR1_RES','NC_FLOAT',[xdimid ydimid zdimid]);
            Pr1_coarse_varid  = netcdf.defVar(ncid,'PR1_COARSE','NC_FLOAT',[xdimid ydimid zdimid]);
            
            
            
            
            dqp_varid = netcdf.defVar(ncid,'DQP','NC_FLOAT',[xdimid ydimid zdimid]);
            dqp_resolved_varid = netcdf.defVar(ncid,'DQP_RESOLVED','NC_FLOAT',[xdimid ydimid zdimid]);
            precip_varid = netcdf.defVar(ncid,'PRECIP','NC_FLOAT',[xdimid ydimid zdimid]);
            precip_resolved_varid = netcdf.defVar(ncid,'PRECIP_RESOLVED','NC_FLOAT',[xdimid ydimid zdimid]);
            precip_energy_varid = netcdf.defVar(ncid,'PRECIP_ENERGY','NC_FLOAT',[xdimid ydimid zdimid]);
            precip_energy_resolved_varid = netcdf.defVar(ncid,'PRECIP_ENERGY_RESOLVED','NC_FLOAT',[xdimid ydimid zdimid]);
            
            
            qp_coarse_init_varid = netcdf.defVar(ncid,'QP_COARSE_INIT','NC_FLOAT',[xdimid ydimid zdimid]);
            qn_coarse_init_varid = netcdf.defVar(ncid,'QN_COARSE_INIT','NC_FLOAT',[xdimid ydimid zdimid]);
            qv_coarse_init_varid = netcdf.defVar(ncid,'QV_COARSE_INIT','NC_FLOAT',[xdimid ydimid zdimid]);
            tabs_resolved_init_varid = netcdf.defVar(ncid,'TABS_RESOLVED_INIT','NC_FLOAT',[xdimid ydimid zdimid]);
            windspeed_coarse_varid = netcdf.defVar(ncid,'WINDSPEED_S_COARSE','NC_FLOAT',[xdimid ydimid]);
            windspeed_resolved_varid = netcdf.defVar(ncid,'WINDSPEED_S_RESOLVED','NC_FLOAT',[xdimid ydimid]);
            
            %momentum vertical advection
            dudt_advect_coarse_varid = netcdf.defVar(ncid,'U_ADV_COARSE','NC_FLOAT',[xdimid ydimid zdimid]);
            dvdt_advect_coarse_varid = netcdf.defVar(ncid,'V_ADV_COARSE','NC_FLOAT',[xdimid ydimid zdimid]);
            dwdt_advect_coarse_varid = netcdf.defVar(ncid,'W_ADV_COARSE','NC_FLOAT',[xdimid ydimid zdimid]);
            
            dudt_advect_resolved_varid = netcdf.defVar(ncid,'U_ADV_RESOLVED','NC_FLOAT',[xdimid ydimid zdimid]);
            dvdt_advect_resolved_varid = netcdf.defVar(ncid,'V_ADV_RESOLVED','NC_FLOAT',[xdimid ydimid zdimid]);
            dwdt_advect_resolved_varid = netcdf.defVar(ncid,'W_ADV_RESOLVED','NC_FLOAT',[xdimid ydimid zdimid]);
            
            
%              %Momentum Surface wind fluxes
              fluxbv_coarse_varid = netcdf.defVar(ncid,'SFLUX_V_COARSE','NC_FLOAT',[xdimid ydimid]);
              fluxbu_coarse_varid = netcdf.defVar(ncid,'SFLUX_U_COARSE','NC_FLOAT',[xdimid ydimid]);
              fluxbv_resolved_varid = netcdf.defVar(ncid,'SFLUX_V_RESOLVED','NC_FLOAT',[xdimid ydimid]);
              fluxbu_resolved_varid = netcdf.defVar(ncid,'SFLUX_U_RESOLVED','NC_FLOAT',[xdimid ydimid]);
              
              fluxbv_samson_coarse_varid = netcdf.defVar(ncid,'SFLUX_V_COARSE_SAM','NC_FLOAT',[xdimid ydimid]);
              fluxbu_samson_coarse_varid = netcdf.defVar(ncid,'SFLUX_U_COARSE_SAM','NC_FLOAT',[xdimid ydimid]);
              fluxbv_samson_resolved_varid = netcdf.defVar(ncid,'SFLUX_V_RESOLVED_SAM','NC_FLOAT',[xdimid ydimid]);
              fluxbu_samson_resolved_varid = netcdf.defVar(ncid,'SFLUX_U_RESOLVED_SAM','NC_FLOAT',[xdimid ydimid]);
%             % change place? 
              % The direct calculation of DQP on T (HL) 
              dqp_t_tend_resolved_varid = netcdf.defVar(ncid,'DQP_T_TEND_RESOLVED','NC_FLOAT',[xdimid ydimid zdimid]);
              dqp_t_tend_coarse_varid = netcdf.defVar(ncid,'DQP_T_TEND_COARSE','NC_FLOAT',[xdimid ydimid zdimid]);
              %Generalization
              rh_init_varid = netcdf.defVar(ncid,'RH_INIT','NC_FLOAT',[xdimid ydimid zdimid]);
              rh_end_resolved_varid = netcdf.defVar(ncid,'RH_END_RESOLVED','NC_FLOAT',[xdimid ydimid zdimid]);
              rh_end_coarse_varid = netcdf.defVar(ncid,'RH_END_COARSE','NC_FLOAT',[xdimid ydimid zdimid]);
              tabs_resolved_end_resolved_varid = netcdf.defVar(ncid,'TABS_RES_END_RESOLVED','NC_FLOAT',[xdimid ydimid zdimid]);
              tabs_resolved_end_coarse_varid = netcdf.defVar(ncid,'TABS_RES_END_COARSE','NC_FLOAT',[xdimid ydimid zdimid]);
              
            % attribute
            netcdf.putAtt(ncid,x_varid,'units','m');
            netcdf.putAtt(ncid,y_varid,'units','m');
            netcdf.putAtt(ncid,z_varid,'units','m');
            netcdf.putAtt(ncid,p_varid,'units','hPa');
            netcdf.putAtt(ncid,rho_varid,'units','kg/m^3');
            
            netcdf.putAtt(ncid,u_varid,'units','m/s');
            netcdf.putAtt(ncid,v_varid,'units','m/s');
            netcdf.putAtt(ncid,w_varid,'units','m/s');
            netcdf.putAtt(ncid,tabs_varid,'units','K');
            netcdf.putAtt(ncid,t_varid,'units','K');
            netcdf.putAtt(ncid,tfull_varid,'units','K');
            netcdf.putAtt(ncid,tfull_init_varid,'units','K');
            netcdf.putAtt(ncid,qv_varid,'units','g/kg');
            netcdf.putAtt(ncid,qn_varid,'units','g/kg');
            netcdf.putAtt(ncid,qp_varid,'units','g/kg');
            netcdf.putAtt(ncid,Qrad_varid,'units','K/day');
            
            netcdf.putAtt(ncid,qv_coarse_init_varid,'units','g/kg');
            netcdf.putAtt(ncid,qn_coarse_init_varid,'units','g/kg');
            netcdf.putAtt(ncid,qp_coarse_init_varid,'units','g/kg');
            netcdf.putAtt(ncid,tabs_resolved_init_varid,'units','K');
            netcdf.putAtt(ncid,windspeed_coarse_varid,'units','m/s');
            netcdf.putAtt(ncid,windspeed_resolved_varid,'units','m/s');
            
           
            
            %advection
            netcdf.putAtt(ncid,tfull_flux_x_resolved_varid,'units','K kg/m^2/s');
            netcdf.putAtt(ncid,tfull_flux_y_resolved_varid,'units','K kg/m^2/s');
            netcdf.putAtt(ncid,tfull_flux_z_resolved_varid,'units','K kg/m^2/s');
            
            netcdf.putAtt(ncid,t_flux_x_resolved_varid,'units','K kg/m^2/s');
            netcdf.putAtt(ncid,t_flux_y_resolved_varid,'units','K kg/m^2/s');
            netcdf.putAtt(ncid,t_flux_z_resolved_varid,'units','K kg/m^2/s');
            
            netcdf.putAtt(ncid,qt_flux_x_resolved_varid,'units','g/m^2/s');
            netcdf.putAtt(ncid,qt_flux_y_resolved_varid,'units','g/m^2/s');
            netcdf.putAtt(ncid,qt_flux_z_resolved_varid,'units','g/m^2/s');
            
            netcdf.putAtt(ncid,qp_flux_x_resolved_varid,'units','g/m^2/s');
            netcdf.putAtt(ncid,qp_flux_y_resolved_varid,'units','g/m^2/s');
            netcdf.putAtt(ncid,qp_flux_z_resolved_varid,'units','g/m^2/s');
            
            %coarse advection
            netcdf.putAtt(ncid,tfull_flux_x_coarse_varid,'units','K kg/m^2/s');
            netcdf.putAtt(ncid,tfull_flux_y_coarse_varid,'units','K kg/m^2/s');
            netcdf.putAtt(ncid,tfull_flux_z_coarse_varid,'units','K kg/m^2/s');
            
            netcdf.putAtt(ncid,t_flux_x_coarse_varid,'units','K kg/m^2/s');
            netcdf.putAtt(ncid,t_flux_y_coarse_varid,'units','K kg/m^2/s');
            netcdf.putAtt(ncid,t_flux_z_coarse_varid,'units','K kg/m^2/s');
            
            netcdf.putAtt(ncid,qt_flux_x_coarse_varid,'units','g/m^2/s');
            netcdf.putAtt(ncid,qt_flux_y_coarse_varid,'units','g/m^2/s');
            netcdf.putAtt(ncid,qt_flux_z_coarse_varid,'units','g/m^2/s');
            
            netcdf.putAtt(ncid,qp_flux_x_coarse_varid,'units','g/m^2/s');
            netcdf.putAtt(ncid,qp_flux_y_coarse_varid,'units','g/m^2/s');
            netcdf.putAtt(ncid,qp_flux_z_coarse_varid,'units','g/m^2/s');
            
            %advection tendency\
            if calc_advection_tend
                netcdf.putAtt(ncid,tfull_flux_x_tend_resolved_varid,'units','K /s');
                netcdf.putAtt(ncid,tfull_flux_y_tend_resolved_varid,'units','K /s');
                netcdf.putAtt(ncid,tfull_flux_z_tend_resolved_varid,'units','K /s');
                
                netcdf.putAtt(ncid,t_flux_x_tend_resolved_varid,'units','K /s');
                netcdf.putAtt(ncid,t_flux_y_tend_resolved_varid,'units','K /s');
                netcdf.putAtt(ncid,t_flux_z_tend_resolved_varid,'units','K /s');
            end
            %diffusion
            netcdf.putAtt(ncid,tfull_diff_flx_x_resolved_varid,'units','K*m /s');
            netcdf.putAtt(ncid,tfull_diff_flx_y_resolved_varid,'units','K*m /s');
            netcdf.putAtt(ncid,tfull_diff_flx_z_resolved_varid,'units','K *kg /s/m^2');
            
            netcdf.putAtt(ncid,t_diff_flx_x_resolved_varid,'units','K*m /s');
            netcdf.putAtt(ncid,t_diff_flx_y_resolved_varid,'units','K*m /s');
            netcdf.putAtt(ncid,t_diff_flx_z_resolved_varid,'units','K *kg /s/m^2');
            
            netcdf.putAtt(ncid,qt_diff_flx_x_resolved_varid,'units','gr*m /s/kg');
            netcdf.putAtt(ncid,qt_diff_flx_y_resolved_varid,'units','gr*m /s/kg');
            netcdf.putAtt(ncid,qt_diff_flx_z_resolved_varid,'units','gr/s/m^2');
            
            netcdf.putAtt(ncid,qp_diff_flx_x_resolved_varid,'units','gr*m /s/kg');
            netcdf.putAtt(ncid,qp_diff_flx_y_resolved_varid,'units','gr*m /s/kg');
            netcdf.putAtt(ncid,qp_diff_flx_z_resolved_varid,'units','gr/s/m^2');
            
            %The coarse diffusive fluxes (in case we want to predict them
            %instead)
            
            %Also using the coarse diffusion - (Undecided if it makes more sense than the coarse-resolved)
            
            netcdf.putAtt(ncid,tfull_diff_flx_x_coarse_varid,'units','K*m /s');
            netcdf.putAtt(ncid,tfull_diff_flx_y_coarse_varid,'units','K*m /s');
            netcdf.putAtt(ncid,tfull_diff_flx_z_coarse_varid,'units','K *kg /s/m^2');
            
            netcdf.putAtt(ncid,t_diff_flx_x_coarse_varid,'units','K*m /s');
            netcdf.putAtt(ncid,t_diff_flx_y_coarse_varid,'units','K*m /s');
            netcdf.putAtt(ncid,t_diff_flx_z_coarse_varid,'units','K *kg /s/m^2');
            
            netcdf.putAtt(ncid,qt_diff_flx_x_coarse_varid,'units','gr*m /s/kg');
            netcdf.putAtt(ncid,qt_diff_flx_y_coarse_varid,'units','gr*m /s/kg');
            netcdf.putAtt(ncid,qt_diff_flx_z_coarse_varid,'units','gr/s/m^2');
            
            netcdf.putAtt(ncid,qp_diff_flx_x_coarse_varid,'units','gr*m /s/kg');
            netcdf.putAtt(ncid,qp_diff_flx_y_coarse_varid,'units','gr*m /s/kg');
            netcdf.putAtt(ncid,qp_diff_flx_z_coarse_varid,'units','gr/s/m^2');
            
            
            
            
            %diffusion tendency
            if calc_diffusive_tend
                netcdf.putAtt(ncid,tfull_diff_x_tend_resolved_varid,'units','K /s');
                netcdf.putAtt(ncid,tfull_diff_y_tend_resolved_varid,'units','K /s');
                netcdf.putAtt(ncid,tfull_diff_z_tend_resolved_varid,'units','K /s');
                
                netcdf.putAtt(ncid,t_diff_x_tend_resolved_varid,'units','K /s');
                netcdf.putAtt(ncid,t_diff_y_tend_resolved_varid,'units','K /s');
                netcdf.putAtt(ncid,t_diff_z_tend_resolved_varid,'units','K /s');
                
                netcdf.putAtt(ncid,qt_diff_x_tend_resolved_varid,'units','gr /s/kg');
                netcdf.putAtt(ncid,qt_diff_y_tend_resolved_varid,'units','gr /s/kg');
                netcdf.putAtt(ncid,qt_diff_z_tend_resolved_varid,'units','gr /s/kg');
                
                netcdf.putAtt(ncid,qp_diff_x_tend_resolved_varid,'units','gr /s/kg');
                netcdf.putAtt(ncid,qp_diff_y_tend_resolved_varid,'units','gr /s/kg');
                netcdf.putAtt(ncid,qp_diff_z_tend_resolved_varid,'units','gr /s/kg');
            end
            
            %precip fall
            netcdf.putAtt(ncid,dqp_fall_resolved_varid,'units','gr /s/kg');
            netcdf.putAtt(ncid,dqp_fall_coarse_varid,'units','gr /s/kg');
            
            netcdf.putAtt(ncid,t_fall_tend_resolved_varid,'units','K /s');
            netcdf.putAtt(ncid,t_fall_tend_coarse_varid,'units','K /s');
            
            %cloud
            netcdf.putAtt(ncid,cloud_qt_tend_resolved_varid,'units','gr /s/kg');
            netcdf.putAtt(ncid,cloud_qt_tend_coarse_varid,'units','gr /s/kg');
            
            netcdf.putAtt(ncid,cloud_lat_heat_resolved_varid,'units','K /s');
            netcdf.putAtt(ncid,cloud_lat_heat_coarse_varid,'units','K /s');
            
            %cloud fluxes
            netcdf.putAtt(ncid,fzt_resolved_varid,'units','K m^2 /s kg');
            netcdf.putAtt(ncid,fz_resolved_varid,'units','gr m^2/s/kg^2');
            
            netcdf.putAtt(ncid,fzt_coarse_varid,'units','K m^2 /s kg');
            netcdf.putAtt(ncid,fz_coarse_varid,'units','gr m^2/s/kg^2');
            %diffusivity
            netcdf.putAtt(ncid,tkz_out_resolved_varid,'units','m^2/s');
            netcdf.putAtt(ncid,tkz_out_coarse_varid,'units','m^2/s');
            
            netcdf.putAtt(ncid,Pr1_resolved_varid,'units','Number -none');
            netcdf.putAtt(ncid,Pr1_coarse_varid,'units','Number -none');
            
            
            netcdf.putAtt(ncid,dqp_varid,'units','g/kg/s');
            netcdf.putAtt(ncid,dqp_resolved_varid,'units','g/kg/s');
            netcdf.putAtt(ncid,precip_varid,'units','kg/m^2/s');
            netcdf.putAtt(ncid,precip_resolved_varid,'units','kg/m^2/s');
            netcdf.putAtt(ncid,precip_energy_varid,'units','K kg/m^2/s');
            netcdf.putAtt(ncid,precip_energy_resolved_varid,'units','K kg/m^2/s');
            
            %momentum vertical advection
            netcdf.putAtt(ncid,dudt_advect_coarse_varid,'units','m/s^2');
            netcdf.putAtt(ncid,dvdt_advect_coarse_varid,'units','m/s^2');
            netcdf.putAtt(ncid,dwdt_advect_coarse_varid,'units','m/s^2');
            
            netcdf.putAtt(ncid,dudt_advect_resolved_varid,'units','m/s^2');
            netcdf.putAtt(ncid,dvdt_advect_resolved_varid,'units','m/s^2');
            netcdf.putAtt(ncid,dwdt_advect_resolved_varid,'units','m/s^2');
            
              %Momentum Surface wind fluxes
              netcdf.putAtt(ncid,fluxbv_coarse_varid,'units','m/s^2 * kg/m^3?');  
              netcdf.putAtt(ncid,fluxbu_coarse_varid,'units','m/s^2 * kg/m^3?');  
              netcdf.putAtt(ncid,fluxbv_resolved_varid,'units','m/s^2 * kg/m^3?');  
              netcdf.putAtt(ncid,fluxbu_resolved_varid,'units','m/s^2 * kg/m^3?');  
              
              netcdf.putAtt(ncid,fluxbv_samson_coarse_varid,'units','m/s^2 * kg/m^3?');  
              netcdf.putAtt(ncid,fluxbu_samson_coarse_varid,'units','m/s^2 * kg/m^3?');  
              netcdf.putAtt(ncid,fluxbv_samson_resolved_varid,'units','m/s^2 * kg/m^3?');  
              netcdf.putAtt(ncid,fluxbu_samson_resolved_varid,'units','m/s^2 * kg/m^3?');  
                            
%             % change place? 
              % The direct calculation of DQP on T (HL) 
              netcdf.putAtt(ncid,dqp_t_tend_resolved_varid,'units','K/s');  
              netcdf.putAtt(ncid,dqp_t_tend_coarse_varid,'units','K/s');  
              %Generalization
              netcdf.putAtt(ncid,rh_init_varid,'units','percent'); 
              netcdf.putAtt(ncid,rh_end_resolved_varid,'units','percent'); 
              netcdf.putAtt(ncid,rh_end_coarse_varid,'units','percent'); 
              netcdf.putAtt(ncid,tabs_resolved_end_resolved_varid,'units','K'); 
              netcdf.putAtt(ncid,tabs_resolved_end_coarse_varid,'units','K'); 
             
            netcdf.endDef(ncid);
            
            % write data
            netcdf.putVar(ncid,x_varid,x_coarse);
            netcdf.putVar(ncid,y_varid,y_coarse);
            netcdf.putVar(ncid,z_varid,z);
            netcdf.putVar(ncid,p_varid,pres);
            netcdf.putVar(ncid,rho_varid,rho);
            
            netcdf.putVar(ncid,u_varid,u_coarse_flip);
            netcdf.putVar(ncid,v_varid,v_coarse_flip);
            netcdf.putVar(ncid,w_varid,w_coarse_flip);
            netcdf.putVar(ncid,tabs_varid,tabs_coarse_flip);
            netcdf.putVar(ncid,t_varid,t_coarse_flip);
            netcdf.putVar(ncid,tfull_varid,tfull_coarse_flip);
            netcdf.putVar(ncid,tfull_init_varid,tfull_coarse_init_flip);
            netcdf.putVar(ncid,qv_varid,qv_coarse_flip*1000.0); % g/kg
            netcdf.putVar(ncid,qn_varid,qn_coarse_flip*1000.0); % g/kg
            netcdf.putVar(ncid,qp_varid,qp_coarse_flip*1000.0); % g/kg
            netcdf.putVar(ncid,Qrad_varid,Qrad_coarse_flip);
            
            
            netcdf.putVar(ncid,tabs_resolved_init_varid,tabs_resolved_init_flip);
            netcdf.putVar(ncid,qv_coarse_init_varid,qv_coarse_init_flip*1000.0); % g/kg
            netcdf.putVar(ncid,qn_coarse_init_varid,qn_coarse_init_flip*1000.0); % g/kg
            netcdf.putVar(ncid,qp_coarse_init_varid,qp_coarse_init_flip*1000.0); % g/kg

            netcdf.putVar(ncid,windspeed_coarse_varid,windspeed_coarse_flip);
            netcdf.putVar(ncid,windspeed_resolved_varid,windspeed_resolved_flip);
            
            %advection flux
            
            netcdf.putVar(ncid,tfull_flux_x_resolved_varid,tfull_flux_x_resolved_flip);
            netcdf.putVar(ncid,tfull_flux_y_resolved_varid,tfull_flux_y_resolved_flip);
            netcdf.putVar(ncid,tfull_flux_z_resolved_varid,tfull_flux_z_resolved_flip);
            
            netcdf.putVar(ncid,t_flux_x_resolved_varid,t_flux_x_resolved_flip);
            netcdf.putVar(ncid,t_flux_y_resolved_varid,t_flux_y_resolved_flip);
            netcdf.putVar(ncid,t_flux_z_resolved_varid,t_flux_z_resolved_flip);
            
            netcdf.putVar(ncid,qt_flux_x_resolved_varid,qt_flux_x_resolved_flip*1000.0);
            netcdf.putVar(ncid,qt_flux_y_resolved_varid,qt_flux_y_resolved_flip*1000.0);
            netcdf.putVar(ncid,qt_flux_z_resolved_varid,qt_flux_z_resolved_flip*1000.0);
            
            netcdf.putVar(ncid,qp_flux_x_resolved_varid,qp_flux_x_resolved_flip*1000.0);
            netcdf.putVar(ncid,qp_flux_y_resolved_varid,qp_flux_y_resolved_flip*1000.0);
            netcdf.putVar(ncid,qp_flux_z_resolved_varid,qp_flux_z_resolved_flip*1000.0);
            
            %coarse advection flux
            
            netcdf.putVar(ncid,tfull_flux_x_coarse_varid,tfull_flux_x_coarse_flip);
            netcdf.putVar(ncid,tfull_flux_y_coarse_varid,tfull_flux_y_coarse_flip);
            netcdf.putVar(ncid,tfull_flux_z_coarse_varid,tfull_flux_z_coarse_flip);
            
            netcdf.putVar(ncid,t_flux_x_coarse_varid,t_flux_x_coarse_flip);
            netcdf.putVar(ncid,t_flux_y_coarse_varid,t_flux_y_coarse_flip);
            netcdf.putVar(ncid,t_flux_z_coarse_varid,t_flux_z_coarse_flip);
            
            netcdf.putVar(ncid,qt_flux_x_coarse_varid,qt_flux_x_coarse_flip*1000.0);
            netcdf.putVar(ncid,qt_flux_y_coarse_varid,qt_flux_y_coarse_flip*1000.0);
            netcdf.putVar(ncid,qt_flux_z_coarse_varid,qt_flux_z_coarse_flip*1000.0);
            
            netcdf.putVar(ncid,qp_flux_x_coarse_varid,qp_flux_x_coarse_flip*1000.0);
            netcdf.putVar(ncid,qp_flux_y_coarse_varid,qp_flux_y_coarse_flip*1000.0);
            netcdf.putVar(ncid,qp_flux_z_coarse_varid,qp_flux_z_coarse_flip*1000.0);
            
            
            %advection
            
            %advection tendency
            if calc_advection_tend
                netcdf.putVar(ncid,tfull_flux_x_tend_resolved_varid,tfull_flux_x_tend_resolved_flip);
                netcdf.putVar(ncid,tfull_flux_y_tend_resolved_varid,tfull_flux_y_tend_resolved_flip);
                netcdf.putVar(ncid,tfull_flux_z_tend_resolved_varid,tfull_flux_z_tend_resolved_flip);
                
                netcdf.putVar(ncid,t_flux_x_tend_resolved_varid,t_flux_x_tend_resolved_flip);
                netcdf.putVar(ncid,t_flux_y_tend_resolved_varid,t_flux_y_tend_resolved_flip);
                netcdf.putVar(ncid,t_flux_z_tend_resolved_varid,t_flux_z_tend_resolved_flip);
                
            end
            %diffusion fluxes and tendencies
            netcdf.putVar(ncid,tfull_diff_flx_x_resolved_varid,tfull_diff_flx_x_resolved_flip);
            netcdf.putVar(ncid,tfull_diff_flx_y_resolved_varid,tfull_diff_flx_y_resolved_flip);
            netcdf.putVar(ncid,tfull_diff_flx_z_resolved_varid,tfull_diff_flx_z_resolved_flip);
            
            netcdf.putVar(ncid,t_diff_flx_x_resolved_varid,t_diff_flx_x_resolved_flip);
            netcdf.putVar(ncid,t_diff_flx_y_resolved_varid,t_diff_flx_y_resolved_flip);
            netcdf.putVar(ncid,t_diff_flx_z_resolved_varid,t_diff_flx_z_resolved_flip);
            
            netcdf.putVar(ncid,qt_diff_flx_x_resolved_varid,qt_diff_flx_x_resolved_flip*1000.0);
            netcdf.putVar(ncid,qt_diff_flx_y_resolved_varid,qt_diff_flx_y_resolved_flip*1000.0);
            netcdf.putVar(ncid,qt_diff_flx_z_resolved_varid,qt_diff_flx_z_resolved_flip*1000.0);
            
            netcdf.putVar(ncid,qp_diff_flx_x_resolved_varid,qp_diff_flx_x_resolved_flip*1000.0);
            netcdf.putVar(ncid,qp_diff_flx_y_resolved_varid,qp_diff_flx_y_resolved_flip*1000.0);
            netcdf.putVar(ncid,qp_diff_flx_z_resolved_varid,qp_diff_flx_z_resolved_flip*1000.0);
            %coarse diffusion fluxes
            netcdf.putVar(ncid,tfull_diff_flx_x_coarse_varid,tfull_diff_flx_x_coarse_flip);
            netcdf.putVar(ncid,tfull_diff_flx_y_coarse_varid,tfull_diff_flx_y_coarse_flip);
            netcdf.putVar(ncid,tfull_diff_flx_z_coarse_varid,tfull_diff_flx_z_coarse_flip);
            
            netcdf.putVar(ncid,t_diff_flx_x_coarse_varid,t_diff_flx_x_coarse_flip);
            netcdf.putVar(ncid,t_diff_flx_y_coarse_varid,t_diff_flx_y_coarse_flip);
            netcdf.putVar(ncid,t_diff_flx_z_coarse_varid,t_diff_flx_z_coarse_flip);
            
            netcdf.putVar(ncid,qt_diff_flx_x_coarse_varid,qt_diff_flx_x_coarse_flip*1000.0);
            netcdf.putVar(ncid,qt_diff_flx_y_coarse_varid,qt_diff_flx_y_coarse_flip*1000.0);
            netcdf.putVar(ncid,qt_diff_flx_z_coarse_varid,qt_diff_flx_z_coarse_flip*1000.0);
            
            netcdf.putVar(ncid,qp_diff_flx_x_coarse_varid,qp_diff_flx_x_coarse_flip*1000.0);
            netcdf.putVar(ncid,qp_diff_flx_y_coarse_varid,qp_diff_flx_y_coarse_flip*1000.0);
            netcdf.putVar(ncid,qp_diff_flx_z_coarse_varid,qp_diff_flx_z_coarse_flip*1000.0);
            
            
            %diffusion tendency
            if calc_diffusive_tend
                netcdf.putVar(ncid,tfull_diff_x_tend_resolved_varid,tfull_diff_x_tend_resolved_flip);
                netcdf.putVar(ncid,tfull_diff_y_tend_resolved_varid,tfull_diff_y_tend_resolved_flip);
                netcdf.putVar(ncid,tfull_diff_z_tend_resolved_varid,tfull_diff_z_tend_resolved_flip);
                
                netcdf.putVar(ncid,t_diff_x_tend_resolved_varid,t_diff_x_tend_resolved_flip);
                netcdf.putVar(ncid,t_diff_y_tend_resolved_varid,t_diff_y_tend_resolved_flip);
                netcdf.putVar(ncid,t_diff_z_tend_resolved_varid,t_diff_z_tend_resolved_flip);
                
                netcdf.putVar(ncid,qt_diff_x_tend_resolved_varid,qt_diff_x_tend_resolved_flip*1000.0);
                netcdf.putVar(ncid,qt_diff_y_tend_resolved_varid,qt_diff_y_tend_resolved_flip*1000.0);
                netcdf.putVar(ncid,qt_diff_z_tend_resolved_varid,qt_diff_z_tend_resolved_flip*1000.0);
                
                netcdf.putVar(ncid,qp_diff_x_tend_resolved_varid,qp_diff_x_tend_resolved_flip*1000.0);
                netcdf.putVar(ncid,qp_diff_y_tend_resolved_varid,qp_diff_y_tend_resolved_flip*1000.0);
                netcdf.putVar(ncid,qp_diff_z_tend_resolved_varid,qp_diff_z_tend_resolved_flip*1000.0);
            end
            % precip fall tedndencies
            netcdf.putVar(ncid,dqp_fall_resolved_varid,dqp_fall_resolved_flip*1000.0);
            netcdf.putVar(ncid,t_fall_tend_resolved_varid,t_fall_tend_resolved_flip);
            
            netcdf.putVar(ncid,dqp_fall_coarse_varid,dqp_fall_coarse_flip*1000.0);
            netcdf.putVar(ncid,t_fall_tend_coarse_varid,t_fall_tend_coarse_flip);
            
            
            %cloud
            netcdf.putVar(ncid,cloud_lat_heat_resolved_varid,cloud_lat_heat_resolved_flip);
            netcdf.putVar(ncid,cloud_qt_tend_resolved_varid,cloud_qt_tend_resolved_flip*1000.0);
            
            netcdf.putVar(ncid,cloud_lat_heat_coarse_varid,cloud_lat_heat_coarse_flip);
            netcdf.putVar(ncid,cloud_qt_tend_coarse_varid,cloud_qt_tend_coarse_flip*1000.0);
            
            %cloud fluxes
            netcdf.putVar(ncid,fzt_resolved_varid,fzt_resolved_flip);
            netcdf.putVar(ncid,fzt_coarse_varid,fzt_coarse_flip);
            
            netcdf.putVar(ncid,fz_resolved_varid,fz_resolved_flip*1000.0);
            netcdf.putVar(ncid,fz_coarse_varid,fz_coarse_flip*1000.0);
            
            % diffusivity
            netcdf.putVar(ncid,tkz_out_resolved_varid,tkz_out_resolved_flip);
            netcdf.putVar(ncid,tkz_out_coarse_varid,tkz_out_coarse_flip);
            
            netcdf.putVar(ncid,Pr1_resolved_varid,Pr1_resolved_flip);
            netcdf.putVar(ncid,Pr1_coarse_varid,Pr1_coarse_flip);
            
            
            netcdf.putVar(ncid,dqp_varid,dqp_coarse_flip*1000.0); % g/kg/s
            netcdf.putVar(ncid,dqp_resolved_varid,dqp_resolved_flip*1000.0); % g/kg/s
            netcdf.putVar(ncid,precip_varid,precip_coarse_flip);
            netcdf.putVar(ncid,precip_resolved_varid,precip_resolved_flip);
            netcdf.putVar(ncid,precip_energy_varid,precip_energy_coarse_flip);
            netcdf.putVar(ncid,precip_energy_resolved_varid,precip_energy_resolved_flip);
            
            %momentum vertical advection 
            netcdf.putVar(ncid,dudt_advect_coarse_varid,dudt_advect_coarse_flip); % g/kg/s
            netcdf.putVar(ncid,dvdt_advect_coarse_varid,dvdt_advect_coarse_flip); % g/kg/s
            netcdf.putVar(ncid,dwdt_advect_coarse_varid,dwdt_advect_coarse_flip); % g/kg/s

            netcdf.putVar(ncid,dudt_advect_resolved_varid,dudt_advect_resolved_flip); % g/kg/s            
            netcdf.putVar(ncid,dvdt_advect_resolved_varid,dvdt_advect_resolved_flip); % g/kg/s            
            netcdf.putVar(ncid,dwdt_advect_resolved_varid,dwdt_advect_resolved_flip); % g/kg/s            
            
            
              %Momentum Surface wind fluxes
              netcdf.putVar(ncid,fluxbu_coarse_varid,fluxbu_coarse_flip);
              netcdf.putVar(ncid,fluxbv_coarse_varid,fluxbv_coarse_flip);      
              netcdf.putVar(ncid,fluxbu_resolved_varid,fluxbu_resolved_flip);
              netcdf.putVar(ncid,fluxbv_resolved_varid,fluxbv_resolved_flip);      
              
              
              netcdf.putVar(ncid,fluxbu_samson_coarse_varid,fluxbu_samson_coarse_flip);
              netcdf.putVar(ncid,fluxbv_samson_coarse_varid,fluxbv_samson_coarse_flip);      
              netcdf.putVar(ncid,fluxbu_samson_resolved_varid,fluxbu_samson_resolved_flip);
              netcdf.putVar(ncid,fluxbv_samson_resolved_varid,fluxbv_samson_resolved_flip); 
%             % change place? 
              % The direct calculation of DQP on T (HL) 
              netcdf.putVar(ncid,dqp_t_tend_resolved_varid,dqp_t_tend_resolved_flip);      
              netcdf.putVar(ncid,dqp_t_tend_coarse_varid,dqp_t_tend_coarse_flip);      
              %Generalization
              netcdf.putVar(ncid,rh_init_varid,rh_init_flip);      
              netcdf.putVar(ncid,rh_end_resolved_varid,rh_end_resolved_flip);      
              netcdf.putVar(ncid,rh_end_coarse_varid,rh_end_coarse_flip);      
              netcdf.putVar(ncid,tabs_resolved_end_resolved_varid,tabs_resolved_end_resolved_flip);      
              netcdf.putVar(ncid,tabs_resolved_end_coarse_varid,tabs_resolved_end_coarse_flip);      
            % close output netcdf
            netcdf.close(ncid);
            
            
            
            if do_show_times, sprintf('whole step 04 (time it took to write vars to netcdf): %f', toc), end
            tic
            
        end
        sprintf('1 iteration took time: %f', toc(tStart))

    end
end


            
    

function [var1,var2] = switch_vars(var1,var2)
tmp = var1;
var1 = var2;
var2 = tmp;
end

