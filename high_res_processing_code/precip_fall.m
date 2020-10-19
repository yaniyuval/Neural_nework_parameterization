function   [dqp_fall,t_fall_tend,precip]= ...
    precip_fall(qp,tabs,rho,rhow,rhofac,num_x,num_y,num_z,dz,adz,dtn,tprmin,a_pr,fac_fus,...
    fac_cond,crain,vrain,tgrmin,a_gr,qp_threshold,vgrau,cgrau,vsnow,csnow)%,do_fall_tend_qp,do_fall_tend_tfull)
                          

%Step 02:07b precip_fall calculation
        
        % calculate precipitation flux
        %%precip_fall
        % from precip.f90
        
        % only kept parts needed to calculate precipitation flux
        % and neglected non-oscillatory option for speed
        t_fall_tend= zeros(size(tabs));
        dqp_fall= zeros(size(tabs));
        precip = zeros(size(tabs));
        
        % qp = dummy111;%NOTE!!!!yANi
        nzm = num_z;
        nz = nzm+1;
        fz = zeros(nz,1);
        tmp_qp = zeros(nz,1);
        mx = zeros(nz,1);
        mn = zeros(nz,1);
        www = zeros(nz,1);
        lfac = zeros(nz,1);
        irhoadz= zeros(nzm,1);
        fz(nz)=0.; %Need to initialize size.
        www(nz)=0.;
        lfac(nz)=0.;
        eps = 1.e-10;
        wp= zeros(nzm,1);
        iwmax = zeros(nzm,1);
        
        for k = 1:num_z
            kb = max(1,k-1);
            wmax       = dz*adz(kb)/dtn; %  ! Velocity equivalent to a cfl of 1.0.
            iwmax(k)   = 1./wmax;
        end
        
        % Compute precipitation velocity and flux column-by-column
        for i=1:num_x
            for j=1:num_y
                prec_cfl = 0.0;
                for k=1:num_z
                    wp(k) = 0.0;
                    omp = max(0.,min(1.,(tabs(k,j,i)-tprmin)*a_pr));
                    lfac(k) = fac_cond+(1.-omp)*fac_fus;
                    
                    if(qp(k,j,i)>qp_threshold)
                        if(omp==1.)
                            wp(k)= rhofac(k)*vrain*(rho(k)*qp(k,j,i))^crain;
                        elseif(omp==0.)
                            omg = max(0.,min(1.,(tabs(k,j,i)-tgrmin)*a_gr));
                            qgg=omg*qp(k,j,i);
                            qss=qp(k,j,i)-qgg;
                            wp(k)= rhofac(k)*(omg*vgrau*(rho(k)*qgg)^cgrau ...
                                +(1.-omg)*vsnow*(rho(k)*qss)^csnow);
                        else
                            omg = max(0.,min(1.,(tabs(k,j,i)-tgrmin)*a_gr));
                            qrr=omp*qp(k,j,i);
                            qss=qp(k,j,i)-qrr;
                            qgg=omg*qss;
                            qss=qss-qgg;
                            wp(k)=rhofac(k)*(omp*vrain*(rho(k)*qrr)^crain ...
                                +(1.-omp)*(omg*vgrau*(rho(k)*qgg)^cgrau ...
                                +(1.-omg)*vsnow*(rho(k)*qss)^csnow));
                        end
                        % note leave out the dtn/dz factor which is removed in write_fields2D.f90
                        % Define upwind precipitation flux
                        prec_cfl = max(prec_cfl,wp(k)*iwmax(k));
                        precip(k,j,i)=qp(k,j,i)*wp(k)*rhow(k);
                        %                    wp_tests(k,j,i) = wp(k); %This was similar (3OOM) to the wp I calculated in Fortran
                        %                 wp_test0(k,j,i) = wp(k);
                        %                 wp_test1(k,j,i) = -wp(k)*rho(k)*dtn/dz;
                        %                 wp_test2(k,j,i) = -wp(k)*rhow(k)*dtn/dz;
                        wp(k) = -wp(k)*rhow(k)*dtn/dz; %more accurate with rhow
                        
                    end % if
                    
                    
                end
                
                
                if (prec_cfl > 0.3) %sub stepping scheme
                    nprec = max(1,ceil(prec_cfl/0.3));
                    for k = 1:nzm
                        wp(k) = wp(k)/nprec;
                    end
                else
                    nprec = 1;
                end
                for lll = 1:nprec
                    %% Added by Yani to take into account precip fall affect (in the calculation of dqp, and maybe later also)
                    for k = 1:nzm
                        tmp_qp(k) = qp(k,j,i); % Temporary array for qp in this column
                        irhoadz(k) = 1./(rho(k)*adz(k)); %! Useful factor - agrees better with the fortran irhoadz var than using rhow.
                    end
                    for k=1:nzm
                        kc=min(nzm,k+1);
                        kb=max(1,k-1);
                        mx(k)=max([tmp_qp(kb),tmp_qp(kc),tmp_qp(k)]);
                        mn(k)=min([tmp_qp(kb),tmp_qp(kc),tmp_qp(k)]);
                        fz(k)=tmp_qp(k)*wp(k);
                    end
                    
                    
                    for k=1:nzm
                        kc=k+1;
                        tmp_qp(k)=tmp_qp(k)-(fz(kc)-fz(k))*irhoadz(k);
                    end
                    
                    for k=1:nzm
                        %             ! Also, compute anti-diffusive correction to previous
                        %             ! (upwind) approximation to the flux
                        kb=max(1,k-1);
                        %             ! The precipitation velocity is a cell-centered quantity,
                        %             ! since it is computed from the cell-centered
                        %             ! precipitation mass fraction.  Therefore, a reformulated
                        %             ! anti-diffusive flux is used here which accounts for
                        %             ! this and results in reduced numerical diffusion.
                        www(k) = 0.5*(1.+wp(k)*irhoadz(k)) ...
                            *(tmp_qp(kb)*wp(kb) - tmp_qp(k)*wp(k)); %! works for wp(k)<0
                    end
                    
                    
                    for k=1:nzm
                        kc=min(nzm,k+1);
                        kb=max(1,k-1);
                        mx(k)=max([tmp_qp(kb),tmp_qp(kc),tmp_qp(k),mx(k)]);
                        mn(k)=min([tmp_qp(kb),tmp_qp(kc),tmp_qp(k),mn(k)]);
                        %                 mn_test(k,j,i) =mn(k); Works well compare to fortran when
                        %                 inputing the exact qp!
                    end
                    
                    for k=1:nzm
                        kc=min(nzm,k+1);
                        mx(k)=rho(k)*adz(k)*(mx(k)-tmp_qp(k)) ...
                            /(pn(www(kc)) + pp(www(k))+eps);
                        mn(k)=rho(k)*adz(k)*(tmp_qp(k)-mn(k)) ...
                            /(pp(www(kc)) + pn(www(k))+eps);
                    end
                    
                    for k=1:nzm
                        kb=max(1,k-1);
                        %                ! Add limited flux correction to fz(k).
                        fz(k) = fz(k) ...                       % ! Upwind flux
                            + pp(www(k))*min([1.,mx(k), mn(kb)]) ...
                            - pn(www(k))*min([1.,mx(kb),mn(k)]); % ! Anti-diffusive flux
                    end
                    
                    for k=1:nzm
                        kc=k+1;
                        %! Update precipitation mass fraction.
                        %! Note that fz is the total flux, including both the
                        %! upwind flux and the anti-diffusive correction.
                        dqp_fall(k,j,i)=dqp_fall(k,j,i)-(fz(kc)-fz(k))*irhoadz(k);
%                         if do_fall_tend_qp
                        qp(k,j,i) = qp(k,j,i)  -(fz(kc)-fz(k))*irhoadz(k);
%                         end
                        %                 negative values?
                        lat_heat = -(lfac(kc)*fz(kc)-lfac(k)*fz(k))*irhoadz(k);
                        t_fall_tend(k,j,i)=t_fall_tend(k,j,i)-lat_heat;
                    end
                    
                    %%
                    
                end
            end
        end
        %Yani:Note this change (It is in order to get the correct tfull for precip_proc
        %calculation. I need to think if it is necessary (and whether I want to
        %model it or not (in the precip fall).
       
        % tfull- makes tfull less accurate for some reason - I need to think of it why and where I have an error.
        % calcu late energy flux associated with precipitation for use in tfull equation (SAM uses something a little different from equation A3 of SAM ref paper)
        %omp  = max(0.,min(1.,(tabs-tprmin)*a_pr)); % need to calculate again as used as scalar in precip_fall

        
        %         precip_energy = precip.*(fac_cond+(1.-omp_3d).*fac_fus); %I need to consider to recalc omp_3d due to changes in tabs
%         
%          if do_fall_tend_tfull
%             tfull = tfull + t_fall_tend; % I found that it is better not to change
%         end
%         dqp_fall = dqp_fall./dtn; %(all the tendencies are devided by dtn - and multiplied in the RF in SAMSON).
%         t_fall_tend = t_fall_tend./dtn;
end
        


