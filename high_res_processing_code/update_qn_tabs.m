function [qn_resolved, tabs_resolved, qt_coarse, qp_coarse,kmin,kmax] = update_qn_tabs(tbgmax, tbgmin, tprmax, tprmin, fac_cond, fac_fus, fac_sub, tgrmax, tgrmin, qp_threshold, num_z, ...
    num_blocks_x, num_blocks_y, qn_coarse, qt_coarse, qp_coarse, tfull_coarse, gamaz, pres)
%A function to calculate qn and tabs (which are not prognostic variables
%   Detailed explanation goes here


tabs_resolved = zeros(size(tfull_coarse));
qn_resolved = zeros(size(tfull_coarse));

an = 1./(tbgmax-tbgmin) ;
bn = tbgmin * an;
ap = 1./(tprmax-tprmin) ;
bp = tprmin * ap;
fac1 = fac_cond+(1+bp)*fac_fus;
fac2 = fac_fus*ap;
%             ag = 1./(tgrmax-tgrmin);

kmax=0;
kmin=num_z+1;

for i=1:num_blocks_x
    for j=1:num_blocks_y
        for k=1:num_z
            qn_resolved(k,j,i) = qn_coarse(k,j,i);
            qt_coarse(k,j,i)=max(0.,qt_coarse(k,j,i));
            %             ! Initail guess for temperature assuming no cloud water/ice:
            tabs_resolved(k,j,i) = tfull_coarse(k,j,i)-gamaz(k); % Yani - modified to tfull.
            %In a newer version of Sam the initial guss makes more sense
            %(It makes a difference which one we use ?!nnn?!
            tabs1=(tabs_resolved(k,j,i)+fac1*qp_coarse(k,j,i))/(1.+fac2*qp_coarse(k,j,i));
%             tabs1=(tabs_resolved(k,j,i)+fac_cond*qp_coarse(k,j,i));
            %             ! Warm cloud:
            if(tabs1 >= tbgmax)
                tabs1=tabs_resolved(k,j,i)+fac_cond*qp_coarse(k,j,i);
                qsat = qsatw(tabs1,pres(k));
                %                 ! Ice cloud:
            elseif(tabs1 <= tbgmin)
                
                tabs1=tabs_resolved(k,j,i)+fac_sub*qp_coarse(k,j,i);
                qsat = qsati(tabs1,pres(k));
                %                 ! Mixed-phase cloud:
            else
                
                %Yani plugged an update to the temperature:
                omp  = max(0.,min(1.,(tabs1-tprmin)*ap));
                tabs1=tabs_resolved(k,j,i)+(fac_cond+(1.-omp).*fac_fus).*qp_coarse(k,j,i);
                %end of update
                om = an*tabs1-bn;
                qsat = om*qsatw(tabs1,pres(k))+(1.-om)*qsati(tabs1,pres(k));
            end
            
            if(qt_coarse(k,j,i) > qsat)
                
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
                        dlstarn=-an*fac_fus;  %corrected by JY yani
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
                        dlstarp=-ap*fac_fus; %corrected by JY yani
                    end
                    fff = tabs_resolved(k,j,i)-tabs1+lstarn*(qt_coarse(k,j,i)-qsat)+lstarp*qp_coarse(k,j,i);
                    dfff=dlstarn*(qt_coarse(k,j,i)-qsat)+dlstarp*qp_coarse(k,j,i)-lstarn*dqsat-1.;
                    dtabs=-fff/dfff;
                    niter=niter+1;
                    tabs1=tabs1+dtabs;
                end
                qsat = qsat + dqsat * dtabs;
                qn_resolved(k,j,i) = max(0.,qt_coarse(k,j,i)-qsat);
            else
                qn_resolved(k,j,i) = 0.;
            end
            tabs_resolved(k,j,i) = tabs1;
            qp_coarse(k,j,i) = max(0.,qp_coarse(k,j,i)); %! just in case
            if(qn_coarse(k,j,i)>qp_threshold)
                kmin = min(kmin,k);
                kmax = max(kmax,k);
            end
        end
    end
end







end

