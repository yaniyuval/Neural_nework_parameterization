function [dudt,dvdt,dwdt] = advect2_mom_z(nx,ny,nz,dz,rho,rhow,adz,adzw,u,v,w)
%Calculate vertical advection of momentum. 

% !       momentum tendency due to the 2nd-order-central vertical advection

dudt= zeros(nz,ny,nx);  
dvdt= zeros(nz,ny,nx);  
dwdt= zeros(nz,ny,nx);  

fuz = zeros(nz+1,ny,nx); %these are one level less compared to the Fortran code (though the top level in the fortran is zero). 
fvz = zeros(nz+1,ny,nx);
fwz = zeros(nz,ny,nx);  %I need to check which levels from w are actually output The level 1 of W is zero in the output. We don't have the top level
% momentum transport then.... 

%The output of W is 1 to nzm (doesn't output the top level)

dz25=1./(4.*dz);
dz2=dz25*2.;

for k=2:nz %Note that I do not have the upper W level- so cannot compute it like in Fortran!
 kb = k-1;
 rhoi = dz25 * rhow(k);
 for j=1:ny
  for i=1:nx
      
      if i == 1
          wxbot = w(k,j,end);
      else
          wxbot = w(k,j,i-1);
      end
          
      if j==1
          wybot = w(k,j,i);
      else
          wybot = w(k,j-1,i);
      end          
   fuz(k,j,i) = rhoi*(w(k,j,i)+wxbot)*(u(k,j,i)+u(kb,j,i));
   fvz(k,j,i) = rhoi*(w(k,j,i)+wybot)*(v(k,j,i)+v(kb,j,i));
  end 
 end 
end 

for k=1:nz - 1  %Note that I do not have the upper W level- so cannot compute it like in Fortran! I inserted nz-1 instead of nz
 kc = k+1;
 for j=1:ny
  for i=1:nx
   fwz(k,j,i)=dz25*(w(kc,j,i)*rhow(kc)+w(k,j,i)*rhow(k))*(w(kc,j,i)+w(k,j,i));   
   end 
 end 
end 

for k=nz  %Take care separately the top index  that I don't have. 
 for j=1:ny
  for i=1:nx
   fwz(k,j,i)=dz25*(2*w(k,j,i)*rhow(k))*(2*w(k,j,i));   
   end 
 end 
end

for k=1:nz  %Note that I do not have the upper W level- so cannot compute it like in Fortran! I inserted nz-1 instead of nz
 kc = k+1;
 rhoi = 1./(rho(k)*adz(k));
 for j=1:ny
  for i=1:nx
   dudt(k,j,i)=dudt(k,j,i)-(fuz(kc,j,i)-fuz(k,j,i))*rhoi;
   dvdt(k,j,i)=dvdt(k,j,i)-(fvz(kc,j,i)-fvz(k,j,i))*rhoi; 
   end 
 end 
end 

for k=2:nz
 kb=k-1;
 rhoi = 1./(rhow(k)*adzw(k));
 for j=1:ny
  for i=1:nx
   dwdt(k,j,i)=dwdt(k,j,i)-(fwz(k,j,i)-fwz(kb,j,i))*rhoi;
  end
 end 
end 


end %matlab function

                                                                  