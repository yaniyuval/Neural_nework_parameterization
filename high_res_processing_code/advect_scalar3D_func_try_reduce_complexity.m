%taken from fortran advect_scalar3D.f90 - when approximation for advection
%didn't work for qp,qt.

function [f_out,flux_x,flux_y,flux_z,tend] = advect_scalar3D_func_try_reduce_complexity(f,u,v,w,rho,rhow,dx,dy,dz,dtn,num_x,num_y, num_z,adz)

eps = 0.0000000001;

% www(:,:,nz)=0.


v_bound = zeros([size(v,1),size(v,2)+2,size(v,3)+2]);
u_bound = zeros([size(v,1),size(v,2)+2,size(v,3)+2]);
w_bound = zeros([size(v,1)+1,size(v,2)+2,size(v,3)+2]);
f_bound = zeros([size(f,1)+1,size(f,2)+2,size(f,3)+2]);
v_bound(:,2:end-1,2:end-1) = v;
u_bound(:,2:end-1,2:end-1) = u;
w_bound(1:end-1,2:end-1,2:end-1) = w;
w_bound(end,:,:) = 0;
f_bound(1:end-1,2:end-1,2:end-1) = f;
f_bound(end,:,:) = 0;
%
% v_bound(:,1,:) = 0;
% v_bound(:,end,:) = 0;
v_bound(:,2:end-1,1) = v(:,:,end-1);
v_bound(:,2:end-1,end) = v(:,:,2);
%         v_bound(:,1,:) = 0;
%         v_bound(:,end,:) = 0;

u_bound(:,2:end-1,1) = u(:,:,end-1);
u_bound(:,2:end-1,end) = u(:,:,2);
%         u_bound(:,1,:) = 0;
%         u_bound(:,end,:) = 0;

w_bound(1:end-1,2:end-1,1) = w(:,:,end-1);
w_bound(1:end-1,2:end-1,end) = w(:,:,2);
%         w_bound(:,1,:) = 0;
%         w_bound(:,end,:) = 0;

f_bound(1:end-1,2:end-1,1) = f(:,:,end-1);
f_bound(1:end-1,2:end-1,end) = f(:,:,2);
%         f_bound(:,1,:) = 0;
%         f_bound(:,end,:) = 0;

% %
%         qp_bound(1:end-1,1,2:end-1) . % I think I want to have qp at the edges to be like to

irho = zeros([num_z,1]);
iadz= zeros([num_z,1]);
irhow= zeros([num_z,1]);


rho_3d_xy = repmat(rho,[1,size(v_bound,2),size(v_bound,3)]);
rhow_z = repmat(rhow(1:end-1),[1,size(v_bound,2),size(v_bound,3)]);
dtdx = rho_3d_xy.*dtn./dx;
dtdy = rho_3d_xy.*dtn./dy;
dtdz = rhow_z.*dtn./dz;

u_bound = u_bound.*dtdx;
v_bound = v_bound.*dtdy;
w_bound = w_bound(1:end-1,:,:).*dtdz;

f_bound = permute(f_bound,[3,2,1]);
v_bound = permute(v_bound,[3,2,1]);
u_bound = permute(u_bound,[3,2,1]);
w_bound = permute(w_bound,[3,2,1]);

uuu =  zeros([size(f_bound,1),size(f_bound,2),size(f_bound,3)]);
www =  zeros([size(f_bound,1),size(f_bound,2),size(f_bound,3)]);
vvv =  zeros([size(f_bound,1),size(f_bound,2),size(f_bound,3)]);
mx = zeros([size(f_bound,1),size(f_bound,2),size(f_bound,3)]);
mn = zeros([size(f_bound,1),size(f_bound,2),size(f_bound,3)]);
tend= zeros([size(f_bound,1),size(f_bound,2),size(f_bound,3)]);
tend_tmp= zeros([size(f_bound,1),size(f_bound,2),size(f_bound,3)]);
flux_x= zeros([size(f_bound,1),size(f_bound,2),size(f_bound,3)]);
flux_y= zeros([size(f_bound,1),size(f_bound,2),size(f_bound,3)]);
flux_z= zeros([size(f_bound,1),size(f_bound,2),size(f_bound,3)]);




for k=1:num_z
    kc=min(num_z,k+1);
    kb=max(1,k-1);
    for j=2:num_y+1
        jb=j-1;
        jc=j+1;
        for i=2:num_x+1
            ib=i-1;
            ic=i+1;
            mx(i,j,k)=max([f_bound(ib,j,k),f_bound(ic,j,k),f_bound(i,jb,k), ...
                f_bound(i,jc,k),f_bound(i,j,kb),f_bound(i,j,kc),f_bound(i,j,k)]);
            mn(i,j,k)=min([f_bound(ib,j,k),f_bound(ic,j,k),f_bound(i,jb,k), ...
                f_bound(i,jc,k),f_bound(i,j,kb),f_bound(i,j,kc),f_bound(i,j,k)]);
        end
    end
end


for k=1:num_z
    for j=2:num_y+1
        for i=2:num_x+1
            uuu(i,j,k)=max(0.,u_bound(i,j,k))*f_bound(i-1,j,k)+min(0.,u_bound(i,j,k))*f_bound(i,j,k);
            
        end
    end
end

for k=1:num_z
    for j=2:num_y+1
        for i=2:num_x+1
            vvv(i,j,k)=max(0.,v_bound(i,j,k))*f_bound(i,j-1,k)+min(0.,v_bound(i,j,k))*f_bound(i,j,k);
        end
    end
end

for k=1:num_z
    kb=max(1,k-1);
    for j=2:num_y+1
        for i=2:num_x+1
            www(i,j,k)=max(0.,w_bound(i,j,k))*f_bound(i,j,kb)+min(0.,w_bound(i,j,k))*f_bound(i,j,k);
        end
    end
end


for k=1:num_z
    irho(k) = 1./rho(k);
    iadz(k) = 1./adz(k);
    for j=2:num_y+1
        for i=2:num_x+1
            
            tend(i,j,k) =  -(uuu(i+1,j,k)-uuu(i,j,k)+vvv(i,j+1,k)-vvv(i,j,k) ...
                +(www(i,j,k+1)-www(i,j,k))*iadz(k))*irho(k);
            flux_x(i,j,k) = uuu(i,j,k);
            flux_y(i,j,k) = vvv(i,j,k);
            flux_z(i,j,k) = www(i,j,k);
            
        end
    end
end

f_bound = f_bound + tend;

for k=1:num_z
    kc=min(num_z,k+1);
    kb=max(1,k-1);
    dd=2./(kc-kb)/adz(k);
    for j=2:num_y+1
        jb=j-1;% YAni???
        jc=j+1;
        for i=2:num_x+1
            ib=i-1;
            uuu(i,j,k)=andiff(f_bound(ib,j,k),f_bound(i,j,k),u_bound(i,j,k),irho(k)) ...
                -(across(f_bound(ib,jc,k)+f_bound(i,jc,k)-f_bound(ib,jb,k)-f_bound(i,jb,k), ...
                u_bound(i,j,k), v_bound(ib,j,k)+v_bound(ib,jc,k)+v_bound(i,jc,k)+v_bound(i,j,k)) ...
                +across(dd*(f_bound(ib,j,kc)+f_bound(i,j,kc)-f_bound(ib,j,kb)-f_bound(i,j,kb)), ...
                u_bound(i,j,k), w_bound(ib,j,k)+w_bound(ib,j,kc)+w_bound(i,j,k)+w_bound(i,j,kc))) *irho(k);
        end
    end
end

for k=1:num_z
    kc=min(num_z,k+1);
    kb=max(1,k-1);
    dd=2./(kc-kb)/adz(k);
    for j=2:num_y+1
        jb=j-1;
        for i=2:num_x+1
            ib=i-1;
            ic=i+1;
            vvv(i,j,k)=andiff(f_bound(i,jb,k),f_bound(i,j,k),v_bound(i,j,k),irho(k))...
                -(across(f_bound(ic,jb,k)+f_bound(ic,j,k)-f_bound(ib,jb,k)-f_bound(ib,j,k), ...
                v_bound(i,j,k), u_bound(i,jb,k)+u_bound(i,j,k)+u_bound(ic,j,k)+u_bound(ic,jb,k)) ...
                +across(dd*(f_bound(i,jb,kc)+f_bound(i,j,kc)-f_bound(i,jb,kb)-f_bound(i,j,kb)), ...
                v_bound(i,j,k), w_bound(i,jb,k)+w_bound(i,j,k)+w_bound(i,j,kc)+w_bound(i,jb,kc))) *irho(k);
        end
    end
end

for k=1:num_z
    kb=max(1,k-1);
    irhow(k)=1./(rhow(k)*adz(k));
    for j=2:num_y+1
        jb=j-1;
        jc=j+1;
        for i=2:num_x+1
            ib=i-1;
            ic=i+1;
            www(i,j,k)=andiff(f_bound(i,j,kb),f_bound(i,j,k),w_bound(i,j,k),irhow(k)) ...
                -(across(f_bound(ic,j,kb)+f_bound(ic,j,k)-f_bound(ib,j,kb)-f_bound(ib,j,k), ...
                w_bound(i,j,k), u_bound(i,j,kb)+u_bound(i,j,k)+u_bound(ic,j,k)+u_bound(ic,j,kb)) ...
                +across(f_bound(i,jc,k)+f_bound(i,jc,kb)-f_bound(i,jb,k)-f_bound(i,jb,kb), ...
                w_bound(i,j,k), v_bound(i,j,kb)+v_bound(i,jc,kb)+v_bound(i,jc,k)+v_bound(i,j,k))) *irho(k);
        end
    end
end

www(:,:,1) = 0.;



for k=1:num_z
    kc=min(num_z,k+1);
    kb=max(1,k-1);
    for j=2:num_y+1
        jb=j-1;
        jc=j+1;
        for i=2:num_x+1
            ib=i-1;
            ic=i+1;
            mx(i,j,k)=max([f_bound(ib,j,k),f_bound(ic,j,k),f_bound(i,jb,k), ...
                f_bound(i,jc,k),f_bound(i,j,kb),f_bound(i,j,kc),f_bound(i,j,k),mx(i,j,k)]);
            mn(i,j,k)=min([f_bound(ib,j,k),f_bound(ic,j,k),f_bound(i,jb,k), ...
                f_bound(i,jc,k),f_bound(i,j,kb),f_bound(i,j,kc),f_bound(i,j,k),mn(i,j,k)]);
        end
    end
end

for k=1:num_z
    kc=min(num_z,k+1);
    for j=2:num_y+1
        jc=j+1;
        for i=2:num_x+1
            ic=i+1;
            mx(i,j,k)=rho(k)*(mx(i,j,k)-f_bound(i,j,k))/ ...
                (pn(uuu(ic,j,k)) + pp(uuu(i,j,k))+ ...
                pn(vvv(i,jc,k)) + pp(vvv(i,j,k))+ ...
                iadz(k)*(pn(www(i,j,kc)) + pp(www(i,j,k)))+eps);
            mn(i,j,k)=rho(k)*(f_bound(i,j,k)-mn(i,j,k))/ ...
                (pp(uuu(ic,j,k)) + pn(uuu(i,j,k))+ ...
                pp(vvv(i,jc,k)) + pn(vvv(i,j,k))+ ...
                iadz(k)*(pp(www(i,j,kc)) + pn(www(i,j,k)))+eps);
        end
    end
end

for k=1:num_z
    for j=2:num_y+1
        for i=2:num_x+1
            ib=i-1;
            uuu(i,j,k)=pp(uuu(i,j,k))*min([1.,mx(i,j,k), mn(ib,j,k)]) ...
                - pn(uuu(i,j,k))*min([1.,mx(ib,j,k),mn(i,j,k)]);
        end
    end
end

for k=1:num_z
    for j=2:num_y+1
        jb=j-1;
        for i=2:num_x+1
            vvv(i,j,k)=pp(vvv(i,j,k))*min([1.,mx(i,j,k), mn(i,jb,k)]) ...
                - pn(vvv(i,j,k))*min([1.,mx(i,jb,k),mn(i,j,k)]);
        end
    end
end

for k=1:num_z
    kb=max(1,k-1);
    for j=2:num_y+1
        for i=2:num_x+1
            www(i,j,k)=pp(www(i,j,k))*min([1.,mx(i,j,k), mn(i,j,kb)]) ...
                - pn(www(i,j,k))*min([1.,mx(i,j,kb),mn(i,j,k)]);
        end
    end
end




for k=1:num_z
    for j=2:num_y+1
        for i=2:num_x+1
            tend_tmp(i,j,k) = -(uuu(i+1,j,k)-uuu(i,j,k)+vvv(i,j+1,k)-vvv(i,j,k) ...
                +(www(i,j,k+1)-www(i,j,k))*iadz(k))*irho(k);
            f_bound(i,j,k) = f_bound(i,j,k)+ tend_tmp(i,j,k);
            flux_x(i,j,k) = flux_x(i,j,k) + uuu(i,j,k);
            flux_y(i,j,k) = flux_y(i,j,k) + vvv(i,j,k);
            flux_z(i,j,k) = flux_z(i,j,k) + www(i,j,k);
            
        end
    end
end

% f = f + tend2_adv;
tend = tend + tend_tmp;
f_bound = permute(f_bound(2:end-1,2:end-1,1:end-1),[3,2,1]);
% v = permute(v,[3,2,1]);
% u = permute(u,[3,2,1]);
% w = permute(w,[3,2,1]);

flux_x = permute(flux_x(2:end-1,2:end-1,1:end-1).*dx,[3,2,1]);
flux_y = permute(flux_y(2:end-1,2:end-1,1:end-1).*dy,[3,2,1]);
flux_z = permute(flux_z(2:end-1,2:end-1,1:end-1).*dz,[3,2,1]);
tend = permute(tend(2:end-1,2:end-1,1:end-1),[3,2,1]);

% u = u./dtdz;
% v = v./dtdz;
% w = w./dtdz;
dtn_inv = 1./dtn;
flux_x = flux_x.*dtn_inv;
flux_y = flux_y.*dtn_inv;
flux_z = flux_z.*dtn_inv;
tend = tend.*dtn_inv;
f_out = f_bound;

% f = f;
% uuu =  [];
% www =  [];
% vvv =  [];
% mx = [];
% mn = [];

end

%helping functions.

function out = andiff(x1,x2,a,b)
out = (abs(a)-a*a*b)*0.5*(x2-x1);
end

function out = across(x1,a1,a2)
out =0.03125*a1*a2*x1;
end
function out = pp(y)
out = max(0.,y);
end
function out = pn(y)
out =-min(0.,y);
end