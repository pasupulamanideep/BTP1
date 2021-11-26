clc;
clear;
%PARAMETER SETTING
omega = 0.2;      %relaxation factor value
tau = 1/omega;
xlen = 100;     %mesh_size in x direction
ylen = 100;    %mesh_size in y direction
Q = 0.5;    %volumetric flow
kvisc = ((1/omega)-0.5)/3; %kinematic viscosity
forcing = (1.5*Q*kvisc/ylen^3);

rho = ones(xlen,ylen);
jx= zeros(xlen,ylen);
jy= zeros(xlen,ylen);
f = zeros(xlen,ylen);
x = 1:xlen;
y = 1:ylen;
ini_iter = 1000;

%initialization
u=jx./rho;
v=jy./rho;
      f(x,y,1) = (rho./9)*(1+3.*u+4.5.*u.*u-1.5.*(u.^2+v.^2));   %for 'plus' corner points of the node
      f(x,y,2) = (rho./9)*(1+3.*v+4.5.*v.*v-1.5*(u.^2+v.^2));
      f(x,y,3) = (rho./9)*(1-3.*u+4.5.*u.*u-1.5*(u.^2+v.^2));
      f(x,y,4) = (rho./9)*(1-3.*v+4.5.*v.*v-1.5*(u.^2+v.^2));

      f(x,y,5) = (rho./36)*(1+3.*(u+v)+4.5*(u+v).^2-1.5*(u.^2+v.^2)); %for 'cross' corner points of the node
      f(x,y,6) = (rho./36)*(1+3.*(-u+v)+4.5*(-u+v).^2-1.5*(u.^2+v.^2));
      f(x,y,7) = (rho./36)*(1-3.*(u+v)+4.5*(u+v).^2-1.5*(u.^2+v.^2));
      f(x,y,8) = (rho./36)*(1+3.*(u-v)+4.5*(u-v).^2-1.5*(u.^2+v.^2));
      
      f(x,y,9) = (4/9)*rho.*(1-1.5*(u.^2+v.^2)); %for 'node centre'


fstar = f;
feq = zeros(xlen,ylen,9);

for i=1:ini_iter
    u=jx./rho;
    v=jy./rho;
      feq(x,2:ylen-1,1) = (rho(x,2:ylen-1)./9).*(1+3.*u(x,2:ylen-1)+4.5.*u(x,2:ylen-1).*u(x,2:ylen-1)-1.5.*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
      feq(x,2:ylen-1,2) = (rho(x,2:ylen-1)./9).*(1+3.*v(x,2:ylen-1)+4.5.*v(x,2:ylen-1).*v(x,2:ylen-1)-1.5.*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
      feq(x,2:ylen-1,3) = (rho(x,2:ylen-1)./9).*(1-3.*u(x,2:ylen-1)+4.5.*u(x,2:ylen-1).*u(x,2:ylen-1)-1.5.*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
      feq(x,2:ylen-1,4) = (rho(x,2:ylen-1)./9).*(1-3.*v(x,2:ylen-1)+4.5.*v(x,2:ylen-1).*v(x,2:ylen-1)-1.5.*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));

      feq(x,2:ylen-1,5) = (rho(x,2:ylen-1)./36).*(1+3.*(u(x,2:ylen-1)+v(x,2:ylen-1))+4.5*(u(x,2:ylen-1)+v(x,2:ylen-1)).^2-1.5*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
      feq(x,2:ylen-1,6) = (rho(x,2:ylen-1)./36).*(1+3.*(-u(x,2:ylen-1)+v(x,2:ylen-1))+4.5*(-u(x,2:ylen-1)+v(x,2:ylen-1)).^2-1.5*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
      feq(x,2:ylen-1,7) = (rho(x,2:ylen-1)./36).*(1-3.*(u(x,2:ylen-1)+v(x,2:ylen-1))+4.5*(u(x,2:ylen-1)+v(x,2:ylen-1)).^2-1.5*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
      feq(x,2:ylen-1,8) = (rho(x,2:ylen-1)./36).*(1+3.*(u(x,2:ylen-1)-v(x,2:ylen-1))+4.5*(u(x,2:ylen-1)-v(x,2:ylen-1)).^2-1.5*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
      feq(x,2:ylen-1,9) = ((4/9)*rho(x,2:ylen-1)).*(1-1.5*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
          
      fstar = (1.0-omega).* f + omega.* feq;
      %collision completed
      
      %forcing from pressure difference in poiseulle flow
      force =forcing/6;
      fstar(x,2:ylen-1,1) = fstar(x,2:ylen-1,1) + force;
      fstar(x,2:ylen-1,3) = fstar(x,2:ylen-1,3) - force;
      fstar(x,2:ylen-1,5) = fstar(x,2:ylen-1,5) + force;
      fstar(x,2:ylen-1,6) = fstar(x,2:ylen-1,6) - force;
      fstar(x,2:ylen-1,7) = fstar(x,2:ylen-1,7) - force;
      fstar(x,2:ylen-1,8) = fstar(x,2:ylen-1,8) + force;
      
      fstar(x,1,:) = f(x,1,:);
      fstar(x,ylen,:) = f(x,ylen,:);
 
      
       
         
%Propagation step... 
       p=2:xlen;
       q=2:ylen;
       r=1:xlen-1;
       s=1:ylen-1;
       f(p,y,1) = fstar(p-1,y,1);   %C1
       f(x,q,2) = fstar(x,q-1,2);   %C2
       f(r,y,3) = fstar(r+1,y,3);   %C3
       f(x,s,4) = fstar(x,s+1,4);   %C4
       f(p,q,5) = fstar(p-1,q-1,5);   %C5
       f(r,q,6) = fstar(r+1,q-1,6);   %C6
       f(r,s,7) = fstar(r+1,s+1,7);   %C7
       f(p,s,8) = fstar(p-1,s+1,8);   %C8
       f(x,y,9) = fstar(x,y,9);   %C9
           
 %Complete Bounce Back Boundary Conditions
   
   %1.Periodic BC.
f(x(1),y,1)=fstar(xlen,y,1);
f(xlen,y,3)=fstar(x(1),y,3);
f(x(1),2:ylen,5)=fstar(xlen,1:ylen-1,5);
f(xlen,2:ylen,6)=fstar(x(1),1:ylen-1,6);
f(xlen,1:ylen-1,7)=fstar(x(1),2:ylen,7);
f(x(1),1:ylen-1,8)=fstar(xlen,2:ylen,8);

   %2.Bounce Back Begins.
temp=f(1:xlen,1,2);
f(1:xlen,1,2)=f(1:xlen,1,4);
f(1:xlen,1,4)=temp;

temp=f(1:xlen,ylen,4);
f(1:xlen,ylen,4)=f(1:xlen,ylen,2);
f(1:xlen,ylen,4)=temp;

temp=f(1:xlen,1,5);
f(1:xlen,1,5)=f(1:xlen,1,7);
f(1:xlen,1,7)=temp;

temp=f(1:xlen,ylen,7);
f(1:xlen,ylen,7)=f(1:xlen,ylen,5);
f(1:xlen,ylen,5)=temp;

temp=f(1:xlen,1,6);
f(1:xlen,1,6)=f(1:xlen,1,8);
f(1:xlen,1,8)=temp;

temp=f(1:xlen,ylen,8);
f(1:xlen,ylen,8)=f(1:xlen,ylen,6);
f(1:xlen,ylen,6)=temp;

rho(x,y) = f(x,y,1)+f(x,y,2)+f(x,y,3)+f(x,y,4)+f(x,y,5)+f(x,y,6)+f(x,y,7)+f(x,y,8)+f(x,y,9);
jx(x,y)=f(x,y,1)-f(x,y,3)+f(x,y,5)-f(x,y,6)-f(x,y,7)+f(x,y,8);%Distributions multiplied by lattice velocities in x-directions
jy(x,y)=f(x,y,2)-f(x,y,4)+f(x,y,5)+f(x,y,6)-f(x,y,7)-f(x,y,8);%Distributions multiplied by lattice velocities in y-directions


end
%initialization completed


diff = 100000;
Upre = zeros(xlen,ylen);
n_iter = ini_iter;
while abs(diff) > 1e-12
      feq(x,2:ylen-1,1) = (rho(x,2:ylen-1)./9).*(1+3.*u(x,2:ylen-1)+4.5.*u(x,2:ylen-1).*u(x,2:ylen-1)-1.5.*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
      feq(x,2:ylen-1,2) = (rho(x,2:ylen-1)./9).*(1+3.*v(x,2:ylen-1)+4.5.*v(x,2:ylen-1).*v(x,2:ylen-1)-1.5.*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
      feq(x,2:ylen-1,3) = (rho(x,2:ylen-1)./9).*(1-3.*u(x,2:ylen-1)+4.5.*u(x,2:ylen-1).*u(x,2:ylen-1)-1.5.*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
      feq(x,2:ylen-1,4) = (rho(x,2:ylen-1)./9).*(1-3.*v(x,2:ylen-1)+4.5.*v(x,2:ylen-1).*v(x,2:ylen-1)-1.5.*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));

      feq(x,2:ylen-1,5) = (rho(x,2:ylen-1)./36).*(1+3.*(u(x,2:ylen-1)+v(x,2:ylen-1))+4.5*(u(x,2:ylen-1)+v(x,2:ylen-1)).^2-1.5*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
      feq(x,2:ylen-1,6) = (rho(x,2:ylen-1)./36).*(1+3.*(-u(x,2:ylen-1)+v(x,2:ylen-1))+4.5*(-u(x,2:ylen-1)+v(x,2:ylen-1)).^2-1.5*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
      feq(x,2:ylen-1,7) = (rho(x,2:ylen-1)./36).*(1-3.*(u(x,2:ylen-1)+v(x,2:ylen-1))+4.5*(u(x,2:ylen-1)+v(x,2:ylen-1)).^2-1.5*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
      feq(x,2:ylen-1,8) = (rho(x,2:ylen-1)./36).*(1+3.*(u(x,2:ylen-1)-v(x,2:ylen-1))+4.5*(u(x,2:ylen-1)-v(x,2:ylen-1)).^2-1.5*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
      feq(x,2:ylen-1,9) = ((4/9)*rho(x,2:ylen-1)).*(1-1.5*(u(x,2:ylen-1).^2+v(x,2:ylen-1).^2));
          
      fstar = (1.0-omega).* f + omega.* feq;
      %collision completed
      
      %forcing from pressure difference in poiseulle flow
      force =forcing/6;
      fstar(x,2:ylen-1,1) = fstar(x,2:ylen-1,1) + force;
      fstar(x,2:ylen-1,3) = fstar(x,2:ylen-1,3) - force;
      fstar(x,2:ylen-1,5) = fstar(x,2:ylen-1,5) + force;
      fstar(x,2:ylen-1,6) = fstar(x,2:ylen-1,6) - force;
      fstar(x,2:ylen-1,7) = fstar(x,2:ylen-1,7) - force;
      fstar(x,2:ylen-1,8) = fstar(x,2:ylen-1,8) + force;
      
      fstar(x,1,:) = f(x,1,:);
      fstar(x,ylen,:) = f(x,ylen,:);
 
      
       
         
%Propagation step... 
       p=2:xlen;
       q=2:ylen;
       r=1:xlen-1;
       s=1:ylen-1;
       f(p,y,1) = fstar(p-1,y,1);   %C1
       f(x,q,2) = fstar(x,q-1,2);   %C2
       f(r,y,3) = fstar(r+1,y,3);   %C3
       f(x,s,4) = fstar(x,s+1,4);   %C4
       f(p,q,5) = fstar(p-1,q-1,5);   %C5
       f(r,q,6) = fstar(r+1,q-1,6);   %C6
       f(r,s,7) = fstar(r+1,s+1,7);   %C7
       f(p,s,8) = fstar(p-1,s+1,8);   %C8
       f(x,y,9) = fstar(x,y,9);   %C9
           
 %Complete Bounce Back Boundary Conditions
   
   %1.Periodic BC.
f(x(1),y,1)=fstar(xlen,y,1);
f(xlen,y,3)=fstar(x(1),y,3);
f(x(1),2:ylen,5)=fstar(xlen,1:ylen-1,5);
f(xlen,2:ylen,6)=fstar(x(1),1:ylen-1,6);
f(xlen,1:ylen-1,7)=fstar(x(1),2:ylen,7);
f(x(1),1:ylen-1,8)=fstar(xlen,2:ylen,8);

   %2.Bounce Back Begins.
temp=f(1:xlen,1,2);
f(1:xlen,1,2)=f(1:xlen,1,4);
f(1:xlen,1,4)=temp;

temp=f(1:xlen,ylen,4);
f(1:xlen,ylen,4)=f(1:xlen,ylen,2);
f(1:xlen,ylen,4)=temp;

temp=f(1:xlen,1,5);
f(1:xlen,1,5)=f(1:xlen,1,7);
f(1:xlen,1,7)=temp;

temp=f(1:xlen,ylen,7);
f(1:xlen,ylen,7)=f(1:xlen,ylen,5);
f(1:xlen,ylen,5)=temp;

temp=f(1:xlen,1,6);
f(1:xlen,1,6)=f(1:xlen,1,8);
f(1:xlen,1,8)=temp;

temp=f(1:xlen,ylen,8);
f(1:xlen,ylen,8)=f(1:xlen,ylen,6);
f(1:xlen,ylen,6)=temp;

rho(x,y) = f(x,y,1)+f(x,y,2)+f(x,y,3)+f(x,y,4)+f(x,y,5)+f(x,y,6)+f(x,y,7)+f(x,y,8)+f(x,y,9);
jx(x,y)=f(x,y,1)-f(x,y,3)+f(x,y,5)-f(x,y,6)-f(x,y,7)+f(x,y,8);%Distributions multiplied by lattice velocities in x-directions
jy(x,y)=f(x,y,2)-f(x,y,4)+f(x,y,5)+f(x,y,6)-f(x,y,7)-f(x,y,8);%Distributions multiplied by lattice velocities in y-directions

u = jx./rho;
v = jy./rho;

%calculating the minimum difference in u velocities for convergence condition
Udiff = u - Upre;
diff = min(Udiff,[],'all');
Upre = u;
n_iter=n_iter+1;

end
VelProfile(1:ylen)=u(xlen/2,y)./u(xlen/2,ylen/2);%Velocity profile in x-direction
h = plot(VelProfile,0:1/ylen:0.99);
set (h,'Color',[0 0 0],'LineStyle','-','LineWidth',2.5);
xlabel('u/umax');
ylabel('y/L');
title('Plot of Timeline');