clc;
clear;
%PARAMETER SETTING STEP
xlen = 100;    %mesh_size in x direction
ylen = 100;    %mesh_size in y direction
uo=0.1;
kvisc=0.01;
Re=uo*ylen/kvisc; 
omega=1.0/(3.*kvisc+0.5); %relaxation factor value
tau = 1/omega;
convg = 1e-8 ; %convergence criteria - u velocity difference
rho = ones(xlen,ylen);
rho_top = ones(xlen);
jx= zeros(xlen,ylen);
jy= zeros(xlen,ylen);
f = zeros(xlen,ylen);
feq = zeros(xlen,ylen,9);
x = 1:xlen;
y = 1:ylen;
d = 1:9;
ini_iter = 100;

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

      %initialization completed



%running some initial iterations without convergence condition
for i=1:ini_iter

      feq(x,y,1) = (rho(x,y)./9).*(1+3.*u(x,y)+4.5.*u(x,y).*u(x,y)-1.5.*(u(x,y).^2+v(x,y).^2));
      feq(x,y,2) = (rho(x,y)./9).*(1+3.*v(x,y)+4.5.*v(x,y).*v(x,y)-1.5.*(u(x,y).^2+v(x,y).^2));
      feq(x,y,3) = (rho(x,y)./9).*(1-3.*u(x,y)+4.5.*u(x,y).*u(x,y)-1.5.*(u(x,y).^2+v(x,y).^2));
      feq(x,y,4) = (rho(x,y)./9).*(1-3.*v(x,y)+4.5.*v(x,y).*v(x,y)-1.5.*(u(x,y).^2+v(x,y).^2));

      feq(x,y,5) = (rho(x,y)./36).*(1+3.*(u(x,y)+v(x,y))+4.5*(u(x,y)+v(x,y)).^2-1.5*(u(x,y).^2+v(x,y).^2));
      feq(x,y,6) = (rho(x,y)./36).*(1+3.*(-u(x,y)+v(x,y))+4.5*(-u(x,y)+v(x,y)).^2-1.5*(u(x,y).^2+v(x,y).^2));
      feq(x,y,7) = (rho(x,y)./36).*(1-3.*(u(x,y)+v(x,y))+4.5*(u(x,y)+v(x,y)).^2-1.5*(u(x,y).^2+v(x,y).^2));
      feq(x,y,8) = (rho(x,y)./36).*(1+3.*(u(x,y)-v(x,y))+4.5*(u(x,y)-v(x,y)).^2-1.5*(u(x,y).^2+v(x,y).^2));
      feq(x,y,9) = ((4/9)*rho(x,y)).*(1-1.5*(u(x,y).^2+v(x,y).^2));
          
      f = (1.0-omega).* f + omega.* feq;
      %collision completed
           
         
%Propagation step... 
       p=xlen:-1:2;
       q=ylen:-1:2;
       r=1:xlen-1;
       s=1:ylen-1;
       f(p,y,1) = f(p-1,y,1);   %C1 right to left update
       f(x,q,2) = f(x,q-1,2);   %C2 top to bottom update
       f(r,y,3) = f(r+1,y,3);   %C3 left to right update
       f(x,s,4) = f(x,s+1,4);   %C4 bottom to top update
       f(p,q,5) = f(p-1,q-1,5);   %C5 top roght to bottom left update
       f(r,q,6) = f(r+1,q-1,6);   %C6 top left to bottom right update
       f(r,s,7) = f(r+1,s+1,7);   %C7 bottom left to top right update
       f(p,s,8) = f(p-1,s+1,8);   %C8 bottom right to top left update
       f(x,y,9) = f(x,y,9);   %C9 node update

%bounce back on west boundary
f(1,y,1) = f(1,y,3);
f(1,y,5) = f(1,y,7);
f(1,y,8) = f(1,y,6);

%bounce back on east boundary
f(xlen,y,3) = f(xlen,y,1);
f(xlen,y,7) = f(xlen,y,5);
f(xlen,y,6) = f(xlen,y,8);

%bounce back on south boundary
f(x,1,1) = f(x,1,3);
f(x,1,5) = f(x,1,7);
f(x,1,8) = f(x,1,6);


%moving lid on north boundary
    for k=2:xlen-1 %top wall except the corners
        rho_top(k)=f(k,ylen,9)+f(k,ylen,1)+f(k,ylen,3)+2*(f(k,ylen,2)+f(k,ylen,6)+f(k,ylen,5));
        f(k,ylen,4)=f(k,ylen,2);
        f(k,ylen,8)=f(k,ylen,6)+rho_top(k)*uo/6;
        f(k,ylen,7)=f(k,ylen,5)-rho_top(k)*uo/6; 
    end


rho(x,y) = f(x,y,1)+f(x,y,2)+f(x,y,3)+f(x,y,4)+f(x,y,5)+f(x,y,6)+f(x,y,7)+f(x,y,8)+f(x,y,9);
jx(x,y)=f(x,y,1)-f(x,y,3)+f(x,y,5)-f(x,y,6)-f(x,y,7)+f(x,y,8);%Distributions multiplied by lattice velocities in x-directions
jy(x,y)=f(x,y,2)-f(x,y,4)+f(x,y,5)+f(x,y,6)-f(x,y,7)-f(x,y,8);%Distributions multiplied by lattice velocities in y-directions
u = jx./rho;
v = jy./rho;
end

diff = 100000; %setting relatively higher value for the parameter of convergence
Upre = zeros(xlen,ylen);
n_iter = ini_iter;

while abs(diff) > convg %convergence condition

    
      feq(x,y,1) = (rho(x,y)./9).*(1+3.*u(x,y)+4.5.*u(x,y).*u(x,y)-1.5.*(u(x,y).^2+v(x,y).^2));
      feq(x,y,2) = (rho(x,y)./9).*(1+3.*v(x,y)+4.5.*v(x,y).*v(x,y)-1.5.*(u(x,y).^2+v(x,y).^2));
      feq(x,y,3) = (rho(x,y)./9).*(1-3.*u(x,y)+4.5.*u(x,y).*u(x,y)-1.5.*(u(x,y).^2+v(x,y).^2));
      feq(x,y,4) = (rho(x,y)./9).*(1-3.*v(x,y)+4.5.*v(x,y).*v(x,y)-1.5.*(u(x,y).^2+v(x,y).^2));

      feq(x,y,5) = (rho(x,y)./36).*(1+3.*(u(x,y)+v(x,y))+4.5*(u(x,y)+v(x,y)).^2-1.5*(u(x,y).^2+v(x,y).^2));
      feq(x,y,6) = (rho(x,y)./36).*(1+3.*(-u(x,y)+v(x,y))+4.5*(-u(x,y)+v(x,y)).^2-1.5*(u(x,y).^2+v(x,y).^2));
      feq(x,y,7) = (rho(x,y)./36).*(1-3.*(u(x,y)+v(x,y))+4.5*(u(x,y)+v(x,y)).^2-1.5*(u(x,y).^2+v(x,y).^2));
      feq(x,y,8) = (rho(x,y)./36).*(1+3.*(u(x,y)-v(x,y))+4.5*(u(x,y)-v(x,y)).^2-1.5*(u(x,y).^2+v(x,y).^2));
      feq(x,y,9) = ((4/9)*rho(x,y)).*(1-1.5*(u(x,y).^2+v(x,y).^2));
          
      f = (1.0-omega).* f + omega.* feq;
      %collision completed

      
%Propagation step... 
       p=xlen:-1:2;
       q=ylen:-1:2;
       r=1:xlen-1;
       s=1:ylen-1;
       f(p,y,1) = f(p-1,y,1);   %C1 right to left update
       f(x,q,2) = f(x,q-1,2);   %C2 top to bottom update
       f(r,y,3) = f(r+1,y,3);   %C3 left to right update
       f(x,s,4) = f(x,s+1,4);   %C4 bottom to top update
       f(p,q,5) = f(p-1,q-1,5);   %C5 top roght to bottom left update
       f(r,q,6) = f(r+1,q-1,6);   %C6 top left to bottom right update
       f(r,s,7) = f(r+1,s+1,7);   %C7 bottom left to top right update
       f(p,s,8) = f(p-1,s+1,8);   %C8 bottom right to top left update
       f(x,y,9) = f(x,y,9);   %C9 node update

%bounce back on west boundary
f(1,y,1) = f(1,y,3);
f(1,y,5) = f(1,y,7);
f(1,y,8) = f(1,y,6);

%bounce back on east boundary
f(xlen,y,3) = f(xlen,y,1);
f(xlen,y,7) = f(xlen,y,5);
f(xlen,y,6) = f(xlen,y,8);

%bounce back on south boundary
f(x,1,1) = f(x,1,3);
f(x,1,5) = f(x,1,7);
f(x,1,8) = f(x,1,6);


%moving lid on north boundary
    for k=2:xlen-1 %top wall except the corners
        rho_top(k)=f(k,ylen,9)+f(k,ylen,1)+f(k,ylen,3)+2*(f(k,ylen,2)+f(k,ylen,6)+f(k,ylen,5));
        f(k,ylen,4)=f(k,ylen,2);
        f(k,ylen,8)=f(k,ylen,6)+rho_top(k)*uo/6;
        f(k,ylen,7)=f(k,ylen,5)-rho_top(k)*uo/6; 
    end


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

fprintf('Total Number of iterations for convergence is %d . \n', n_iter);

%rotating velocity arrays for plotting
UU=rot90(u);
VV=rot90(v);
U=flipud(UU);
V=flipud(VV);
figure(1)
contourf(x,y,hypot(U,V),30)
xlabel('x')
ylabel('y')
title('Plot of velocity contours')

figure(2)
quiver(x,y,U,V)
xlabel('x')
ylabel('y')
title('Plot of velocity vectors')

figure(3)
[startx,starty] = meshgrid(0:5:xlen,0:5:ylen);
verts = stream2(x,y,U,V,startx,starty);
h = streamline(verts);
axis tight
set( h, 'Color', [0 0 0] )
xlabel('x')
ylabel('y')
title('Plot of streamlines')


figure(4)
VelProfile(1:ylen)=u(xlen/2,y)./uo;%Velocity profile in x-direction
VP = transpose(VelProfile);
plot(VelProfile,0.01:1/xlen:1)
xlabel('u/uo')
ylabel('y/L')
title('u velocity Profile at cavity centre')

figure(5)
VelProfile_v(1:xlen)=v(x,ylen/2)./uo;%Velocity profile in x-direction
VP_v = transpose(VelProfile_v);
plot(0.01:1/xlen:1,VelProfile_v)
xlabel('x/L')
ylabel('v/uo')
title('v velocity Profile at cavity centre')

figure(6)
g = streamslice(x,y,U,V,2);
set( g, 'Color', [0 0 0], 'LineWidth',2.0)
xlabel('x')
ylabel('y')
title('Plot of streamlines')