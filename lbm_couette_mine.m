clc;
clear;
omega=1;
tau =1/omega;
xlen=101;
ylen=101;
kvisc = ((tau)-0.5)/3;
Q = 0.5;    %volumetric flow
uo=0.25; %top wall velocity

u=zeros(xlen,ylen);
v=zeros(xlen,ylen);
rho=ones(xlen,ylen);
rho_top = ones(xlen);
f=zeros(xlen,ylen);
x=1:xlen;
y=1:ylen;
n=1;
ini_iter=14000;
%Initializataion of distributions
      f(x,y,1) = (1/9).*(rho+3.*u+4.5.*u.*u-1.5.*(u.^2+v.^2));
      f(x,y,2) = (1/9).*(rho+3.*v+4.5.*v.*v-1.5*(u.^2+v.^2));
      f(x,y,3) = (1/9).*(rho-3.*u+4.5.*u.*u-1.5*(u.^2+v.^2));
      f(x,y,4) = (1/9).*(rho-3.*v+4.5.*v.*v-1.5*(u.^2+v.^2));

      f(x,y,5) = (1/36).*(rho+3.*(u+v)+4.5*(u+v).^2-1.5*(u.^2+v.^2));
      f(x,y,6) = (1/36).*(rho+3.*(-u+v)+4.5*(-u+v).^2-1.5*(u.^2+v.^2));
      f(x,y,7) = (1/36).*(rho-3.*(u+v)+4.5*(u+v).^2-1.5*(u.^2+v.^2));
      f(x,y,8) = (1/36).*(rho+3.*(u-v)+4.5*(u-v).^2-1.5*(u.^2+v.^2));

      f(x,y,9) = (4/9).*(rho-1.5*(u.^2+v.^2));
%Initialization completed
%Assign first fprop to equilibrium distributions.
      fstar=f;  
      feq=zeros(xlen,ylen,9);
      u(x,ylen-1)=uo;
      feq(x,ylen-1,1) = (1/9).*(rho(x,ylen-1)+3.*u(x,ylen-1)+4.5.*u(x,ylen-1).*u(x,ylen-1)-1.5.*(u(x,ylen-1).^2+v(x,ylen-1).^2));
      feq(x,ylen-1,2) = (1/9).*(rho(x,ylen-1)+3.*v(x,ylen-1)+4.5.*v(x,ylen-1).*v(x,ylen-1)-1.5.*(u(x,ylen-1).^2+v(x,ylen-1).^2));
      feq(x,ylen-1,3) = (1/9).*(rho(x,ylen-1)-3.*u(x,ylen-1)+4.5.*u(x,ylen-1).*u(x,ylen-1)-1.5.*(u(x,ylen-1).^2+v(x,ylen-1).^2));
      feq(x,ylen-1,4) = (1/9).*(rho(x,ylen-1)-3.*v(x,ylen-1)+4.5.*v(x,ylen-1).*v(x,ylen-1)-1.5.*(u(x,ylen-1).^2+v(x,ylen-1).^2));

      feq(x,ylen-1,5) = (1/36).*(rho(x,ylen-1)+3.*(u(x,ylen-1)+v(x,ylen-1))+4.5*(u(x,ylen-1)+v(x,ylen-1)).^2-1.5*(u(x,ylen-1).^2+v(x,ylen-1).^2));
      feq(x,ylen-1,6) = (1/36).*(rho(x,ylen-1)+3.*(-u(x,ylen-1)+v(x,ylen-1))+4.5*(-u(x,ylen-1)+v(x,ylen-1)).^2-1.5*(u(x,ylen-1).^2+v(x,ylen-1).^2));
      feq(x,ylen-1,7) = (1/36).*(rho(x,ylen-1)-3.*(u(x,ylen-1)+v(x,ylen-1))+4.5*(u(x,ylen-1)+v(x,ylen-1)).^2-1.5*(u(x,ylen-1).^2+v(x,ylen-1).^2));
      feq(x,ylen-1,8) = (1/36).*(rho(x,ylen-1)+3.*(u(x,ylen-1)-v(x,ylen-1))+4.5*(u(x,ylen-1)-v(x,ylen-1)).^2-1.5*(u(x,ylen-1).^2+v(x,ylen-1).^2));
      feq(x,ylen-1,9) = (4/9).*(rho(x,ylen-1)-1.5*(u(x,ylen-1).^2+v(x,ylen-1).^2));
for iter=1:ini_iter

      feq(x,ylen-1,7) = (1/36).*(rho(x,ylen-1)-3.*(u(x,ylen-1)+v(x,ylen-1))+4.5*(u(x,ylen-1)+v(x,ylen-1)).^2-1.5*(u(x,ylen-1).^2+v(x,ylen-1).^2));     
      feq(x,ylen-1,4) = (1/9).*(rho(x,ylen-1)-3.*v(x,ylen-1)+4.5.*v(x,ylen-1).*v(x,ylen-1)-1.5.*(u(x,ylen-1).^2+v(x,ylen-1).^2));
      feq(x,ylen-1,8) = (1/36).*(rho(x,ylen-1)+3.*(u(x,ylen-1)-v(x,ylen-1))+4.5*(u(x,ylen-1)-v(x,ylen-1)).^2-1.5*(u(x,ylen-1).^2+v(x,ylen-1).^2));
      
      feq(x,1:ylen-2,1) = (1/9).*(rho(x,1:ylen-2)+3.*u(x,1:ylen-2)+4.5.*u(x,1:ylen-2).*u(x,1:ylen-2)-1.5.*(u(x,1:ylen-2).^2+v(x,1:ylen-2).^2));
      feq(x,1:ylen-2,2) = (1/9).*(rho(x,1:ylen-2)+3.*v(x,1:ylen-2)+4.5.*v(x,1:ylen-2).*v(x,1:ylen-2)-1.5.*(u(x,1:ylen-2).^2+v(x,1:ylen-2).^2));
      feq(x,1:ylen-2,3) = (1/9).*(rho(x,1:ylen-2)-3.*u(x,1:ylen-2)+4.5.*u(x,1:ylen-2).*u(x,1:ylen-2)-1.5.*(u(x,1:ylen-2).^2+v(x,1:ylen-2).^2));
      feq(x,1:ylen-2,4) = (1/9).*(rho(x,1:ylen-2)-3.*v(x,1:ylen-2)+4.5.*v(x,1:ylen-2).*v(x,1:ylen-2)-1.5.*(u(x,1:ylen-2).^2+v(x,1:ylen-2).^2));

      feq(x,1:ylen-2,5) = (1/36).*(rho(x,1:ylen-2)+3.*(u(x,1:ylen-2)+v(x,1:ylen-2))+4.5*(u(x,1:ylen-2)+v(x,1:ylen-2)).^2-1.5*(u(x,1:ylen-2).^2+v(x,1:ylen-2).^2));
      feq(x,1:ylen-2,6) = (1/36).*(rho(x,1:ylen-2)+3.*(-u(x,1:ylen-2)+v(x,1:ylen-2))+4.5*(-u(x,1:ylen-2)+v(x,1:ylen-2)).^2-1.5*(u(x,1:ylen-2).^2+v(x,1:ylen-2).^2));
      feq(x,1:ylen-2,7) = (1/36).*(rho(x,1:ylen-2)-3.*(u(x,1:ylen-2)+v(x,1:ylen-2))+4.5*(u(x,1:ylen-2)+v(x,1:ylen-2)).^2-1.5*(u(x,1:ylen-2).^2+v(x,1:ylen-2).^2));
      feq(x,1:ylen-2,8) = (1/36).*(rho(x,1:ylen-2)+3.*(u(x,1:ylen-2)-v(x,1:ylen-2))+4.5*(u(x,1:ylen-2)-v(x,1:ylen-2)).^2-1.5*(u(x,1:ylen-2).^2+v(x,1:ylen-2).^2));
      feq(x,1:ylen-2,9) = (4/9)*(rho(x,1:ylen-2)-1.5*(u(x,1:ylen-2).^2+v(x,1:ylen-2).^2));
      
%propagating distributions fprop by applying kinetic equation.

      fstar = (1.0-omega).* f + omega.* feq;
      
            fstar(x,1,:) = f(x,1,:);
            fstar(x,ylen,:) = f(x,ylen,:);
      
            
 p=2:xlen;
 q=2:ylen;
 r=1:xlen-1;
 s=1:ylen-1;
 g1=f(x,y,1);
 g2=f(x,y,2);
 g3=f(x,y,3);
 g4=f(x,y,4);
 g5=f(x,y,5);
 g6=f(x,y,6);
 g7=f(x,y,7);
 g8=f(x,y,8);
 g9=f(x,y,9);
 m1=g1';
 m2=g2';
 m3=g3';
 m4=g4';
 m5=g5';
 m6=g6';
 m7=g7';
 m8=g8';
 m9=g9';
 h1=fstar(x,y,1);
 h2=fstar(x,y,2);
 h3=fstar(x,y,3);
 h4=fstar(x,y,4);
 h5=fstar(x,y,5);
 h6=fstar(x,y,6);
 h7=fstar(x,y,7);
 h8=fstar(x,y,8);
 h9=fstar(x,y,9);
 n1=h1';
 n2=h2';
 n3=h3';
 n4=h4';
 n5=h5';
 n6=h6';
 n7=h7';
 n8=h8';
 n9=h9';
 
%streaming
      m1(y,p) = n1(y,p-1);%c1
      m2(q,x) = n2(q-1,x);%c2
      m3(y,r) = n3(y,r+1);%c3  
      m4(s,x) = n4(s+1,x);%c4
      m5(q,p) = n5(q-1,p-1);%c5
      m6(q,r) = n6(q-1,r+1);%c6
      m7(s,r) = n7(s+1,r+1);%c7
      m8(s,p) = n8(s+1,p-1);%c8
      m9(y,x) = n9(y,x);%c9
      
%Complete Bounce Back Boundary Conditions

%1.Implementing Periodic BC.
m1(y,x(1))=n1(y,xlen);
m3(y,xlen)=n3(y,x(1));
m5(2:ylen,x(1))=n5(1:ylen-1,xlen);
m6(2:ylen,xlen)=n6(1:ylen-1,x(1));
m7(1:ylen-1,xlen)=n7(2:ylen,x(1));
m8(1:ylen-1,x(1))=n8(2:ylen,xlen);

%Reassigning values to f
f(x,y,1)=m1';
f(x,y,2)=m2'; 
f(x,y,3)=m3';
f(x,y,4)=m4';
f(x,y,5)=m5';
f(x,y,6)=m6';
f(x,y,7)=m7';
f(x,y,8)=m8';
f(x,y,9)=m9';

%2.Bounce Back Begins.

%South Boundary
f(x,1,1) = f(x,1,3);
f(x,1,5) = f(x,1,7);
f(x,1,8) = f(x,1,6);

for k=x %North boundary (Zou He moving wall boundary condition)
        rho_top(k)=f(k,ylen,9)+f(k,ylen,1)+f(k,ylen,3)+2*(f(k,ylen,2)+f(k,ylen,6)+f(k,ylen,5));
        f(k,ylen,4)=f(k,ylen,2);
        f(k,ylen,8)=f(k,ylen,6)+rho_top(k)*uo/6;
        f(k,ylen,7)=f(k,ylen,5)-rho_top(k)*uo/6;
 end


rho(x,y) = f(x,y,1)+f(x,y,2)+f(x,y,3)+f(x,y,4)+f(x,y,5)+f(x,y,6)+f(x,y,7)+f(x,y,8)+f(x,y,9);
u(x,1:ylen)=f(x,1:ylen,1)-f(x,1:ylen,3)+f(x,1:ylen,5)-f(x,1:ylen,6)-f(x,1:ylen,7)+f(x,1:ylen,8);
v(x,1:ylen)=f(x,1:ylen,2)-f(x,1:ylen,4)+f(x,1:ylen,5)+f(x,1:ylen,6)-f(x,1:ylen,7)-f(x,1:ylen,8);

uprofile(1:ylen)=sum(u(1:xlen,1:ylen))/xlen;

up = uprofile(1:ylen-1)/uo;
if mod(n,1000)==0
h = plot(up,0.01:1/(ylen-1):1);
set( h, 'Color', [0 0 0], 'LineWidth',2.0);
ylabel('y/L')
xlabel('u/uo')
title('Plot of TIme-lines')
hold on
end
n=n+1;
end
