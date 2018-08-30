%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Gauss-Laguerre Acoustic Beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% Physical Parametes
   dx = 1.5;    % [mm] voxel size
    f = 1E6;    % 1E6[Hz] = [MHz] freq.
   c0 = 1500;   % [m/s] speed of sound in water
lambda= c0/f;   % [m] wave length
omega = 2*pi*f; % [rads] radial freq.
   mm = 1E-3;   % [m] conversion factor
   t0 = 1E-6;   % [s] period of an acoustical signal

% Model Parameters
w0 = 10*lambda; % width of the beam at z=0;
 k = 100;       % wave number;
zR = k*w0^2/2;  % Rayleigh distance;
 l =-1;         % topological charge;
 n = abs(l)+0;	% radial index; n=|l|,|l|+2,|l|+4 ...
 D = 1000;      % is a constant for normalization;
   
% Discrete domain
xc=-40.5:dx:40.5; %[mm]
yc=-40.5:dx:40.5; %[mm]
zc = 100; %[mm] distance from the transduce device
[X,Y,Z] = meshgrid(xc*mm,yc*mm,zc*mm); %[m]
[TH,R,Z] = cart2pol(X,Y,Z);

% Discrete planes for debuging (do not delete!)
% only RZ-plane
%[R,Z] = meshgrid(-40.5*mm:dx:40.5*mm,0*mm:1*mm:100*mm);
% only YX-plane
%[X,Y] = meshgrid(-40.5*mm:dx:40.5*mm,-40.5*mm:dx:40.5*mm);
%R=hypot(X,Y); X(R>11*mm)=NaN; Y(R>11*mm)=NaN; % Trimp experimental data
%[TH,R] = cart2pol(X,Y);

% Analytical functions
w = @(z) w0*sqrt(1+(z/zR).^2);
A = @(r,z) (sqrt(2)*r./w(z)).^(abs(l)).*LaguerreL((n-abs(l))/2,abs(l),2*r.^2./w(z).^2);
G = @(r,z) D./sqrt(1+z.^2/zR^2).*exp(r.^2./w(z).^2).*exp(-0.5i*k*(z.*r.^2)./(z.^2+zR^2));
PHI = @(th) exp(1i*l*th);
PSI = @(z) exp(-1i*(n+1)*atan(z/zR));
P = @(th,r,z,t) G(r,z).*A(r,z).*PHI(th).*exp(1i*(k*z-omega*t)).*PSI(z);

% Deguging routines (do not delete!)
%surf(angle(PHI(TH))); % ok
%surf(R,Z,A(R,Z)); % ok
%surf(R,Z,abs(G(R,Z))); %ok
%surf(R,Z,abs(G(R,Z).*A(R,Z))) %ok

% Compute profile for a seleted time 't':
p1=P(TH,R,Z,0);
p2=P(TH,R,Z,t0/4);
p3=P(TH,R,Z,t0/2);
p4=P(TH,R,Z,3*t0/4);

% Plot a single slice of the presure profile
figure(1); fontsize=12;
set(gcf,'position',[100,100,600,800])
subplot(421), imagesc(xc,yc,flipud(abs(p1))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('magnitude','interpreter','latex','fontsize',fontsize);
subplot(422), imagesc(xc,yc,flipud(angle(p1))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('phase at $t_0$','interpreter','latex','fontsize',fontsize);
subplot(423), imagesc(xc,yc,flipud(abs(p2))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('magnitude','interpreter','latex','fontsize',fontsize);
subplot(424), imagesc(xc,yc,flipud(angle(p2))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('phase at $t_0$/4','interpreter','latex','fontsize',fontsize);
subplot(425), imagesc(xc,yc,flipud(abs(p3))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('magnitude','interpreter','latex','fontsize',fontsize);
subplot(426), imagesc(xc,yc,flipud(angle(p3))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('phase at $t_0$/2','interpreter','latex','fontsize',fontsize);
subplot(427), imagesc(xc,yc,flipud(abs(p4))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('magnitude','interpreter','latex','fontsize',fontsize);
subplot(428), imagesc(xc,yc,flipud(angle(p4))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('phase at 3/4 $t_0$','interpreter','latex','fontsize',fontsize);
