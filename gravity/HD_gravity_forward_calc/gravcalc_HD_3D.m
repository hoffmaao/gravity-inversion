function gravcalc_HD_3D(LN)

% Gravity anomaly is referenced to the geoid whereas disturbance is 
% referenced to the ellipsoid. This code uses the observed gravity 
% disturbance because the elevation data are all referenced to WGS1984
% ellipsoid. 
% 
% This code takes about two minutes to run.
% 
% Input
%   LN = a character string specifying a PolarGAP flight line, P15, P16,
%        P17 or P18.
% 
% Output
%   A figure of the ice-surface and the bed topography and gravity disturbances

% Written by Atsuhiro Muto
% Dept. of Earth & Environ. Sci., Temple Univ.
% amuto@temple.edu
% Last updated March 4, 2021

%% Load data 
load('polargap_HD')
data = polargap_HD.(LN);

%%

load BM_5km
A = BM_5km;

dx = 5000;                          % prism dimension in x-direction
dy = 5000;                          % prism dimension in y-direction
xx = -850000:dx:150000;
yy = -450000:dy:150000;

mx = length(xx);                   % number of prisms in x-direction
my = length(yy);                   % number of prisms in y-direction

%%

% define density
rhoi = 900;
rhoc = 2670;
rho1 = [rhoi;rhoc];                 % density contrast between bedrock (2670 kg/m^3)
rho2 = rhoi-rhoc;                   % density contrast between bedrock (2670 kg/m^3)

dn = 150000;                        % distance for neighbouring prisms for gravity calculation

Pp = makep(mx,my,dx,dy);
Pp = [Pp(:,1) Pp(:,1)+dx Pp(:,2) Pp(:,2)+dy];
Pobs = [data(:,8)-xx(1)+dx/2 data(:,9)-yy(1)+dy/2];

% flight level
flevel = 4000;

[m01,m02] = msplit(-A(:,3:end),flevel);
m01 = [m01 ones(my*mx,1)*flevel+1e-3];

%%

% clear memory carried over from previous mex calls
clear mex

% tic
g = plouff3(Pobs,Pp,m01,m02,rho1,rho2,dx,dy,dn);
% toc

%% low-pass filter

hmax = 20000;
hmin = 10000;
glpf = fftfilt_1d(data(:,10),g,hmax,hmin)';
glpf = glpf*1e5;

%% Long-wavelength signal from GOCE

load boug_goce5000

S = size(boug_goce5000);
XI = reshape(boug_goce5000(:,:,1)',S(1)*S(2),1);
YI = reshape(boug_goce5000(:,:,2)',S(1)*S(2),1);
ZI = reshape(boug_goce5000(:,:,3)',S(1)*S(2),1);

F = scatteredInterpolant(XI,YI,ZI,'natural');
boug_goce = F(data(:,8),data(:,9));

gmodel = glpf+boug_goce;

%% Plotting

% Plot density model

% 2D Polygon for plotting
mm = mmodtwodpoly([data(:,10) -data(:,13) -data(:,11)]);
id = isnan(mm(:,2));
mm(id,:) = [];

col = ['c','c','c','c','c'];
figure
subplot(211);
title(LN)

npolygon = size(mm,2);
for ip=1:npolygon-1
    zp = mm(:,ip+1);          % z-position
    fill(mm(:,1)/1000,-1*zp,col(ip)); 
    hold on
end
title(LN)

xlim([min(data(:,10))/1000 max(data(:,10))/1000])
ylabel('Elevation (m)')
grid on

id = 1:3:length(data);

% Gravity anomaly
subplot(212)
plot(data(id,10)/1000,data(id,5),'ko');
hold all
plot(data(:,10)/1000,glpf,'-','LineWidth',2,'Color',[0 0.447 0.741]);
plot(data(:,10)/1000,boug_goce,'-','LineWidth',2,'Color',[0.466 0.674 0.188]);
plot(data(:,10)/1000,gmodel,'-','LineWidth',2,'Color',[0.850 0.325 0.098]);
title(LN)

xlim([min(data(:,10))/1000 max(data(:,10))/1000])
xlabel('Distance along the profile (km)')
ylabel('Gravity anomaly (miligals)')
grid on
legend('Observed gravity','Gravity effect of ice and topography',...
       'Gravity effect of the Moho','Sum of modelled gravity')